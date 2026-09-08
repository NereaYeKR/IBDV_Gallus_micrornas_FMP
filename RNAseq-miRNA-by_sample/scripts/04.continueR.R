#!/usr/bin/env Rscript
# ==============================================================================
# miRNA Differential Expression, Target Annotation and GSEA Dotplot Pipeline
#
# Description:
#   End-to-end pipeline that:
#     1. Runs DESeq2 differential expression analysis on miRNA count data.
#     2. Filters significant miRNAs and predicts their target genes using
#        miRDB, annotating targets with gene symbols/Entrez IDs via
#        org.Gg.eg.db (Gallus gallus).
#     3. Integrates mRNA GO/KEGG functional enrichment results
#        (clusterProfiler GSEA objects) with the miRNA target predictions.
#     4. Produces summarized top GO/KEGG tables enriched with miRNA info.
#     5. Generates GSEA dotplots (Activated vs Suppressed) for every
#        GO/KEGG comparison and saves them as PNG files.
#
# Requirements:
#   - A "config.sh" that exports the DESEQ2_CONTRASTS environment variable,
#     e.g.: export DESEQ2_CONTRASTS="contrastName,factor,levelA,levelB;..."
#   - Input files (relative to the project root):
#       data/SampleInfo.tab
#       results/count_featurecounts/mat_count_fc_clean.txt
#       data/miRDB_v6.0_prediction_result.txt
#       data/mrna_input/GO&KEGG_cProf/*.RData (or .RDS)
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 0: LOAD REQUIRED LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)          # Differential expression analysis
  library(dplyr)           # Data manipulation
  library(tibble)          # Data frame utilities (rownames_to_column, tribble)
  library(tidyr)           # Data tidying utilities (separate_rows)
  library(AnnotationDbi)   # Gene annotation queries
  library(org.Gg.eg.db)    # Gallus gallus (chicken) annotation database
  library(clusterProfiler) # Functional enrichment analysis
  library(enrichplot)      # Visualization of enrichment results
  library(pathview)        # Pathway visualization
  library(ggplot2)         # Plotting
  library(writexl)
})

# ------------------------------------------------------------------------------
# STEP 1: DEFINE PROJECT PATHS AND LOAD INPUT DATA
# ------------------------------------------------------------------------------
cat("\n[1] LOADING DATA...\n")

base_project <- getwd()
results_dir  <- file.path(base_project, "results")
data_dir     <- file.path(base_project, "data")

colData_path    <- file.path(data_dir, "SampleInfo.tab")
mat_counts_path <- file.path(results_dir, "count_featurecounts", "mat_count_fc_clean.txt")
mirdb_path      <- file.path(data_dir, "miRDB_v6.0_prediction_result.txt")

for (required_file in c(colData_path, mat_counts_path, mirdb_path)) {
  if (!file.exists(required_file)) {
    stop("Error: required input file not found: ", required_file)
  }
}

colData <- read.table(
  colData_path,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = TRUE,
  check.names      = FALSE
)

mat_counts <- read.table(
  mat_counts_path,
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

rownames(colData) <- gsub("\\.", "-", colData$file)
mat_counts <- mat_counts[, rownames(colData)]

if (!all(colnames(mat_counts) == rownames(colData))) {
  stop("Error: column names in mat_counts do not match colData.")
}

# ------------------------------------------------------------------------------
# OPTIONAL STEP 1.5: CUSTOM GROUPING (status) AND FILTERING MATRICES
# ------------------------------------------------------------------------------
# NOTE: If you only want to run DESeq2 by grouping samples into 'Control' and 
# 'Treatment' (status), uncomment the block below before running STEP 2.
# 
# cat("\n[1.5] CREATING 'status' VARIABLE AND FILTERING SAMPLES...\n")
# 
# colData <- colData %>%
#   mutate(status = case_when(
#     sample %in% c("DF.1_uninfected", "DF.1PC_uninfected") ~ "Control",
#     sample %in% c("DF.1_IBDV", "DF.1PC_IBDV")             ~ "Treatment",
#     TRUE ~ NA_character_ # Discard samples that are not part of this analysis
#   ))
# 
# colData <- colData %>% filter(!is.na(status))
# colData$status <- factor(colData$status, levels = c("Control", "Treatment"))
# mat_counts <- mat_counts[, rownames(colData)]
# 
# if (ncol(mat_counts) == 0) {
#   stop("Error: No samples matched the Control/Treatment criteria in colData.")
# }
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# STEP 2: BUILD THE DESeq2 DATASET AND RUN THE MODEL
# ------------------------------------------------------------------------------
cat("\n[2] PREPARING DESeq2 OBJECT...\n")

# NOTE: If you uncommented STEP 1.5 above to use 'status', you MUST change 
# `design = ~ sample` to `design = ~ status` in the function below.
dds <- DESeqDataSetFromMatrix(
  countData = round(mat_counts),
  colData   = colData,
  design    = ~ sample 
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# OPTIONAL NOTE FOR DESEQ2 'STATUS' ANALYSIS:
# If you only wanted to obtain the differential expression results grouped by 
# 'status', you could technically stop the script after STEP 4. The rest of the 
# pipeline integrates with functional enrichment downstream.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# STEP 3: READ CONTRAST DEFINITIONS FROM ENVIRONMENT CONFIGURATION
# ------------------------------------------------------------------------------
cat("\n[3] READING CONTRASTS FROM config.sh...\n")

contrast_str <- Sys.getenv("DESEQ2_CONTRASTS")

if (contrast_str == "") {
  stop("Error: DESEQ2_CONTRASTS variable is not defined in config.sh")
}

contrasts_list <- list()
pairs <- strsplit(contrast_str, ";")[[1]]

for (p in pairs) {
  parts <- strsplit(p, ",")[[1]]
  if (length(parts) != 4) {
    stop("Error: malformed contrast definition: '", p,
         "'. Expected format 'name,factor,levelA,levelB'.")
  }
  contrast_name <- parts[1]
  contrasts_list[[contrast_name]] <- c(parts[2], parts[3], parts[4])
}

cat("The following contrasts will be analyzed:\n")
print(names(contrasts_list))

# ------------------------------------------------------------------------------
# STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS AND SIGNIFICANCE FILTERING
#         (|log2FC| >= 1, padj < 0.05)
# ------------------------------------------------------------------------------
cat("\n[4] DESeq2 ANALYSIS AND FILTERING (|log2FC| >= 1, padj < 0.05)...\n")

deseq_output_dir <- file.path(results_dir, "deseq2_results")
dir.create(deseq_output_dir, recursive = TRUE, showWarnings = FALSE)

for (contrast_name in names(contrasts_list)) {
  
  res <- results(dds, contrast = contrasts_list[[contrast_name]])
  
  df_filt <- as.data.frame(res) %>%
    rownames_to_column("miRNA") %>%
    filter(
      padj < 0.05 &
        abs(log2FoldChange) >= 1
    )
  
  write.csv(
    df_filt,
    file.path(deseq_output_dir, paste0(contrast_name, "_L2FC_1.csv")),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# STEP 5: miRNA TARGET PREDICTION AND GENE ANNOTATION
# ------------------------------------------------------------------------------
cat("\n[5] miRNA TARGET PREDICTION AND ANNOTATION...\n")

chicken_targets <- read.table(mirdb_path, sep = "\t", header = FALSE)
colnames(chicken_targets)[1:3] <- c("miRNA", "gene", "score")
chicken_targets <- chicken_targets[grepl("^gga-", chicken_targets$miRNA), ]
chicken_targets$score <- as.numeric(as.character(chicken_targets$score))

targets_output_dir <- file.path(results_dir, "targets_annotated")
dir.create(targets_output_dir, recursive = TRUE, showWarnings = FALSE)

for (contrast_name in names(contrasts_list)) {
  
  de_results_file <- file.path(deseq_output_dir, paste0(contrast_name, "_L2FC_1.csv"))
  df_dea <- read.csv(de_results_file)
  
  if (nrow(df_dea) == 0) next
  
  rdb <- merge(df_dea, chicken_targets, by = "miRNA") %>%
    filter(score >= 80)
  
  if (nrow(rdb) == 0) next
  
  anno <- AnnotationDbi::select(
    org.Gg.eg.db,
    keys    = as.character(unique(rdb$gene)),
    columns = c("SYMBOL", "ENTREZID", "GENENAME"),
    keytype = "REFSEQ"
  ) %>%
    distinct(REFSEQ, .keep_all = TRUE)
  
  targets_annotated <- rdb %>%
    left_join(anno, by = c("gene" = "REFSEQ")) %>%
    rename(description = GENENAME)
  
  write.csv(
    targets_annotated,
    file.path(targets_output_dir, paste0(contrast_name, "_targets_final80.csv")),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# STEP 6: CLEAN ANNOTATED TARGETS (REMOVE ENTRIES WITHOUT ENTREZ ID)
# ------------------------------------------------------------------------------
cat("\n[6] CLEANING ANNOTATED TARGETS (REMOVING NA ENTREZID)...\n")

clean_targets <- function(df) {
  target_cols <- c(
    "miRNA",
    "log2FoldChange",
    "score",
    "SYMBOL",
    "ENTREZID",
    "description"
  )
  
  df %>%
    select(any_of(target_cols)) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(miRNA, ENTREZID, .keep_all = TRUE) %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    rename(
      L2FC_miRNA = log2FoldChange,
      core_enrichment = ENTREZID
    )
}

targets_noNA_dir <- file.path(results_dir, "targets_annotated_no_NA")
dir.create(targets_noNA_dir, recursive = TRUE, showWarnings = FALSE)

targets_files <- list.files(
  targets_output_dir,
  pattern    = "_targets_final80\\.csv$",
  full.names = TRUE
)

for (file in targets_files) {
  df  <- clean_targets(read.csv(file))
  out <- sub("_targets_final80", "_noNA", basename(file))
  
  write.csv(
    df,
    file.path(targets_noNA_dir, out),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# STEP 7: EXPAND CORE_ENRICHMENT IN ORIGINAL GO/KEGG OBJECTS
# ------------------------------------------------------------------------------
cat("\n[7] EXPANDING CORE_ENRICHMENT IN GO/KEGG OBJECTS...\n")

go_kegg_input_dir    <- file.path(data_dir, "mrna_input", "GO&KEGG_cProf")
go_kegg_expanded_dir <- file.path(data_dir, "mrna_input", "GO&KEGG_coreEnrichmentExpanded")

dir.create(go_kegg_expanded_dir, recursive = TRUE, showWarnings = FALSE)

enrich_files <- list.files(
  go_kegg_input_dir,
  pattern     = "\\.(RData|RDS|rds)$",
  full.names  = TRUE,
  ignore.case = TRUE
)

for (input_file in enrich_files) {
  
  object_name <- tools::file_path_sans_ext(basename(input_file))
  obj <- readRDS(input_file)
  
  obj@result <- obj@result %>%
    separate_rows(
      core_enrichment,
      sep = "/"
    )
  
  saveRDS(
    obj,
    file.path(go_kegg_expanded_dir, paste0(object_name, ".RDS"))
  )
}

# ------------------------------------------------------------------------------
# STEP 8: INTEGRATE mRNA ENRICHMENT WITH miRNA TARGETS
# ------------------------------------------------------------------------------
cat("\n[8] MERGING mRNA ENRICHMENT WITH miRNA TARGETS...\n")

rnam_mirna_all_dir <- file.path(results_dir, "RNAm&miRNA_all")
dir.create(rnam_mirna_all_dir, recursive = TRUE, showWarnings = FALSE)

targets_list <- list()
target_csv_files <- list.files(
  targets_noNA_dir,
  pattern    = "\\.csv$",
  full.names = TRUE
)

for (input_file in target_csv_files) {
  object_name <- tools::file_path_sans_ext(basename(input_file))
  targets_list[[object_name]] <- read.csv(input_file)
}

comparisons <- do.call(rbind, lapply(names(contrasts_list), function(contrast_name) {
  data.frame(
    csv  = paste0(contrast_name, "_noNA"),
    go   = paste0("GO_sample_", contrast_name),
    kegg = paste0("KEGG_sample_", contrast_name),
    stringsAsFactors = FALSE
  )
  
}))

for (i in seq_len(nrow(comparisons))) {
  
  csv_name <- comparisons$csv[i]
  
  if (!csv_name %in% names(targets_list)) {
    cat(paste0("  - Warning: targets file '", csv_name, "' not found. Skipping.\n"))
    next
  }
  
  targets_table <- targets_list[[csv_name]] %>%
    mutate(core_enrichment = as.character(core_enrichment))
  
  for (go_kegg_name in c(comparisons$go[i], comparisons$kegg[i])) {
    
    expanded_file <- file.path(go_kegg_expanded_dir, paste0(go_kegg_name, ".RDS"))
    
    if (!file.exists(expanded_file)) {
      alt_name <- sub("_sample_", "_status_", go_kegg_name)
      alt_file <- file.path(go_kegg_expanded_dir, paste0(alt_name, ".RDS"))
      if (file.exists(alt_file)) expanded_file <- alt_file
    }
    
    if (file.exists(expanded_file)) {
      obj <- readRDS(expanded_file)
      
      obj@result <- obj@result %>%
        mutate(core_enrichment = as.character(core_enrichment)) %>%
        left_join(
          targets_table,
          by = "core_enrichment",
          relationship = "many-to-many"
        )
      
      saveRDS(
        obj,
        file.path(rnam_mirna_all_dir, paste0(go_kegg_name, "_RNAm_miRNA.RDS"))
      )
    } else {
      cat(paste0("  - Warning: expanded GO/KEGG file for '", go_kegg_name, "' not found. Skipping.\n"))
    }
  }
}

# ------------------------------------------------------------------------------
# STEP 9: FILTER MERGED RESULTS (REMOVE ROWS WITH NA SCORE)
# ------------------------------------------------------------------------------
cat("\n[9] FILTERING MERGED RESULTS (REMOVING NA SCORES)...\n")

rnam_mirna_nona_dir <- file.path(results_dir, "RNAm&miRNA_no_NA")
dir.create(rnam_mirna_nona_dir, recursive = TRUE, showWarnings = FALSE)

merged_files <- list.files(
  rnam_mirna_all_dir,
  pattern    = "\\.RDS$",
  full.names = TRUE
)

for (input_file in merged_files) {
  
  obj <- readRDS(input_file)
  
  obj@result <- obj@result %>%
    filter(!is.na(score))
  
  saveRDS(
    obj,
    file.path(rnam_mirna_nona_dir, basename(input_file))
  )
}

# ------------------------------------------------------------------------------
# STEP 9B: FILTER NO-NA RESULTS BY ADJUSTED P-VALUE (padj < 0.05) & EXPORT TO EXCEL
# ------------------------------------------------------------------------------
cat("\n[9B] FILTERING NO-NA RESULTS BY padj < 0.05 AND EXPORTING TO EXCEL...\n")

rnam_mirna_nona_padj_dir <- file.path(results_dir, "RNAm&miRNA_no_NA_padj0.05")
dir.create(rnam_mirna_nona_padj_dir, recursive = TRUE, showWarnings = FALSE)

nona_files <- list.files(
  rnam_mirna_nona_dir,
  pattern     = "\\.RDS$",
  full.names  = TRUE,
  ignore.case = TRUE
)

for (input_file in nona_files) {
  
  obj <- readRDS(input_file)
  
  if (!all(c("p.adjust", "ID") %in% colnames(obj@result))) {
    cat(paste0(
      "  - Warning: 'p.adjust'/'ID' columns not found in '",
      basename(input_file), "'. Skipping.\n"
    ))
    next
  }
  
  # Filter by padj < 0.05
  obj@result <- obj@result %>%
    filter(p.adjust < 0.05)
  
  # Save the filtered RDS (original script behavior)
  saveRDS(
    obj,
    file.path(rnam_mirna_nona_padj_dir, basename(input_file))
  )
  df_res <- as.data.frame(obj@result)
  df_res[] <- lapply(df_res, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = ", "))
    } else {
      col
    }
  })
  
  # Define the Excel file name keeping the same base name but with .xlsx extension
  excel_name <- sub("\\.RDS$", ".xlsx", basename(input_file), ignore.case = TRUE)
  
  # Save the Excel file in the same rnam_mirna_nona_padj_dir folder
  write_xlsx(df_res, path = file.path(rnam_mirna_nona_padj_dir, excel_name))
}

# ------------------------------------------------------------------------------
# STEP 10: SUMMARIZE FINAL GO/KEGG RESULTS WITH miRNA INFORMATION
#          (Top 10 up / 10 down by NES+padj for GO; all pathways for KEGG)
# ------------------------------------------------------------------------------
cat("\n[10] SUMMARIZING FINAL GO/KEGG RESULTS WITH miRNA INFORMATION...\n")

rnam_mirna_top_padj_dir <- file.path(results_dir, "RNAm&miRNA_top_padj0.05")

dir.create(rnam_mirna_top_padj_dir, recursive = TRUE, showWarnings = FALSE)

make_GO_summary_S4 <- function(gsea_obj, n = 10) {
  obj <- gsea_obj
  is_kegg <- !"ONTOLOGY" %in% colnames(obj@result)
  
  group_cols <- if (is_kegg) {
    c("ID", "Description", "setSize", "enrichmentScore", "NES", "p.adjust")
  } else {
    c("ONTOLOGY", "ID", "Description", "setSize", "enrichmentScore", "NES", "p.adjust")
  }
  
  df <- obj@result
  
  # --- Target Genes: keep the highest score per gene within each pathway ---
  genes_summary <- df %>%
    filter(!is.na(SYMBOL)) %>%
    group_by(across(all_of(group_cols)), SYMBOL) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      `Target Genes`   = paste(sort(unique(SYMBOL)), collapse = ", "),
      `Target_Details` = paste(
        paste0(SYMBOL, " (score=", round(score, 1), ")"),
        collapse = "<br>"
      ),
      .groups = "drop"
    )
  
  # --- miRNAs: the L2FC is constant per miRNA, so we keep a single value ---
  mirna_summary <- df %>%
    filter(!is.na(miRNA)) %>%
    distinct(across(all_of(group_cols)), miRNA, L2FC_miRNA) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      miRNAs         = paste(sort(unique(miRNA)), collapse = ", "),
      `miRNA_Details` = paste(
        paste0(miRNA, " (L2FC=", round(L2FC_miRNA, 2), ")"),
        collapse = "<br>"
      ),
      .groups = "drop"
    )
  
  # --- Merge both summaries: a single unique row per ID/pathway ---
  summary_df <- genes_summary %>%
    left_join(mirna_summary, by = group_cols) %>%
    distinct(across(all_of(group_cols)), .keep_all = TRUE)
  
  # --- Apply Top 10 filtering ONLY to GO ---
  if (!is_kegg) {
    
    top_up <- summary_df %>%
      filter(NES > 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = n)
    
    top_down <- summary_df %>%
      filter(NES < 0) %>%
      arrange(p.adjust) %>%
      slice_head(n = n)
    
    summary_df <- bind_rows(top_up, top_down)
  }
  
  obj@result <- summary_df
  return(obj)
}

generate_top_tables <- function(input_dir, output_dir, n = 10) {
  
  input_files <- list.files(input_dir, pattern = "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
  
  for (input_file in input_files) {
    obj <- readRDS(input_file)
    obj_summary <- make_GO_summary_S4(obj, n = n)
    
    output_name <- sub("_RNAm_miRNA\\.RDS$", "_RNAm_miRNA_Top.RDS", basename(input_file))
    saveRDS(obj_summary, file.path(output_dir, output_name))
    
    df_to_save <- obj_summary@result
    
    output_name_csv <- sub("_RNAm_miRNA\\.RDS$", "_RNAm_miRNA_Top.csv", basename(input_file))
    write.table(
      df_to_save,
      file      = file.path(output_dir, output_name_csv),
      sep       = "\t",
      row.names = FALSE,
      quote     = TRUE,
      na        = "NA"
    )
  }
}

# Only the Top folder filtered by padj < 0.05 is generated
generate_top_tables(rnam_mirna_nona_padj_dir, rnam_mirna_top_padj_dir, n = 10)

# ------------------------------------------------------------------------------
# STEP 11: GENERATE GSEA DOTPLOTS (ACTIVATED vs SUPPRESSED)
# ------------------------------------------------------------------------------
cat("\n[11] GENERATING GSEA DOTPLOTS...\n")

plot_gsea_dotplot <- function(gsea_object, padj_filter = NULL) {
  
  df <- gsea_object@result
  
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  if (!is.null(padj_filter)) {
    df <- df %>% filter(p.adjust < padj_filter)
  }
  
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  df_plot <- df %>%
    mutate(Status = ifelse(NES > 0, "Activated", "Suppressed")) %>%
    arrange(NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot(
    df_plot,
    aes(x = NES, y = Description, size = setSize, colour = p.adjust)
  ) +
    geom_point(alpha = 0.9) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    facet_grid(. ~ Status, scales = "free_x", space = "free_x") +
    scale_size_continuous(name = "setSize") +
    scale_colour_gradient(
      name   = "p.adjust",
      low    = "red",
      high   = "blue",
      limits = c(max(df_plot$p.adjust), min(df_plot$p.adjust)),
      guide  = guide_colourbar(reverse = TRUE)
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 12),
      axis.text.y      = element_text(size = 9),
      legend.title     = element_text(face = "bold"),
      plot.title       = element_blank(),
      plot.margin      = margin(t = 10, r = 15, b = 10, l = 10, unit = "mm")
    ) +
    xlab("Normalized Enrichment Score (NES)") +
    ylab("")
  
  return(p)
}

plots_output_padj_dir <- file.path(results_dir, "Activated_Suppressed_Graphics_padj0.05")

dir.create(plots_output_padj_dir, recursive = TRUE, showWarnings = FALSE)

analyses_to_plot <- c(comparisons$go, comparisons$kegg)

generate_gsea_dotplots <- function(top_dir, output_dir, analyses) {
  
  for (go_kegg_name in analyses) {
    
    top_file <- file.path(top_dir, paste0(go_kegg_name, "_RNAm_miRNA_Top.RDS"))
    
    if (!file.exists(top_file)) {
      cat(paste0("  - Warning: summarized results for '", go_kegg_name, "' not found. Skipping plot.\n"))
      next
    }
    
    gsea_obj <- readRDS(top_file)
    p <- plot_gsea_dotplot(gsea_obj)
    
    if (is.null(p)) {
      cat(paste0("  - Warning: no data to plot for '", go_kegg_name, "'. Skipping.\n"))
      next
    }
    
    ggsave(
      filename = file.path(output_dir, paste0(go_kegg_name, "_dotplot.png")),
      plot     = p,
      width    = 13,
      height   = 7.5,
      dpi      = 300
    )
  }
}

# Dotplots are only generated from the Top folder filtered by padj < 0.05
generate_gsea_dotplots(rnam_mirna_top_padj_dir, plots_output_padj_dir, analyses_to_plot)

cat("\n[R script completed]\n")
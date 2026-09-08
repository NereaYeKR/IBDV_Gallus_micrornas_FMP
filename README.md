# RNAseq-miRNA-IBDV-DF1

miRNA-seq pipeline that runs differential expression analysis (DESeq2), predicts target genes (miRDB), and integrates the results with mRNA GO/KEGG enrichment (clusterProfiler/GSEA), generating an HTML report per comparison.

Part of the Master's Thesis *"Identification of Differentially Expressed microRNAs in Gallus gallus Infected with IBDV"*

## Context

IBDV (infectious bursal disease virus) infects chicken B lymphocytes and causes severe immunosuppression. This project studies the miRNA-mediated host response using three cell models derived from DF-1 (chicken embryo fibroblasts):

- **DF-1** — wild type, acute infection
- **DF-1P** — persistently infected
- **DF-1PC** — "recovered" cells

## Experimental design

- 18 miRNA-seq samples (single-end, Illumina), 3 biological replicates per condition
- 9 infected vs. 9 uninfected, across the three cell models
- Main comparisons:
  - Infected vs. control (overall effect of infection)
  - Infected DF-1P vs infected DF-1.
  - Uninfected DF-1P vs uninfected DF-1.
  - Uninfected DF-1P vs uninfected DF-1PC.
 
<img width="960" height="489" alt="Captura de pantalla 2026-09-08 a las 16 47 28" src="https://github.com/user-attachments/assets/62dfeb47-6969-4a52-8adc-50a72333bd13" />


## Pipeline

```
miRBase reference → Subread index/alignment → featureCounts
    → DESeq2 + miRDB targets + GO/KEGG integration → HTML reports
```

## Structure

```
.
├── run_pipeline.sh
├── config.sh.example
└── scripts/
    ├── 00.prepare_mirbase_reference.sh
    ├── 01.create_subread_index.sh
    ├── 02.fasta_to_gff.sh
    ├── 02.align_subread.sh
    ├── 03.counts_fc_subread.sh
    ├── 04.continueR.sh
    ├── 04.continueR.R
    └── 05.generate_report.py
```

`data/` and `results/` aren't included (large files) — set the paths in `config.sh`.

## Requirements

- `bash`, `wget`, `seqkit`
- `Subread` (index, alignment, featureCounts)
- `R` ≥ 4.2 — `DESeq2`, `clusterProfiler`, `org.Gg.eg.db`, `enrichplot`, `pathview`, `ggplot2`
- `python3` — `pandas`

## Results

- `deseq2_results/` — differential miRNAs (padj < 0.05, |log2FC| ≥ 1)
- `targets_annotated_no_NA/` — predicted targets (miRDB, score ≥ 80)
- `RNAm&miRNA_no_NA_padj0.05/`, `RNAm&miRNA_top_padj0.05/` — miRNA/mRNA integration
- `Activated_Suppressed_Graphics_padj0.05/` — GSEA dotplots
- `reports/` — HTML reports per comparison + `index.html`

## Author

Nerea Ye Keituqwa Rebollo

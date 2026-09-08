#!/bin/bash
# ==============================================================================
# Main Pipeline - miRNA-seq Analysis
#
# Description:
#   Orchestrates the complete end-to-end pipeline:
#     0. Preparation of the miRBase reference
#     1. Genome indexing (Subread)
#     2. Annotation (FASTA -> GFF) and alignment
#     3. Quantification (featureCounts)
#     4. Downstream analysis in R (DESeq2, miRNA targets, GO/KEGG)
#     5. HTML report generation
#
# Usage:
#   bash run_pipeline.sh [start_step]
#
#   start_step allows the pipeline to be resumed from a specific point
#   (e.g. './run_pipeline.sh 3'), skipping previous steps if their results
#   already exist. It defaults to step 0.
#
# Example:
#   bash run_pipeline.sh        # runs the full pipeline
#   bash run_pipeline.sh 3      # resumes from the quantification step
# ==============================================================================

source config.sh
set -e # stop if any command fails

START_STEP=${1:-0}

mkdir -p "$RESULTS" "$SUBREAD_INDEX_DIR" "$ALIGN_DIR" "$COUNTS_DIR"

# --- STEP 0: PREPARE miRBASE REFERENCE ---
if [ "$START_STEP" -le 0 ]; then
echo "STEP 0 - Preparing miRBase reference"

bash scripts/00.prepare_mirbase_reference.sh \
"$SPECIE" \
"data"
fi

# --- STEP 1: INDEXING ---
if [ "$START_STEP" -le 1 ]; then
echo "STEP 1 - Creating genome index"

# Option A: Subread (default)
bash scripts/01.create_subread_index.sh "$GENOME" "$SUBREAD_INDEX_DIR"

# Option B: Bowtie2
# bash scripts/01.create_bowtie2_index.sh "$GENOME" "$RESULTS/index_bowtie2"
fi

# --- STEP 2: ANNOTATION + ALIGNMENT ---
if [ "$START_STEP" -le 2 ]; then
echo "STEP 2 - miRNA annotation and alignment"

# 2.1 Create miRNA annotation (GFF)
bash scripts/02.fasta_to_gff.sh "$MATURE" "$SPECIE" "$MATURE_GFF"

# 2.2 Align reads with Subread (default)
bash scripts/02.align_subread.sh \
"$SUBREAD_INDEX_DIR/$(basename ${GENOME%.*})" \
"$MATURE_GFF" \
"$FASTQ_DIR" \
"$ALIGN_DIR"

# Option B: Bowtie2 alignment
# bash scripts/02.align_bowtie2.sh "$RESULTS/index_bowtie2" "$FASTQ_DIR" "$ALIGN_DIR"
fi

# --- STEP 3: QUANTIFICATION / COUNTS ---
if [ "$START_STEP" -le 3 ]; then
echo "STEP 3 - Quantification with FeatureCounts"

bash scripts/03.counts_fc_subread.sh \
"$ALIGN_DIR" \
"$MATURE_GFF" \
"$COUNTS_DIR/counts.txt"

# Optional alternative: Corset
# bash scripts/03.counts_corset_subread.sh \
# "$ALIGN_DIR" \
# "$CORSET_DIR"
#
# echo "STEP 3c - Rename clusters to miRNAs"
# bash scripts/03.rename_clusters_corset.sh \
# "$CORSET_DIR/mature-clusters.txt" \
# "$CORSET_DIR/mature-counts.txt" \
# "$COUNTS_DIR/mature-counts-gga.txt"
fi

# --- STEP 4: R ANALYSIS ---
if [ "$START_STEP" -le 4 ]; then
echo "STEP 4 - R downstream analysis"
bash scripts/04.continueR.sh
fi

# --- STEP 5: HTML REPORT GENERATION ---
if [ "$START_STEP" -le 5 ]; then
echo "STEP 5 - Generating HTML Reports"
python3 scripts/05.generate_report.py
fi

echo "PIPELINE COMPLETED SUCCESSFULLY"

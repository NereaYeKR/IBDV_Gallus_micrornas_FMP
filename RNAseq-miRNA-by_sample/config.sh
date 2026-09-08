#!/bin/bash

############################################################
# CONFIGURATION FILE - miRNA Subread Pipeline
############################################################

# Species code according to miRBase
SPECIE="gga"

############################################################
# INPUT REFERENCES
############################################################

# Path to the mature miRNA sequences (FASTA)
MATURE="data/mature.fa"

# Path to the genome reference file specific to the selected species
GENOME="data/mature-${SPECIE}-DNA.fa"

# GFF annotation file for mature miRNA coordinates
MATURE_GFF="data/mature_${SPECIE}.gff"

############################################################
# INPUT READS & METADATA
############################################################

# Directory containing processed/QC-filtered FASTQ files
FASTQ_DIR="data/fastq-qc-named"

# Path to tab-separated metadata file (required for R downstream processing)
SAMPLE_INFO="data/SampleInfo.tab"

# External database predictions for miRNA target validation
MIRDB_PREDICTIONS="data/miRDB_v6.0_prediction_result.txt"

############################################################
# OUTPUT DIRECTORIES
############################################################

# Root directory for all analysis outputs
RESULTS="results"

# Subread-specific subdirectories
SUBREAD_INDEX_DIR="$RESULTS/index_subread"
ALIGN_DIR="$RESULTS/alig_subread"
COUNTS_DIR="$RESULTS/count_featurecounts"

# Directories for downstream R analysis and generated reports
R_RESULTS_DIR="$RESULTS"
REPORTS_DIR="$RESULTS/reports"

############################################################
# VARIABLES EXCLUSIVAS PARA BOWTIE2 (Para evitar sobrescribir)
############################################################
# Directorio y prefijo del índice de Bowtie2
BOWTIE2_INDEX_DIR="$RESULTS/index_bowtie2"
BOWTIE2_INDEX="$BOWTIE2_INDEX_DIR/mature-${SPECIE}-DNA"

# Directorios de salida aislados para el alineamiento y el conteo
ALIGN_DIR_BOWTIE2="$RESULTS/alig_bowtie2"
COUNTS_DIR_BOWTIE2="$RESULTS/count_fc_bowtie2"

# Hilos de procesamiento para Bowtie2
THREADS=32

############################################################
# OPTIONAL: Corset analysis
############################################################
#CORSET_DIR="$RESULTS/corset_mature"

############################################################
# DESEQ2 CONTRASTS (Dynamic 2 vs 2)
############################################################
# Añade, edita o elimina las líneas entre los paréntesis libremente.
#
# Formato esperado: "NombreContraste,Columna,Numerador,Denominador"
#   - Numerador   = LFC > 0 (Condición a evaluar / Tratada)
#   - Denominador = LFC < 0 (Control / Baseline)

CONTRASTES=(
    "DF.1P_uninfected_vs_DF.1_uninfected,sample,DF.1P_uninfected,DF.1_uninfected"
    "DF.1P_IBDV_vs_DF.1_IBDV,sample,DF.1P_IBDV,DF.1_IBDV"
    "DF.1P_uninfected_vs_DF.1PC_uninfected,sample,DF.1P_uninfected,DF.1PC_uninfected"
)

# Empaqueta la lista en una sola variable exportable para que R pueda leerla
export DESEQ2_CONTRASTS=$(IFS=";"; echo "${CONTRASTES[*]}")

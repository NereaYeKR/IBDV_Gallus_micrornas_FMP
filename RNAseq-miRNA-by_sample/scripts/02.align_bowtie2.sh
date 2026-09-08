#!/bin/bash

# Bowtie2 alignment (miRNA-seq)
# Usage:
#   bash 02.align_bowtie2.sh index_dir fastq_dir output_dir
#
# Example (individual use):
#   bash 02.align_bowtie2.sh index_bowtie2 fastq-qc results/alig_bowtie2
#
# NOTE: Requires bowtie2 and samtools installed


source config.sh 2>/dev/null

INDEX=${1:-$BOWTIE2_INDEX}
FASTQ_DIR=${2:-$FASTQ_DIR}
OUTPUT_DIR=${3:-$ALIGN_DIR}

BOWTIE2_OPTS="-a --very-sensitive-local -p $THREADS"

mkdir -p "$OUTPUT_DIR"

for fq in ${FASTQ_DIR}/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    echo "Processing $base..."

    bowtie2 $BOWTIE2_OPTS -x "$INDEX" -U "$fq" | \
    samtools view -bS - | \
    samtools sort -o "${OUTPUT_DIR}/${base}_sorted.bam"

done

echo "Done!"

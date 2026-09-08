#!/bin/bash
# ==============================================================================
# 02 - Subread alignment (miRNA-seq), POSIX-compatible
#
# Description:
#   Aligns each miRNA-seq sample (single-end or paired-end) against the
#   Subread index generated in step 01, using the mature miRNA GFF
#   annotation generated in step 02 (fasta_to_gff) to report annotated
#   alignments. The script automatically iterates over every FASTQ file
#   containing "R1" in its name within the given directory, detecting its
#   "R2" counterpart if present (paired-end) or treating it as single-end
#   otherwise.
#
# Usage:
#   bash scripts/02.align_subread.sh [index] [gff] [fastq_dir] [output_dir]
#
# Example (standalone):
#   bash scripts/02.align_subread.sh \
#       results/index_subread/genome data/mature_gga.gff fastq-qc results/alig_subread
# ==============================================================================

# Load config.sh to get the environment variables (GENOME, SUBREAD_INDEX_DIR, etc.)
[ -f config.sh ] && source config.sh

# Assign arguments, with defaults taken from config.sh
ref_index=${1:-$SUBREAD_INDEX_DIR/$(basename ${GENOME%.*})}
ref_gtf=${2:-$MATURE_GFF}
fastq_dir=${3:-$FASTQ_DIR}
output_dir=${4:-$ALIGN_DIR}

# Subread constants (alignment mode)
RNAseq=0
DNAseq=1

# Prepare the GFF annotation argument, if available
if [ -n "$ref_gtf" ] && [ -f "$ref_gtf" ]; then
    ref_gtf_arg="-a $ref_gtf --gtfFeature miRNA --gtfAttr Name"
else
    ref_gtf_arg=''
fi

mkdir -p "$output_dir"

# Iterate over every sample (R1), detecting its R2 counterpart if present
for R1 in "$fastq_dir"/*R1*.{fastq,fq,fastq.gz,fq.gz} ; do
    [ -e "$R1" ] || continue

    base=$(basename "$R1" .fastq.gz | sed 's/.fq.gz//')
    R2=$(echo "$R1" | sed -e 's/R1/R2/g')
    bam="$output_dir/${base}_sorted.bam"

    echo "Processing $base..."

    if [ ! -e "$R2" ] ; then
        # Single-end sample
        subread-align -t $DNAseq \
                      -n 35 \
                      -m 4 \
                      -M 3 \
                      -T 10 \
                      -I 0 \
                      --multiMapping \
                      -B 10 \
                      -i "$ref_index" \
                      -r "$R1" \
                      -o "$bam" \
                      $ref_gtf_arg \
                      --sv \
                      --sortReadsByCoordinates
    else
        # Paired-end sample
        subread-align -t $DNAseq \
                      -n 35 \
                      -m 4 \
                      -M 3 \
                      -T 16 \
                      -I 0 \
                      --multiMapping \
                      -B 10 \
                      -i "$ref_index" \
                      -r "$R1" -R "$R2" \
                      -o "$bam" \
                      $ref_gtf_arg \
                      --sv \
                      --sortReadsByCoordinates
    fi
done

#!/bin/bash
#PBS -N HiC_Trimming
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=100gb
#PBS -l walltime=04:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 0: PRE-PROCESSING HI-C READS
# ==============================================================================
# Description:
# Quality filters and trims adapters from raw Hi-C Illumina reads using fastp.
# This generates the 'trimmed_R1.fq.gz' files expected by Hifiasm and YAHS.
# ==============================================================================

module add fastp

RAW_DIR="/path/to/cluster/project_dir/raw_hic"
OUT_DIR="/path/to/cluster/project_dir/hic/trimmed_data"
mkdir -p "$OUT_DIR"

echo "Trimming Hi-C Reads..."
fastp -i "$RAW_DIR/raw_hic_R1.fastq.gz" -I "$RAW_DIR/raw_hic_R2.fastq.gz" \
      -o "$OUT_DIR/trimmed_R1.fq.gz" -O "$OUT_DIR/trimmed_R2.fq.gz" \
      --cut_front --cut_tail \
      --length_required 30 \
      --thread 16 \
      --json "$OUT_DIR/hic_fastp.json" --html "$OUT_DIR/hic_fastp.html"

echo "Hi-C Trimming Complete."

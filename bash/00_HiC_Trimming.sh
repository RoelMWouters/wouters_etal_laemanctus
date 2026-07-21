#!/bin/bash
#PBS -N HiC_Trimming
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=100gb
#PBS -l walltime=04:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 0: PRE-PROCESSING ARIMA HI-C READS
# ==============================================================================
# Description:
# Quality filters and trims adapters from raw Arima Hi-C Illumina reads using fastp.
# Applies Arima-specific 5bp 5' trimming to remove restriction enzyme motifs,
# ensuring reads do not contain artificial overhangs that cause soft-clipping.
# Enforces strict paired-end synchronization for downstream tools (Hifiasm/YAHS).
# ==============================================================================

module add fastp

# Define paths (Users should update these variables)
RAW_DIR="/path/to/cluster/project_dir/raw_hic"
OUT_DIR="/path/to/cluster/project_dir/hic/trimmed_data"

mkdir -p "$OUT_DIR"

# Move into local scratch space to handle temporary files
cd $SCRATCHDIR

echo "=== GENERATING ARIMA ADAPTER FASTA ==="
cat << 'EOF' > adapters.fasta
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
EOF

echo "=== TRIMMING HI-C READS (ARIMA OPTIMIZED) ==="
# Notes on parameters:
# --unpaired1/2: Routes singletons to scratch to keep the R1/R2 outputs perfectly synced.
# --trim_front1/2: Removes the 5bp Arima restriction motif overhangs.
# --length_required: Permissive 30bp filter to retain short but informative chimeric reads.

fastp -i "$RAW_DIR/raw_hic_R1.fastq.gz" \
      -I "$RAW_DIR/raw_hic_R2.fastq.gz" \
      -o "$OUT_DIR/trimmed_R1.fq.gz" \
      -O "$OUT_DIR/trimmed_R2.fq.gz" \
      --unpaired1 "unpaired_R1.fq.gz" \
      --unpaired2 "unpaired_R2.fq.gz" \
      --adapter_fasta adapters.fasta \
      --trim_front1 5 \
      --trim_front2 5 \
      --length_required 30 \
      --thread 16 \
      --json "$OUT_DIR/hic_fastp.json" \
      --html "$OUT_DIR/hic_fastp.html"

echo "=== HI-C TRIMMING COMPLETE ==="

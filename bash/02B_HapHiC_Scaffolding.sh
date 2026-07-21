#!/bin/bash
#PBS -N HapHiC_Scaffolding
#PBS -l select=1:ncpus=64:mem=256gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 2B: ALLELE-AWARE GENOME SCAFFOLDING USING HAPHIC
# ==============================================================================
# Description:
# Scaffolds a draft assembly using Hi-C data utilizing the HapHiC pipeline. 
# Optimized for highly heterozygous or phased diploid assemblies.
#
# Dependencies: bwa, samtools, HapHiC
# ==============================================================================

# --- Configuration Block ---
PREFIX="species_haphic"
PROJECT_DIR="/path/to/cluster/project_dir/haphic_assembly"
DRAFT_ASSEMBLY="/path/to/cluster/project_dir/genome_dir/genome.fasta"
HIC_R1="/path/to/cluster/project_dir/hic/trimmed_data/trimmed_R1.fq.gz"
HIC_R2="/path/to/cluster/project_dir/hic/trimmed_data/trimmed_R2.fq.gz"

THREADS=64
BWA_THREADS=48
SAMTOOLS_THREADS=16

FINAL_DIR="$PROJECT_DIR/final_results"
mkdir -p "$FINAL_DIR/Alignment" "$FINAL_DIR/HapHiC_Output"

rescue_data() {
    echo "=== JOB ENDED OR CRASHED: RESCUING DATA FROM SCRATCH ==="
    cp $SCRATCHDIR/Alignment/*.bam "$FINAL_DIR/Alignment/" 2>/dev/null || true
    cp -r $SCRATCHDIR/HapHiC_Output/* "$FINAL_DIR/HapHiC_Output/" 2>/dev/null || true
}
trap rescue_data EXIT

cd "$SCRATCHDIR" 
set -e

# --- Environment Setup ---
module add mambaforge
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /path/to/cluster/envs/haphic_env
HAPHIC_TOOL="/path/to/cluster/tools/HapHiC/haphic"

mkdir -p "$SCRATCHDIR/Alignment"
cd "$SCRATCHDIR/Alignment"

echo "=== COPYING DRAFT ASSEMBLY TO SCRATCH ==="
cp "$DRAFT_ASSEMBLY" draft.fasta

echo "=== INDEXING DRAFT ASSEMBLY FOR BWA ==="
bwa index draft.fasta
samtools faidx draft.fasta

# --- Mapping ---
if [ -f "$FINAL_DIR/Alignment/hic_mapped.bam" ]; then
    echo "=== CHECKPOINT 1 FOUND: Skipping Hi-C alignment. ==="
    cp "$FINAL_DIR/Alignment/hic_mapped.bam" "$SCRATCHDIR/Alignment/"
else
    echo "=== CHECKPOINT 1 NOT FOUND: Mapping Hi-C reads ==="
    bwa mem -5SP -t $BWA_THREADS draft.fasta "$HIC_R1" "$HIC_R2" | \
        samtools view -@ $SAMTOOLS_THREADS -bh - > hic_mapped.bam
    cp hic_mapped.bam "$FINAL_DIR/Alignment/"
fi

# --- Scaffolding ---
mkdir -p "$SCRATCHDIR/HapHiC_Output"
cd "$SCRATCHDIR/HapHiC_Output"

echo "=== RUNNING HAPHIC ALLELE-AWARE SCAFFOLDING ==="
$HAPHIC_TOOL pipeline "$SCRATCHDIR/Alignment/draft.fasta" \
                "$SCRATCHDIR/Alignment/hic_mapped.bam" \
                36 \
                --threads $THREADS

echo "=== PIPELINE COMPLETE ==="
trap - EXIT
rescue_data

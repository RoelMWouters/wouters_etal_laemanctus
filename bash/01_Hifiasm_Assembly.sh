#!/bin/bash
#PBS -N Hifiasm_ONT_Only
#PBS -l select=1:ncpus=64:mem=256gb:scratch_local=500gb
#PBS -l walltime=48:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 1: DE NOVO GENOME ASSEMBLY (HIFIASM - ONT ONLY)
# ==============================================================================
# Description:
# This script performs a de novo genome assembly using Oxford Nanopore (ONT) 
# long reads ONLY. It runs NanoPlot for raw read quality control, Hifiasm 
# for the assembly (generating primary and alternate contigs), QUAST for 
# assembly metrics, and BUSCO for completeness assessment.
#
# Dependencies: mambaforge/conda, hifiasm, nanoplot, busco, gfatools, quast, samtools
# ==============================================================================

# --- Configuration Block ---
PREFIX="species_name"
PROJECT_DIR="/path/to/cluster/project_dir/hifiasm_assembly"
ONT_RAW="/path/to/cluster/project_dir/raw_ont/reads.fastq.gz"

# BUSCO Lineage (e.g., sauropsida_odb10, vertebrata_odb10)
BUSCO_LINEAGE="sauropsida_odb10"

FINAL_DIR="$PROJECT_DIR/final_results"
mkdir -p "$FINAL_DIR"

# --- Emergency Rescue Trap ---
rescue_data() {
    echo "=== JOB ENDED OR CRASHED: RESCUING DATA FROM SCRATCH ==="
    cp -r $SCRATCHDIR/QC_Reports "$FINAL_DIR/" 2>/dev/null || true
    cp -r $SCRATCHDIR/Assembly_Output "$FINAL_DIR/" 2>/dev/null || true
    cp -r $SCRATCHDIR/BUSCO_Output "$FINAL_DIR/" 2>/dev/null || true
    cp -r $SCRATCHDIR/Assembly_Output/QUAST_Output "$FINAL_DIR/" 2>/dev/null || true
}
trap rescue_data EXIT

cd "$SCRATCHDIR" 

# --- Environment Setup ---
module add mambaforge
export CONDA_PKGS_DIRS=$SCRATCHDIR/conda_pkgs
export CONDA_ENVS_DIRS=$SCRATCHDIR/conda_envs

echo "=== INITIALIZING BIOINFORMATICS ENVIRONMENT ==="
mamba create -n hifiasm_env -c conda-forge -c bioconda hifiasm nanoplot busco fastp gfatools quast samtools -y

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hifiasm_env
set -e

# --- 1. Raw Read QC ---
mkdir -p "$SCRATCHDIR/QC_Reports/ONT_Raw"
echo "=== RUNNING DORADO SUP QUALITY CHECK (NanoPlot) ==="
NanoPlot -t 63 --fastq "$ONT_RAW" -o "$SCRATCHDIR/QC_Reports/ONT_Raw" --plots dot --legacy hex

# --- 2. Assembly ---
mkdir -p "$SCRATCHDIR/Assembly_Output"
cd "$SCRATCHDIR/Assembly_Output"

echo "=== RUNNING HIFIASM (ONT ONLY) ==="
hifiasm -o ${PREFIX}_assembly -t 63 \
        --ont "$ONT_RAW" \
        --telo-m CCCTAA \
        --max-od-ec 0.05

echo "=== EXTRACTING PRIMARY AND ALTERNATE CONTIGS TO FASTA ==="
gfatools gfa2fa ${PREFIX}_assembly.p_ctg.gfa > ${PREFIX}_primary.fasta
gfatools gfa2fa ${PREFIX}_assembly.a_ctg.gfa > ${PREFIX}_alternate.fasta

# --- 3. Quality Assessment ---
echo "=== GENERATING ASSEMBLY STATISTICS (QUAST) ==="
quast.py ${PREFIX}_primary.fasta ${PREFIX}_alternate.fasta -o QUAST_Output -t 31 --large

echo "=== CREATING FASTA INDEXES ==="
samtools faidx ${PREFIX}_primary.fasta
samtools faidx ${PREFIX}_alternate.fasta

mkdir -p "$SCRATCHDIR/BUSCO_Output"
mkdir -p "$SCRATCHDIR/busco_downloads"
cd "$SCRATCHDIR/BUSCO_Output"

echo "=== PRE-DOWNLOADING BUSCO LINEAGE TO PREVENT RACE CONDITIONS ==="
busco --download "$BUSCO_LINEAGE" --download_path "$SCRATCHDIR/busco_downloads"

echo "=== RUNNING BUSCO ON PRIMARY AND ALTERNATE ASSEMBLIES ==="
busco -i "$SCRATCHDIR/Assembly_Output/${PREFIX}_primary.fasta" -o BUSCO_Primary -l "$BUSCO_LINEAGE" -m geno -c 31 --download_path "$SCRATCHDIR/busco_downloads" --offline &
busco -i "$SCRATCHDIR/Assembly_Output/${PREFIX}_alternate.fasta" -o BUSCO_Alternate -l "$BUSCO_LINEAGE" -m geno -c 31 --download_path "$SCRATCHDIR/busco_downloads" --offline
wait

echo "=== PIPELINE COMPLETE ==="

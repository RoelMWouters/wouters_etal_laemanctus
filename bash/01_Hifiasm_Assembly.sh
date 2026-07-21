#!/bin/bash
#PBS -N Hifiasm_Assembly
#PBS -l select=1:ncpus=64:mem=256gb:scratch_local=500gb
#PBS -l walltime=48:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 1: DE NOVO GENOME ASSEMBLY (HIFIASM)
# ==============================================================================
# Description:
# This script performs a de novo genome assembly using Oxford Nanopore (ONT) 
# long reads, phased with pre-trimmed Hi-C data. It runs NanoPlot for raw 
# read quality control, Hifiasm for the actual assembly, QUAST for assembly 
# metrics, and BUSCO for completeness assessment on both haplotypes.
#
# Dependencies: mambaforge/conda, hifiasm, nanoplot, busco, gfatools, quast, samtools
# ==============================================================================

# --- Configuration Block ---
PREFIX="species_name"
PROJECT_DIR="/path/to/cluster/project_dir/hifiasm_assembly"
ONT_RAW="/path/to/cluster/project_dir/raw_ont/reads.fastq.gz"
TRIMMED_HIC_DIR="/path/to/cluster/project_dir/hic/trimmed_data"

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

# --- 2. Stage Hi-C Data ---
mkdir -p "$SCRATCHDIR/HIC_INPUT"
echo "=== CHECKING FOR PRE-TRIMMED HI-C DATA ==="
if [ -f "$TRIMMED_HIC_DIR/trimmed_R1.fq.gz" ]; then
    echo "=== Copying trimmed reads from network ==="
    cp "$TRIMMED_HIC_DIR"/trimmed_R*.fq.gz "$SCRATCHDIR/HIC_INPUT/"
else
    echo "=== ERROR: Trimmed files not found! Run the interactive prep first. ==="
    exit 1
fi

# --- 3. Assembly ---
mkdir -p "$SCRATCHDIR/Assembly_Output"
cd "$SCRATCHDIR/Assembly_Output"

echo "=== RUNNING HIFIASM (ONT + HI-C PHASING) ==="
hifiasm -o ${PREFIX}_phased -t 63 \
        --ont "$ONT_RAW" \
        --h1 "$SCRATCHDIR/HIC_INPUT/trimmed_R1.fq.gz" \
        --h2 "$SCRATCHDIR/HIC_INPUT/trimmed_R2.fq.gz" \
        --telo-m CCCTAA \
        --max-od-ec 0.05

echo "=== EXTRACTING HAPLOTYPES TO FASTA ==="
gfatools gfa2fa ${PREFIX}_phased.hic.hap1.p_ctg.gfa > ${PREFIX}_hap1.fasta
gfatools gfa2fa ${PREFIX}_phased.hic.hap2.p_ctg.gfa > ${PREFIX}_hap2.fasta

# --- 4. Quality Assessment ---
echo "=== GENERATING ASSEMBLY STATISTICS (QUAST) ==="
quast.py ${PREFIX}_hap1.fasta ${PREFIX}_hap2.fasta -o QUAST_Output -t 31 --large

echo "=== CREATING FASTA INDEXES ==="
samtools faidx ${PREFIX}_hap1.fasta
samtools faidx ${PREFIX}_hap2.fasta

mkdir -p "$SCRATCHDIR/BUSCO_Output"
mkdir -p "$SCRATCHDIR/busco_downloads"
cd "$SCRATCHDIR/BUSCO_Output"

echo "=== PRE-DOWNLOADING BUSCO LINEAGE TO PREVENT RACE CONDITIONS ==="
busco --download "$BUSCO_LINEAGE" --download_path "$SCRATCHDIR/busco_downloads"

echo "=== RUNNING BUSCO ON BOTH HAPLOTYPES SIMULTANEOUSLY ==="
busco -i "$SCRATCHDIR/Assembly_Output/${PREFIX}_hap1.fasta" -o BUSCO_Hap1 -l "$BUSCO_LINEAGE" -m geno -c 31 --download_path "$SCRATCHDIR/busco_downloads" --offline &
busco -i "$SCRATCHDIR/Assembly_Output/${PREFIX}_hap2.fasta" -o BUSCO_Hap2 -l "$BUSCO_LINEAGE" -m geno -c 31 --download_path "$SCRATCHDIR/busco_downloads" --offline
wait

echo "=== PIPELINE COMPLETE ==="

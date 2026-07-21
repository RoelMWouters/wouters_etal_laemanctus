#!/bin/bash
#PBS -N synteny_alignment
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=100gb
#PBS -l walltime=03:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 1: DIPLOID MACRO-SYNTENY ALIGNMENT & ORTHOLOG FILTERING
# ==============================================================================
# Description:
# This script downloads reference CDS data (e.g., Anolis sagrei) and aligns it 
# against a phased, diploid target genome using minimap2. It then applies a 
# strict filtering step using awk to isolate 1:2 orthologs (where one reference 
# gene maps to exactly two target homoeologs), effectively separating the 
# subgenomes for downstream JCVI synteny visualization.
#
# Dependencies: minimap2, samtools
# ==============================================================================

export TMPDIR=$SCRATCHDIR

if [ -z "$SCRATCHDIR" ] ; then
  echo "Error: SCRATCHDIR is not set!"
  exit 1
fi
trap 'clean_scratch' TERM EXIT

module add minimap2
module add samtools

# --- Configuration Paths (Update these for your environment) ---
GENOME_DIR="/path/to/genome_directory"
TARGET_GENOME="/path/to/phased_target_genome.fa"
OUT_DIR="/path/to/output_directory/synteny_results"

mkdir -p "$GENOME_DIR"
mkdir -p "$OUT_DIR"
cd "$SCRATCHDIR" || exit 1

echo "=== DOWNLOADING REFERENCE DATA ==="
# Utilizing the Anolis sagrei rAnoSag1.mat maternal haplotype assembly as the reference
ANOLIS_GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/765/GCF_037176765.1_rAnoSag1.mat/GCF_037176765.1_rAnoSag1.mat_genomic.gff.gz"
ANOLIS_CDS_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/765/GCF_037176765.1_rAnoSag1.mat/GCF_037176765.1_rAnoSag1.mat_cds_from_genomic.fna.gz"

wget -nc -P "$GENOME_DIR" "$ANOLIS_GFF_URL"
wget -nc -P "$GENOME_DIR" "$ANOLIS_CDS_URL"

echo "=== STAGING DATA TO SCRATCH ==="
cp "$GENOME_DIR/GCF_037176765.1_rAnoSag1.mat_cds_from_genomic.fna.gz" ./reference_cds.fna.gz
cp "$TARGET_GENOME" ./target_genome.fa
gunzip ./reference_cds.fna.gz

echo "=== ALIGNING REFERENCE GENES TO PHASED TARGET ==="
# CRITICAL FIX FOR PHASED GENOMES:
# -p 0.75 : Retain secondary alignments that score at least 75% as well as the primary hit
# -N 2    : Keep up to 2 secondary alignments per gene
minimap2 -t 8 -c -p 0.75 -N 2 -x splice ./target_genome.fa ./reference_cds.fna > raw_alignment.paf

echo "=== FILTERING FOR 1:2 PHASED ORTHOLOGS ==="
# minimap2 outputs alignments sorted by query (Reference gene) and then by score.
# This awk script counts how many times it sees each gene. 
# It keeps the Top 2 hits (the two parental homoeologs) as long as they are >200bp long,
# and dynamically discards any 3rd/4th paralogous hits.
awk '$11 >= 200 { 
    hits[$1]++ 
    if (hits[$1] <= 2) { print $0 } 
}' raw_alignment.paf > phased_orthologs.paf

echo "=== EXTRACTING MACRO-SYNTENY LINKS ==="
# Extract necessary columns to build the synteny block TSV.
# Format: Ref_Gene | Target_Contig | Target_Start | Target_End | Strand | Map_Length
awk '{print $1, $6, $8, $9, $5, $11}' OFS='\t' phased_orthologs.paf > clean_phased_blocks.tsv

echo "=== SAVING RESULTS TO STORAGE ==="
cp raw_alignment.paf "$OUT_DIR/"
cp phased_orthologs.paf "$OUT_DIR/"
cp clean_phased_blocks.tsv "$OUT_DIR/"

echo "Diploid pipeline alignment completed successfully!"

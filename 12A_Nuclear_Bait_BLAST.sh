#!/bin/bash
#PBS -N Nuclear_Bait_BLAST
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=20gb
#PBS -l walltime=04:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 12A: NUCLEAR BAIT BLAST SEARCH (PHASED HITS)
# ==============================================================================
# Description:
# Blasts nuclear bait sequences (target loci) against the phased subgenomes 
# to generate the '*_phased_hits.txt' files required by the Python parser in 
# Script 12. 
# ==============================================================================

module add blast-plus

WORK_DIR="/path/to/cluster/project_dir"
GENOME_DIR="$WORK_DIR/genome_dir/haploid_genome"
BAIT_DIR="$WORK_DIR/output_dir/nuclear_baits"
OUT_DIR="$WORK_DIR/output_dir/blast_results"

mkdir -p "$OUT_DIR"
cd "$SCRATCHDIR"

# Standardized Subgenome names from Script 6
cp "$GENOME_DIR/Subgenome_A.fasta" .
cp "$GENOME_DIR/Subgenome_B.fasta" .
cat Subgenome_A.fasta Subgenome_B.fasta > combined_phased.fasta

makeblastdb -in combined_phased.fasta -dbtype nucl -out phased_db

echo "Blasting Nuclear Baits..."
for BAIT in "$BAIT_DIR"/*.fasta; do
    GENE_NAME=$(basename "$BAIT" .fasta)
    
    tblastn -query "$BAIT" \
            -db phased_db \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -out "$OUT_DIR/${GENE_NAME}_phased_hits.txt"
done

echo "Nuclear BLAST complete."

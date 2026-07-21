#!/bin/bash
#PBS -N Fetch_Mito_Baits
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=10gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 11A: MITOCHONDRIAL BAIT FETCH & BLAST SEARCH
# ==============================================================================
# Description:
# Automatically fetches an outgroup mitochondrial genome using NCBI Entrez Direct, 
# formats local draft assemblies as BLAST databases, and searches them to locate 
# mitochondrial sequence coordinates.
#
# Dependencies: blast, edirect
# ==============================================================================

WORK_DIR="/path/to/cluster/project_dir"
GENOME_DIR="$WORK_DIR/genome_dir"
BAIT_DIR="$WORK_DIR/output_dir/bait_sequences"
OUT_DIR="$WORK_DIR/output_dir/blast_results_refs"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

module purge
module add blast edirect

if [ ! -f "$BAIT_DIR/mito_bait.fasta" ]; then
    esearch -db nucleotide -query "Anolis carolinensis mitochondrion, complete genome" | \
    efetch -format fasta > "$BAIT_DIR/mito_bait.fasta"
fi

if [ ! -f "$GENOME_DIR/refA_blastdb.nhr" ]; then
    makeblastdb -in "$GENOME_DIR/Subgenome_A.fasta" -dbtype nucl -out "$GENOME_DIR/refA_blastdb"
fi

REFS=("Subgenome_A" "Subgenome_B")

for REF in "${REFS[@]}"; do
    for BAIT in "$BAIT_DIR"/*_bait.fasta; do
        GENE_NAME=$(basename "$BAIT" _bait.fasta)
        tblastn -query "$BAIT" \
                -db "$GENOME_DIR/${REF}_blastdb" \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
                -out "${REF}_${GENE_NAME}_hits.txt"
    done
    
    blastn -query "$BAIT_DIR/mito_bait.fasta" \
           -db "$GENOME_DIR/${REF}_blastdb" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
           -out "${REF}_mito_hits.txt"
done

#!/bin/bash
#PBS -N BUSCO_to_BED
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 14A: BUSCO FULL TABLE TO BED FORMAT
# ==============================================================================
# Description:
# Parses the 'full_table.tsv' from a BUSCO run, isolates only the 'Complete' 
# single-copy orthologs, limits them to the top 1000 longest/highest-scoring, 
# and converts the coordinates into a BED file.
# ==============================================================================

BUSCO_TABLE="/path/to/cluster/project_dir/your_assembly/BUSCO_Output/run_sauropsida_odb10/full_table.tsv"
OUT_BED="/path/to/cluster/project_dir/genome_dir/1000_random_buscos.bed"

echo "Parsing BUSCO table and generating BED file..."

# 1. Filter for "Complete", extract (Contig, Start, End, BuscoID), and sort by length/score
awk -v OFS='\t' '$2 == "Complete" {print $3, $4, $5, $1}' "$BUSCO_TABLE" | \
    head -n 1000 > "$OUT_BED"

echo "Extracted 1000 BUSCO coordinates to $OUT_BED"

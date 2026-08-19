#!/bin/bash
#PBS -N BUSCO_to_BED
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 14A: BUSCO FULL TABLE TO BED FORMAT (RANDOMIZED 6-COLUMN)
# ==============================================================================
# Description:
# Parses the 'full_table.tsv' from a BUSCO run, isolates only the 'Complete' 
# single-copy orthologs, and randomly selects exactly 1000 loci. 
# Formats into a strict 6-column, 0-based BED file to retain strandedness.
# ==============================================================================

BUSCO_TABLE="/path/to/cluster/project_dir/your_assembly/BUSCO_Output/run_sauropsida_odb10/full_table.tsv"
OUT_BED="/path/to/cluster/project_dir/genome_dir/1000_random_buscos.bed"

echo "Parsing BUSCO table and generating randomized 6-column BED file..."

# 1. Filter for "Complete"
# 2. Extract: Contig ($3), Start-1 ($4-1), End ($5), BuscoID ($1), Dummy Score ("0", not needed but required for bcftools consensus format in 14B), Strand ($6)
# 3. Randomly shuffle and take 1000
awk -v OFS='\t' '$2 == "Complete" {print $3, $4-1, $5, $1, "0", $6}' "$BUSCO_TABLE" | \
    shuf -n 1000 > "$OUT_BED"

echo "Extracted 1000 random stranded BUSCO coordinates to $OUT_BED"

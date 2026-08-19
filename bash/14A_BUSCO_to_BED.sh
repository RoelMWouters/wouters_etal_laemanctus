#!/bin/bash
#PBS -N BUSCO_to_BED_Duplicated
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 14A: BUSCO FULL TABLE TO BED FORMAT (1000 DUPLICATED PAIRS)
# ==============================================================================
# Description:
# Parses the 'full_table.tsv' from a BUSCO run, isolates only the 'Duplicated' 
# orthologs (present in exactly two copies), pairs them, and randomly selects 
# 1000 gene pairs (2000 total sequences). Formats into a strict 6-column, 
# 0-based BED file with exact _Copy1 and _Copy2 nomenclature.
# ==============================================================================

# --- Configuration Paths ---
BUSCO_TABLE="/path/to/your/project/BUSCO_Output/run_sauropsida_odb10/full_table.tsv"
OUT_BED="/path/to/your/project/genome_dir/1000_random_buscos.bed"

echo "Parsing BUSCO table for strict duplicated pairs..."

awk -v OFS='\t' '
  # 1. Only look at Duplicated genes
  $2 == "Duplicated" {
    busco_id = $1
    contig = $3
    start = $4 - 1  # 0-based conversion for BED format
    end = $5
    strand = $6
    score = $7
    
    count[busco_id]++
    
    # 2. Tag the first instance as Copy1 and the second as Copy2
    if (count[busco_id] == 1) {
      copy1[busco_id] = contig OFS start OFS end OFS busco_id "_Copy1" OFS score OFS strand
    } else if (count[busco_id] == 2) {
      copy2[busco_id] = contig OFS start OFS end OFS busco_id "_Copy2" OFS score OFS strand
      valid[busco_id] = 1 # Mark as a valid pair
    } else {
      valid[busco_id] = 0 # Invalidate if there are 3 or more copies
    }
  }
  
  # 3. Output the valid pairs glued together by a pipe character
  END {
    for (id in valid) {
      if (valid[id] == 1) {
        print copy1[id] "|" copy2[id]
      }
    }
  }
' "$BUSCO_TABLE" | \
shuf -n 1000 | \
awk -F'|' -v OFS='\n' '{print $1, $2}' > "$OUT_BED"

echo "Extracted 1000 random duplicated BUSCO pairs (2000 sequences) to $OUT_BED"

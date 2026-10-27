#!/bin/bash
#PBS -N IBS_Distances
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=50gb
#PBS -l walltime=01:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# IBS & GENETIC DISTANCE
# Description: Uses PLINK to calculate pairwise absolute SNP differences and 
# percentage divergence (1 - IBS) across haploid mapped clones. Enforces strict 
# 100% missing-data thresholds to prevent missingness from inflating distance.
# ==============================================================================

set -euo pipefail
export TMPDIR=$SCRATCHDIR

module add bcftools
module add vcftools

# ==============================================================================
# 1. GLOBAL PATHS & CONFIGURATION
# ==============================================================================
PROJECT_DIR="/path/to/project"
RAW_VCF="${PROJECT_DIR}/bam_output/haploid_bam/loh_ljul/vcfs/ljul.heterozygosity.vcf.gz"
STATS_DIR="${PROJECT_DIR}/output/clonal_inheritance/stats_absolute_divergence"

SPECIES_ID="ljul_clones_haploid"
FILTER_STATS="$STATS_DIR/filtering_stats.txt"

mkdir -p "$STATS_DIR" "$SCRATCHDIR"

echo "=== HAPLOID FILTERING STATISTICS ===" > "$FILTER_STATS"
echo "Target VCF: $RAW_VCF" >> "$FILTER_STATS"
echo "-----------------------------------------------" >> "$FILTER_STATS"

# ==============================================================================
# 2. PRE-FILTERING DATASET CHECK
# ==============================================================================
echo "Samples included in this run:" >> "$FILTER_STATS"
bcftools query -l "$RAW_VCF" | sed 's/^/  - /' >> "$FILTER_STATS"

RAW_SNPS=$(bcftools view -H "$RAW_VCF" | wc -l)
echo -e "\nRaw SNPs before filtering: $RAW_SNPS" >> "$FILTER_STATS"

# ==============================================================================
# 3. DOWNLOAD PLINK STANDALONE (TEMPORARY SCRATCH)
# ==============================================================================
wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O "$SCRATCHDIR/plink.zip"
unzip -jo "$SCRATCHDIR/plink.zip" plink -d "$SCRATCHDIR/" > /dev/null
chmod +x "$SCRATCHDIR/plink"
PLINK_BIN="$SCRATCHDIR/plink"

# ==============================================================================
# 4. STRICT CLONE-SAFE FILTERING 
# ==============================================================================
# --max-missing 1.0 forces 100% data presence across all samples
vcftools --gzvcf "$RAW_VCF" \
    --max-missing 1.0 \
    --minQ 30 \
    --minDP 10 \
    --maxDP 60 \
    --recode --out "$SCRATCHDIR/${SPECIES_ID}_filtered"

FILTERED_SNPS=$(grep "After filtering, kept" "$SCRATCHDIR/${SPECIES_ID}_filtered.log" | awk '{print $5}')
echo "Strictly filtered SNPs retained: $FILTERED_SNPS" >> "$FILTER_STATS"

RETENTION_PCT=$(awk -v filt="$FILTERED_SNPS" -v raw="$RAW_SNPS" 'BEGIN { printf "%.2f", (filt/raw)*100 }')
echo "Percentage of SNPs retained: ${RETENTION_PCT}%" >> "$FILTER_STATS"

# ==============================================================================
# 5. PLINK DISTANCE MATRICES
# ==============================================================================
echo "=== CALCULATING ABSOLUTE SNP DIFFERENCES ==="
$PLINK_BIN --vcf "$SCRATCHDIR/${SPECIES_ID}_filtered.recode.vcf" \
      --double-id --allow-extra-chr --keep-allele-order \
      --distance square allele-ct \
      --out "$STATS_DIR/${SPECIES_ID}_absolute_diffs"

echo "=== CALCULATING PERCENTAGE DIVERGENCE (1 - IBS) ==="
$PLINK_BIN --vcf "$SCRATCHDIR/${SPECIES_ID}_filtered.recode.vcf" \
      --double-id --allow-extra-chr --keep-allele-order \
      --distance square 1-ibs \
      --out "$STATS_DIR/${SPECIES_ID}_percent_divergence"

# Cleanup artifacts
rm -f "$STATS_DIR"/*.log "$STATS_DIR"/*.nosex
echo "Part PLINK pipeline finished."

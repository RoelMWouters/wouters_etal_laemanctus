#!/bin/bash
#PBS -N PLINK_Stats
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=10gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 8: POPULATION GENOMICS STATISTICS (PLINK)
# ==============================================================================
# Description:
# Uses PLINK to perform linkage disequilibrium (LD) pruning, followed by 
# calculating Identity-By-Descent (IBD), Principal Component Analysis (PCA), 
# and Identity-By-State (IBS) distances from a strictly filtered VCF.
#
# Dependencies: bcftools, unzip, wget, plink (downloads internally)
# ==============================================================================

export TMPDIR=$SCRATCHDIR
module add bcftools

FINAL_VCF="/path/to/cluster/project_dir/output_dir/clonality/vcfs/population_final_filtered.vcf.gz"
STATS_DIR="/path/to/cluster/project_dir/output_dir/clonality/stats"
SPECIES_ID="species"

mkdir -p "$STATS_DIR" "$SCRATCHDIR"
TEMP_VCF="$SCRATCHDIR/${SPECIES_ID}.temp.vcf.gz"

wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O "$SCRATCHDIR/plink.zip"
unzip -jo "$SCRATCHDIR/plink.zip" plink -d "$SCRATCHDIR/" > /dev/null
chmod +x "$SCRATCHDIR/plink"
PLINK_BIN="$SCRATCHDIR/plink"

# Thin VCF to speed up PLINK parsing (optional, depending on marker density)
bcftools view "$FINAL_VCF" | awk '/^#/ {print; next} {if (NR % 5 == 0) print}' | bgzip -c > "$TEMP_VCF"

$PLINK_BIN --vcf "$TEMP_VCF" \
      --double-id \
      --allow-extra-chr \
      --indep-pairwise 50 10 0.1 \
      --out "$SCRATCHDIR/${SPECIES_ID}_pruned"

$PLINK_BIN --vcf "$TEMP_VCF" \
      --double-id \
      --allow-extra-chr \
      --extract "$SCRATCHDIR/${SPECIES_ID}_pruned.prune.in" \
      --genome \
      --out "$STATS_DIR/${SPECIES_ID}_IBD"

$PLINK_BIN --vcf "$TEMP_VCF" \
      --double-id \
      --allow-extra-chr \
      --extract "$SCRATCHDIR/${SPECIES_ID}_pruned.prune.in" \
      --pca \
      --out "$STATS_DIR/${SPECIES_ID}_PCA"

$PLINK_BIN --vcf "$TEMP_VCF" \
      --double-id \
      --allow-extra-chr \
      --extract "$SCRATCHDIR/${SPECIES_ID}_pruned.prune.in" \
      --distance 1-ibs square \
      --out "$STATS_DIR/${SPECIES_ID}_IBS"

rm -f "$STATS_DIR"/*.log "$STATS_DIR"/*.nosex

#!/bin/bash
#PBS -N Calc_BAF
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 10: B-ALLELE FREQUENCY EXTRACTION
# ==============================================================================
# Description:
# Extracts raw Allelic Depths (AD) from high-quality variants in a VCF. Used 
# subsequently to calculate B-Allele Frequency for hybridization/ploidy checks.
#
# Dependencies: bcftools
# ==============================================================================

export TMPDIR=$SCRATCHDIR
module add bcftools

VCF="/path/to/cluster/project_dir/bam_output/haploid_bam/vcfs/ljul.heterozygosity.vcf.gz"
OUT_DIR="/path/to/cluster/project_dir/output_dir/baf_analysis"
mkdir -p "$OUT_DIR"

bcftools query -i "QUAL>=30 && MQ>=40" \
               -f '%CHROM\t%POS[\t%AD]\n' "$VCF" > "$SCRATCHDIR/raw_baf_depths.txt"

echo -e "CHROM\tPOS\tS01_AD\tS02_AD\tS03_AD\tS04_AD" > "$OUT_DIR/Final_BAF_Data.tsv"
cat "$SCRATCHDIR/raw_baf_depths.txt" >> "$OUT_DIR/Final_BAF_Data.tsv"

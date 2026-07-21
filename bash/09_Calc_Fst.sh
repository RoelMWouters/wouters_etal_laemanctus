#!/bin/bash
#PBS -N Calc_Fst
#PBS -l select=1:ncpus=4:mem=8gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 9: PAIRWISE FIXATION INDEX (FST)
# ==============================================================================
# Description:
# Calculates sliding-window Weir & Cockerham Fst between two defined 
# subpopulations to detect genomic regions of high divergence.
#
# Dependencies: bcftools, vcftools
# ==============================================================================

export TMPDIR=$SCRATCHDIR
module add bcftools vcftools

FINAL_VCF="/path/to/cluster/project_dir/output_dir/clonality/vcfs/population_final_filtered.vcf.gz"
STATS_DIR="/path/to/cluster/project_dir/output_dir/clonality/stats"

mkdir -p "$STATS_DIR" "$SCRATCHDIR"

# Define Populations
echo -e "S04\nS05\nS06" > "$SCRATCHDIR/popA.txt"
echo -e "S01\nS02\nS03" > "$SCRATCHDIR/popB.txt"
echo -e "S07\nS08" > "$SCRATCHDIR/popC.txt"

vcftools --gzvcf "$FINAL_VCF" \
    --weir-fst-pop "$SCRATCHDIR/popA.txt" \
    --weir-fst-pop "$SCRATCHDIR/popB.txt" \
    --fst-window-size 100000 \
    --fst-window-step 50000 \
    --out "$STATS_DIR/fst_A_vs_B"

vcftools --gzvcf "$FINAL_VCF" \
    --weir-fst-pop "$SCRATCHDIR/popA.txt" \
    --weir-fst-pop "$SCRATCHDIR/popC.txt" \
    --fst-window-size 100000 \
    --fst-window-step 50000 \
    --out "$STATS_DIR/fst_A_vs_C"

vcftools --gzvcf "$FINAL_VCF" \
    --weir-fst-pop "$SCRATCHDIR/popB.txt" \
    --weir-fst-pop "$SCRATCHDIR/popC.txt" \
    --fst-window-size 100000 \
    --fst-window-step 50000 \
    --out "$STATS_DIR/fst_B_vs_C"



    

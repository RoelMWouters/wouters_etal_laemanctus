#!/bin/bash

# ==============================================================================
# Laemanctus Hi-C Figure Pipeline
# Requirements: Execute within the 'hic_figures' conda environment.
# ==============================================================================
## Create the environment with the required tools
module add mambaforge
mamba create -n hic_figures python=3.9 hictk hic2cool hicexplorer -c bioconda -c conda-forge -y

# Navigate to the working directory
cd your/path/to/hic/file

# Activate the environment
source activate hic_figures

# ------------------------------------------------------------------------------
# 1. FORMAT CONVERSION (.hic v9 to .cool)
# ------------------------------------------------------------------------------
# Using hictk to convert the Juicer/YaHS .hic output into cooler format.
# We extract both 50kb (for insets) and 250kb (to prevent OOM errors on the full genome).

hictk convert Ljul_curated_map.hic laemanctus_50kb.cool --resolutions 50000
hictk convert Ljul_curated_map.hic laemanctus_250kb.cool --resolutions 250000


# ------------------------------------------------------------------------------
# 2. WHOLE GENOME DIPLOID HEATMAP (250kb)
# ------------------------------------------------------------------------------
# Plotted sequentially (All HapA -> All HapB). 
# Resolution lowered to 250kb to successfully render within standard cluster memory limits.

hicPlotMatrix \
  --matrix laemanctus_250kb.cool \
  --outFileName diploid_whole_genome_interleaved_250kb.pdf \
  --chromosomeOrder Chr1_hapA Chr1_hapB Chr2_hapA Chr2_hapB Chr3_hapA Chr3_hapB Chr4_hapA Chr4_hapB Chr5_hapA Chr5_hapB Chr6_hapA Chr6_hapB Chr7_hapA Chr7_hapB Chr8_hapA Chr8_hapB Chr9_hapA Chr9_hapB Chr10_hapA Chr10_hapB Chr11_hapA Chr11_hapB Chr12_hapA Chr12_hapB Chr13_hapA Chr13_hapB Chr14_hapA Chr14_hapB Chr15_hapA Chr15_hapB Chr16_hapA Chr16_hapB Chr17_hapA Chr17_hapB ChrX_hapA ChrX_hapB \
  --vMax 50 \
  --colorMap Reds \
  --clearMaskedBins \
  --rotationX 45 \
  --fontsize 4 \
  --increaseFigureWidth 8 \
  --increaseFigureHeight 8 \
  --dpi 600 \
  --title ""


# ------------------------------------------------------------------------------
# 3. INTERLEAVED MICROCHROMOSOMES & X INSET (50kb)
# ------------------------------------------------------------------------------
# Plotted at high resolution (50kb).
# Arranged alternating (HapA, HapB, HapA, HapB...) to visualize homologous trans-contacts.
# Note: --log1p is removed and replaced with a hard ceiling (--vMax 50) to mimic 
# Juicebox's color saturation slider and reveal the faint "ghost line" diagonals.

hicPlotMatrix \
  --matrix laemanctus_50kb.cool \
  --outFileName diploid_microchromosomes_interleaved_vmax50.pdf \
  --chromosomeOrder Chr7_hapA Chr7_hapB Chr8_hapA Chr8_hapB Chr9_hapA Chr9_hapB Chr10_hapA Chr10_hapB Chr11_hapA Chr11_hapB Chr12_hapA Chr12_hapB Chr13_hapA Chr13_hapB Chr14_hapA Chr14_hapB Chr15_hapA Chr15_hapB Chr16_hapA Chr16_hapB Chr17_hapA Chr17_hapB ChrX_hapA ChrX_hapB \
  --colorMap Reds \
  --clearMaskedBins \
  --vMax 15 \
  --rotationX 45 \
  --fontsize 4 \
  --increaseFigureWidth 8 \
  --increaseFigureHeight 8 \
  --dpi 600 \
  --title ""

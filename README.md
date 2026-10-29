# wouters_etal_laemanctus
Bioinformatic pipelines and scripts accompanying **Wouters et al. Origin of the largest obligately parthenogenetic vertebrate**.

# Laemanctus julioi genome assembly and phylogenomics pipeline #

**Author:** Roel Wouters | Faculty of Science, Charles University, Prague, Czech Republic 
**Description:** A comprehensive bioinformatics pipeline for *de novo* genome assembly, Hi-C scaffolding, population genomics, and species tree inference in *Laemanctus* lizards.

This repository contains the bash and R scripts used to process raw sequencing reads (ONT/Illumina/Hi-C) through to final figures and phylogenomic topologies.

## Pipeline Overview

The pipeline is designed for execution on a PBS-scheduled High-Performance Computing (HPC) cluster (e.g., MetaCentrum), utilising strict environment isolation and modular job arrays. It is divided into three primary phases:

### Phase 1: Genome Assembly & Scaffolding
*   **`00_HiC_Trimming.sh`**: Pre-processing and adapter trimming of raw Hi-C reads.
*   **`01_Hifiasm_Assembly.sh`**: *De novo* assembly of ONT long reads.
*   **`02A_YAHS_Scaffolding.sh`**: Scaffolding of uncollapsed assemblies for Juicebox curation.
*   **`02B_HapHiC_Scaffolding.sh`**: Allele-aware scaffolding for highly heterozygous genomes.
*   **`02C_HiCExplorer_figures.sh`**: HiCExplorer for generating HiC plots.
*   **`06_Subgenome_Phasing.sh`**: Segregation of draft assemblies into distinct parental subgenomes based on mapping bias.

### Phase 2: Population Genomics & Variant Calling
*   **`03_FastP_BWA_Alignment.sh`**: Job array for trimming Illumina short reads and mapping to the curated reference.
*   **`04_GVCF_Generation.sh`**: Highly parallelized GATK HaplotypeCaller scatter-gather execution.
*   **`05_Joint_Calling.sh`**: Population joint genotyping and strict paralog purging based on depth thresholds.
*   **`07_Kmer_Mash.sh`**: Fast distance matrix estimation from raw reads.
*   **`08_IBS_pairwise_snp_divergence.sh`**: Calculates pairwise absolute SNP difference between samples.
*   **`09_AB_Histograms_Stats.sh`**: Calculates Allele Balance (AB) in offspring compared to mother and shows the distribution of AB in histograms with a modelled curve.
*   **`10A_LOH_AB_ZooSuyapa.sh`** & **`10B_LOH_AB_UNAH.sh`**: Generate Manhattan plots from normalized heterozygosity (AB) and concordance statistics for all samples.
*   **`11_Calc_BAF.sh`**: Extraction of allelic depths for B-Allele Frequency analysis.

### Phase 3: Phylogenomics & Synteny
*   **`12A_Mitochondrial_RAxML_tree.sh`** & **`12B_Nuclear_RAxML_tree.sh`**: Pipeline for Mitochondrial and Nuclear gene extraction, alignments and computation of ML trees with RAxML.
*   **`13A_TargetCapture_BUSCO.sh`** & **`13B_MAFFT_IQTREE_ML_GCF_trees.sh`** & `**13C_Filtered_loci_check.sh`**: Target capture 2000 random BUSCO loci (CDS), align, filter and generate ML and concordance trees and check filtered genes.
*   **`14A_Synteny_Alignment.sh` & `14B_Synteny_JCVI.sh`**: 1:2 ortholog filtering and monkey-patched 3-track JCVI karyotype visualization.

---

## Data Visualisation (R Scripts)

Following the upstream bash pipeline, data visualisation is handled locally using a suite of custom R scripts. These scripts generate the following publication figures:

*   **BAF Manhattan Plot:** Maps the whole-genome B-Allele Frequency landscape across chromosomes using `ggpointdensity`.
*   **Haplotype Networks:** Builds mitochondrial haplotype networks using `ape` and `pegas`, mapping taxonomic sample frequencies and capping extreme mutational steps.
*   **ML and eGCF plots:** Make plots for the uncollapsed Maximum Likelihood tree and collapsed effective Gene Concordance tree.
*   **Parental Coverage:** Creates a Manhattan plot visualising the parental mapping coverage  across the 36 chromosomes.
*   **Mitochondrial Topology:** Plots the RAxML topology for the mitochondrial markers
*   **Nuclear Topology**:** Plots the RAxML topology for the nuclear markers
*   **Pairwise Matrix:** Plots the pairwise SNP differences (absolute and 1-IBS) into a matrix.
*   **Distribution map for the genus Laemanctus** Generate the distribution map for the genus Laemanctus based on iNaturalist and own observations.
---

## Data Availability
All sequencing data for this project can be found under NCBI BioProject/SRA PRJNA1507630 and the genome assembly can be found under NCBI BioProject/SRA PRJNA1514599. 

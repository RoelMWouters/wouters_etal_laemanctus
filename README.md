# wouters_etal_laemanctus
Bioinformatic pipelines and scripts accompanying Wouters et al., (2026). Field discovery and hybrid origin of the largest all-female clonal vertebrate, the casquehead lizard Laemanctus julioi

# Laemanctus Phylogenomic & Hybrid Assembly Pipeline

**Author:** Roel Wouters | Charles University 
**Description:** A comprehensive bioinformatics pipeline for *de novo* genome assembly, Hi-C scaffolding, population genomics, and species tree inference in *Laemanctus* lizards.

This repository contains the exact bash and R scripts used to process raw sequencing reads (ONT/Illumina/Hi-C) through to final publication-quality figures and phylogenomic topologies.

## Pipeline Overview

The pipeline is designed for execution on a PBS-scheduled High-Performance Computing (HPC) cluster (e.g., MetaCentrum), utilising strict environment isolation and modular job arrays. It is divided into three primary phases:

### Phase 1: Genome Assembly & Scaffolding
*   **`00_HiC_Trimming.sh`**: Pre-processing and adapter trimming of raw Hi-C reads.
*   **`01_Hifiasm_Assembly.sh`**: *De novo* assembly of ONT long reads phased with Hi-C data.
*   **`02A_YAHS_Scaffolding.sh`**: Scaffolding of uncollapsed assemblies for Juicebox curation.
*   **`02B_HapHiC_Scaffolding.sh`**: Allele-aware scaffolding for highly heterozygous genomes.
*   **`06_Subgenome_Phasing.sh`**: Segregation of draft assemblies into distinct parental subgenomes based on mapping bias.

### Phase 2: Population Genomics & Variant Calling
*   **`03_FastP_BWA_Alignment.sh`**: Job array for trimming Illumina short reads and mapping to the curated reference.
*   **`04_GVCF_Generation.sh`**: Highly parallelized GATK HaplotypeCaller scatter-gather execution.
*   **`05_Joint_Calling.sh`**: Population joint genotyping and strict paralog purging based on depth thresholds.
*   **`07_Kmer_Mash.sh`**: Fast distance matrix estimation from raw reads.
*   **`08_PLINK_Stats.sh`**: LD pruning, PCA, and Identity-By-Descent calculations.
*   **`09_Calc_Fst.sh`**: Sliding-window Weir & Cockerham Fst calculations.
*   **`10_Calc_BAF.sh`**: Extraction of allelic depths for B-Allele Frequency analysis.

### Phase 3: Phylogenomics & Synteny
*   **`11A_Mito_Bait_BLAST.sh` & `11B_Extract_Mito.sh`**: Fetching, blasting, and extracting mitochondrial markers from BAM files.
*   **`12A_Nuclear_Bait_Blast.sh` & `12B_Nuclear_Marker_Extract.sh`**: Hybrid Python/Bash coordinate calculation and allele extraction via `bcftools consensus`.
*   **`13A_Supermatrix_Builder.sh`**: MAFFT alignment and native Python supermatrix concatenation.
*   **`13B_RAxML_Phylogeny.sh`**: Maximum Likelihood tree inference with Rapid Bootstrapping.
*   **`14_BUSCO_Gene_Trees.sh`**: Extraction, ghost-allele purging, and ML gene tree inference for 300 target BUSCO loci.
*   **`15_Phylogenomic_Pipeline.sh`**: Unified ASTRAL coalescent species tree inference and majority-rule consensus generation.
*   **`16A_Synteny_Alignment.sh` & `16B_Synteny_JCVI.sh`**: 1:2 ortholog filtering and monkey-patched 3-track JCVI karyotype visualization.

---

## Data Visualisation (R Scripts)

Following the upstream bash pipeline, data visualisation is handled locally using a suite of custom R scripts[cite: 1]. These scripts generate the following publication figures:

*   **Clonality & Genetic Distance Plots:** Generates a PI_HAT Identity-by-Descent heatmap and an Identity-by-State (IBS) distance matrix using `ggplot2`[cite: 1]. It also produces an upgraded PCA scatterplot to distinguish lineages[cite: 1].
*   **Fst Density:** Visualises the distribution of genomic divergence across 100kb windows between populations[cite: 1].
*   **BAF Manhattan Plot:** Maps the whole-genome B-Allele Frequency landscape across chromosomes using `ggpointdensity`[cite: 1].
*   **Parental Coverage:** Creates a bar chart visualising the parental mapping coverage percentages across the top 50 contigs[cite: 1].
*   **Phylogenetic Topologies:** Visualises RAxML trees with bootstrap support using `ggtree`, dynamically truncating extremely long outgroup branches[cite: 1].
*   **Haplotype Networks:** Builds mitochondrial haplotype networks using `ape` and `pegas`, mapping taxonomic sample frequencies and capping extreme mutational steps[cite: 1].

---

## Data Availability
*(Insert Zenodo DOI or NCBI SRA accession numbers here prior to publication).*

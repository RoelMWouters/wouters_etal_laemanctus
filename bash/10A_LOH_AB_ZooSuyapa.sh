#!/bin/bash
#PBS -N Concordance_LOH_Manhattan
#PBS -l select=1:ncpus=4:mem=64gb:scratch_local=100gb
#PBS -l walltime=02:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SUYAPA/ZOO CONCORDANCE & LOH
# Description: Extracts VCF data for the Zoo and Suyapa clonal lineages, computes a 
# Concordance matrix (Clonal/LOH/Discordant) from bi-allelic read depths normalized 
# to the respective mothers and maps structural LOH via 1Mb-windowed Manhattan plots.
# ==============================================================================

set -euo pipefail
export TMPDIR=$SCRATCHDIR
cd "$SCRATCHDIR" || exit 1

# ==============================================================================
# 1. GLOBAL PATHS & CONFIGURATION
# ==============================================================================
PROJECT_DIR="/path/to/project"
VCF_FILE="${PROJECT_DIR}/bam_output/haploid_bam/loh_ljul/vcfs/ljul.heterozygosity.vcf.gz" 
FAI_FILE="${PROJECT_DIR}/genomes/Genome_C09.fasta.fai"
OUTPUT_DIR="${PROJECT_DIR}/output/clonal_inheritance/suyapa_zoo_loh"

mkdir -p "$OUTPUT_DIR"
STATS_FILE="filtering_stats.txt"

# ==============================================================================
# 2. ISOLATED CONDA ENVIRONMENT
# ==============================================================================
echo "=== INITIALIZING ISOLATED ENVIRONMENT ==="
module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
export CONDA_ENVS_DIRS="$SCRATCHDIR/conda_envs"

mamba create -n loh_env -c conda-forge -c bioconda r-base r-ggplot2 r-data.table bcftools -y
eval "$(conda shell.bash hook)"
conda activate loh_env

# ==============================================================================
# 3. WRITE R SCRIPT 1: CONCORDANCE MATRIX
# ==============================================================================
cat << 'EOF' > "$SCRATCHDIR/calculate_concordance.R"
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]

parse_ad <- function(ad_vec) {
    parsed <- tstrsplit(ad_vec, ",", type.convert=TRUE, fill=0)
    if(length(parsed) < 2) parsed[[2]] <- rep(0, length(ad_vec))
    ref <- parsed[[1]]; alt1 <- parsed[[2]]
    if(length(parsed) > 2) { other <- rowSums(as.data.table(parsed)[, 3:length(parsed), with=FALSE], na.rm=TRUE)
    } else { other <- rep(0, length(ad_vec)) }
    return(list(ref, alt1, other))
}

process_lineage <- function(tsv_file, mother_name, offspring_names) {
    num_offspring <- length(offspring_names)
    cc <- c("character", "integer", "character", rep("character", num_offspring))
    df <- fread(tsv_file, sep="\t", header=FALSE, colClasses=cc)
    colnames(df) <- c("CHROM", "POS", "MOTHER_AD", paste0("O_AD_", 1:num_offspring))
    
    df <- df[!(CHROM %in% c("ptg000008l", "ptg000012l", "ptg000024l"))]
    df <- df[!(CHROM == "ptg000005l" & POS > 278731350)]
    
    df[, c("M_REF", "M_ALT1", "M_OTHER") := parse_ad(MOTHER_AD)]
    df[, M_DP := M_REF + M_ALT1 + M_OTHER]
    df[, M_AB_REF := M_REF / M_DP]
    df[, M_AB_ALT1 := M_ALT1 / M_DP]
    df[, M_AB_OTHER := M_OTHER / M_DP]
    
    df <- df[!is.na(M_DP) & M_DP > 0]
    df <- df[M_AB_REF >= 0.3 & M_AB_REF <= 0.7 & M_AB_ALT1 >= 0.3 & M_AB_ALT1 <= 0.7 & M_AB_OTHER < 0.05]
    
    concordance_list <- list()
    for(i in 1:num_offspring) {
        ad_col <- paste0("O_AD_", i)
        df[, c("O_REF", "O_ALT1", "O_OTHER") := parse_ad(get(ad_col))]
        df[, O_DP := O_REF + O_ALT1 + O_OTHER]
        df[O_DP == 0, O_DP := NA] 
        
        df[, Fate := fcase(
            is.na(O_DP), "No_Coverage",
            O_OTHER > 0, "Biological_Impossibility_Novel_Allele",
            O_REF > 0 & O_ALT1 > 0, "Concordant_Retention",
            O_REF > 0 & O_ALT1 == 0, "Putative_LOH_to_Reference",
            O_REF == 0 & O_ALT1 > 0, "Putative_LOH_to_Alternate",
            default = "Ambiguous"
        )]
        
        stats <- df[, .N, by=Fate]
        stats[, Offspring := offspring_names[i]]
        stats[, Maternal_Lineage := mother_name]
        concordance_list[[i]] <- stats
    }
    return(rbindlist(concordance_list))
}

master_stats <- rbindlist(list(
    process_lineage("F_lineage.tsv", "F01", c("F03")),
    process_lineage("C_lineage.tsv", "C04", c("C05", "C06"))
))

concordance_matrix <- dcast(master_stats, Maternal_Lineage + Offspring ~ Fate, value.var="N", fill=0)
expected_fates <- c("Concordant_Retention", "Putative_LOH_to_Reference", "Putative_LOH_to_Alternate", "Biological_Impossibility_Novel_Allele", "No_Coverage")
for (f in expected_fates) { if (!(f %in% names(concordance_matrix))) concordance_matrix[[f]] <- 0 }

concordance_matrix[, Total_Validated_Maternal_SNPs := rowSums(.SD), .SDcols = expected_fates]
concordance_matrix[, Percent_Clonal := round((Concordant_Retention / Total_Validated_Maternal_SNPs) * 100, 3)]
concordance_matrix[, Percent_LOH := round(((Putative_LOH_to_Reference + Putative_LOH_to_Alternate) / Total_Validated_Maternal_SNPs) * 100, 3)]
concordance_matrix[, Percent_Error_Baseline := round((Biological_Impossibility_Novel_Allele / Total_Validated_Maternal_SNPs) * 100, 3)]

write.csv(concordance_matrix, file.path(out_dir, "Offspring_TriState_Error_Matrix.csv"), row.names=FALSE)
EOF

# ==============================================================================
# 4. WRITE R SCRIPT 2: LOH MANHATTAN PLOTS
# ==============================================================================
cat << 'EOF' > "$SCRATCHDIR/plot_manhattan_only.R"
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
fai_path <- args[2]

window_size <- 1000000
col_F_prim <- "#bc3933"; col_F_shad <- "#7a2420"
col_C_prim <- "#6444a3"; col_C_shad <- "#3a2760"

fai <- fread(fai_path, col.names=c("CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"))
fai <- fai[!(CHROM %in% c("ptg000008l", "ptg000012l", "ptg000024l"))]
fai[CHROM == "ptg000005l", LENGTH := 278731350]
fai <- fai[LENGTH >= window_size]
fai <- fai[order(-LENGTH)]
fai[, cum_len := cumsum(as.numeric(LENGTH)) - LENGTH]

process_manhattan <- function(tsv_file, mother_name, offspring_names, lineage_tag) {
    num_offspring <- length(offspring_names)
    cc <- c("character", "integer", "character", rep("character", num_offspring))
    df <- fread(tsv_file, sep="\t", header=FALSE, colClasses=cc)
    colnames(df) <- c("CHROM", "POS", "MOTHER_AD", paste0("O_AD_", 1:num_offspring))
    
    df <- df[!(CHROM %in% c("ptg000008l", "ptg000012l", "ptg000024l"))]
    df <- df[!(CHROM == "ptg000005l" & POS > 278731350)]
    df <- df[CHROM %in% fai$CHROM]
    
    df[, c("M_REF", "M_ALT") := tstrsplit(MOTHER_AD, ",", type.convert=TRUE)]
    df[, M_AB_raw := M_ALT / (M_REF + M_ALT)]
    df <- df[!is.na(M_AB_raw) & M_AB_raw >= 0.3 & M_AB_raw <= 0.7]
    df[, M_Het_State := 0.5] 
    
    for(i in 1:num_offspring) {
        ad_col <- paste0("O_AD_", i)
        df[, c("O_REF", "O_ALT") := tstrsplit(get(ad_col), ",", type.convert=TRUE)]
        df[, paste0("Het_State_", offspring_names[i]) := ifelse(O_REF > 0 & O_ALT > 0, 0.5, 0.0)]
    }
    
    df[, Bin := floor(POS / window_size) * window_size + (window_size / 2)]
    het_cols <- c("M_Het_State", paste0("Het_State_", offspring_names))
    windowed <- df[, lapply(.SD, mean, na.rm=TRUE), by=.(CHROM, Bin), .SDcols=het_cols]
    setnames(windowed, het_cols, c(paste("Mother", mother_name), paste("Offspring", offspring_names)))
    
    man_melt <- melt(windowed, id.vars=c("CHROM", "Bin"), measure.vars=c(paste("Mother", mother_name), paste("Offspring", offspring_names)), variable.name="Individual", value.name="Mean_AB")
    man_melt[, Lineage := lineage_tag]
    return(man_melt)
}

F_man <- process_manhattan("F_lineage.tsv", "F01", c("F03"), "F")
C_man <- process_manhattan("C_lineage.tsv", "C04", c("C05", "C06"), "C")
master_man <- rbindlist(list(F_man, C_man))

master_man <- merge(master_man, fai[, .(CHROM, cum_len)], by="CHROM")
master_man[, abs_pos := Bin + cum_len]
master_man$CHROM <- factor(master_man$CHROM, levels = fai$CHROM)
master_man[, color_group := as.numeric(CHROM) %% 2]

master_man[, ColorKey := paste0(Lineage, "_", color_group)]
man_colors <- c("F_0" = col_F_prim, "F_1" = col_F_shad, "C_0" = col_C_prim, "C_1" = col_C_shad)
master_man$Individual <- factor(master_man$Individual, levels=c("Mother F01", "Offspring F03", "Mother C04", "Offspring C05", "Offspring C06"))

p2 <- ggplot(master_man, aes(x=abs_pos, y=Mean_AB, color=ColorKey, group=CHROM)) +
  geom_point(size=1.2, alpha=0.8) +
  geom_line(linewidth=0.6, alpha=0.5) +
  facet_grid(Individual ~ .) +
  scale_color_manual(values=man_colors) +
  scale_y_continuous(limits=c(-0.02, 0.52), breaks=seq(0, 0.5, 0.1)) +
  theme_classic(base_size = 14) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(r = 15)),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        strip.text = element_text(face="bold", size=14, angle=270), strip.background = element_blank(),
        panel.spacing = unit(1, "lines")) +
  labs(x="Scaffolds (Ordered by Length)", y="Heterozygosity Retention\n(0 = LOH, 0.5 = Clonal)")

ggsave(file.path(out_dir, "Master_LOH_Manhattan.pdf"), plot=p2, width=14, height=7)
ggsave(file.path(out_dir, "Master_LOH_Manhattan.png"), plot=p2, width=14, height=7, dpi=300)
EOF

# ==============================================================================
# 5. EXECUTE EXTRACTION & R SCRIPTS
# ==============================================================================
echo "=== SNP FILTERING STATISTICS (MULTI-ALLELIC ENABLED) ===" > "$STATS_FILE"

echo "=== EXTRACTING LINEAGE F (F01 -> F03) ==="
RAW_F=$(bcftools view -s "F01,F03" -m2 -M2 -v snps "$VCF_FILE" -H | wc -l)
bcftools view -s "F01,F03" -m2 -M2 -v snps "$VCF_FILE" -O u | \
bcftools filter -i "QUAL >= 30 && MQ >= 40 && MIN(FORMAT/DP) >= 10 && MAX(FORMAT/DP) <= 60" -O u | \
bcftools query -f '%CHROM\t%POS[\t%AD]\n' > "F_lineage.tsv"
FILT_F=$(wc -l < "F_lineage.tsv")

echo "Lineage F (F01, F03):" >> "$STATS_FILE"
echo "  Raw biallelic SNPs: $RAW_F" >> "$STATS_FILE"
echo "  Retained after quality filters: $FILT_F" >> "$STATS_FILE"

echo "=== EXTRACTING LINEAGE C (C04 -> C05, C06) ==="
RAW_C=$(bcftools view -s "C04,C05,C06" -m2 -M2 -v snps "$VCF_FILE" -H | wc -l)
bcftools view -s "C04,C05,C06" -m2 -M2 -v snps "$VCF_FILE" -O u | \
bcftools filter -i "QUAL >= 30 && MQ >= 40 && MIN(FORMAT/DP) >= 10 && MAX(FORMAT/DP) <= 60" -O u | \
bcftools query -f '%CHROM\t%POS[\t%AD]\n' > "C_lineage.tsv"
FILT_C=$(wc -l < "C_lineage.tsv")

echo -e "\nLineage C (C04, C05, C06):" >> "$STATS_FILE"
echo "  Raw biallelic SNPs: $RAW_C" >> "$STATS_FILE"
echo "  Retained after quality filters: $FILT_C" >> "$STATS_FILE"

echo "=== GENERATING CONCORDANCE MATRIX ==="
Rscript calculate_concordance.R "$OUTPUT_DIR"

echo "=== GENERATING LOH MANHATTAN PLOTS ==="
Rscript plot_manhattan_only.R "$OUTPUT_DIR" "$FAI_FILE"

# ==============================================================================
# 6. EXPORT
# ==============================================================================
cp "F_lineage.tsv" "C_lineage.tsv" "$STATS_FILE" "$OUTPUT_DIR/"
echo "Part 5A pipeline completed successfully!"

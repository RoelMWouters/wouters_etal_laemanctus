#!/bin/bash
#PBS -N UNAH_Mock_Pedigree_LOH
#PBS -l select=1:ncpus=4:mem=64gb:scratch_local=100gb
#PBS -l walltime=00:30:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# UNAH MOCK PEDIGREE
# Description: Evaluates the UNAH siblings by establishing each clone as a "mock" 
# mother to track theoretical inheritance/LOH across the other two. Generates 
# a 6-way concordance matrix and 3 LOH Manhattan plots.
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
OUTPUT_DIR="${PROJECT_DIR}/output/clonal_inheritance/unah_loh"

mkdir -p "$OUTPUT_DIR"

# ==============================================================================
# 2. ISOLATED CONDA ENVIRONMENT
# ==============================================================================
module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
export CONDA_ENVS_DIRS="$SCRATCHDIR/conda_envs"

mamba create -n loh_env -c conda-forge -c bioconda r-base r-ggplot2 r-data.table bcftools -y
eval "$(conda shell.bash hook)"
conda activate loh_env

# ==============================================================================
# 3. WRITE THE R SCRIPT (CONCORDANCE + MANHATTAN)
# ==============================================================================
cat << 'EOF' > "$SCRATCHDIR/unah_tri_directional_no_hist.R"
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
fai_path <- args[2]

window_size <- 1000000 
col_U_prim <- "#E69F00"; col_U_shad <- "#a87400" 

fai <- fread(fai_path, col.names=c("CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"))
fai <- fai[!(CHROM %in% c("ptg000008l", "ptg000012l", "ptg000024l"))]
fai[CHROM == "ptg000005l", LENGTH := 278731350]
fai <- fai[LENGTH >= window_size]
fai <- fai[order(-LENGTH)]
fai[, cum_len := cumsum(as.numeric(LENGTH)) - LENGTH]

parse_ad <- function(ad_vec) {
  parsed <- tstrsplit(ad_vec, ",", type.convert=TRUE, fill=0)
  if(length(parsed) < 2) parsed[[2]] <- rep(0, length(ad_vec))
  ref <- parsed[[1]]; alt1 <- parsed[[2]]
  if(length(parsed) > 2) { other <- rowSums(as.data.table(parsed)[, 3:length(parsed), with=FALSE], na.rm=TRUE)
  } else { other <- rep(0, length(ad_vec)) }
  return(list(ref, alt1, other))
}

run_mock_pedigree <- function(master_df, mother_id, offspring_ids, exp_name) {
  cat("\n--- Running Experiment:", exp_name, "---\n")
  cols_to_keep <- c("CHROM", "POS", paste0(mother_id, "_AD"), paste0(offspring_ids, "_AD"))
  df <- master_df[, ..cols_to_keep]
  setnames(df, cols_to_keep, c("CHROM", "POS", "MOTHER_AD", "O1_AD", "O2_AD"))
  
  df[, c("M_REF", "M_ALT1", "M_OTHER") := parse_ad(MOTHER_AD)]
  df[, M_DP := M_REF + M_ALT1 + M_OTHER]
  df[, M_AB_REF := M_REF / M_DP]
  df[, M_AB_ALT1 := M_ALT1 / M_DP]
  df[, M_AB_OTHER := M_OTHER / M_DP]
  
  df <- df[!is.na(M_DP) & M_DP > 0]
  df <- df[M_AB_REF >= 0.3 & M_AB_REF <= 0.7 & M_AB_ALT1 >= 0.3 & M_AB_ALT1 <= 0.7 & M_AB_OTHER < 0.05]
  df[, M_Het_State := 0.5]
  
  concordance_list <- list()
  for(i in 1:2) {
    ad_col <- paste0("O", i, "_AD")
    off_name <- offspring_ids[i]
    
    df[, c("O_REF", "O_ALT1", "O_OTHER") := parse_ad(get(ad_col))]
    df[, paste0("O_DP_", off_name) := O_REF + O_ALT1 + O_OTHER]
    df[, paste0("Het_State_", off_name) := ifelse(O_REF > 0 & O_ALT1 > 0, 0.5, 0.0)]
    
    df[, Fate := fcase(
      get(paste0("O_DP_", off_name)) == 0, "No_Coverage",
      O_OTHER > 0, "Biological_Impossibility_Novel_Allele",
      O_REF > 0 & O_ALT1 > 0, "Concordant_Retention",
      O_REF > 0 & O_ALT1 == 0, "Putative_LOH_to_Reference",
      O_REF == 0 & O_ALT1 > 0, "Putative_LOH_to_Alternate",
      default = "Ambiguous"
    )]
    
    stats <- df[, .N, by=Fate]
    stats[, Mock_Offspring := off_name]
    stats[, Mock_Mother := mother_id]
    concordance_list[[i]] <- stats
  }
  
  # MANHATTAN PLOT GENERATION
  df[, Bin := floor(POS / window_size) * window_size + (window_size / 2)]
  het_cols <- c("M_Het_State", paste0("Het_State_", offspring_ids))
  
  windowed <- df[, c(lapply(.SD, mean, na.rm=TRUE), list(SNP_Count = .N)), by=.(CHROM, Bin), .SDcols=het_cols]
  median_snps <- median(windowed$SNP_Count, na.rm=TRUE)
  windowed[SNP_Count < (median_snps * 0.10), (het_cols) := NA]
  
  setnames(windowed, het_cols, c(paste("Mother", mother_id), paste("Offspring", offspring_ids)))
  man_melt <- melt(windowed, id.vars=c("CHROM", "Bin"), measure.vars=c(paste("Mother", mother_id), paste("Offspring", offspring_ids)), variable.name="Individual", value.name="Mean_AB")
  
  master_man <- merge(man_melt, fai[, .(CHROM, cum_len)], by="CHROM")
  master_man[, abs_pos := Bin + cum_len]
  master_man$CHROM <- factor(master_man$CHROM, levels = fai$CHROM)
  master_man[, color_group := as.numeric(CHROM) %% 2]
  master_man[, ColorKey := paste0("UNAH_", color_group)]
  
  p_man <- ggplot(master_man, aes(x=abs_pos, y=Mean_AB, color=ColorKey, group=CHROM)) +
    geom_point(size=1.2, alpha=0.8) +
    geom_line(linewidth=0.6, alpha=0.5) +
    facet_grid(Individual ~ .) +
    scale_color_manual(values=c("UNAH_0"=col_U_prim, "UNAH_1"=col_U_shad)) +
    scale_y_continuous(limits=c(-0.02, 0.52), breaks=seq(0, 0.5, 0.1)) +
    theme_classic(base_size = 14) +
    theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), strip.text=element_text(face="bold", size=14, angle=270), strip.background=element_blank()) +
    labs(x="Scaffolds (Ordered by Length)", y="Heterozygosity Retention", title=paste("LOH Architecture (Baseline:", mother_id, ")"))
  
  ggsave(file.path(out_dir, paste0(exp_name, "_Manhattan.png")), plot=p_man, width=14, height=7, dpi=300)
  ggsave(file.path(out_dir, paste0(exp_name, "_Manhattan.pdf")), plot=p_man, width=14, height=7, dpi=300)
  
  return(rbindlist(concordance_list))
}

cc <- c("character", "integer", "character", "character", "character")
master_df <- fread("UNAH_lineage.tsv", sep="\t", header=FALSE, colClasses=cc)
colnames(master_df) <- c("CHROM", "POS", "UNAH1_AD", "UNAH2_AD", "UNAH3_AD")

master_df <- master_df[CHROM %in% fai$CHROM]
master_df <- master_df[!(CHROM == "ptg000005l" & POS > 278731350)]

stats_exp1 <- run_mock_pedigree(master_df, "UNAH1", c("UNAH2", "UNAH3"), "Exp1_U1_Mother")
stats_exp2 <- run_mock_pedigree(master_df, "UNAH2", c("UNAH1", "UNAH3"), "Exp2_U2_Mother")
stats_exp3 <- run_mock_pedigree(master_df, "UNAH3", c("UNAH1", "UNAH2"), "Exp3_U3_Mother")

all_stats <- rbindlist(list(stats_exp1, stats_exp2, stats_exp3))
concordance_matrix <- dcast(all_stats, Mock_Mother + Mock_Offspring ~ Fate, value.var="N", fill=0)

expected_fates <- c("Concordant_Retention", "Putative_LOH_to_Reference", "Putative_LOH_to_Alternate", "Biological_Impossibility_Novel_Allele", "No_Coverage")
for (f in expected_fates) { if (!(f %in% names(concordance_matrix))) concordance_matrix[[f]] <- 0 }

concordance_matrix[, Total_Validated_Maternal_SNPs := rowSums(.SD), .SDcols = expected_fates]
concordance_matrix[, Percent_Clonal := round((Concordant_Retention / Total_Validated_Maternal_SNPs) * 100, 3)]
concordance_matrix[, Percent_LOH := round(((Putative_LOH_to_Reference + Putative_LOH_to_Alternate) / Total_Validated_Maternal_SNPs) * 100, 3)]
concordance_matrix[, Percent_Discordant := round((Biological_Impossibility_Novel_Allele / Total_Validated_Maternal_SNPs) * 100, 3)]

write.csv(concordance_matrix, file.path(out_dir, "UNAH_6Way_Concordance_Matrix.csv"), row.names=FALSE)
EOF

# ==============================================================================
# 4. EXECUTE PIPELINE & CAPTURE STATS WITH CHECKPOINT
# ==============================================================================
STATS_FILE="unah_tri_directional_stats.txt"
TSV_FILE="UNAH_lineage.tsv"
PERMANENT_TSV="$OUTPUT_DIR/$TSV_FILE"

echo "=== UNAH LINEAGE PRE-FILTERING (C01, C02, C03) ===" > "$STATS_FILE"

if [ -f "$PERMANENT_TSV" ]; then
    echo "Found existing $TSV_FILE in $OUTPUT_DIR. Copying to scratch..."
    cp "$PERMANENT_TSV" "$SCRATCHDIR/"
else
    echo "Extracting UNAH lineage from VCF..."
    RAW_U=$(bcftools view -s "C01,C02,C03" -v snps "$VCF_FILE" -H | wc -l)
    
    bcftools view -s "C01,C02,C03" -v snps "$VCF_FILE" -O u | \
    bcftools filter -i "QUAL >= 30 && MQ >= 40 && MIN(FORMAT/DP) >= 10 && MAX(FORMAT/DP) <= 60" -O u | \
    bcftools query -f '%CHROM\t%POS[\t%AD]\n' > "$TSV_FILE"

    FILT_U=$(wc -l < "$TSV_FILE")
    echo "Raw SNPs: $RAW_U" >> "$STATS_FILE"
    echo "Retained after quality filters: $FILT_U" >> "$STATS_FILE"
    
    cp "$TSV_FILE" "$OUTPUT_DIR/"
fi

echo "=== EXECUTING 3-WAY MOCK PEDIGREE ==="
Rscript unah_tri_directional_no_hist.R "$OUTPUT_DIR" "$FAI_FILE"

# ==============================================================================
# 5. EXPORT
# ==============================================================================
cp "$STATS_FILE" "$OUTPUT_DIR/"
echo "pipeline completed successfully!"

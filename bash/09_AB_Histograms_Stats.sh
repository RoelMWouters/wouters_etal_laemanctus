#!/bin/bash
#PBS -N BetaBinomial_Mixture_Histograms
#PBS -l select=1:ncpus=2:mem=128gb:scratch_local=100gb
#PBS -l walltime=10:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# BETA-BINOMIAL HISTOGRAMS
# Description: Takes the exact lineage TSVs as in script 10A and 10B, extracts the 
# bi-allelic loci normalized to the respecitve (mock) mothers, fits a Beta-Binomial mixture model to track 
# allele balance, and plots the resulting density histograms.
# ==============================================================================

set -euo pipefail
export TMPDIR=$SCRATCHDIR
cd "$SCRATCHDIR" || exit 1

# ==============================================================================
# 1. GLOBAL PATHS & CONFIGURATION
# ==============================================================================
PROJECT_DIR="/path/to/project"
DIR_SUYAPA_ZOO="${PROJECT_DIR}/output/clonal_inheritance/suyapa_zoo_loh"
DIR_UNAH="${PROJECT_DIR}/output/clonal_inheritance/unah_loh"
OUTPUT_DIR="${PROJECT_DIR}/output/clonal_inheritance/beta_binomial_histograms"

mkdir -p "$OUTPUT_DIR"

echo "=== COPYING VCF TSV DATA TO SCRATCH ==="
cp "$DIR_SUYAPA_ZOO/F_lineage.tsv" "$DIR_SUYAPA_ZOO/C_lineage.tsv" "$DIR_UNAH/UNAH_lineage.tsv" ./

# ==============================================================================
# 2. ISOLATED CONDA ENVIRONMENT
# ==============================================================================
echo "=== INITIALIZING ISOLATED ENVIRONMENT ==="
module add mambaforge
export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
export CONDA_ENVS_DIRS="$SCRATCHDIR/conda_envs"

mamba create -n mixture_env -c conda-forge r-base r-data.table r-ggplot2 -y
eval "$(conda shell.bash hook)"
conda activate mixture_env

# ==============================================================================
# 3. WRITE FULL PAIRWISE EXTRACTION SCRIPT
# ==============================================================================
cat << 'EOF' > "$SCRATCHDIR/extract_unbiased_ALL.R"
library(data.table)

extract_offspring_binomial_full <- function(file_path, m_col, o_col, m_name, o_name) {
  cat("Extracting Target:", o_name, "(using Baseline:", m_name, ")...\n")
  
  df <- fread(file_path, sep="\t", header=FALSE, select=c(1, 2, m_col, o_col), colClasses="character")
  colnames(df) <- c("CHROM", "POS", "M_AD", "O_AD")
  df[, POS := as.integer(POS)]
  
  m_parsed <- tstrsplit(df$M_AD, ",", type.convert=FALSE, fill="0")
  df[, M_A := as.numeric(m_parsed[[1]])]
  if(length(m_parsed) >= 2) { df[, M_B := as.numeric(m_parsed[[2]])] } else { df[, M_B := 0] }
  df[, M_DP := M_A + M_B]
  df[, M_AB := M_B / M_DP]
  
  o_parsed <- tstrsplit(df$O_AD, ",", type.convert=FALSE, fill="0")
  df[, O_A := as.numeric(o_parsed[[1]])]
  if(length(o_parsed) >= 2) { df[, O_B := as.numeric(o_parsed[[2]])] } else { df[, O_B := 0] }
  df[, O_DP := O_A + O_B]
  df[, O_AB := O_B / O_DP]
  
  valid_loci <- df[!is.na(M_DP) & M_DP >= 10 & M_DP <= 60 & M_AB >= 0.3 & M_AB <= 0.7 & !is.na(O_DP) & O_DP >= 1]
  
  rel_name <- paste0(o_name, " (ref. ", m_name, ")")
  result <- valid_loci[, .(Individual = rel_name, DP = O_DP, A_Alleles = O_A, B_Alleles = O_B, Allele_Balance = O_AB)]
  
  return(result)
}

f03 <- extract_offspring_binomial_full("F_lineage.tsv", 3, 4, "F01", "F03")
c05 <- extract_offspring_binomial_full("C_lineage.tsv", 3, 4, "C04", "C05")
c06 <- extract_offspring_binomial_full("C_lineage.tsv", 3, 5, "C04", "C06")

u1_2 <- extract_offspring_binomial_full("UNAH_lineage.tsv", 3, 4, "UNAH1", "UNAH2")
u1_3 <- extract_offspring_binomial_full("UNAH_lineage.tsv", 3, 5, "UNAH1", "UNAH3")
u2_1 <- extract_offspring_binomial_full("UNAH_lineage.tsv", 4, 3, "UNAH2", "UNAH1")
u2_3 <- extract_offspring_binomial_full("UNAH_lineage.tsv", 4, 5, "UNAH2", "UNAH3")
u3_1 <- extract_offspring_binomial_full("UNAH_lineage.tsv", 5, 3, "UNAH3", "UNAH1")
u3_2 <- extract_offspring_binomial_full("UNAH_lineage.tsv", 5, 4, "UNAH3", "UNAH2")

master_binomial <- rbindlist(list(f03, c05, c06, u1_2, u1_3, u2_1, u2_3, u3_1, u3_2))
fwrite(master_binomial, "Master_Unbiased_Binomial_ALL.csv")
EOF

# ==============================================================================
# 4. WRITE OPTIMIZED STATISTICAL MODEL SCRIPT
# ==============================================================================
cat << 'EOF' > "$SCRATCHDIR/colleague_mixture_model_optimized.R"
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]

dat <- fread("Master_Unbiased_Binomial_ALL.csv")

dat$Individual <- as.character(dat$Individual)
dat$DP <- as.integer(dat$DP)
dat$A_Alleles <- as.integer(dat$A_Alleles)
dat$B_Alleles <- as.integer(dat$B_Alleles)
dat$Allele_Balance <- as.numeric(dat$Allele_Balance)

inv_scaled <- function(x, lo, hi) {
  p <- (x-lo)/(hi-lo)
  p <- pmin(pmax(p,1e-8),1-1e-8)
  qlogis(p)
}

logspace_add <- function(a,b) {
  m <- pmax(a,b)
  m + log(exp(a-m)+exp(b-m))
}

dbetabinom_log <- function(k,n,alpha,beta) {
  lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) + lbeta(k+alpha,n-k+beta)-lbeta(alpha,beta)
}

unpack_mix <- function(par) {
  list(
    mu_minor  = 0.001 + 0.2490*plogis(par[1]),
    phi_minor = exp(par[2]),
    mu_main   = 0.30  + 0.40*plogis(par[3]),
    phi_main  = exp(par[4]),
    w_minor   = 0.0001 + 0.2999*plogis(par[5])
  )
}

nll_mix_opt <- function(par,k,n,counts) {
  p <- unpack_mix(par)
  a1 <- p$mu_minor*p$phi_minor; b1 <- (1-p$mu_minor)*p$phi_minor
  a2 <- p$mu_main*p$phi_main; b2 <- (1-p$mu_main)*p$phi_main
  l1 <- log(p$w_minor) + dbetabinom_log(k,n,a1,b1)
  l2 <- log(1-p$w_minor) + dbetabinom_log(k,n,a2,b2)
  -sum(counts * logspace_add(l1,l2))
}

unpack_single <- function(par) { list(mu=0.05+0.90*plogis(par[1]), phi=exp(par[2])) }

nll_single_opt <- function(par,k,n,counts) {
  p <- unpack_single(par)
  a <- p$mu*p$phi; b <- (1-p$mu)*p$phi
  -sum(counts * dbetabinom_log(k,n,a,b))
}

unpack_sym <- function(par) {
  list(
    mu_minor  = 0.001 + 0.2490*plogis(par[1]),
    phi_minor = exp(par[2]),
    phi_main  = exp(par[3]),
    w_minor   = 0.0001 + 0.2999*plogis(par[4])
  )
}

nll_sym_opt <- function(par,k,n,counts) {
  p <- unpack_sym(par)
  a1 <- p$mu_minor*p$phi_minor; b1 <- (1-p$mu_minor)*p$phi_minor
  a2 <- 0.5*p$phi_main; b2 <- 0.5*p$phi_main
  l1 <- log(p$w_minor) + dbetabinom_log(k,n,a1,b1)
  l2 <- log(1-p$w_minor) + dbetabinom_log(k,n,a2,b2)
  -sum(counts * logspace_add(l1,l2))
}

fit_one <- function(d) {
  N <- nrow(d)
  agg <- d[, .(Count = .N), by=.(B_Alleles, DP)]
  k <- agg$B_Alleles; n <- agg$DP; counts <- agg$Count
  
  starts <- expand.grid(mu_minor=c(0.02,0.05,0.10,0.15), mu_main=c(0.47,0.49,0.50), w_minor=c(0.005,0.01,0.02,0.05))
  best <- NULL
  for(i in seq_len(nrow(starts))) {
    s <- starts[i,]
    par0 <- c(inv_scaled(s$mu_minor,0.001,0.25), log(20), inv_scaled(s$mu_main,0.30,0.70), log(40), inv_scaled(s$w_minor,0.0001,0.30))
    fit <- optim(par0, nll_mix_opt, k=k, n=n, counts=counts, method="L-BFGS-B", lower=c(-15,-7,-15,-7,-15), upper=c(15,12,15,12,15), control=list(maxit=1000))
    if(is.null(best) || fit$value < best$value) best <- fit
  }
  p <- unpack_mix(best$par)
  logLik_mix <- -best$value
  
  fit_single <- optim(c(inv_scaled(0.48,0.05,0.95),log(30)), nll_single_opt, k=k, n=n, counts=counts, method="L-BFGS-B", lower=c(-15,-7),upper=c(15,12))
  logLik_single <- -fit_single$value
  
  par_sym0 <- c(inv_scaled(p$mu_minor,0.001,0.25), log(p$phi_minor), log(p$phi_main), inv_scaled(p$w_minor,0.0001,0.30))
  fit_sym <- optim(par_sym0, nll_sym_opt, k=k, n=n, counts=counts, method="L-BFGS-B", lower=c(-15,-7,-7,-15), upper=c(15,12,12,15))
  logLik_sym <- -fit_sym$value
  
  AIC_mix <- 2*5 - 2*logLik_mix
  BIC_mix <- 5*log(N) - 2*logLik_mix
  AIC_single <- 2*2 - 2*logLik_single
  BIC_single <- 2*log(N) - 2*logLik_single
  AIC_sym <- 2*4 - 2*logLik_sym
  LR_sym <- 2*(logLik_mix-logLik_sym)
  p_sym <- pchisq(max(LR_sym,0),df=1,lower.tail=FALSE)
  
  list(parameters=p, summary=data.frame(N=N, main_p=p$mu_main, main_phi=p$phi_main, main_rho=1/(p$phi_main+1), minor_weight=p$w_minor, minor_p=p$mu_minor, minor_phi=p$phi_minor, minor_rho=1/(p$phi_minor+1), delta_AIC_single_minus_mix=AIC_single-AIC_mix, delta_BIC_single_minus_mix=BIC_single-BIC_mix, delta_AIC_main_p_0.5_minus_free=AIC_sym-AIC_mix, LR_p_main_p_eq_0.5=p_sym, logLik=logLik_mix))
}

individuals <- unique(dat$Individual)
fit_objects <- list(); summary_list <- list()

for(id in individuals) {
  cat("\nFitting statistical models for:", id, "\n")
  d <- dat[dat$Individual==id,]
  ff <- fit_one(d)
  fit_objects[[id]] <- ff
  tmp <- ff$summary; tmp$Individual <- id
  summary_list[[id]] <- tmp
}

summary_table <- do.call(rbind,summary_list)
summary_table <- summary_table[,c("Individual",setdiff(names(summary_table),"Individual"))]
write.csv(summary_table, file.path(out_dir, "beta_binomial_mixture_summary.csv"), row.names=FALSE)

predicted_support <- function(d,p) {
  depth_tab <- table(d$DP)
  a_main <- p$mu_main*p$phi_main; b_main <- (1-p$mu_main)*p$phi_main
  a_minor <- p$mu_minor*p$phi_minor; b_minor <- (1-p$mu_minor)*p$phi_minor
  pieces <- list()
  z <- 1
  for(nn_chr in names(depth_tab)) {
    nn <- as.integer(nn_chr)
    cnt <- as.numeric(depth_tab[nn_chr])
    kk <- 0:nn; rr <- kk/nn
    main_mass <- exp(dbetabinom_log(kk,nn,a_main,b_main))*cnt*(1-p$w_minor)
    minor_mass <- exp(dbetabinom_log(kk,nn,a_minor,b_minor))*cnt*p$w_minor
    pieces[[z]] <- data.frame(r=rr,main=main_mass,minor=minor_mass)
    z <- z+1
  }
  all <- do.call(rbind,pieces)
  main_agg <- aggregate(main~r,data=all,FUN=sum)
  minor_agg <- aggregate(minor~r,data=all,FUN=sum)
  out <- merge(main_agg,minor_agg,by="r",all=TRUE)
  out[is.na(out)] <- 0
  out
}

smooth_support <- function(r,mass,grid,h=0.04) {
  y <- numeric(length(grid))
  for(j in seq_along(r)) {
    y <- y + mass[j]*(dnorm(grid,mean=r[j],sd=h) + dnorm(grid,mean=-r[j],sd=h) + dnorm(grid,mean=2-r[j],sd=h))
  }
  y
}

breaks <- seq(0,1,length.out=31); binwidth <- diff(breaks)[1]
grid <- seq(0,1,length.out=1000); smooth_bandwidth <- 0.025

curve_list <- list()
for(id in individuals) {
  cat("Generating smoothed curves for:", id, "\n")
  d <- dat[dat$Individual==id,]
  p <- fit_objects[[id]]$parameters
  supp <- predicted_support(d,p)
  sm_main <- smooth_support(supp$r,supp$main,grid,h=smooth_bandwidth)
  sm_minor <- smooth_support(supp$r,supp$minor,grid,h=smooth_bandwidth)
  curve_list[[id]] <- data.frame(Individual=id, AB=grid, Total=(sm_main+sm_minor)*binwidth, Main=sm_main*binwidth, Minor=sm_minor*binwidth)
}
curves <- do.call(rbind,curve_list)

plot_data <- function(subset_data, title, filename) {
  p <- ggplot() +
    geom_histogram(data=dat[dat$Individual %in% subset_data, ], aes(x=Allele_Balance), breaks=breaks, fill="grey75", colour="black") +
    geom_line(data=curves[curves$Individual %in% subset_data, ], aes(x=AB,y=Total), linewidth=0.9) +
    geom_vline(xintercept=0.5, linetype="dashed", alpha=0.4) +
    facet_wrap(~Individual, ncol=2, scales="free_y") +
    labs(x="Allele balance (B / [A + B])", y="Count", title=title) +
    theme_bw(base_size=11)
  
  ggsave(file.path(out_dir, paste0(filename, ".pdf")), p, width=9, height=10)
  ggsave(file.path(out_dir, paste0(filename, ".png")), p, width=9, height=10, dpi=300)
}

plot_data(c("F03 (ref. F01)", "C05 (ref. C04)", "C06 (ref. C04)"), "Pedigree Models (Full Genome)", "beta_binomial_pedigree")
plot_data(c("UNAH2 (ref. UNAH1)", "UNAH3 (ref. UNAH1)", "UNAH1 (ref. UNAH2)", "UNAH3 (ref. UNAH2)", "UNAH1 (ref. UNAH3)", "UNAH2 (ref. UNAH3)"), "UNAH Pairwise Models (Full Genome)", "beta_binomial_unah")

write.csv(curves, file.path(out_dir, "beta_binomial_smoothed_curves.csv"), row.names=FALSE)
cat("\nOutputs successfully written.\n")
EOF

# ==============================================================================
# 5. EXECUTE R SCRIPTS
# ==============================================================================
echo "=== EXECUTING R SCRIPT 1: EXTRACTING ALL UNBIASED LOCI ==="
Rscript extract_unbiased_ALL.R

echo "=== EXECUTING R SCRIPT 2: RUNNING OPTIMIZED MIXTURE MODEL ==="
Rscript colleague_mixture_model_optimized.R "$OUTPUT_DIR"

echo "pipeline completed successfully!"

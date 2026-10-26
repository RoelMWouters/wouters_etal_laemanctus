#!/bin/bash
#PBS -N busco_parental_audit
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=100gb
#PBS -l walltime=02:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT OVERVIEW
# Generates a diagnostic audit of discarded BUSCO loci by checking mapping 
# depth, mapping quality, and clipped read proportions in 100-bp windows in BAM files. 
# Compares F1 hybrid alignments against two parental baselines to differentiate 
# true biological recombination/LOH from assembly artefacts and dropouts.
# ==============================================================================

set -euo pipefail

# Clean scratch directory on exit or failure
cleanup() { rm -rf "${SCRATCHDIR:?}"/* 2>/dev/null || true; }
trap cleanup EXIT TERM

cd "$SCRATCHDIR" || exit 1
echo "Job started on $(hostname) at $(date)"

# ==============================================================================
# 1. ENVIRONMENT SETUP
# ==============================================================================
module purge
module add mambaforge

export CONDA_PKGS_DIRS="${SCRATCHDIR}/conda_pkgs"
export CONDA_ENVS_DIRS="${SCRATCHDIR}/conda_envs"
mkdir -p "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS"

ENV_PREFIX="${SCRATCHDIR}/conda_envs/audit_env"
echo "=== CREATING AUDIT ENVIRONMENT ==="
mamba create -p "${ENV_PREFIX}" -y -c conda-forge -c bioconda python=3.11 pysam pandas

eval "$(conda shell.bash hook)"
conda activate "${ENV_PREFIX}"

# ==============================================================================
# 2. PATHS & DATA STAGING (USER CONFIGURABLE)
# ==============================================================================
# Replace these placeholder paths with your actual project directories
QC_FILE="/path/to/project/output_dir/GENE_QC_CLASSIFICATION.tsv"
BUSCO_TABLE="/path/to/project/busco_run/full_table.tsv"
BAM_DIR="/path/to/project/bam_alignments"
BASE_OUT="/path/to/project/audit_results"

mkdir -p "$BASE_OUT"
cp "$QC_FILE" .
cp "$BUSCO_TABLE" .

# ==============================================================================
# 3. PYTHON DIAGNOSTIC (PARENTAL CONTROLS + F1)
# ==============================================================================
cat << 'EOF' > audit_with_parents.py
import pysam
import pandas as pd
import sys
import os

# --- COMMAND LINE ARGUMENTS ---
qc_file = sys.argv[1]
busco_table = sys.argv[2]
bam_dir = sys.argv[3]

# --- USER CONFIGURABLE COHORT VARIABLES ---
# Define the prefixes for your BAM files (excluding the _final_markdup.bam suffix)
f1_samples = ["F1_Sample_1", "F1_Sample_2", "F1_Sample_3"] 
parent1_id = "Parent1_Sample_ID"
parent2_id = "Parent2_Sample_ID"

avg_genome_depth = 25.0
window_size = 100

print("=== DIAGNOSTIC CHECKPOINTS ===")
qc_df = pd.read_csv(qc_file, sep="\t")
qc_df['BUSCO_ID'] = qc_df['BUSCO_ID'].astype(str).str.strip()

# Isolate loci that failed the initial ortholog QC pipeline
discarded = qc_df[qc_df['QC_Status'] != "0_PASSED_CLEAN_ORTHOLOG"]

busco_df = pd.read_csv(busco_table, sep="\t", comment="#", header=None, 
                       names=["Busco_id", "Status", "Contig", "Start", "End", "Strand", "Score", "Length", "URL", "Desc"])
busco_df['Busco_id'] = busco_df['Busco_id'].astype(str).str.strip()

# Target physically duplicated BUSCO regions for depth/artifact auditing
busco_groups = busco_df[busco_df['Status'] == 'Duplicated'].groupby('Busco_id')

overlap = set(discarded['BUSCO_ID']).intersection(set(busco_groups.groups.keys()))
print(f"Discarded loci matched to Duplicated BUSCO coordinates: {len(overlap)}")

# ---------------------------------------------------------
# PARENTAL EVALUATION
# ---------------------------------------------------------
def get_parental_consensus_and_depth(bam, contig, start, end):
    """
    Evaluates parental bam files across the target window.
    Extracts consensus SNPs and minimum/maximum read depth to detect
    reference dropouts or chimeric assemblies in the baseline.
    """
    consensus = {}
    length = end - start
    if length <= 0: return consensus, 0, 0
    
    depths = [0] * length
    
    # Pileup for consensus SNPs (Min MAPQ 20 to avoid repeats)
    for col in bam.pileup(contig, start, end, truncate=True, min_mapping_quality=20):
        idx = col.pos - start
        if 0 <= idx < length:
            depths[idx] = col.nsegments
            counts = {}
            for read in col.pileups:
                if not read.is_del and not read.is_refskip and read.query_position is not None:
                    b = read.alignment.query_sequence[read.query_position]
                    counts[b] = counts.get(b, 0) + 1
            if counts:
                top_base = max(counts, key=counts.get)
                consensus[idx] = top_base
                
    # Sliding window for min/max depth
    min_w_d = 9999.0
    max_w_d = 0.0
    for w in range(0, length, window_size):
        w_end = min(w + window_size, length)
        w_avg = sum(depths[w:w_end]) / (w_end - w)
        if w_avg < min_w_d: min_w_d = w_avg
        if w_avg > max_w_d: max_w_d = w_avg
        
    if min_w_d == 9999.0: min_w_d = 0.0
    return consensus, min_w_d, max_w_d

# ---------------------------------------------------------
# F1 HYBRID SLIDING WINDOW
# ---------------------------------------------------------
def get_f1_worst_window(bam, contig, start, end):
    """
    Scans the F1 BAM file in sliding windows to find structural
    alignment garbage (Depth spikes, MAPQ 0 pileups, discordant pairs, 
    and soft-clipping) indicative of collapsed paralogs or repeat misalignment.
    """
    length = end - start
    if length <= 0: return 0.0, 0.0, 0.0, 0.0
    
    depth = [0] * length
    mapq0 = [0] * length
    discord = [0] * length
    clipped = [0] * length
    aligned = [0] * length
    
    for read in bam.fetch(contig, start, end):
        r_start = max(0, read.reference_start - start)
        r_end = min(length, read.reference_end - start) if read.reference_end else r_start + 1
        
        is_mapq0 = 1 if read.mapping_quality < 10 else 0
        is_disc = 1 if not read.is_proper_pair else 0
        
        for i in range(r_start, r_end):
            depth[i] += 1
            mapq0[i] += is_mapq0
            discord[i] += is_disc
            
        if read.cigarstring:
            clip_len = sum(l for op, l in read.cigartuples if op == 4)
            aln_len = sum(l for op, l in read.cigartuples if op in (0, 1, 2, 7, 8))
            if clip_len > 0:
                for i in range(r_start, r_end):
                    clipped[i] += clip_len
                    aligned[i] += aln_len

    max_d = 0; max_mq = 0.0; max_dc = 0.0; max_cl = 0.0
    
    for w in range(0, length, window_size):
        w_end = min(w + window_size, length)
        w_len = w_end - w
        
        sum_d = sum(depth[w:w_end])
        w_depth = sum_d / w_len
        if w_depth > max_d: max_d = w_depth
        
        if sum_d > 0:
            mq_pct = (sum(mapq0[w:w_end]) / sum_d) * 100
            if mq_pct > max_mq: max_mq = mq_pct
            
            dc_pct = (sum(discord[w:w_end]) / sum_d) * 100
            if dc_pct > max_dc: max_dc = dc_pct
            
        sum_aln = sum(aligned[w:w_end])
        sum_clip = sum(clipped[w:w_end])
        if (sum_aln + sum_clip) > 0:
            cl_pct = (sum_clip / (sum_aln + sum_clip)) * 100
            if cl_pct > max_cl: max_cl = cl_pct

    return max_d, max_mq, max_dc, max_cl

# ---------------------------------------------------------
# EXECUTE AUDIT
# ---------------------------------------------------------
p1_bam = pysam.AlignmentFile(os.path.join(bam_dir, f"{parent1_id}_final_markdup.bam"), "rb")
p2_bam = pysam.AlignmentFile(os.path.join(bam_dir, f"{parent2_id}_final_markdup.bam"), "rb")

f1_bams = {}
for s in f1_samples:
    p = os.path.join(bam_dir, f"{s}_final_markdup.bam")
    if os.path.exists(p):
        f1_bams[s] = pysam.AlignmentFile(p, "rb")

results = []

for _, row in discarded.iterrows():
    busco_id = row['BUSCO_ID']
    qc_status = row['QC_Status']
    if busco_id not in busco_groups.groups: continue
    
    loci = busco_groups.get_group(busco_id)
    
    # --- 1. PARENTAL BASELINE ---
    parental_snps = 0
    p1_min_d_overall = 9999.0
    p2_min_d_overall = 9999.0
    total_length = 0
    
    for _, locus in loci.iterrows():
        s, e = int(locus['Start']), int(locus['End'])
        if s > e: s, e = e, s
        total_length += (e - s)
        
        p1_cons, p1_min, _ = get_parental_consensus_and_depth(p1_bam, locus['Contig'], s, e)
        p2_cons, p2_min, _ = get_parental_consensus_and_depth(p2_bam, locus['Contig'], s, e)
        
        if p1_min < p1_min_d_overall: p1_min_d_overall = p1_min
        if p2_min < p2_min_d_overall: p2_min_d_overall = p2_min
        
        shared_pos = set(p1_cons.keys()).intersection(set(p2_cons.keys()))
        parental_snps += sum(1 for p in shared_pos if p1_cons[p] != p2_cons[p])
    
    # Establish Parental Status
    if p1_min_d_overall < (avg_genome_depth * 0.2) or p2_min_d_overall < (avg_genome_depth * 0.2):
        parental_status = "Reference Chimera (Parental Dropout)"
    elif parental_snps < 3:
        parental_status = "Ultra-Conserved (Parents Identical)"
    else:
        parental_status = "Valid Biological Anchor"

    # --- 2. F1 HYBRID AUDIT ---
    for f1_name, f1_bam in f1_bams.items():
        gene_max_d = 0; gene_max_mq0 = 0.0; gene_max_disc = 0.0; gene_max_clip = 0.0
        
        for _, locus in loci.iterrows():
            s, e = int(locus['Start']), int(locus['End'])
            if s > e: s, e = e, s
            
            d, mq, dc, cl = get_f1_worst_window(f1_bam, locus['Contig'], s, e)
            if d > gene_max_d: gene_max_d = d
            if mq > gene_max_mq0: gene_max_mq0 = mq
            if dc > gene_max_disc: gene_max_disc = dc
            if cl > gene_max_clip: gene_max_clip = cl
            
        # --- 3. FINAL LOGIC GATE ---
        verdict = "Unknown"
        
        # Check F1 BAM Garbage Metrics First
        if gene_max_d > (avg_genome_depth * 2.0):
            verdict = "Artifact: Collapsed Paralog (F1 Depth Spike)"
        elif gene_max_mq0 > 30.0 or gene_max_disc > 25.0 or gene_max_clip > 15.0:
            verdict = "Artifact: Assembly/Repeat Misalignment (F1 Garbage Metrics)"
        elif gene_max_d < (avg_genome_depth * 0.3):
            verdict = "Artifact: F1 Allelic Dropout"
        else:
            # F1 is structurally sound. Check Parental Reality.
            if parental_status == "Ultra-Conserved (Parents Identical)":
                verdict = "Artifact: Ultra-Conserved False Positive"
            elif parental_status == "Reference Chimera (Parental Dropout)":
                verdict = "Artifact: Assembly Chimera (Reference Error)"
            else:
                # Both F1 and Parents are perfectly sound and divergent.
                if "Het<0.5%" in qc_status:
                    verdict = "BIOLOGICAL: Confirmed LOH"
                else:
                    verdict = "BIOLOGICAL: Confirmed Recombinant"

        results.append({
            "Sample": f1_name,
            "BUSCO_ID": busco_id,
            "Gene_Length_bp": total_length,
            "Parental_Fixed_SNPs": parental_snps,
            "Parent1_MinDepth": round(p1_min_d_overall, 1),
            "Parent2_MinDepth": round(p2_min_d_overall, 1),
            "Parental_Status": parental_status,
            "F1_Worst_Depth": round(gene_max_d, 1),
            "F1_Worst_Clipped_Pct": round(gene_max_clip, 1),
            "Final_Verdict": verdict
        })

p1_bam.close()
p2_bam.close()
for b in f1_bams.values(): b.close()

columns = [
    "Sample", "BUSCO_ID", "Gene_Length_bp", 
    "Parental_Fixed_SNPs", "Parent1_MinDepth", "Parent2_MinDepth", "Parental_Status", 
    "F1_Worst_Depth", "F1_Worst_Clipped_Pct", "Final_Verdict"
]
output_df = pd.DataFrame(results, columns=columns)
output_df.to_csv("Master_Cohort_Audit_Report.tsv", sep="\t", index=False)

print("\n=== FINAL VERDICT SUMMARY ===")
if not output_df.empty:
    print(output_df['Final_Verdict'].value_counts())
    print("\n=== CONFIRMED BIOLOGICAL EVENTS ===")
    recombinants = output_df[output_df['Final_Verdict'].str.contains("BIOLOGICAL")]
    if not recombinants.empty:
        print(recombinants[['Sample', 'BUSCO_ID', 'Parental_Fixed_SNPs', 'Final_Verdict']])
    else:
        print("Zero true biological exceptions detected. 100% of discarded loci confirmed as artifacts.")
else:
    print("No loci were processed.")
EOF

# ==============================================================================
# 4. EXECUTE & EXPORT
# ==============================================================================
echo "=== EXECUTING MULTI-SAMPLE AUDIT WITH PARENTAL CONTROLS ==="
python3 audit_with_parents.py "GENE_QC_CLASSIFICATION.tsv" "full_table.tsv" "$BAM_DIR" > audit_log.txt

cp Master_Cohort_Audit_Report.tsv "$BASE_OUT/"
cp audit_log.txt "$BASE_OUT/"

echo "=== PIPELINE COMPLETE. Results in $BASE_OUT ==="

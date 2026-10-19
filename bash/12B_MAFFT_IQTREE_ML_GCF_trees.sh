#!/bin/bash
#PBS -N Part2_Phylogenomics_Merged
#PBS -l select=1:ncpus=32:mem=64gb:scratch_local=100gb
#PBS -l walltime=12:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# PHYLOGENOMIC PIPELINE - PART 2: ALIGN, PHASE, QC & TREES
# Takes raw extracted CDS loci from Part 1, compiles them, aligns 
# with MAFFT, trims with trimAl, and uses a distance-based engine to phase alleles.
# Phased alignments then undergo rigorous SNP QC (LOH & phase-switch detection)
# before calculating unphased NJ distances, single-gene trees, supermatrices, 
# and gene/site concordance factors (gCF/sCF).
# ==============================================================================

set -euo pipefail
cleanup() {
    echo "=== Cleaning up SCRATCHDIR ==="
    rm -rf "${SCRATCHDIR:?}"/* 2>/dev/null || true
}
trap cleanup EXIT TERM

cd "$SCRATCHDIR" || { echo "No SCRATCHDIR – abort"; exit 1; }
echo "Job started on $(hostname) at $(date)"

# ==============================================================================
# 1. GLOBAL PATHS & WORKSPACE SETUP
# ==============================================================================
PROJECT_DIR="/path/to/project"
PART1_DIR="${PROJECT_DIR}/output/part1"
INPUT_META="${PART1_DIR}/metadata"
INPUT_SEQS="${PART1_DIR}/1_extracted_sequences"

BASE_OUT="${PROJECT_DIR}/output/part2_phylogenomics"

# Internal Scratch Architecture
mkdir -p unaligned aligned trimmed phased 
mkdir -p clean_alns_19tip clean_4tax_mat clean_4tax_pat iqtree_out report
mkdir -p "$BASE_OUT"

# ==============================================================================
# 2. ISOLATED CONDA ENVIRONMENT
# ==============================================================================
module purge
module add mambaforge

export CONDA_PKGS_DIRS="${SCRATCHDIR}/conda_pkgs"
export CONDA_ENVS_DIRS="${SCRATCHDIR}/conda_envs"
ENV_PREFIX="${SCRATCHDIR}/conda_envs/phylo_env"

echo "=== Creating Conda Environment ==="
mamba create -p "${ENV_PREFIX}" -y -c conda-forge -c bioconda \
    python=3.10 mafft trimal iqtree parallel biopython r-base r-ape

eval "$(conda shell.bash hook)"
conda activate "${ENV_PREFIX}"

# ==============================================================================
# 3. COMPILATION, ALIGNMENT & TRIMMING
# ==============================================================================
echo "=== COMPILING UNALIGNED FASTA MATRICES ==="
tr -d '\r' < "${INPUT_META}/target_busco_ids.txt" > clean_busco_ids.txt
rm -f alignment_commands.txt

while read -r BUSCO; do
    [ -z "$BUSCO" ] && continue 
    OUT_FASTA="unaligned/${BUSCO}.fasta"
    
    # Merge references and exonerated CDS
    cat "${INPUT_SEQS}/ref_and_ont_seqs/${BUSCO}_"*.fasta \
        "${INPUT_SEQS}/exonerated_cds/"*"_${BUSCO}_"*.fasta > "$OUT_FASTA" 2>/dev/null

    # Strictly omit any sequence that Exonerate fragmented into multiple pieces
    awk '/^>/ {id=$0; if (!seen[id]++) order[++n]=id; count[id]++; next} {seqs[id]=seqs[id]""$0} END {for(i=1;i<=n;i++) if(count[order[i]]==1) print order[i]"\n"seqs[order[i]]}' "$OUT_FASTA" > "${OUT_FASTA}.tmp" && mv "${OUT_FASTA}.tmp" "$OUT_FASTA"

    echo "mafft --quiet --auto unaligned/${BUSCO}.fasta > aligned/${BUSCO}.fasta 2>/dev/null && trimal -in aligned/${BUSCO}.fasta -out trimmed/${BUSCO}.fasta -automated1 > /dev/null 2>&1" >> alignment_commands.txt
done < clean_busco_ids.txt

find unaligned -type f -empty -delete

if [ -s alignment_commands.txt ]; then
    echo "=== RUNNING ALIGNMENTS ON 32 CORES ==="
    parallel -j 32 < alignment_commands.txt
fi

# ==============================================================================
# 4. PYTHON ALLELE PHASING
# ==============================================================================
echo "=== RUNNING HIERARCHICAL ALLELE PHASING ==="
cat << 'EOF' > phase_alignments.py
import os
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

input_dir = 'trimmed'
output_dir = 'phased'
hybrid_samples = ['C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'F01', 'F03']

for filename in os.listdir(input_dir):
    if not filename.endswith(('.fasta', '.fa', '.phy')): continue
    filepath = os.path.join(input_dir, filename)
    try: aln = AlignIO.read(filepath, 'fasta')
    except: continue
        
    length = aln.get_alignment_length()
    if length < 150: continue
        
    # Grab References
    ref_records = {}
    sample_records = {s: [] for s in hybrid_samples + ['C10', 'C09']}
    
    for rec in aln:
        if rec.id.startswith('Parent1_Ref'):
            rec.id = 'Parent1_Ref'; rec.name = ''; rec.description = ''
            ref_records['Parent1_Ref'] = rec
        elif rec.id.startswith('Outgroup'):
            rec.id = 'Outgroup'; rec.name = ''; rec.description = ''
            ref_records['Outgroup'] = rec
        else:
            for s in sample_records.keys():
                if rec.id.startswith(s) or (s in rec.id):
                    sample_records[s].append(rec)
                    break
                    
    if 'Parent1_Ref' not in ref_records: continue
    pat_anchor = str(ref_records['Parent1_Ref'].seq).upper()
    
    def get_dist(seq_a, seq_b):
        valid = sum(1 for a, b in zip(seq_a, seq_b) if a in 'ACGT' and b in 'ACGT')
        m = sum(1 for a, b in zip(seq_a, seq_b) if a in 'ACGT' and b in 'ACGT' and a != b)
        return m / valid if valid > 0 else 999.0

    new_records = list(ref_records.values())

    # Pure Parent Reference Selection
    if sample_records['C10']:
        best_c10 = min(sample_records['C10'], key=lambda r: get_dist(str(r.seq).upper(), pat_anchor))
        best_c10.id = 'C10'; best_c10.name = ''; best_c10.description = ''
        new_records.append(best_c10)
        pat_anchor = str(best_c10.seq).upper()

    mat_anchor = None
    if sample_records['C09']:
        best_c09 = max(sample_records['C09'], key=lambda r: get_dist(str(r.seq).upper(), pat_anchor))
        best_c09.id = 'C09'; best_c09.name = ''; best_c09.description = ''
        new_records.append(best_c09)
        mat_anchor = str(best_c09.seq).upper()

    # Phase Hybrids
    for s in hybrid_samples:
        recs = sample_records[s]
        if not recs: continue
        
        if mat_anchor:
            scored = [{'rec': r, 'd_pat': get_dist(str(r.seq).upper(), pat_anchor), 'd_mat': get_dist(str(r.seq).upper(), mat_anchor)} for r in recs]
                
            if len(scored) == 1:
                best = scored[0]
                best['rec'].id = f"{s}_Pat" if best['d_pat'] < best['d_mat'] else f"{s}_Mat"
                best['rec'].name = ''; best['rec'].description = ''
                new_records.append(best['rec'])
            else:
                scored.sort(key=lambda x: x['d_pat'])
                best_pat = scored[0]
                scored.sort(key=lambda x: x['d_mat'])
                best_mat = scored[0]
                if best_pat['rec'] is best_mat['rec']: best_mat = scored[1]
                
                best_pat['rec'].id = f"{s}_Pat"; best_mat['rec'].id = f"{s}_Mat"
                best_pat['rec'].name = best_mat['rec'].name = ''
                best_pat['rec'].description = best_mat['rec'].description = ''
                new_records.extend([best_pat['rec'], best_mat['rec']])

    phased_aln = MultipleSeqAlignment(new_records)
    if len(phased_aln) >= 4: 
        AlignIO.write(phased_aln, os.path.join(output_dir, filename), 'fasta')
EOF
python3 phase_alignments.py

# ==============================================================================
# 5. PYTHON QC, DISTANCES & SUPERMATRIX BUILDER
# ==============================================================================
echo "=== RUNNING QC, UNPHASED DISTANCES & MATRIX BUILDER ==="
cat << 'EOF' > step2_qc_and_unphased.py
import os, re
from collections import Counter, defaultdict

input_dir     = "phased"
clean_19_dir  = "clean_alns_19tip"
clean_mat_dir = "clean_4tax_mat"
clean_pat_dir = "clean_4tax_pat"

unah_ids = ["C01", "C02", "C03"]
suya_ids = ["C04", "C05", "C06", "F01", "F03"]
all_hybrids = unah_ids + suya_ids

unah_mat = [f"{x}_Mat" for x in unah_ids]
suya_mat = [f"{x}_Mat" for x in suya_ids]
unah_pat = [f"{x}_Pat" for x in unah_ids]
suya_pat = [f"{x}_Pat" for x in suya_ids]

diploid_taxa = ["Outgroup", "Parent1_Ref", "C10", "C09"] + all_hybrids
all_19_taxa  = ["Outgroup", "Parent1_Ref", "C10", "C09"] + unah_mat + suya_mat + unah_pat + suya_pat
mat_4_taxa   = ["Parent1_Ref", "C09", "UNAH_Mat", "SuyapaZoo_Mat"]
pat_4_taxa   = ["C09", "C10", "UNAH_Pat", "SuyapaZoo_Pat"]

def read_fasta(fp):
    seqs = {}
    cur = None
    with open(fp) as FH:
        for line in FH:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}

def majority_consensus(seqs, taxa_list, min_votes):
    present = [seqs[t] for t in taxa_list if t in seqs]
    if len(present) < min_votes: return None
    L = len(present[0])
    out = []
    for i in range(L):
        bases = [s[i] for s in present if s[i] in "ACGT"]
        if len(bases) < min_votes: out.append("N")
        else:
            top, cnt = Counter(bases).most_common(1)[0]
            out.append(top if cnt > len(bases) / 2.0 else "N")
    return "".join(out)

qc_counts  = Counter()
qc_details = []
stats = {
    "RAW_MAT_SNPS": Counter(), "RAW_MAT_GENES": Counter(),
    "CLEAN_MAT_SNPS": Counter(), "CLEAN_MAT_GENES": Counter(),
    "RAW_PAT_SNPS": Counter(), "RAW_PAT_GENES": Counter(),
    "CLEAN_PAT_SNPS": Counter(), "CLEAN_PAT_GENES": Counter(),
    "CLEAN_MAT_UNINF_REASONS": Counter(),
    "CLEAN_PAT_UNINF_REASONS": Counter()
}
pair_diff_all, pair_tot_all     = defaultdict(float), defaultdict(int)
pair_diff_clean, pair_tot_clean = defaultdict(float), defaultdict(int)

super_19  = {t: [] for t in all_19_taxa}
super_mat = {t: [] for t in mat_4_taxa}
super_pat = {t: [] for t in pat_4_taxa}
part_19, part_mat, part_pat = [], [], []
pos_19 = pos_mat = pos_pat = 1

for fname in sorted(os.listdir(input_dir)):
    if not fname.endswith(('.fasta', '.fa', '.phy')): continue
    busco_id = re.sub(r"\.fasta$", "", fname)
    fpath = os.path.join(input_dir, fname)
    seqs = read_fasta(fpath)

    u_mat = majority_consensus(seqs, unah_mat, 2)
    s_mat = majority_consensus(seqs, suya_mat, 3)
    u_pat = majority_consensus(seqs, unah_pat, 2)
    s_pat = majority_consensus(seqs, suya_pat, 3)

    core_ok = all(k in seqs for k in ("C09", "C10", "Parent1_Ref", "Outgroup"))
    if not (u_mat and s_mat and u_pat and s_pat and core_ok):
        qc_counts["1_Excluded_Missing_Core_Taxa_or_Consensus"] += 1
        qc_details.append((busco_id, "Excluded_Missing_Core", 0.0, 0.0, 0, 0))
        continue

    c09, c10, p1, out = seqs["C09"], seqs["C10"], seqs["Parent1_Ref"], seqs["Outgroup"]
    L = len(c09)
    u_het = s_het = 0.0
    valid = 0
    for i in range(L):
        if all(x[i] in "ACGT" for x in (u_mat, u_pat, s_mat, s_pat)):
            valid += 1
            if u_mat[i] != u_pat[i]: u_het += 1
            if s_mat[i] != s_pat[i]: s_het += 1
    if valid < 100:
        qc_counts["1_Excluded_Missing_Core_Taxa_or_Consensus"] += 1
        continue
    u_het /= valid
    s_het /= valid

    mat_g = Counter(); pat_g = Counter()
    mat_sw = pat_sw = 0
    mat_cln = mat_c09 = pat_cln = pat_c10 = 0

    for i in range(L):
        bu_m, bs_m = u_mat[i], s_mat[i]
        bu_p, bs_p = u_pat[i], s_pat[i]
        bc09, bc10, bp1, bout = c09[i], c10[i], p1[i], out[i]

        if bu_m != bs_m and all(b in "ACGT" for b in (bu_m, bs_m)): mat_cln += 1
        if bu_p != bs_p and all(b in "ACGT" for b in (bu_p, bs_p)): pat_cln += 1

        if bout == bp1 and all(b in "ACGT" for b in (bu_m, bs_m, bc09, bp1, bout)):
            if bc09 == bu_m and bc09 != bout and bs_m == bout:
                mat_g["Hyp1_DualOrigin_(Parent+UNAH)"] += 1
                if bs_p == bc09: mat_sw += 1
            elif bu_m == bs_m and bu_m != bout and bc09 == bout:
                mat_g["Hyp2_SingleOrigin_(UNAH+SuyapaZoo)"] += 1
            elif bc09 == bs_m and bc09 != bout and bu_m == bout:
                mat_g["Hyp3_AltDual_(Parent+SuyapaZoo)"] += 1
                if bu_p == bc09: mat_sw += 1
            elif bc09 != bout and bu_m == bout and bs_m == bout:
                mat_c09 += 1

        if bout == bc09 and all(b in "ACGT" for b in (bu_p, bs_p, bc10, bc09, bout)):
            if bc10 == bu_p and bc10 != bout and bs_p == bout:
                pat_g["Hyp1_DualOrigin_(Parent+UNAH)"] += 1
                if bs_m == bc10: pat_sw += 1
            elif bu_p == bs_p and bu_p != bout and bc10 == bout:
                pat_g["Hyp2_SingleOrigin_(UNAH+SuyapaZoo)"] += 1
            elif bc10 == bs_p and bc10 != bout and bu_p == bout:
                pat_g["Hyp3_AltDual_(Parent+SuyapaZoo)"] += 1
                if bu_m == bc10: pat_sw += 1
            elif bc10 != bout and bu_p == bout and bs_p == bout:
                pat_c10 += 1

    def winner(cnt):
        if not cnt: return "0_Uninformative"
        top = cnt.most_common(2)
        if len(top) > 1 and top[0][1] == top[1][1]: return "4_Conflicted_Tie"
        return top[0][0]

    for k, v in mat_g.items(): stats["RAW_MAT_SNPS"][k] += v
    stats["RAW_MAT_GENES"][winner(mat_g)] += 1
    for k, v in pat_g.items(): stats["RAW_PAT_SNPS"][k] += v
    stats["RAW_PAT_GENES"][winner(pat_g)] += 1

    het_ratio = max(u_het, s_het) / max(min(u_het, s_het), 1e-4)

    if u_het > 0.10 or s_het > 0.10:
        status = "2_Excluded_Paralog_or_Frameshift_(Het>10%)"
    elif u_het < 0.005 or s_het < 0.005 or het_ratio > 3.0:
        status = "3_Excluded_LOH_or_AllelicDropout_(Het<0.5%_or_Asymmetric)"
    elif mat_sw > 0 or pat_sw > 0:
        status = "4_Excluded_Subgenome_PhaseSwitch_Chimera"
    else:
        status = "0_PASSED_CLEAN_ORTHOLOG"

    qc_counts[status] += 1
    qc_details.append((busco_id, status, u_het, s_het, mat_sw, pat_sw))

    # Unphased distances
    gene_alleles = {}
    for t in ["Outgroup", "Parent1_Ref", "C10", "C09"]:
        if t in seqs:
            gene_alleles[t] = [(b, b) if b in "ACGT" else ("N", "N") for b in seqs[t]]
    for h in all_hybrids:
        mk, pk = f"{h}_Mat", f"{h}_Pat"
        if mk in seqs and pk in seqs:
            gene_alleles[h] = [
                tuple(sorted((seqs[mk][i], seqs[pk][i]))) if (seqs[mk][i] in "ACGT" and seqs[pk][i] in "ACGT") else ("N", "N")
                for i in range(L)
            ]

    if len(gene_alleles) == len(diploid_taxa):
        for i_idx, t1 in enumerate(diploid_taxa):
            g1 = gene_alleles[t1]
            for t2 in diploid_taxa[i_idx+1:]:
                g2 = gene_alleles[t2]
                diffs = valid_d = 0
                for a, b in zip(g1, g2):
                    if a[0] != "N" and b[0] != "N":
                        valid_d += 1
                        if a == b: continue
                        elif a[0] in b or a[1] in b: diffs += 0.5
                        else: diffs += 1.0
                pair_diff_all[(t1, t2)] += diffs
                pair_tot_all[(t1, t2)]  += valid_d
                if status == "0_PASSED_CLEAN_ORTHOLOG":
                    pair_diff_clean[(t1, t2)] += diffs
                    pair_tot_clean[(t1, t2)]  += valid_d

    if status == "0_PASSED_CLEAN_ORTHOLOG":
        with open(os.path.join(clean_19_dir, fname), "w") as FH19:
            for t in all_19_taxa:
                if t in seqs: FH19.write(f">{t}\n{seqs[t]}\n")

        with open(os.path.join(clean_mat_dir, fname), "w") as FH:
            FH.write(f">UNAH_Mat\n{u_mat}\n>SuyapaZoo_Mat\n{s_mat}\n>C09\n{c09}\n>Parent1_Ref\n{p1}\n")
        with open(os.path.join(clean_pat_dir, fname), "w") as FH:
            FH.write(f">UNAH_Pat\n{u_pat}\n>SuyapaZoo_Pat\n{s_pat}\n>C10\n{c10}\n>C09\n{c09}\n")

        for t in all_19_taxa: super_19[t].append(seqs.get(t, "-" * L))
        part_19.append(f"DNA, {busco_id} = {pos_19}-{pos_19 + L - 1};")
        pos_19 += L

        mat_map = {"UNAH_Mat": u_mat, "SuyapaZoo_Mat": s_mat, "C09": c09, "Parent1_Ref": p1}
        for t in mat_4_taxa: super_mat[t].append(mat_map[t])
        part_mat.append(f"DNA, {busco_id} = {pos_mat}-{pos_mat + L - 1};")
        pos_mat += L

        pat_map = {"UNAH_Pat": u_pat, "SuyapaZoo_Pat": s_pat, "C10": c10, "C09": c09}
        for t in pat_4_taxa: super_pat[t].append(pat_map[t])
        part_pat.append(f"DNA, {busco_id} = {pos_pat}-{pos_pat + L - 1};")
        pos_pat += L

        for k, v in mat_g.items(): stats["CLEAN_MAT_SNPS"][k] += v
        stats["CLEAN_MAT_GENES"][winner(mat_g)] += 1
        if winner(mat_g) == "0_Uninformative":
            if mat_cln == 0 and mat_c09 == 0: stats["CLEAN_MAT_UNINF_REASONS"]["1_100%_Identical_C09_UNAH_Suyapa"] += 1
            elif mat_cln == 0 and mat_c09 > 0: stats["CLEAN_MAT_UNINF_REASONS"]["2_Clones_100%_Identical_Only_C09_Singletons"] += 1
            else: stats["CLEAN_MAT_UNINF_REASONS"]["3_Only_Clone_Autapomorphy_or_Unpolarized"] += 1

        for k, v in pat_g.items(): stats["CLEAN_PAT_SNPS"][k] += v
        stats["CLEAN_PAT_GENES"][winner(pat_g)] += 1
        if winner(pat_g) == "0_Uninformative":
            if pat_cln == 0 and pat_c10 == 0: stats["CLEAN_PAT_UNINF_REASONS"]["1_100%_Identical_C10_UNAH_Suyapa"] += 1
            elif pat_cln == 0 and pat_c10 > 0: stats["CLEAN_PAT_UNINF_REASONS"]["2_Clones_100%_Identical_Only_C10_Singletons"] += 1
            else: stats["CLEAN_PAT_UNINF_REASONS"]["3_Only_Clone_Autapomorphy_or_Unpolarized"] += 1

def write_supermatrix(s_dict, p_list, fasta_out, nex_out):
    with open(fasta_out, "w") as F:
        for t, chunks in s_dict.items(): F.write(f">{t}\n{''.join(chunks)}\n")
    with open(nex_out, "w") as N:
        N.write("#nexus\nbegin sets;\n")
        for line in p_list: N.write(f"  charset {line}\n")
        N.write("end;\n")

write_supermatrix(super_19, part_19, "supermatrix_19tip.fasta", "supermatrix_19tip.nex")
write_supermatrix(super_mat, part_mat, "supermatrix_4tax_mat.fasta", "supermatrix_4tax_mat.nex")
write_supermatrix(super_pat, part_pat, "supermatrix_4tax_pat.fasta", "supermatrix_4tax_pat.nex")

with open("GENE_QC_CLASSIFICATION.tsv", "w") as O:
    O.write("BUSCO_ID\tQC_Status\tUNAH_Het\tSuyapa_Het\tMat_PhaseSwitches\tPat_PhaseSwitches\n")
    for r in qc_details: O.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.4f}\t{r[4]}\t{r[5]}\n")

with open("FILTER_SUMMARY.tsv", "w") as O:
    O.write("FILTER_CATEGORY\tLoci_Count\n")
    for k in sorted(qc_counts.keys()): O.write(f"{k}\t{qc_counts[k]}\n")

def write_matrix(diff, tot, tsv, nex):
    taxa = diploid_taxa
    def d(a, b):
        if (a, b) in tot and tot[(a, b)] > 0: return diff[(a, b)] / tot[(a, b)]
        if (b, a) in tot and tot[(b, a)] > 0: return diff[(b, a)] / tot[(b, a)]
        return 0.0
    with open(tsv, "w") as F:
        F.write("\t" + "\t".join(taxa) + "\n")
        for a in taxa: F.write(a + "\t" + "\t".join(f"{d(a, b):.6f}" for b in taxa) + "\n")
    with open(nex, "w") as N:
        N.write("#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=%d;\n  TAXLABELS %s;\nEND;\n\n" % (len(taxa), " ".join(taxa)))
        N.write("BEGIN DISTANCES;\n  DIMENSIONS NTAX=%d;\n  FORMAT TRIANGLE=BOTH DIAGONAL LABELS;\n  MATRIX\n" % len(taxa))
        for a in taxa: N.write("    %-15s %s\n" % (a, " ".join(f"{d(a, b):.6f}" for b in taxa)))
        N.write("  ;\nEND;\n")

write_matrix(pair_diff_all, pair_tot_all, "unphased_distances_ALL_LOCI.tsv", "unphased_network_ALL_LOCI.nex")
write_matrix(pair_diff_clean, pair_tot_clean, "unphased_distances_CLEAN_LOCI.tsv", "unphased_network_CLEAN_LOCI.nex")

with open("STEP2_QC_SUMMARY.txt", "w") as sf:
    def log(msg): print(msg); sf.write(msg + "\n")
    log("==========================================================================")
    log("1. GENE-BY-GENE QUALITY CONTROL (QC) FILTERING ACCOUNTING")
    log("==========================================================================")
    for k, v in sorted(qc_counts.items()): log(f"  {k:60s}: {v:5d} loci")
    log("\n==========================================================================")
    log("2. MATERNAL SUBGENOME: BEFORE vs. AFTER QC")
    log("==========================================================================")
    for k in ["Hyp1_DualOrigin_(Parent+UNAH)", "Hyp2_SingleOrigin_(UNAH+SuyapaZoo)", "Hyp3_AltDual_(Parent+SuyapaZoo)"]:
        log(f"  RAW   {k:36s}: {stats['RAW_MAT_SNPS'][k]:5d} SNPs | {stats['RAW_MAT_GENES'][k]:4d} Genes")
    for k in ["Hyp1_DualOrigin_(Parent+UNAH)", "Hyp2_SingleOrigin_(UNAH+SuyapaZoo)", "Hyp3_AltDual_(Parent+SuyapaZoo)"]:
        log(f"  CLEAN {k:36s}: {stats['CLEAN_MAT_SNPS'][k]:5d} SNPs | {stats['CLEAN_MAT_GENES'][k]:4d} Genes")
    log(f"  CLEAN {'0_Uninformative':36s}:       ---  | {stats['CLEAN_MAT_GENES']['0_Uninformative']:4d} Genes")
    for k, v in sorted(stats["CLEAN_MAT_UNINF_REASONS"].items()): log(f"    -> {k:50s}: {v:4d} Genes")

    log("\n==========================================================================")
    log("3. PATERNAL SUBGENOME: BEFORE vs. AFTER QC")
    log("==========================================================================")
    for k in ["Hyp1_DualOrigin_(Parent+UNAH)", "Hyp2_SingleOrigin_(UNAH+SuyapaZoo)", "Hyp3_AltDual_(Parent+SuyapaZoo)"]:
        log(f"  RAW   {k:36s}: {stats['RAW_PAT_SNPS'][k]:5d} SNPs | {stats['RAW_PAT_GENES'][k]:4d} Genes")
    for k in ["Hyp1_DualOrigin_(Parent+UNAH)", "Hyp2_SingleOrigin_(UNAH+SuyapaZoo)", "Hyp3_AltDual_(Parent+SuyapaZoo)"]:
        log(f"  CLEAN {k:36s}: {stats['CLEAN_PAT_SNPS'][k]:5d} SNPs | {stats['CLEAN_PAT_GENES'][k]:4d} Genes")
    log(f"  CLEAN {'0_Uninformative':36s}:       ---  | {stats['CLEAN_PAT_GENES']['0_Uninformative']:4d} Genes")
    for k, v in sorted(stats["CLEAN_PAT_UNINF_REASONS"].items()): log(f"    -> {k:50s}: {v:4d} Genes")
EOF
python3 step2_qc_and_unphased.py

# ==============================================================================
# 6. IQ-TREE – SINGLE-GENE TREES & SUPERMATRIX
# ==============================================================================
mkdir -p iqtree_out/single_19tip iqtree_out/single_4tax_mat iqtree_out/single_4tax_pat

for f in clean_alns_19tip/*.fasta; do
    echo "iqtree -s $f -st DNA -m HKY+G -czb -T 1 --quiet --prefix iqtree_out/single_19tip/$(basename "$f" .fasta)"
done > cmd_trees.txt
for f in clean_4tax_mat/*.fasta; do
    echo "iqtree -s $f -st DNA -m HKY+G -czb -T 1 --quiet --prefix iqtree_out/single_4tax_mat/$(basename "$f" .fasta)"
done >> cmd_trees.txt
for f in clean_4tax_pat/*.fasta; do
    echo "iqtree -s $f -st DNA -m HKY+G -czb -T 1 --quiet --prefix iqtree_out/single_4tax_pat/$(basename "$f" .fasta)"
done >> cmd_trees.txt

echo "=== Running single-gene IQ-TREE jobs in parallel ==="
parallel -j 32 < cmd_trees.txt

cat iqtree_out/single_19tip/*.treefile > iqtree_out/clean_19tip_all.trees
cat iqtree_out/single_4tax_mat/*.treefile > iqtree_out/clean_4tax_mat_all.trees
cat iqtree_out/single_4tax_pat/*.treefile > iqtree_out/clean_4tax_pat_all.trees

echo "=== Running Supermatrix & Concordance (19-tip, Mat 4-tax, Pat 4-tax) ==="
iqtree -s supermatrix_19tip.fasta -st DNA -m GTR+G -bb 1000 -o Outgroup -T 16 --quiet --prefix iqtree_out/clean_19tip_supermatrix
iqtree -t iqtree_out/clean_19tip_supermatrix.treefile --gcf iqtree_out/clean_19tip_all.trees -s supermatrix_19tip.fasta -st DNA --scf 1000 -T 16 --quiet --prefix iqtree_out/clean_19tip_concordance

iqtree -s supermatrix_4tax_mat.fasta -st DNA -m GTR+G -o Parent1_Ref -T 16 --quiet --prefix iqtree_out/clean_4tax_mat_supermatrix
iqtree -t iqtree_out/clean_4tax_mat_supermatrix.treefile --gcf iqtree_out/clean_4tax_mat_all.trees -s supermatrix_4tax_mat.fasta -st DNA --scf 1000 -T 16 --quiet --prefix iqtree_out/clean_4tax_mat_concordance

iqtree -s supermatrix_4tax_pat.fasta -st DNA -m GTR+G -o C09 -T 16 --quiet --prefix iqtree_out/clean_4tax_pat_supermatrix
iqtree -t iqtree_out/clean_4tax_pat_supermatrix.treefile --gcf iqtree_out/clean_4tax_pat_all.trees -s supermatrix_4tax_pat.fasta -st DNA --scf 1000 -T 16 --quiet --prefix iqtree_out/clean_4tax_pat_concordance

# ==============================================================================
# 7. R – UNPHASED NJ TREES & CONCORDANCE SUMMARY TABLES (gCF/eGCF/sCF)
# ==============================================================================
cat << 'EOF' > step3_summarize_and_trees.R
library(ape)

d_all   <- as.matrix(read.table("unphased_distances_ALL_LOCI.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE))
d_clean <- as.matrix(read.table("unphased_distances_CLEAN_LOCI.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE))

nj_all   <- root(nj(as.dist(d_all)), outgroup="Outgroup", resolve.root=TRUE)
nj_clean <- root(nj(as.dist(d_clean)), outgroup="Outgroup", resolve.root=TRUE)

write.tree(nj_all, file="report/UNPHASED_DIPLOID_NJ_ALL_LOCI.treefile")
write.tree(nj_clean, file="report/UNPHASED_DIPLOID_NJ_CLEAN_LOCI.treefile")

process_cf <- function(stat_file, out_csv) {
  cf <- read.table(stat_file, header=TRUE, sep="\t", comment.char="#")
  cf$Total_Genes    <- round(cf$gCF_N + cf$gDF1_N + cf$gDF2_N + cf$gDFP_N)
  cf$Resolved_Genes <- cf$gCF_N + cf$gDF1_N + cf$gDF2_N
  cf$eGCF           <- ifelse(cf$Resolved_Genes > 0, round(100 * cf$gCF_N / cf$Resolved_Genes, 2), 0.0)
  cf$eDF1           <- ifelse(cf$Resolved_Genes > 0, round(100 * cf$gDF1_N / cf$Resolved_Genes, 2), 0.0)
  cf$eDF2           <- ifelse(cf$Resolved_Genes > 0, round(100 * cf$gDF2_N / cf$Resolved_Genes, 2), 0.0)
  cf$Binom_P_vs_1_3 <- sapply(1:nrow(cf), function(i) {
    if (cf$Resolved_Genes[i] == 0) return(1.0)
    binom.test(cf$gCF_N[i], cf$Resolved_Genes[i], p = 1/3, alternative = "greater")$p.value
  })
  write.csv(cf, out_csv, row.names=FALSE)
  return(cf)
}

cf19  <- process_cf("iqtree_out/clean_19tip_concordance.cf.stat", "report/CLEAN_19TIP_CONCORDANCE_SUMMARY.csv")
cfmat <- process_cf("iqtree_out/clean_4tax_mat_concordance.cf.stat", "report/CLEAN_4TAX_MAT_CONCORDANCE_SUMMARY.csv")
cfpat <- process_cf("iqtree_out/clean_4tax_pat_concordance.cf.stat", "report/CLEAN_4TAX_PAT_CONCORDANCE_SUMMARY.csv")

cat("\n==========================================================================\n")
cat("4. UNPHASED DIPLOID NJ TREES (ALL LOCI vs CLEAN LOCI)\n")
cat("==========================================================================\n")
cat("NJ Tree (All Loci):  ", write.tree(nj_all), "\n")
cat("NJ Tree (Clean Loci):", write.tree(nj_clean), "\n")

cat("\n==========================================================================\n")
cat("5. MATERNAL & PATERNAL 4-TAXON CONCORDANCE (gCF / eGCF / sCF)\n")
cat("==========================================================================\n")
cat("Maternal 4-Taxon Tree:", readLines("iqtree_out/clean_4tax_mat_concordance.cf.tree"), "\n")
print(cfmat[, c("ID", "gCF_N", "gDF1_N", "gDF2_N", "gDFP_N", "Resolved_Genes", "gCF", "eGCF", "Binom_P_vs_1_3", "sCF", "sCF_N", "sDF1_N", "sDF2_N")])

cat("\nPaternal 4-Taxon Tree:", readLines("iqtree_out/clean_4tax_pat_concordance.cf.tree"), "\n")
print(cfpat[, c("ID", "gCF_N", "gDF1_N", "gDF2_N", "gDFP_N", "Resolved_Genes", "gCF", "eGCF", "Binom_P_vs_1_3", "sCF", "sCF_N", "sDF1_N", "sDF2_N")])

cat("\n==========================================================================\n")
cat("6. CLEAN 19-TIP SUPERMATRIX TREE & CONCORDANCE\n")
cat("==========================================================================\n")
cat("19-Tip Branch IDs: ", readLines("iqtree_out/clean_19tip_concordance.cf.branch"), "\n")
cat("19-Tip Concordance:", readLines("iqtree_out/clean_19tip_concordance.cf.tree"), "\n")
print(cf19[, c("ID", "gCF_N", "gDF1_N", "gDF2_N", "gDFP_N", "Resolved_Genes", "gCF", "eGCF", "Binom_P_vs_1_3", "sCF", "Label")])
EOF

Rscript step3_summarize_and_trees.R | tee STEP3_TREES_AND_CONCORDANCE.txt
cat STEP2_QC_SUMMARY.txt STEP3_TREES_AND_CONCORDANCE.txt > report/MASTER_PIPELINE_REPORT.txt

# ==============================================================================
# 8. COPY FINAL PRODUCTS BACK TO PERMANENT STORAGE
# ==============================================================================
echo "=== Copying final results to permanent storage ==="
cp -r iqtree_out clean_alns_19tip clean_4tax_mat clean_4tax_pat phased report "${BASE_OUT}/"
cp GENE_QC_CLASSIFICATION.tsv FILTER_SUMMARY.tsv *.tsv *.nex *.fasta "${BASE_OUT}/"

echo "=== PART 2 PIPELINE FINISHED SUCCESSFULLY at $(date) ==="

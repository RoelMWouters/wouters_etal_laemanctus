#!/bin/bash
#PBS -N jcvi_3track_synteny
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=50gb
#PBS -l walltime=04:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 2: 3-TRACK JCVI KARYOTYPE GENERATION & VISUALIZATION
# ==============================================================================
# Description:
# This script generates publication-quality, 3-track synteny karyotypes. It takes 
# the filtered ortholog mapping from Script 1, formats it into BED/BLAST mock 
# inputs required by JCVI, and runs the native MCScan algorithm. 
# 
# Key Features:
# - Strict HPC environment isolation to prevent library clashes (e.g., NumPy/Brotli).
# - Dynamic Python "monkey-patching" of JCVI's core graphics modules to remove 
#   obstructive bounding boxes and adjust synteny ribbon transparency.
# - Custom sorting to maintain biological subgenome phasing nomenclature (_lo / _se).
#
# Dependencies: mambaforge/conda, python=3.10, pandas, jcvi, ete4
# ==============================================================================

export TMPDIR=$SCRATCHDIR

if [ -z "$SCRATCHDIR" ] ; then
  echo "Error: SCRATCHDIR is not set!"
  exit 1
fi
trap 'clean_scratch' TERM EXIT

cd "$SCRATCHDIR" || exit 1

# Force headless backend for Matplotlib to prevent UI crashes on HPC nodes
echo "backend: Agg" > matplotlibrc
echo "font.family: sans-serif" >> matplotlibrc
echo "font.sans-serif: DejaVu Sans" >> matplotlibrc
echo "font.size: 14" >> matplotlibrc

echo "=== INITIALIZING ISOLATED MAMBA ENVIRONMENT ==="
module add mambaforge

export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
export CONDA_ENVS_PATH="$SCRATCHDIR/conda_envs"
mkdir -p "$CONDA_PKGS_DIRS" "$CONDA_ENVS_PATH"

eval "$(conda shell.bash hook)"

echo "=== BUILDING JCVI ENVIRONMENT ==="
mamba create -y -p "$SCRATCHDIR/jcvi_env" -c conda-forge -c bioconda python=3.10 pandas jcvi
mamba activate "$SCRATCHDIR/jcvi_env"

# Prevent Python from loading corrupted/conflicting packages from ~/.local
export PYTHONNOUSERSITE=1
unset PYTHONPATH

# Pin numpy to < 2.0 to prevent Pandas compatibility crashes
pip install ete4 "numpy<2"

echo "=== TESTING MODULE IMPORTS ==="
if ! python -c "import pandas; import jcvi; import ete4"; then
    echo "CRITICAL ERROR: Failed to load modules. Check the Python traceback above."
    exit 1
fi

# ==============================================================================
# PYTHON MONKEY-PATCH (FIXING JCVI DEFAULTS GLOBALLY)
# ==============================================================================
echo "=== APPLYING AESTHETIC PATCHES TO JCVI ==="
python3 - << 'EOF'
import os
import glob
import re

graphics_dirs = glob.glob(os.path.join(os.environ.get("SCRATCHDIR", ""), "jcvi_env/lib/python*/site-packages/jcvi/graphics"))
if graphics_dirs:
    # Loop through ALL python files in the jcvi graphics folder to enforce aesthetic changes
    for file_path in glob.glob(os.path.join(graphics_dirs[0], "*.py")):
        with open(file_path, "r") as f:
            code = f.read()

        # 1. Stop JCVI from truncating custom subgenome suffixes (e.g., '1_lo' into 'L')
        code = code.replace("def make_circle_name(sid, rev):", "def make_circle_name(sid, rev):\n    return sid\n\ndef old_make_circle_name(sid, rev):")

        # 2. Obliterate the white backgrounds on text bounding boxes across ALL scripts
        code = code.replace('fc="w"', 'fc="none"').replace("fc='w'", "fc='none'")
        code = code.replace('ec="w"', 'ec="none"').replace("ec='w'", "ec='none'")
        
        # Failsafe regex to strip the circle bounding box entirely while keeping the text label
        code = re.sub(r'bbox\s*=\s*dict\([^)]+boxstyle\s*=\s*["\']circle["\'][^)]*\)', 'bbox=dict(fc="none", ec="none", pad=0)', code)

        # 3. Dim the synteny ribbons to make chromosomes stand out better.
        code = code.replace("alpha=0.2", "alpha=0.08")

        with open(file_path, "w") as f:
            f.write(code)
    print("JCVI successfully monkey-patched across all graphics files!")
else:
    print("Error: Could not find jcvi/graphics directory to patch.")
EOF

# ==============================================================================
# PATHS & DATA DOWNLOAD
# ==============================================================================
# Update these paths to point to your Script 1 output and desired final output directory
TSV_INPUT="/path/to/synteny_results/clean_phased_blocks.tsv"
OUT_DIR="/path/to/output_directory/jcvi_karyotype_plots"

mkdir -p "$OUT_DIR"

echo "=== DOWNLOADING REFERENCE GFF ==="
wget -nc -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/765/GCF_037176765.1_rAnoSag1.mat/GCF_037176765.1_rAnoSag1.mat_genomic.gff.gz
gunzip -f GCF_037176765.1_rAnoSag1.mat_genomic.gff.gz

echo "=== STAGING TSV DATA ==="
cp "$TSV_INPUT" ./clean_phased_blocks.tsv

# ==============================================================================
# PREPARE BED FILES & MOCK BLAST FORMAT
# ==============================================================================
echo "=== FORMATTING DATA & BUILDING MOCK BLAST ==="
python3 - << 'EOF'
import pandas as pd
import re

def clean_bed(df):
    df = df.iloc[:, :6].copy()
    df.columns = ["chrom", "start", "end", "name", "score", "strand"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df[df["strand"].isin(["+", "-"])]
    return df

# Load the filtered 1:2 orthologs from Script 1
df = pd.read_csv("clean_phased_blocks.tsv", sep='\t', header=None,
                 names=["Anolis_Raw", "Contig", "Start", "End", "Strand", "Map_Length"])

def extract_prot(x):
    try:
        return x.split("_cds_")[1].rsplit("_", 1)[0]
    except:
        return x
df["Anolis_Gene"] = df["Anolis_Raw"].apply(extract_prot)

# Parse the reference GFF to establish coordinate space
gff_data = []
with open("GCF_037176765.1_rAnoSag1.mat_genomic.gff", "r") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.split("\t")
        if len(parts) < 9 or parts[2] != "CDS": continue
        
        match = re.search(r"protein_id=([^;\r\n]+)", parts[8])
        if match:
            prot_id = match.group(1).replace('"', '').strip()
            chrom_clean = parts[0].strip()
            gff_data.append({
                "chrom": chrom_clean, "start": parts[3], "end": parts[4],
                "name": prot_id, "score": "0", "strand": parts[6].strip()
            })

anolis_df = pd.DataFrame(gff_data)
anolis_df["start"] = pd.to_numeric(anolis_df["start"], errors="coerce")
anolis_df["end"] = pd.to_numeric(anolis_df["end"], errors="coerce")
anolis_df = anolis_df.dropna()

anolis_bed_raw = anolis_df.groupby(["name", "chrom", "strand"]).agg({"start": "min", "end": "max"}).reset_index()
anolis_bed_raw["score"] = "0"
anolis_bed = clean_bed(anolis_bed_raw[["chrom", "start", "end", "name", "score", "strand"]])

# Define raw scaffold groupings for Subgenome 1 and Subgenome 2
# UPDATE THESE LISTS with your specific scaffold/chromosome IDs
raw_subgenome_1 = ["scaffold_104", "scaffold_105", "scaffold_106", "scaffold_107", "scaffold_108", "scaffold_109", "scaffold_11", "scaffold_110", "scaffold_13", "scaffold_15", "scaffold_17", "scaffold_18", "scaffold_2", "scaffold_21", "scaffold_23", "scaffold_25", "scaffold_27", "scaffold_29", "scaffold_31", "scaffold_32", "scaffold_35", "scaffold_37", "scaffold_40", "scaffold_42", "scaffold_43", "scaffold_44", "scaffold_46", "scaffold_49", "scaffold_5", "scaffold_54", "scaffold_55", "scaffold_56", "scaffold_61", "scaffold_7", "scaffold_9"]
raw_subgenome_2 = ["scaffold_1", "scaffold_10", "scaffold_12", "scaffold_14", "scaffold_16", "scaffold_19", "scaffold_20", "scaffold_22", "scaffold_24", "scaffold_26", "scaffold_28", "scaffold_3", "scaffold_30", "scaffold_33", "scaffold_34", "scaffold_36", "scaffold_38", "scaffold_39", "scaffold_4", "scaffold_41", "scaffold_45", "scaffold_47", "scaffold_48", "scaffold_50", "scaffold_51", "scaffold_52", "scaffold_53", "scaffold_57", "scaffold_58", "scaffold_6", "scaffold_60", "scaffold_62", "scaffold_66", "scaffold_68", "scaffold_69", "scaffold_70", "scaffold_73", "scaffold_8", "scaffold_80", "scaffold_81", "scaffold_82", "scaffold_90", "scaffold_92", "scaffold_93", "scaffold_94"]

# Explicitly assign exact pairs based on biological layout (e.g., Pair 18 is XX)
# UPDATE THESE TUPLES based on your assembly's paired scaffolds
pair_list = [
    ("scaffold_1", "scaffold_2"),
    ("scaffold_3", "scaffold_5"),
    ("scaffold_6", "scaffold_7"),
    ("scaffold_8", "scaffold_9"),
    ("scaffold_10", "scaffold_11"),
    ("scaffold_12", "scaffold_13"),
    ("scaffold_14", "scaffold_15"),
    ("scaffold_16", "scaffold_17"),
    ("scaffold_18", "scaffold_19"),
    ("scaffold_20", "scaffold_21"),
    ("scaffold_22", "scaffold_23"),
    ("scaffold_24", "scaffold_25"),
    ("scaffold_26", "scaffold_27"),
    ("scaffold_30", "scaffold_31"),
    ("scaffold_32", "scaffold_33"),
    ("scaffold_34", "scaffold_35"),
    ("scaffold_36", "scaffold_37"),
    ("scaffold_28", "scaffold_29")  # Sex Chromosome Pair (XX)
]

# Map original scaffold names to biological nomenclature (e.g., 1_lo, 1_se)
target_mapping = {}
for idx, (s1, s2) in enumerate(pair_list):
    pair_num = "X" if idx == 17 else str(idx + 1)
    
    if s1 in raw_subgenome_1:
        target_mapping[s1] = f"{pair_num}_lo"
        target_mapping[s2] = f"{pair_num}_se"
    elif s2 in raw_subgenome_1:
        target_mapping[s2] = f"{pair_num}_lo"
        target_mapping[s1] = f"{pair_num}_se"
    else:
        target_mapping[s1] = f"{pair_num}_lo"
        target_mapping[s2] = f"{pair_num}_se"

df_sub1_raw = df[df["Contig"].isin(raw_subgenome_1)].copy()
df_sub2_raw = df[df["Contig"].isin(raw_subgenome_2)].copy()

sub1_bed = clean_bed(df_sub1_raw[df_sub1_raw["Contig"].isin(target_mapping.keys())][["Contig", "Start", "End", "Anolis_Gene", "Map_Length", "Strand"]])
sub2_bed = clean_bed(df_sub2_raw[df_sub2_raw["Contig"].isin(target_mapping.keys())][["Contig", "Start", "End", "Anolis_Gene", "Map_Length", "Strand"]])

mapped_anolis_genes = set(sub1_bed["name"]).union(set(sub2_bed["name"]))
mapped_anolis_bed = anolis_bed[anolis_bed["name"].isin(mapped_anolis_genes)]
top_anolis_contigs = sorted(mapped_anolis_bed["chrom"].value_counts().head(15).index.tolist())

# Rename Reference chromosomes for clean visualization (e.g., mapping NC_ identifiers to integers)
anolis_rename_map = {}
for i, old in enumerate(top_anolis_contigs):
    idx = i + 1
    if idx == 14:
        anolis_rename_map[old] = "X"
    elif idx == 15:
        anolis_rename_map[old] = "Y"
    else:
        anolis_rename_map[old] = str(idx)

anolis_bed = anolis_bed[anolis_bed["chrom"].isin(top_anolis_contigs)].copy()
anolis_bed["chrom"] = anolis_bed["chrom"].map(anolis_rename_map)

sub1_bed["chrom"] = sub1_bed["chrom"].map(target_mapping)
sub2_bed["chrom"] = sub2_bed["chrom"].map(target_mapping)

anolis_bed.to_csv("reference.bed", sep='\t', header=False, index=False)
sub1_bed.to_csv("subgenome1.bed", sep='\t', header=False, index=False)
sub2_bed.to_csv("subgenome2.bed", sep='\t', header=False, index=False)

# GENERATE MOCK BLAST OUTPUTS FOR JCVI NATIVE SCAN
common_rs1 = set(anolis_bed["name"]).intersection(set(sub1_bed["name"]))
with open("reference.subgenome1.blast", "w") as f:
    for g in common_rs1:
        f.write(f"{g}\t{g}\t100.0\t1000\t0\t0\t1\t1000\t1\t1000\t1e-50\t1000\n")

common_s1s2 = set(sub1_bed["name"]).intersection(set(sub2_bed["name"]))
with open("subgenome1.subgenome2.blast", "w") as f:
    for g in common_s1s2:
        f.write(f"{g}\t{g}\t100.0\t1000\t0\t0\t1\t1000\t1\t1000\t1e-50\t1000\n")

# Lock in strictly explicit order for the tracks
ref_seqs = [str(i) for i in range(1, 14)] + ["X", "Y"]
sub1_seqs = [f"{i}_lo" for i in range(1, 18)] + ["X_lo"]
sub2_seqs = [f"{i}_se" for i in range(1, 18)] + ["X_se"]

with open("seqids", "w") as f:
    f.write(",".join(ref_seqs) + "\n")
    f.write(",".join(sub1_seqs) + "\n")
    f.write(",".join(sub2_seqs) + "\n")
EOF

# ==============================================================================
# NATIVE MCSCAN ALGORITHM
# ==============================================================================
echo "=== RUNNING NATIVE MCSCAN ORTHOLOG PIPELINE ==="

python -m jcvi.compara.synteny scan reference.subgenome1.blast reference.subgenome1.anchors --qbed=reference.bed --sbed=subgenome1.bed
python -m jcvi.compara.synteny screen reference.subgenome1.anchors reference.subgenome1.filtered.anchors --minspan=5 --simple --qbed=reference.bed --sbed=subgenome1.bed

python -m jcvi.compara.synteny scan subgenome1.subgenome2.blast subgenome1.subgenome2.anchors --qbed=subgenome1.bed --sbed=subgenome2.bed
python -m jcvi.compara.synteny screen subgenome1.subgenome2.anchors subgenome1.subgenome2.filtered.anchors --minspan=15 --simple --qbed=subgenome1.bed --sbed=subgenome2.bed

# ==============================================================================
# PLOT CONFIGURATION & RENDERING
# ==============================================================================
echo "=== GENERATING KARYOTYPE ==="

cat << 'EOF' > layout
0.85,    0.22,    0.98,      0,      #FF4500, Reference, top, reference.bed
0.50,    0.22,    0.98,      0,      #78D134, Subgenome 1, top, subgenome1.bed
0.15,    0.22,    0.98,      0,      #12E0AC, Subgenome 2, top, subgenome2.bed
# edges mapping to the .simple output files generated by the screen command
e, 0, 1, reference.subgenome1.filtered.simple
e, 1, 2, subgenome1.subgenome2.filtered.simple
EOF

python -m jcvi.graphics.karyotype seqids layout --basepair --figsize=36x15 --dpi=300 -o 3Track_Synteny_Plot.png
python -m jcvi.graphics.karyotype seqids layout --basepair --figsize=36x15 -o 3Track_Synteny_Plot.svg

echo "=== SAVING RESULTS TO STORAGE ==="
cp 3Track_Synteny_Plot.* "$OUT_DIR/"
cp *.bed *.blast *.anchors* seqids layout "$OUT_DIR/"

echo "Pipeline fully completed! Native MCScan karyotypes saved."

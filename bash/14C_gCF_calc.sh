#!/bin/bash
#PBS -N gcf_busco_pipeline
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=20gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 14C: SUPERMATRIX AND GENE CONCORDANCE FACTOR (gCF) PIPELINE
# ==============================================================================
# Description:
# 1. Sets up an isolated Python environment inside local fast scratch storage.
# 2. Concatenates previously cleaned alignments into a supermatrix.
# 3. Infers individual gene trees using IQ-TREE.
# 4. Builds a rooted reference tree from the supermatrix.
# 5. Calculates Gene Concordance Factors (gCF) with verbose discordance tables.
# ==============================================================================

# Load required modules
module add python
module add iqtree/3.0.1-mpi

# --- Configuration Paths ---
STORAGE_DIR="/path/to/your/project/output_dir/busco_phylogeny"
ALIGN_DIR="${STORAGE_DIR}/alignments_clean"
FINAL_OUT="${STORAGE_DIR}/gcf_final_output"

mkdir -p "${FINAL_OUT}"

# ------------------------------------------------------------------------------
# HPC SCRATCHDIR SETUP
# ------------------------------------------------------------------------------
if [ -z "$SCRATCHDIR" ]; then
    echo "ERROR: SCRATCHDIR is not defined!"
    exit 1
fi

echo "[$(date +'%T')] >>> Initializing scratch workspace in $SCRATCHDIR"
trap 'rm -rf "$SCRATCHDIR"/*' TERM EXIT

SCRATCH_ALIGN="$SCRATCHDIR/alignments_clean"
SCRATCH_OUT="$SCRATCHDIR/gcf_final"
SCRATCH_FILTERED_ALIGN="${SCRATCH_OUT}/alignments_strictly_filtered"

mkdir -p "${SCRATCH_OUT}"
mkdir -p "${SCRATCH_FILTERED_ALIGN}"

# Checkpoint restoration
echo "[$(date +'%T')] >>> Checking for previous progress..."
if [ "$(ls -A ${FINAL_OUT} 2>/dev/null)" ]; then
    echo "[$(date +'%T')]     -> Found existing files. Restoring checkpoints..."
    cp -r "${FINAL_OUT}/"* "${SCRATCH_OUT}/"
fi

cp -r "${ALIGN_DIR}" "$SCRATCHDIR/" || { echo "Failed to copy input alignments"; exit 1; }

cd "$SCRATCHDIR" || exit 1

# ------------------------------------------------------------------------------
# 1. STANDALONE PYTHON ENVIRONMENT
# ------------------------------------------------------------------------------
if [ ! -s "${SCRATCH_OUT}/filtered_supermatrix.phy" ]; then
    echo "[$(date +'%T')] >>> [1/5] Setting up isolated Python environment..."
    python3 -m venv python_env_gcf
    source python_env_gcf/bin/activate
    pip install --upgrade pip setuptools > /dev/null 2>&1
    pip install biopython numpy > /dev/null 2>&1
else
    echo "[$(date +'%T')] >>> [1/5] Supermatrix already exists. Skipping setup."
fi

# ------------------------------------------------------------------------------
# 2. BUILD SUPERMATRIX
# ------------------------------------------------------------------------------
if [ -s "${SCRATCH_OUT}/filtered_supermatrix.phy" ]; then
    echo "[$(date +'%T')] >>> [2/5] CHECKPOINT: Supermatrix found. Skipping."
else
    echo "[$(date +'%T')] >>> [2/5] Filtering and building supermatrix..."
    python3 -c "
import os
import re
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

align_dir = '${SCRATCH_ALIGN}'
filtered_align_dir = '${SCRATCH_FILTERED_ALIGN}'
out_supermatrix = os.path.join('${SCRATCH_OUT}', 'filtered_supermatrix.phy')

alignment_files = [f for f in os.listdir(align_dir) if f.endswith(('.fasta', '.fa', '.phy', '.phylip'))]
surviving_files = []
dropped_excessive = 0
dropped_too_small = 0

exclude_patterns = []

for filename in alignment_files:
    in_path = os.path.join(align_dir, filename)
    fmt = 'fasta' if filename.endswith(('.fasta', '.fa')) else 'phylip'
    try:
        align = AlignIO.read(in_path, fmt)
    except:
        try:
            align = AlignIO.read(in_path, 'fasta')
        except:
            continue
            
    cleaned_records = [rec for rec in align if not any(excl in rec.id for excl in exclude_patterns)]
            
    has_excessive = False
    sample_counts = {}
    
    for record in cleaned_records:
        name = record.id
        if re.search(r'_Copy[3-9]', name):
            has_excessive = True
            break
        
        base = re.sub(r'_Copy\d+', '', name)
        sample_counts[base] = sample_counts.get(base, 0) + 1
        if sample_counts[base] > 2:
            has_excessive = True
            break

    if has_excessive:
        dropped_excessive += 1
        continue
        
    if len(cleaned_records) < 4:
        dropped_too_small += 1
        continue

    cleaned_align = MultipleSeqAlignment(cleaned_records)
    out_path = os.path.join(filtered_align_dir, filename)
    AlignIO.write(cleaned_align, out_path, fmt)
    surviving_files.append(out_path)

print(f'\n>>> ALIGNMENTS: Dropped {dropped_excessive} containing >2 copies.')
print(f'>>> ALIGNMENTS: Dropped {dropped_too_small} for having <4 taxa.')
print(f'>>> ALIGNMENTS: Retained {len(surviving_files)} completely clean alignments.')

all_taxa = set()
locus_data = []

for filepath in surviving_files:
    fmt = 'fasta' if filepath.endswith(('.fasta', '.fa')) else 'phylip'
    align = AlignIO.read(filepath, fmt)
    length = align.get_alignment_length()
    seq_dict = {r.id: str(r.seq) for r in align}
    all_taxa.update(seq_dict.keys())
    locus_data.append((length, seq_dict))

all_taxa = sorted(list(all_taxa))
concat_seqs = {t: '' for t in all_taxa}

for length, seq_dict in locus_data:
    for t in all_taxa:
        concat_seqs[t] += seq_dict.get(t, '-' * length)

final_align = MultipleSeqAlignment([SeqRecord(Seq(concat_seqs[t]), id=t, description='') for t in all_taxa])
AlignIO.write(final_align, out_supermatrix, 'phylip-relaxed')
print(f'>>> SUPERMATRIX: Built with {len(all_taxa)} taxa and {final_align.get_alignment_length()} bp.\n')
"
    cp -r "${SCRATCH_FILTERED_ALIGN}" "${FINAL_OUT}/"
    cp "${SCRATCH_OUT}/filtered_supermatrix.phy" "${FINAL_OUT}/"
fi

# ------------------------------------------------------------------------------
# 3. REGENERATE PRISTINE GENE TREES 
# ------------------------------------------------------------------------------
if [ -s "${SCRATCH_OUT}/filtered_genetrees.trees" ]; then
    echo "[$(date +'%T')] >>> [3/5] CHECKPOINT: Gene trees found. Skipping."
else
    echo "[$(date +'%T')] >>> [3/5] Inferring pristine gene trees..."
    cd "${SCRATCH_FILTERED_ALIGN}" || exit 1

    GENE_TREES_OUT="${SCRATCH_OUT}/filtered_genetrees.trees"
    > "${GENE_TREES_OUT}"

    count=0
    for aln in *.{phy,phylip,fasta,fa}; do
        [ -e "$aln" ] || continue
        
        iqtree3-mpi -s "$aln" -m MFP -fast -T 4 --quiet > /dev/null 2>&1
        
        if [ -f "${aln}.treefile" ]; then
            cat "${aln}.treefile" >> "${GENE_TREES_OUT}"
            rm -f "${aln}."*
            ((count++))
        fi
    done
    echo "[$(date +'%T')]     -> Successfully generated ${count} pristine gene trees."
    
    cp "${GENE_TREES_OUT}" "${FINAL_OUT}/"
fi

# ------------------------------------------------------------------------------
# 4. REFERENCE TOPOLOGY
# ------------------------------------------------------------------------------
cd "${SCRATCH_OUT}" || exit 1

if [ -s "reference_tree.treefile" ]; then
    echo "[$(date +'%T')] >>> [4/5] CHECKPOINT: Reference tree found. Skipping."
else
    echo "[$(date +'%T')] >>> [4/5] Building ROOTED reference topology..."
    iqtree3-mpi -s filtered_supermatrix.phy -m MFP -B 1000 -T 8 -o OUTGROUP --prefix reference_tree

    if [ ! -f "reference_tree.treefile" ]; then
        echo "ERROR: IQ-TREE failed to build reference_tree.treefile"
        exit 1
    fi
    
    cp reference_tree* "${FINAL_OUT}/"
fi

# ------------------------------------------------------------------------------
# 5. CONCORDANCE FACTOR ANALYSIS
# ------------------------------------------------------------------------------
echo "[$(date +'%T')] >>> [5/5] Running IQ-TREE Concordance Factor (gCF) analysis..."
iqtree3-mpi -t reference_tree.treefile --gcf filtered_genetrees.trees --prefix subgenome_concordance -T 8 --cf-verbose

# ------------------------------------------------------------------------------
# FINAL TRANSFER
# ------------------------------------------------------------------------------
echo "[$(date +'%T')] >>> Copying final results back to permanent storage..."
cp -r "${SCRATCH_OUT}/"* "${FINAL_OUT}/" || { echo "Failed to copy results back to storage"; exit 1; }

echo "[$(date +'%T')] >>> Pipeline completely finished! Your files are ready for R in ${FINAL_OUT}"

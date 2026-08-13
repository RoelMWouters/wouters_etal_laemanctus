#!/bin/bash
#PBS -N gcf_busco_pipeline
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=20gb
#PBS -l walltime=04:00:00
#PBS -M your_email@institution.cz
#PBS -m abe
#PBS -j oe

# ==============================================================================
# PIPELINE OVERVIEW:
# 1. Sets up an isolated Python environment inside local fast scratch storage.
# 2. Filters individual gene alignments (drops loci with >2 copies per sample,
#    removes specific unwanted taxa/references, and discards alignments <4 taxa).
# 3. Concatenates surviving alignments into a supermatrix.
# 4. Infers individual gene trees using IQ-TREE.
# 5. Builds a rooted reference tree from the supermatrix.
# 6. Calculates Gene Concordance Factors (gCF) and outputs summary statistics.
# ==============================================================================

# Load required cluster modules (adjust module names according to your HPC environment)
module add python
module add iqtree/3.0.1-mpi

# Define paths on permanent storage (MODIFY THESE TO MATCH YOUR DIRECTORY STRUCTURE)
PROJECT_DIR="/path/to/your/project/directory"
STORAGE_DIR="${PROJECT_DIR}/output_dir/busco_phylogeny"
ALIGN_DIR="${STORAGE_DIR}/alignments_clean"
FINAL_OUT_EXCLUDED="${STORAGE_DIR}/gcf_300_final_excluded"

# ------------------------------------------------------------------------------
# META CENTRUM / HPC SCRATCHDIR SETUP
# ------------------------------------------------------------------------------
if [ -z "$SCRATCHDIR" ]; then
    echo "ERROR: SCRATCHDIR is not defined!"
    exit 1
fi

echo ">>> Initializing scratch workspace in $SCRATCHDIR"
# Bulletproof scratch cleanup trap (removes scratch contents on exit or termination)
trap 'rm -rf "$SCRATCHDIR"/*' TERM EXIT

# Copy raw alignment files into fast local scratch storage
cp -r "${ALIGN_DIR}" "$SCRATCHDIR/" || { echo "Failed to copy input alignments"; exit 1; }

# Define internal scratch working paths
SCRATCH_ALIGN="$SCRATCHDIR/alignments_clean"
SCRATCH_OUT="$SCRATCHDIR/gcf_300_final"
SCRATCH_FILTERED_ALIGN="${SCRATCH_OUT}/alignments_strictly_filtered"

mkdir -p "${SCRATCH_OUT}"
mkdir -p "${SCRATCH_FILTERED_ALIGN}"

# Move entirely into scratch space for high-performance file operations
cd "$SCRATCHDIR" || exit 1

# ------------------------------------------------------------------------------
# 1. STANDALONE PYTHON ENVIRONMENT IN SCRATCH
# ------------------------------------------------------------------------------
echo ">>> [1/5] Setting up isolated Python environment in SCRATCHDIR..."
python3 -m venv python_env_gcf
source python_env_gcf/bin/activate
pip install --upgrade pip setuptools > /dev/null 2>&1
pip install biopython numpy > /dev/null 2>&1

# ------------------------------------------------------------------------------
# 2. FILTER ALIGNMENTS & BUILD SUPERMATRIX
# ------------------------------------------------------------------------------
echo ">>> [2/5] Running strict filtering and building supermatrix..."
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

# Define specific sample IDs or patterns to completely exclude from the analysis
exclude_patterns = ['C09', 'C10', 'PATERNAL_REF']

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
            
    # Purge any records matching the exclusion patterns
    cleaned_records = [rec for rec in align if not any(excl in rec.id for excl in exclude_patterns)]
            
    has_excessive = False
    sample_counts = {}
    
    # Check for paralogs / copy-number violations (>2 copies per base taxon)
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
        
    # Ensure the alignment maintains a minimum taxon count for phylogenetic inference
    if len(cleaned_records) < 4:
        dropped_too_small += 1
        continue

    cleaned_align = MultipleSeqAlignment(cleaned_records)
    out_path = os.path.join(filtered_align_dir, filename)
    AlignIO.write(cleaned_align, out_path, fmt)
    surviving_files.append(out_path)

print(f'\n>>> ALIGNMENTS: Dropped {dropped_excessive} containing excessive copies.')
print(f'>>> ALIGNMENTS: Dropped {dropped_too_small} for having <4 taxa after exclusions.')
print(f'>>> ALIGNMENTS: Retained {len(surviving_files)} completely clean alignments.')

# Concatenate surviving alignments into a supermatrix
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

# ------------------------------------------------------------------------------
# 3. REGENERATE INDIVIDUAL GENE TREES 
# ------------------------------------------------------------------------------
echo ">>> [3/5] Inferring individual gene trees from the filtered alignments..."
cd "${SCRATCH_FILTERED_ALIGN}" || exit 1

GENE_TREES_OUT="${SCRATCH_OUT}/filtered_genetrees.trees"
> "${GENE_TREES_OUT}"

count=0
for aln in *.{phy,phylip,fasta,fa}; do
    [ -e "$aln" ] || continue
    
    # Fast ML tree search per locus using IQ-TREE
    iqtree3-mpi -s "$aln" -m MFP -fast -T 4 --quiet > /dev/null 2>&1
    
    if [ -f "${aln}.treefile" ]; then
        cat "${aln}.treefile" >> "${GENE_TREES_OUT}"
        rm -f "${aln}."*  # Clean temporary output files to conserve scratch space
        ((count++))
    fi
done
echo "   -> Successfully generated ${count} gene trees."

# ------------------------------------------------------------------------------
# 4. REFERENCE TOPOLOGY INFERENCE
# ------------------------------------------------------------------------------
echo ">>> [4/5] Building rooted reference topology from the supermatrix..."
cd "${SCRATCH_OUT}" || exit 1

# Model selection and tree inference with ultrafast bootstraps, rooted explicitly on OUTGROUP
iqtree3-mpi -s filtered_supermatrix.phy -m MFP -B 1000 -T 8 -o OUTGROUP --prefix reference_tree

if [ ! -f "reference_tree.treefile" ]; then
    echo "ERROR: IQ-TREE failed to build reference_tree.treefile"
    exit 1
fi

# ------------------------------------------------------------------------------
# 5. CONCORDANCE FACTOR ANALYSIS (gCF)
# ------------------------------------------------------------------------------
echo ">>> [5/5] Running IQ-TREE Gene Concordance Factor (gCF) analysis..."
iqtree3-mpi -t reference_tree.treefile --gcf filtered_genetrees.trees --prefix subgenome_concordance300 -T 8

# ------------------------------------------------------------------------------
# TRANSFER RESULTS BACK TO PERMANENT STORAGE
# ------------------------------------------------------------------------------
echo ">>> Copying results back to permanent storage..."
mkdir -p "${FINAL_OUT_EXCLUDED}"
cp -r "${SCRATCH_OUT}/"* "${FINAL_OUT_EXCLUDED}/" || { echo "Failed to copy results back to storage"; exit 1; }

echo ">>> Pipeline completely finished! Results saved in ${FINAL_OUT_EXCLUDED}"

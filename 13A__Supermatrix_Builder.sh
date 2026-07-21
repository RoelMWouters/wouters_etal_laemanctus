#!/bin/bash
#PBS -N Supermatrix_Builder
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=50gb
#PBS -l walltime=12:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 13A: ALIGNMENT & SUPERMATRIX CONCATENATION
# ==============================================================================
# Description:
# Takes the individual sample/allele FASTA files generated in Script 12, groups
# them by gene, aligns them with MAFFT, and concatenates the aligned genes into
# a single 'RAxML_Supermatrix.fasta' required by Script 13.
# ==============================================================================

export TMPDIR=$SCRATCHDIR
module add mafft python

IN_DIR="/path/to/cluster/project_dir/output_dir/final_markers"
OUT_DIR="/path/to/cluster/project_dir/output_dir/supermatrix"

mkdir -p "$OUT_DIR"
cd "$SCRATCHDIR"

echo "=== 1. GROUPING & ALIGNING GENES ==="
cp "$IN_DIR"/*.fasta .

# Find all unique gene names based on the file naming convention from Script 12
ls *_Allele*.fasta | cut -d'_' -f2 | sort | uniq > gene_list.txt

while read GENE; do
    echo "Aligning $GENE..."
    cat *_${GENE}_Allele*.fasta > ${GENE}_pool.fasta
    mafft --auto --quiet ${GENE}_pool.fasta > ${GENE}_aligned.fasta
done < gene_list.txt

echo "=== 2. CONCATENATING SUPERMATRIX ==="
# A pure-Python native concatenator (no external AMAS dependency needed)
python3 - << 'EOF'
import glob
from collections import defaultdict

alignments = glob.glob("*_aligned.fasta")
super_dict = defaultdict(str)

for aln in alignments:
    with open(aln, 'r') as f:
        header = ""
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    super_dict[header] += "".join(seq)
                header = line[1:] # Remove '>'
                seq = []
            else:
                seq.append(line)
        if header:
            super_dict[header] += "".join(seq)

# Ensure all sequences are the same length by padding gaps if necessary
max_len = max(len(s) for s in super_dict.values())
for k, v in super_dict.items():
    if len(v) < max_len:
        super_dict[k] = v + ("-" * (max_len - len(v)))

with open("RAxML_Supermatrix.fasta", "w") as out:
    for header, sequence in super_dict.items():
        out.write(f">{header}\n{sequence}\n")
EOF

cp RAxML_Supermatrix.fasta "$OUT_DIR/"
echo "Supermatrix built successfully!"

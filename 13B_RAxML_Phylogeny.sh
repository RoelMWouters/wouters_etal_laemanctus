#!/bin/bash
#PBS -N RAxML_Phylogeny
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=10gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 13B: RAxML TREE INFERENCE
# ==============================================================================
# Description:
# Converts an aligned FASTA supermatrix into PHYLIP format using a short inline 
# Python script, then runs RAxML Maximum Likelihood tree inference with Rapid 
# Bootstrapping (-N 100).
#
# Dependencies: python3, raxml
# ==============================================================================

WORK_DIR="/path/to/cluster/project_dir/output_dir/supermatrix"
cd "$WORK_DIR"

module purge
module add raxml python

# --- Format Conversion (FASTA to PHYLIP) ---
cat << 'EOF' > fasta_to_phylip.py
import sys
fasta_in = "RAxML_Supermatrix.fasta"
phy_out = "RAxML_Supermatrix.phy"

with open(fasta_in, "r") as f:
    raw_data = f.read().split(">")[1:]

seqs = []
for block in raw_data:
    lines = block.strip().split("\n")
    name = lines[0].strip()
    seq = "".join(lines[1:]).strip()
    seqs.append((name, seq))

with open(phy_out, "w") as f:
    f.write(f"{len(seqs)} {len(seqs[0][1])}\n")
    for name, seq in seqs:
        f.write(f"{name:<30} {seq}\n")
EOF

python3 fasta_to_phylip.py

# --- RAxML Execution ---
raxmlHPC-PTHREADS-SSE3 -T 8 \
    -f a \
    -x 12345 -p 12345 \
    -N 100 \
    -m GTRGAMMA \
    -o outgroup_name \
    -s RAxML_Supermatrix.phy \
    -n Topology_Final

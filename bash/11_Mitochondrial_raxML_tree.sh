#!/bin/bash
#PBS -N Mito_Phylogeny
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=100gb
#PBS -l walltime=12:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# PHYLOGENOMIC PIPELINE - PART 3: MITOCHONDRIAL PHYLOGENY (GENE-BY-GENE)
# Description: Uses Biopython to auto-fetch specific mitochondrial coding genes 
# (COI, ND4, CYTB) as baits. Identifies mapping coordinates across reference and 
# outgroup genomes, extracts BAM consensus locus-by-locus, aligns each gene 
# independently with MAFFT, and concatenates them into a RAxML supermatrix.
# ==============================================================================

set -euo pipefail
cleanup() { rm -rf "${SCRATCHDIR:?}"/* 2>/dev/null || true; }
trap cleanup EXIT TERM

cd "$SCRATCHDIR" || { echo "No SCRATCHDIR – abort"; exit 1; }
echo "Job started on $(hostname) at $(date)"

# ==============================================================================
# 1. GLOBAL PATHS & WORKSPACE SETUP
# ==============================================================================
PROJECT_DIR="/path/to/project"
GENOME_DIR="${PROJECT_DIR}/genomes"
BAM_DIR="${PROJECT_DIR}/bam_output"
BASE_OUT="${PROJECT_DIR}/output/mitochondrial"

# EXPLICIT REFERENCE DEFINITIONS
OUTGROUP_REF="${GENOME_DIR}/basvit_genome.fna"
PARENT1_REF="${GENOME_DIR}/lser_ref.fna"
BAM_MAPPING_REF="${GENOME_DIR}/ljul_genome.fasta"  
# TARGET GENES
TARGET_MT_GENES=("COI" "ND4" "CYTB")

mkdir -p "$BASE_OUT"/{1_baits,2_extracted_mt,3_alignments,4_raxml_tree}
mkdir -p work_blast work_bams
cd work_blast

# ==============================================================================
# 2. ISOLATED CONDA ENVIRONMENT
# ==============================================================================
module purge
module add mambaforge

export CONDA_PKGS_DIRS="${SCRATCHDIR}/conda_pkgs"
export CONDA_ENVS_DIRS="${SCRATCHDIR}/conda_envs"
ENV_PREFIX="${SCRATCHDIR}/conda_envs/mito_env"

echo "=== Creating Conda Environment ==="
mamba create -p "${ENV_PREFIX}" -y -c conda-forge -c bioconda \
    blast bedtools samtools bcftools seqtk mafft raxml python biopython entrez-direct

eval "$(conda shell.bash hook)"
conda activate "${ENV_PREFIX}"

# ==============================================================================
# 3. AUTOMATED BAIT FETCHING (BIOPYTHON)
# ==============================================================================
echo "=== FETCHING SPECIFIC MT BAITS (COI, ND4, CYTB) ==="
cd "${BASE_OUT}/1_baits"

cat << 'EOF' > fetch_mt_baits.py
import os
from Bio import Entrez, SeqIO

Entrez.email = "user@institute.edu"
print("Downloading NC_010966.1 (Anolis carolinensis mitochondrion)...")
try:
    handle = Entrez.efetch(db="nucleotide", id="NC_010966.1", rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    # Map GenBank annotations to our desired file names
    targets = {"COX1": "COI", "ND4": "ND4", "CYTB": "CYTB"}
    extracted = 0
    
    for feature in record.features:
        if feature.type == "CDS":
            gene = feature.qualifiers.get("gene", [""])[0]
            if gene in targets:
                out_name = targets[gene]
                seq = feature.extract(record.seq)
                with open(f"{out_name}_bait.fasta", "w") as f:
                    f.write(f">{out_name}_Anolis_carolinensis\n{seq}\n")
                print(f"  -> Successfully extracted {out_name}")
                extracted += 1
                
    if extracted < 3: print("WARNING: Not all target genes were found.")
except Exception as e:
    print("Error fetching baits:", e)
EOF

python3 fetch_mt_baits.py
cd "$SCRATCHDIR/work_blast"

# ==============================================================================
# 4. GENE-BY-GENE BLAST SEARCH
# ==============================================================================
echo "=== CREATING BLAST DATABASES ==="
makeblastdb -in "$BAM_MAPPING_REF" -dbtype nucl -out mapping_ref_db >/dev/null
makeblastdb -in "$PARENT1_REF" -dbtype nucl -out parent_db >/dev/null
makeblastdb -in "$OUTGROUP_REF" -dbtype nucl -out outgroup_db >/dev/null

echo "=== LOCATING MT COORDINATES FOR EACH GENE ==="
> coords_mapping_all.bed
> coords_parent_all.bed
> coords_outgroup_all.bed

for GENE in "${TARGET_MT_GENES[@]}"; do
    BAIT="${BASE_OUT}/1_baits/${GENE}_bait.fasta"
    
    # 4A. Find gene in the BAM mapping reference
    blastn -task blastn -query "$BAIT" -db mapping_ref_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_mapping_${GENE}.txt
    awk -v g="$GENE" '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, g; else print $1, $3-1, $2, g}' blast_mapping_${GENE}.txt >> coords_mapping_all.bed

    # 4B. Find gene in Pure Parent 1 (L. serratus)
    blastn -task blastn -query "$BAIT" -db parent_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_parent_${GENE}.txt
    awk -v g="$GENE" '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, g; else print $1, $3-1, $2, g}' blast_parent_${GENE}.txt >> coords_parent_all.bed

    # 4C. Find gene in Outgroup
    blastn -task blastn -query "$BAIT" -db outgroup_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_outgroup_${GENE}.txt
    awk -v g="$GENE" '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, g; else print $1, $3-1, $2, g}' blast_outgroup_${GENE}.txt >> coords_outgroup_all.bed
done

# ==============================================================================
# 5. PURE REFERENCE MT EXTRACTION
# ==============================================================================
echo "=== EXTRACTING OUTGROUP & PARENT REFERENCE GENES ==="
cd ../work_bams

# Extract from Parent and dynamically name files based on the BED gene column
bedtools getfasta -fi "$PARENT1_REF" -bed ../work_blast/coords_parent_all.bed -nameOnly > parent1_all_mt.fasta
awk '/^>/ {gene=substr($1,2); file="Parent1_Ref_" gene ".fasta"; print ">Parent1_Ref" > file; next} {print >> file}' parent1_all_mt.fasta

# Extract from Outgroup
bedtools getfasta -fi "$OUTGROUP_REF" -bed ../work_blast/coords_outgroup_all.bed -nameOnly > outgroup_all_mt.fasta
awk '/^>/ {gene=substr($1,2); file="Outgroup_" gene ".fasta"; print ">Outgroup" > file; next} {print >> file}' outgroup_all_mt.fasta

# ==============================================================================
# 6. SAMPLE BAM CONSENSUS EXTRACTION
# ==============================================================================
echo "=== EXTRACTING CONSENSUS GENES FROM BAM FILES ==="
for BAM in "${BAM_DIR}"/*.bam; do
    SAMPLE=$(basename "${BAM}" _final_markdup.bam)
    echo "Processing $SAMPLE..."
    
    while read -r contig start end gene; do
        REGION="${contig}:${start}-${end}"
        
        # Generates a consensus FASTA for this specific gene
        samtools view -b "${BAM}" "${REGION}" | \
          bcftools mpileup -Ou -f "${BAM_MAPPING_REF}" - 2>/dev/null | \
          bcftools call -c - 2>/dev/null | \
          vcfutils.pl vcf2fq 2>/dev/null | \
          seqtk seq -A > "${SAMPLE}_${gene}.fasta"
          
        # Rename header to the sample name
        sed -i "1s/.*/\>${SAMPLE}/" "${SAMPLE}_${gene}.fasta"
    done < ../work_blast/coords_mapping_all.bed
done

# Compile all sequences by gene
for GENE in "${TARGET_MT_GENES[@]}"; do
    cat Parent1_Ref_${GENE}.fasta Outgroup_${GENE}.fasta *_${GENE}.fasta > "${BASE_OUT}/2_extracted_mt/${GENE}_compiled.fasta"
done

# ==============================================================================
# 7. ALIGNMENT & SUPERMATRIX CONCATENATION
# ==============================================================================
echo "=== ALIGNING MT GENES LOCUS-BY-LOCUS ==="
cd "${BASE_OUT}/3_alignments"

for GENE in "${TARGET_MT_GENES[@]}"; do
    echo "Aligning $GENE..."
    mafft --auto --quiet "${BASE_OUT}/2_extracted_mt/${GENE}_compiled.fasta" > "${GENE}_aligned.fasta"
done

echo "=== BUILDING CONCATENATED SUPERMATRIX ==="
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
            if not line: continue
            if line.startswith(">"):
                if header: super_dict[header] += "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header: super_dict[header] += "".join(seq)

# Pad missing data to ensure a perfectly rectangular alignment matrix
max_len = max(len(s) for s in super_dict.values())
for k, v in super_dict.items():
    if len(v) < max_len: super_dict[k] = v + ("-" * (max_len - len(v)))

with open("RAxML_MT_Supermatrix.fasta", "w") as out:
    for header, sequence in super_dict.items():
        out.write(f">{header}\n{sequence}\n")
EOF

# ==============================================================================
# 8. RAxML TREE INFERENCE
# ==============================================================================
echo "=== RUNNING RAxML INFERENCE ==="
cd "${BASE_OUT}/4_raxml_tree"
cp ../3_alignments/RAxML_MT_Supermatrix.fasta .

# Convert FASTA to PHYLIP strictly required by older RAxML binaries
python3 - << 'EOF'
fasta_in = "RAxML_MT_Supermatrix.fasta"
phy_out = "RAxML_MT_Supermatrix.phy"

with open(fasta_in, "r") as f:
    raw_data = f.read().strip().split(">")[1:]

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

# Run RAxML Rapid Bootstrapping (-N 100)
raxmlHPC-PTHREADS-SSE3 -T 8 \
    -f a \
    -x 12345 -p 12345 \
    -N 100 \
    -m GTRGAMMA \
    -o Outgroup \
    -s RAxML_MT_Supermatrix.phy \
    -n mt_markers

echo "=== PART 3 PIPELINE FINISHED SUCCESSFULLY at $(date) ==="

#!/bin/bash
#PBS -N Mito_Phylogeny
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=100gb
#PBS -l walltime=12:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# PHYLOGENOMIC PIPELINE MITOCHONDRIAL PHYLOGENY
# Automatically fetches mitochondrial baits, locates MT coordinates 
# in the reference genomes, extracts consensus MT sequences from sample BAM files, 
# aligns them with MAFFT, and infers a Maximum Likelihood tree using RAxML.
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
BASE_OUT="${PROJECT_DIR}/output/part3_mitochondrial"

# EXPLICIT REFERENCE DEFINITIONS
OUTGROUP_REF="${GENOME_DIR}/basvit_genome.fna"
PARENT1_REF="${GENOME_DIR}/lser_ref.fna"

# CRITICAL: This MUST be the genome the BAMs were mapped against
BAM_MAPPING_REF="${GENOME_DIR}/ljul_genome.fasta"  

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
# 3. BAIT FETCHING & BLAST SEARCH
# ==============================================================================
echo "=== FETCHING MT BAIT & CREATING BLAST DATABASES ==="
BAIT_FASTA="${BASE_OUT}/1_baits/mito_bait.fasta"

if [ ! -f "$BAIT_FASTA" ]; then
    esearch -db nucleotide -query "Anolis carolinensis mitochondrion, complete genome" | \
    efetch -format fasta > "$BAIT_FASTA"
fi

# Create BLAST DBs for the all relevant genomes
makeblastdb -in "$BAM_MAPPING_REF" -dbtype nucl -out mapping_ref_db >/dev/null
makeblastdb -in "$PARENT1_REF" -dbtype nucl -out parent_db >/dev/null
makeblastdb -in "$OUTGROUP_REF" -dbtype nucl -out outgroup_db >/dev/null

echo "=== LOCATING MT COORDINATES ==="
# 3A. Find MT in the BAM mapping reference (ONT Genome)
blastn -task blastn -query "$BAIT_FASTA" -db mapping_ref_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_mapping.txt
awk '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, $4; else print $1, $3-1, $2, $4}' blast_mapping.txt > coords_mapping.bed

# 3B. Find MT in the Pure Parent 1 (L. serratus)
blastn -task blastn -query "$BAIT_FASTA" -db parent_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_parent.txt
awk '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, $4; else print $1, $3-1, $2, $4}' blast_parent.txt > coords_parent.bed

# 3C. Find MT in the Outgroup reference
blastn -task blastn -query "$BAIT_FASTA" -db outgroup_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_outgroup.txt
awk '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, $4; else print $1, $3-1, $2, $4}' blast_outgroup.txt > coords_outgroup.bed

# ==============================================================================
# 4. REFERENCE MT EXTRACTION
# ==============================================================================
echo "=== EXTRACTING PURE OUTGROUP & PARENT REFERENCE MT SEQUENCES ==="
cd ../work_bams

bedtools getfasta -fi "$PARENT1_REF" -bed ../work_blast/coords_parent.bed -nameOnly | \
    sed 's/^>.*/>Parent1_Ref/' > Parent1_Ref_mt.fasta

bedtools getfasta -fi "$OUTGROUP_REF" -bed ../work_blast/coords_outgroup.bed -nameOnly | \
    sed 's/^>.*/>Outgroup/' > Outgroup_mt.fasta

# ==============================================================================
# 5. SAMPLE BAM CONSENSUS EXTRACTION
# ==============================================================================
echo "=== EXTRACTING CONSENSUS MT FROM BAM FILES ==="
for BAM in "${BAM_DIR}"/*.bam; do
    SAMPLE=$(basename "${BAM}" _final_markdup.bam)
    echo "Processing $SAMPLE..."
    
    while read -r contig start end gene; do
        REGION="${contig}:${start}-${end}"
        
        # Generates a consensus FASTA for this specific region using the ONT reference
        samtools view -b "${BAM}" "${REGION}" | \
          bcftools mpileup -Ou -f "${BAM_MAPPING_REF}" - 2>/dev/null | \
          bcftools call -c - 2>/dev/null | \
          vcfutils.pl vcf2fq 2>/dev/null | \
          seqtk seq -A > "${SAMPLE}_${gene}.fasta"
          
        # Rename header to the sample name
        sed -i "1s/.*/\>${SAMPLE}/" "${SAMPLE}_${gene}.fasta"
    done < ../work_blast/coords_mapping.bed
done

# Compile all sequences by bait target
grep "^>" "$BAIT_FASTA" | awk '{print $1}' | sed 's/>//' > gene_list.txt
while read -r gene; do
    cat Parent1_Ref_mt.fasta Outgroup_mt.fasta *_${gene}.fasta > "${BASE_OUT}/2_extracted_mt/${gene}_compiled.fasta"
done < gene_list.txt

# ==============================================================================
# 6. ALIGNMENT & SUPERMATRIX CONCATENATION
# ==============================================================================
echo "=== ALIGNING MT SEQUENCES ==="
cd "${BASE_OUT}/3_alignments"

while read -r gene; do
    echo "Aligning $gene..."
    mafft --auto --quiet "${BASE_OUT}/2_extracted_mt/${gene}_compiled.fasta" > "${gene}_aligned.fasta"
done < ../../work_bams/gene_list.txt

echo "=== BUILDING SUPERMATRIX ==="
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

# Pad missing data to ensure perfectly rectangular matrix
max_len = max(len(s) for s in super_dict.values())
for k, v in super_dict.items():
    if len(v) < max_len: super_dict[k] = v + ("-" * (max_len - len(v)))

with open("RAxML_MT_Supermatrix.fasta", "w") as out:
    for header, sequence in super_dict.items():
        out.write(f">{header}\n{sequence}\n")
EOF

# ==============================================================================
# 7. RAxML TREE INFERENCE
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
# Outgroup name (-o) explicitly matches the extracted reference
raxmlHPC-PTHREADS-SSE3 -T 8 \
    -f a \
    -x 12345 -p 12345 \
    -N 100 \
    -m GTRGAMMA \
    -o Outgroup \
    -s RAxML_MT_Supermatrix.phy \
    -n Mito_Topology_Final

echo "=== PART 3 PIPELINE FINISHED SUCCESSFULLY at $(date) ==="

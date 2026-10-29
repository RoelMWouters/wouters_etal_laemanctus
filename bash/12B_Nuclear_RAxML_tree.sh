#!/bin/bash
#PBS -N master_phylogenomics_pipeline
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=50gb
#PBS -l walltime=15:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# MASTER PHYLOGENOMICS PIPELINE 
# Description: Maps nuclear baits against reference genomes and phased hybrid 
# assemblies, extracts allelic coordinates, aligns with MAFFT, and builds a 
# concatenated supermatrix for RAxML phylogenetic inference.
# ==============================================================================

set -euo pipefail

# ==============================================================================
# 0. GLOBAL DIRECTORIES & MODULES
# ==============================================================================
# USER CONFIGURABLE: Set this to your project's root directory
WORK_DIR="/path/to/your/project/directory"

GENOME_DIR="$WORK_DIR/genome_dir"
BAM_DIR="$WORK_DIR/bam_output"
BAIT_DIR="$WORK_DIR/output_dir/bait_sequences"

BLAST_PHASED="$WORK_DIR/output_dir/blast_results"
BLAST_REFS="$WORK_DIR/output_dir/blast_results_refs"
MARKERS_DIR="$WORK_DIR/output_dir/final_markers"
ALIGN_DIR="$WORK_DIR/output_dir/alignments"
SUPERMATRIX_DIR="$WORK_DIR/output_dir/supermatrix"

mkdir -p "$BLAST_PHASED" "$BLAST_REFS" "$MARKERS_DIR" "$ALIGN_DIR" "$SUPERMATRIX_DIR"

cd "$SCRATCHDIR" || { echo "No SCRATCHDIR available – aborting"; exit 1; }

# Load required modules (Adjust these names to match your cluster's module system)
module purge
module add blast edirect bcftools samtools python mafft raxml

# ==============================================================================
# 1. BLAST REFERENCES & HYBRID BACKBONE
# ==============================================================================
echo "=== STEP 1: BLASTING NUCLEAR BAITS ==="

# Build BLAST databases if they don't exist
if [ ! -f "$GENOME_DIR/lser_blastdb.nhr" ]; then
    makeblastdb -in "$GENOME_DIR/lser_ref.fna" -dbtype nucl -out "$GENOME_DIR/lser_blastdb"
fi
if [ ! -f "$GENOME_DIR/basvit_blastdb.nhr" ]; then
    makeblastdb -in "$GENOME_DIR/basvit_genome.fna" -dbtype nucl -out "$GENOME_DIR/basvit_blastdb"
fi
# Build database for the hybrid mapping reference
if [ ! -f "$GENOME_DIR/ljul_blastdb.nhr" ]; then
    makeblastdb -in "$GENOME_DIR/ljul_genome.fasta" -dbtype nucl -out "$GENOME_DIR/ljul_blastdb"
fi

for BAIT in "$BAIT_DIR"/*_anolis_bait.fasta; do
    GENE_NAME=$(basename "$BAIT" _anolis_bait.fasta)
    echo " -> Mapping $GENE_NAME..."
    
    # 1A. Map against pure unphased references (Outgroup & Parent)
    for REF in "lser" "basvit"; do
        tblastn -query "$BAIT" -db "$GENOME_DIR/${REF}_blastdb" \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
                -out "$BLAST_REFS/${REF}_${GENE_NAME}_hits.txt"
    done
    
    # 1B. Map against hybrid reference to generate phased BAM coordinates
    tblastn -query "$BAIT" -db "$GENOME_DIR/ljul_blastdb" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -out "$BLAST_PHASED/${GENE_NAME}_phased_hits.txt"
done

# ==============================================================================
# 2. MASTER EXTRACTOR
# ==============================================================================
echo "=== STEP 2: EXTRACTING MARKERS ==="
export GENOME_DIR OUT_DIR="$MARKERS_DIR" BLAST_REFS BLAST_PHASED

cat << 'EOF' > calculate_coordinates.py
import os, glob

def get_gene_bounds(file_path, top_n_alleles=1):
    hits = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12: continue
            contig = parts[1]
            sstart, send = int(parts[8]), int(parts[9])
            evalue, bitscore = float(parts[10]), float(parts[11])
            if evalue < 1e-10:
                hits.append({'contig': contig, 'min_c': min(sstart, send), 'max_c': max(sstart, send), 'bitscore': bitscore})
    
    if not hits: return []
    contig_hits = {}
    for h in hits:
        contig_hits.setdefault(h['contig'], []).append(h)

    results = []
    sorted_contigs = sorted(contig_hits.keys(), key=lambda c: max(h['bitscore'] for h in contig_hits[c]), reverse=True)
    
    for c in sorted_contigs[:top_n_alleles]:
        top_hit = max(contig_hits[c], key=lambda h: h['bitscore'])
        golden_center = (top_hit['min_c'] + top_hit['max_c']) / 2
        valid_exons = [h for h in contig_hits[c] if abs(((h['min_c'] + h['max_c']) / 2) - golden_center) <= 100000]
        min_val = max(1, min(h['min_c'] for h in valid_exons) - 500)
        max_val = max(h['max_c'] for h in valid_exons) + 500
        results.append((c, min_val, max_val))
    return results

# Write Reference Extraction Commands
with open("extract_refs.sh", "w") as f_out:
    for hit_file in glob.glob(os.path.join(os.environ['BLAST_REFS'], "*_hits.txt")):
        filename = os.path.basename(hit_file)
        if filename.startswith("lser_"):
            genome_fasta, prefix = os.path.join(os.environ['GENOME_DIR'], "lser_ref.fna"), "lser"
        elif filename.startswith("basvit_"):
            genome_fasta, prefix = os.path.join(os.environ['GENOME_DIR'], "basvit_genome.fna"), "basvit"
        else: continue
            
        gene_name = filename.replace(f"{prefix}_", "").replace("_hits.txt", "")
        bounds = get_gene_bounds(hit_file, top_n_alleles=1)
        if bounds:
            contig, start, end = bounds[0]
            out_name = os.path.join(os.environ['OUT_DIR'], f"{prefix}_{gene_name}.fasta")
            f_out.write(f"samtools faidx {genome_fasta} {contig}:{start}-{end} | sed 's/>.*/>{prefix}/' > {out_name}\n")

# Write BAM Regions (from Hybrid Mapping)
with open("bam_regions.txt", "w") as f_out:
    for hit_file in glob.glob(os.path.join(os.environ['BLAST_PHASED'], "*_phased_hits.txt")):
        gene_name = os.path.basename(hit_file).replace("_phased_hits.txt", "")
        for idx, (contig, start, end) in enumerate(get_gene_bounds(hit_file, top_n_alleles=2)):
            f_out.write(f"{contig}:{start}-{end}|{gene_name}_Allele{idx+1}\n")
EOF

python3 calculate_coordinates.py
bash extract_refs.sh

REF_FASTA="$GENOME_DIR/ljul_genome.fasta"

# USER CONFIGURABLE: List your BAM file prefixes here
BAMS=(
    "C01_final_markdup.bam" "C02_final_markdup.bam" "C03_final_markdup.bam"
    "C04_final_markdup.bam" "C05_final_markdup.bam" "C06_final_markdup.bam"
    "F01_final_markdup.bam" "F03_final_markdup.bam"
    "C09_final_markdup.bam" "C10_final_markdup.bam"
)

if [ -s bam_regions.txt ]; then
    cut -d'|' -f2 bam_regions.txt | cut -d'_' -f1 | sort | uniq > unique_genes.txt
    for BAM in "${BAMS[@]}"; do
        BAM_PATH="$BAM_DIR/$BAM"
        if [ ! -f "$BAM_PATH" ]; then continue; fi
        SAMPLE_NAME=$(basename "$BAM" _final_markdup.bam)
        
        while read GENE; do
            COORD_A1=$(grep "${GENE}_Allele1" bam_regions.txt | cut -d'|' -f1 || true)
            COORD_A2=$(grep "${GENE}_Allele2" bam_regions.txt | cut -d'|' -f1 || true)
            
            if [ -n "$COORD_A1" ]; then
                samtools faidx "$REF_FASTA" "$COORD_A1" > "tmp_ref.fasta"
                bcftools mpileup -Ou -f "$REF_FASTA" -r "$COORD_A1" "$BAM_PATH" 2>/dev/null | bcftools call -mv -Oz -o "tmp_vars.vcf.gz" 2>/dev/null
                tabix -f -p vcf "tmp_vars.vcf.gz" 2>/dev/null
                cat "tmp_ref.fasta" | bcftools consensus "tmp_vars.vcf.gz" 2>/dev/null | sed "s/>.*/>${SAMPLE_NAME}_Allele1/" > "$MARKERS_DIR/${SAMPLE_NAME}_${GENE}_Allele1.fasta"
            fi
            if [ -n "$COORD_A2" ]; then
                samtools faidx "$REF_FASTA" "$COORD_A2" > "tmp_ref.fasta"
                bcftools mpileup -Ou -f "$REF_FASTA" -r "$COORD_A2" "$BAM_PATH" 2>/dev/null | bcftools call -mv -Oz -o "tmp_vars.vcf.gz" 2>/dev/null
                tabix -f -p vcf "tmp_vars.vcf.gz" 2>/dev/null
                cat "tmp_ref.fasta" | bcftools consensus "tmp_vars.vcf.gz" 2>/dev/null | sed "s/>.*/>${SAMPLE_NAME}_Allele2/" > "$MARKERS_DIR/${SAMPLE_NAME}_${GENE}_Allele2.fasta"
            fi
        done < unique_genes.txt
    done
fi

# ==============================================================================
# 3. MAFFT ALIGNMENTS
# ==============================================================================
echo "=== STEP 3: MAFFT ALIGNMENTS ==="
GENES=("CHRM4" "FGB" "MC1R" "PDC" "RELN")
for GENE in "${GENES[@]}"; do
    cat "$MARKERS_DIR"/*_"${GENE}"*.fasta > "$ALIGN_DIR/${GENE}_unaligned.fasta" 2>/dev/null || true
    if [ -s "$ALIGN_DIR/${GENE}_unaligned.fasta" ]; then
        mafft --adjustdirection --auto "$ALIGN_DIR/${GENE}_unaligned.fasta" > "$ALIGN_DIR/${GENE}_aligned.fasta" 2> /dev/null
    fi
done

# ==============================================================================
# 4. SUPERMATRIX CONCATENATION 
# ==============================================================================
echo "=== STEP 4: CONCATENATING SUPERMATRIX ==="
cat << 'EOF' > concat_supermatrix.py
import os, glob
alignments = glob.glob(os.environ['ALIGN_DIR'] + "/*_aligned.fasta")
taxa_seqs = {}

for aln in alignments:
    with open(aln, 'r') as f:
        blocks = f.read().split(">")[1:]
        seq_len = len("".join(blocks[0].strip().split("\n")[1:])) if blocks else 0
        
        current_taxa = []
        for block in blocks:
            lines = block.strip().split("\n")
            taxon = lines[0].strip().replace("_R_", "") # Clean MAFFT reverse indicators
            seq = "".join(lines[1:]).strip()
            taxa_seqs.setdefault(taxon, []).append(seq)
            current_taxa.append(taxon)
            
        # Pad missing taxa with gaps to keep matrix rectangular
        for t in taxa_seqs.keys():
            if t not in current_taxa:
                taxa_seqs[t].append("-" * seq_len)

with open(os.path.join(os.environ['SUPERMATRIX_DIR'], "RAxML_Supermatrix.fasta"), "w") as f_out:
    for taxon, seq_list in taxa_seqs.items():
        f_out.write(f">{taxon}\n{''.join(seq_list)}\n")
EOF
export ALIGN_DIR SUPERMATRIX_DIR
python3 concat_supermatrix.py

# ==============================================================================
# 5. RAXML PHYLOGENETICS
# ==============================================================================
echo "=== STEP 5: RAXML PHYLOGENETICS ==="
cd "$SUPERMATRIX_DIR"

cat << 'EOF' > fasta_to_phylip.py
with open("RAxML_Supermatrix.fasta", "r") as f:
    raw_data = f.read().split(">")[1:]
seqs = [(b.strip().split("\n")[0].strip(), "".join(b.strip().split("\n")[1:]).strip()) for b in raw_data]

with open("RAxML_Supermatrix.phy", "w") as f:
    f.write(f"{len(seqs)} {len(seqs[0][1])}\n")
    for name, seq in seqs:
        f.write(f"{name:<30} {seq}\n")
EOF
python3 fasta_to_phylip.py

raxmlHPC-PTHREADS-SSE3 -T 8 -f a -x 12345 -p 12345 -N 1000 -m GTRGAMMA -o basvit -s RAxML_Supermatrix.phy -n F1_Hybrid_Topology
echo "=== PIPELINE COMPLETE ==="

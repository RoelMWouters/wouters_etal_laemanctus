#!/bin/bash
#PBS -N Extract_Nuclear_Markers
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=50gb
#PBS -l walltime=06:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 12: NUCLEAR MARKER CONSENSUS EXTRACTION
# ==============================================================================
# Description:
# Uses a hybrid Python/Bash approach to parse phased BLAST results, calculate 
# accurate genomic coordinates for top allelic hits (managing exon bounds), 
# and extracts allelic variants per sample from BAM files using BCFtools consensus.
#
# Dependencies: python3, samtools, bcftools, tabix
# ==============================================================================

WORK_DIR="/path/to/cluster/project_dir"
GENOME_DIR="$WORK_DIR/genome_dir"
BAM_DIR="$WORK_DIR/bam_output"

BLAST_PHASED="$WORK_DIR/output_dir/blast_results"
BLAST_REFS="$WORK_DIR/output_dir/blast_results_refs"
OUT_DIR="$WORK_DIR/output_dir/final_markers"

mkdir -p "$OUT_DIR"

module purge
module add bcftools samtools python
cd "$SCRATCHDIR"

# --- Hybrid Python Parsing Script ---
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
        c = h['contig']
        if c not in contig_hits: contig_hits[c] = []
        contig_hits[c].append(h)

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

with open("extract_refs.sh", "w") as f_out:
    for hit_file in glob.glob(os.path.join(os.environ['BLAST_REFS'], "*_hits.txt")):
        filename = os.path.basename(hit_file)
        if filename.startswith("Subgenome_A"):
            genome_fasta = os.path.join(os.environ['GENOME_DIR'], "Subgenome_A.fasta")
            gene_name = filename.replace("Subgenome_A", "").replace("_hits.txt", "")
            prefix = "refA"
        else: continue
            
        bounds = get_gene_bounds(hit_file, top_n_alleles=1)
        if bounds:
            contig, start, end = bounds[0]
            out_name = os.path.join(os.environ['OUT_DIR'], f"{prefix}_{gene_name}.fasta")
            f_out.write(f"samtools faidx {genome_fasta} {contig}:{start}-{end} > {out_name}\n")

with open("bam_regions.txt", "w") as f_out:
    for hit_file in glob.glob(os.path.join(os.environ['BLAST_PHASED'], "*_phased_hits.txt")):
        filename = os.path.basename(hit_file)
        gene_name = filename.replace("_phased_hits.txt", "")
        bounds = get_gene_bounds(hit_file, top_n_alleles=2)
        for idx, (contig, start, end) in enumerate(bounds):
            f_out.write(f"{contig}:{start}-{end}|{gene_name}_Allele{idx+1}\n")
EOF

export GENOME_DIR OUT_DIR BLAST_REFS BLAST_PHASED
python3 calculate_coordinates.py

chmod +x extract_refs.sh
./extract_refs.sh

REF_FASTA="$GENOME_DIR/genome.fasta"
BAMS=("S01_final_markdup.bam" "S02_final_markdup.bam")

cut -d'|' -f2 bam_regions.txt | cut -d'_' -f1 | sort | uniq > unique_genes.txt

for BAM in "${BAMS[@]}"; do
    BAM_PATH="$BAM_DIR/$BAM"
    if [ ! -f "$BAM_PATH" ]; then continue; fi
    SAMPLE_NAME=$(basename "$BAM" _final_markdup.bam)
    
    while read GENE; do
        COORD_A1=$(grep "${GENE}_Allele1" bam_regions.txt | cut -d'|' -f1)
        COORD_A2=$(grep "${GENE}_Allele2" bam_regions.txt | cut -d'|' -f1)
        
        OUT_A1="$OUT_DIR/${SAMPLE_NAME}_${GENE}_Allele1.fasta"
        OUT_A2="$OUT_DIR/${SAMPLE_NAME}_${GENE}_Allele2.fasta"

        READS_A1=$(samtools view -c "$BAM_PATH" "$COORD_A1")
        READS_A2=$(samtools view -c "$BAM_PATH" "$COORD_A2")
        
        EXTRACT_A1=true; EXTRACT_A2=true
        
        if [ "$READS_A1" -lt 10 ]; then EXTRACT_A1=false; fi
        if [ "$READS_A2" -lt 10 ]; then EXTRACT_A2=false; fi

        if [ "$EXTRACT_A1" = true ]; then
            samtools faidx "$REF_FASTA" "$COORD_A1" > "$SCRATCHDIR/tmp_ref.fasta"
            bcftools mpileup -Ou -f "$REF_FASTA" -r "$COORD_A1" "$BAM_PATH" 2>/dev/null | \
            bcftools call -mv -Oz -o "$SCRATCHDIR/tmp_vars.vcf.gz" 2>/dev/null
            tabix -f -p vcf "$SCRATCHDIR/tmp_vars.vcf.gz" 2>/dev/null
            cat "$SCRATCHDIR/tmp_ref.fasta" | bcftools consensus "$SCRATCHDIR/tmp_vars.vcf.gz" 2>/dev/null | \
            sed "s/>.*/>${SAMPLE_NAME}_${GENE}_Allele1/" > "$OUT_A1"
        else
            rm -f "$OUT_A1"
        fi
        
        if [ "$EXTRACT_A2" = true ]; then
            samtools faidx "$REF_FASTA" "$COORD_A2" > "$SCRATCHDIR/tmp_ref.fasta"
            bcftools mpileup -Ou -f "$REF_FASTA" -r "$COORD_A2" "$BAM_PATH" 2>/dev/null | \
            bcftools call -mv -Oz -o "$SCRATCHDIR/tmp_vars.vcf.gz" 2>/dev/null
            tabix -f -p vcf "$SCRATCHDIR/tmp_vars.vcf.gz" 2>/dev/null
            cat "$SCRATCHDIR/tmp_ref.fasta" | bcftools consensus "$SCRATCHDIR/tmp_vars.vcf.gz" 2>/dev/null | \
            sed "s/>.*/>${SAMPLE_NAME}_${GENE}_Allele2/" > "$OUT_A2"
        else
            rm -f "$OUT_A2"
        fi
    done < unique_genes.txt
done

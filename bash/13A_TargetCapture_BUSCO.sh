#!/bin/bash
#PBS -N Part1_TargetCapture
#PBS -l select=1:ncpus=32:mem=256gb:scratch_local=250gb
#PBS -l walltime=24:00:00
#PBS -M user@institute.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# PHYLOGENOMIC PIPELINE - PART 1: LOCUS ASSEMBLY & EXTRACTION
# Uses BUSCO to identify orthologs. Extracts these regions from 
# Illumina BAM files, performs local SPAdes assembly, and uses Exonerate to 
# pull pure CDS sequences. Exonerate is run on the parental 
# and outgroup reference genomes as well to ensure perfectly matching exon boundaries.
# ==============================================================================

set -euo pipefail
export TMPDIR=$SCRATCHDIR
cd "$SCRATCHDIR" || exit 1

# ==============================================================================
# 1. CONFIGURATION & LOGICAL PATHS
# ==============================================================================
PROJECT_DIR="/path/to/project"
GENOME_DIR="${PROJECT_DIR}/genomes"
BAM_DIR="${PROJECT_DIR}/bam_output"

OUTGROUP="${GENOME_DIR}/basvit_genome.fna"
PARENT1_REF="${GENOME_DIR}/lser_ref.fna"
HYBRID_ONT="${GENOME_DIR}/ljul_genome.fasta"

ILLUMINA_SAMPLES=("C10" "C09" "C01" "C02" "C03" "C04" "C05" "C06" "F01" "F03")

BASE_OUT="${PROJECT_DIR}/output/part1"
BUSCO_BACKUP="${BASE_OUT}/0_busco_checkpoints"
SEQ_OUT="${BASE_OUT}/1_extracted_sequences"
META_OUT="${BASE_OUT}/metadata"

mkdir -p "$BUSCO_BACKUP" "$SEQ_OUT/exonerated_cds" "$SEQ_OUT/ref_seqs" "$META_OUT"
SUMMARY_FILE="${META_OUT}/pipeline_summary_spades.txt"

echo "=== PART 1: BAM-TO-SPADES + EXONERATE TARGET CAPTURE ===" > "$SUMMARY_FILE"
echo "Pipeline started at: $(date)" >> "$SUMMARY_FILE"

# ==============================================================================
# 2. ISOLATED ENVIRONMENT SETUP
# ==============================================================================
module purge
module add mambaforge

unset PYTHONPATH
unset PYTHONHOME

export CONDA_PKGS_DIRS="$SCRATCHDIR/conda_pkgs"
export CONDA_ENVS_DIRS="$SCRATCHDIR/conda_envs"
mkdir -p "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS"

mamba create -n bam_spades_env -c conda-forge -c bioconda \
    python=3.11 spades=3.15.5 samtools=1.19 busco=5.7.1 biopython=1.80 seqtk parallel exonerate -y
eval "$(conda shell.bash hook)"
conda activate bam_spades_env

# ==============================================================================
# 3. INDEPENDENT BUSCO RUNS
# ==============================================================================
echo "=== CHECKING FOR EXISTING BUSCO RUNS ==="

if [ -d "${BUSCO_BACKUP}/BUSCO_Lser" ]; then
    cp -r "${BUSCO_BACKUP}/BUSCO_Lser" .
else
    seqtk seq -L 10000 "$PARENT1_REF" > parent1_filtered.fna
    busco -i parent1_filtered.fna -o BUSCO_Lser -l sauropsida_odb10 -m genome -c 16 -f --metaeuk
    cp -r BUSCO_Lser "$BUSCO_BACKUP/"
fi

if [ -d "${BUSCO_BACKUP}/BUSCO_Basvit" ]; then
    cp -r "${BUSCO_BACKUP}/BUSCO_Basvit" .
else
    seqtk seq -L 10000 "$OUTGROUP" > outgroup_filtered.fna
    busco -i outgroup_filtered.fna -o BUSCO_Basvit -l sauropsida_odb10 -m genome -c 16 -f --metaeuk
    cp -r BUSCO_Basvit "$BUSCO_BACKUP/"
fi

if [ -d "${BUSCO_BACKUP}/BUSCO_Hybrid_ONT" ]; then
    cp -r "${BUSCO_BACKUP}/BUSCO_Hybrid_ONT" .
else
    seqtk seq -L 10000 "$HYBRID_ONT" > hybrid_ont_filtered.fna
    busco -i hybrid_ont_filtered.fna -o BUSCO_Hybrid_ONT -l sauropsida_odb10 -m genome -c 16 -f --metaeuk
    cp -r BUSCO_Hybrid_ONT "$BUSCO_BACKUP/"
fi

# ==============================================================================
# 4. GENE SELECTION (FILTERING)
# ==============================================================================
echo "=== SELECTING GENES ==="
awk '$2=="Complete" {print $1}' BUSCO_Lser/run_sauropsida_odb10/full_table.tsv | sort > lser_sc.txt
awk '$2=="Complete" {print $1}' BUSCO_Basvit/run_sauropsida_odb10/full_table.tsv | sort > basvit_sc.txt
awk '$2=="Duplicated" {print $1}' BUSCO_Hybrid_ONT/run_sauropsida_odb10/full_table.tsv | sort | uniq -c | awk '$1==2 {print $2}' > ont_dup.txt

comm -12 lser_sc.txt basvit_sc.txt > common_single.txt
comm -12 common_single.txt ont_dup.txt > valid_targets.txt

TOTAL_VALID=$(wc -l < valid_targets.txt)
if [ "$TOTAL_VALID" -gt 2000 ]; then
    shuf -n 2000 valid_targets.txt > target_busco_ids.txt
else
    cp valid_targets.txt target_busco_ids.txt
fi

echo "Genes Selected for BAM Extraction: $(wc -l < target_busco_ids.txt)" >> "$SUMMARY_FILE"

# ==============================================================================
# 5. EXTRACT HYBRID GENOME COORDINATES AND BUILD HYBRID TASK LIST
# ==============================================================================
awk 'NR==FNR{a[$1]; next} $1 in a && $2=="Duplicated" {print $1, $3, $4, $5}' \
    target_busco_ids.txt BUSCO_Hybrid_ONT/run_sauropsida_odb10/full_table.tsv | sort > ont_coords_raw.txt

awk '{
    if ($1 == prev_id) {copy="copy2"} else {copy="copy1"}
    prev_id=$1;
    print $1, $2, $3, $4, copy
}' ont_coords_raw.txt > ont_coords.txt

> parallel_tasks.txt
for SAMPLE in "${ILLUMINA_SAMPLES[@]}"; do
    awk -v s="$SAMPLE" '{print s, $1, $2, $3, $4, $5}' ont_coords.txt >> parallel_tasks.txt
done

# ==============================================================================
# 6. PARALLEL HYBRID SPAdes & EXONERATE
# ==============================================================================
process_gene_assembly() {
    SAMPLE=$1; BUSCO=$2; CHROM=$3; START=$4; END=$5; COPY=$6
    BAM_DIR=$7; OUT_DIR=$8; SCRATCH=$9
    
    PAD_START=$(( START > 1000 ? START - 1000 : 1 ))
    PAD_END=$(( END + 1000 ))
    REGION="${CHROM}:${PAD_START}-${PAD_END}"
    
    WORK_DIR="${SCRATCH}/work_${SAMPLE}_${BUSCO}_${COPY}"
    mkdir -p "$WORK_DIR"
    
    samtools view -b "${BAM_DIR}/${SAMPLE}_final_markdup.bam" "$REGION" > "${WORK_DIR}/region.bam" 2>/dev/null
    samtools fastq -F 256 "${WORK_DIR}/region.bam" > "${WORK_DIR}/reads.fq" 2>/dev/null
    
    if [ $(wc -l < "${WORK_DIR}/reads.fq") -gt 80 ]; then
        spades.py --careful --only-assembler -s "${WORK_DIR}/reads.fq" -o "${WORK_DIR}/spades" -m 4 -t 2 >/dev/null 2>&1
        
        if [ -f "${WORK_DIR}/spades/contigs.fasta" ]; then
            REF_SEQ="${SCRATCH}/BUSCO_Lser/run_sauropsida_odb10/busco_sequences/single_copy_busco_sequences/${BUSCO}.fna"
            
            if [ -f "$REF_SEQ" ]; then
                exonerate --model coding2genome --query "$REF_SEQ" --target "${WORK_DIR}/spades/contigs.fasta" \
                    --bestn 1 --showvulgar no --showalignment no \
                    --ryo "\n@@@START@@@\n%tcs\n@@@END@@@\n" > "${WORK_DIR}/exo.out" 2>/dev/null
                
                awk -v name=">${SAMPLE}_${BUSCO}_${COPY}" '
                    /^@@@START@@@$/ {flag=1; print name; next}
                    /^@@@END@@@$/ {flag=0; next}
                    flag { gsub(/-/, ""); if (length($0) > 0) print toupper($0) }
                ' "${WORK_DIR}/exo.out" > "${OUT_DIR}/${SAMPLE}_${BUSCO}_${COPY}.fasta"
                
                if [ -f "${OUT_DIR}/${SAMPLE}_${BUSCO}_${COPY}.fasta" ] && [ $(wc -l < "${OUT_DIR}/${SAMPLE}_${BUSCO}_${COPY}.fasta") -lt 2 ]; then
                    rm -f "${OUT_DIR}/${SAMPLE}_${BUSCO}_${COPY}.fasta"
                fi
            fi
        fi
    fi
    rm -rf "$WORK_DIR"
}
export -f process_gene_assembly

# ==============================================================================
# 7. PARALLEL REFERENCE EXONERATE
# ==============================================================================
process_reference_exonerate() {
    TAXON=$1; BUSCO=$2; GENOME=$3; TABLE=$4; OUT_DIR=$5; SCRATCH=$6
    
    COORD_LINE=$(awk -v b="$BUSCO" '$1==b && $2=="Complete" {print $3, $4, $5; exit}' "$TABLE")
    if [ -z "$COORD_LINE" ]; then return; fi
    read -r CHROM START END <<< "$COORD_LINE"
    
    if [ "$START" -lt "$END" ]; then
        TRUE_START=$START; TRUE_END=$END
    else
        TRUE_START=$END; TRUE_END=$START
    fi
    
    PAD_START=$(( TRUE_START > 1000 ? TRUE_START - 1000 : 1 ))
    PAD_END=$(( TRUE_END + 1000 ))
    REGION="${CHROM}:${PAD_START}-${PAD_END}"
    
    WORK_DIR="${SCRATCH}/work_${TAXON}_${BUSCO}"
    mkdir -p "$WORK_DIR"
    
    samtools faidx "$GENOME" "$REGION" > "${WORK_DIR}/target.fasta" 2>/dev/null
    
    # Query is ALWAYS the Parent1 BUSCO CDS to ensure identical models
    REF_SEQ="${SCRATCH}/BUSCO_Lser/run_sauropsida_odb10/busco_sequences/single_copy_busco_sequences/${BUSCO}.fna"
    
    if [ -f "$REF_SEQ" ] && [ -s "${WORK_DIR}/target.fasta" ]; then
        exonerate --model coding2genome --query "$REF_SEQ" --target "${WORK_DIR}/target.fasta" \
            --bestn 1 --showvulgar no --showalignment no \
            --ryo "\n@@@START@@@\n%tcs\n@@@END@@@\n" > "${WORK_DIR}/exo.out" 2>/dev/null
            
        awk -v name=">${TAXON}" '
            /^@@@START@@@$/ {flag=1; print name; next}
            /^@@@END@@@$/ {flag=0; next}
            flag { gsub(/-/, ""); if (length($0) > 0) print toupper($0) }
        ' "${WORK_DIR}/exo.out" > "${OUT_DIR}/${BUSCO}_${TAXON}.fasta"
        
        if [ -f "${OUT_DIR}/${BUSCO}_${TAXON}.fasta" ] && [ $(wc -l < "${OUT_DIR}/${BUSCO}_${TAXON}.fasta") -lt 2 ]; then
            rm -f "${OUT_DIR}/${BUSCO}_${TAXON}.fasta"
        fi
    fi
    rm -rf "$WORK_DIR"
}
export -f process_reference_exonerate

# ==============================================================================
# 8. RUN BOTH PARALLEL EXTRACTIONS
# ==============================================================================
echo "=== RUNNING LOCAL ASSEMBLIES & HYBRID EXONERATE ==="
mkdir -p scratch_out scratch_ref

parallel --colsep ' ' --joblog "${META_OUT}/parallel_spades.log" -j 16 \
    "process_gene_assembly {1} {2} {3} {4} {5} {6} $BAM_DIR scratch_out $SCRATCHDIR" :::: parallel_tasks.txt

echo "=== RUNNING REFERENCE GENOME EXONERATE ==="
# Index references for faidx
samtools faidx "$PARENT1_REF"
samtools faidx "$OUTGROUP"

> ref_tasks.txt
while read -r BUSCO; do
    [ -z "$BUSCO" ] && continue
    echo "Parent1_Ref $BUSCO $PARENT1_REF BUSCO_Lser/run_sauropsida_odb10/full_table.tsv" >> ref_tasks.txt
    echo "Outgroup $BUSCO $OUTGROUP BUSCO_Basvit/run_sauropsida_odb10/full_table.tsv" >> ref_tasks.txt
done < target_busco_ids.txt

parallel --colsep ' ' --joblog "${META_OUT}/parallel_references.log" -j 16 \
    "process_reference_exonerate {1} {2} {3} {4} scratch_ref $SCRATCHDIR" :::: ref_tasks.txt

# ==============================================================================
# 9. EXPORT TO PERMANENT STORAGE
# ==============================================================================
rsync -a scratch_out/ "$SEQ_OUT/exonerated_cds/"
rsync -a scratch_ref/ "$SEQ_OUT/ref_seqs/"
cp target_busco_ids.txt parallel_tasks.txt ref_tasks.txt "${META_OUT}/"

echo "Total Hybrid SPAdes Contigs Exonerated: $(ls scratch_out 2>/dev/null | wc -l)" >> "$SUMMARY_FILE"
echo "Total References Exonerated: $(ls scratch_ref 2>/dev/null | wc -l)" >> "$SUMMARY_FILE"
echo "Part 1 completed successfully at: $(date)" >> "$SUMMARY_FILE"

#!/bin/bash
#PBS -N FastP_BWA_Alignment
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=100gb
#PBS -l walltime=12:00:00
#PBS -J 1-50
#PBS -j oe

# ==============================================================================
# SCRIPT 3: READ TRIMMING & GENOME ALIGNMENT (JOB ARRAY)
# ==============================================================================
# Description:
# This script processes Illumina raw reads. It trims adapters using fastp, 
# aligns the reads to a reference genome using bwa-mem, merges lanes (if any), 
# and marks PCR duplicates. It operates as a Job Array utilizing a sample list.
#
# Dependencies: fastp, bwa, samtools
# ==============================================================================

# --- Configuration Variables ---
LIST_FILE="/path/to/cluster/project_dir/alignment_list.txt"
SAMPLE_NAME=$(sed -n "${PBS_ARRAY_INDEX}p" "$LIST_FILE" | tr -d '\r')

if [ -z "$SAMPLE_NAME" ]; then
    echo "Error: Could not retrieve sample name from line $PBS_ARRAY_INDEX"
    exit 1
fi

export TMPDIR=$SCRATCHDIR
mkdir -p "$SCRATCHDIR/tmp"

INPUT_DIR="/path/to/cluster/project_dir/raw_fastqs/${SAMPLE_NAME}"
OUTPUT_DIR="/path/to/cluster/project_dir/bam_output"
GENOME_DIR="/path/to/cluster/project_dir/genome_dir"
REF_GENOME="Genome.fasta"

R1_SUFFIX="_1.fq.gz"
R2_SUFFIX="_2.fq.gz"

P5_ADAPTER=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
P7_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

THREADS_TOTAL=16
THREADS_BWA=12
THREADS_SORT=4

# --- Setup & Staging ---
set -e
module purge
module add samtools bwa fastp

mkdir -p "$SCRATCHDIR/raw" "$SCRATCHDIR/trimmed" "$SCRATCHDIR/qc_reports" "$SCRATCHDIR/bams" "$SCRATCHDIR/ref"

cp "$GENOME_DIR/$REF_GENOME"* "$SCRATCHDIR/ref/" || { echo "Failed to copy reference"; exit 1; }
REF="$SCRATCHDIR/ref/$REF_GENOME"

# Indexing Reference
if [ ! -f "${REF}.bwt" ]; then bwa index "$REF"; fi
if [ ! -f "${REF}.fai" ]; then samtools faidx "$REF"; fi
if [ ! -f "${REF%.*}.dict" ]; then samtools dict "$REF" > "${REF%.*}.dict"; fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "CRITICAL ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if ls "$INPUT_DIR"/*"$R1_SUFFIX" 1> /dev/null 2>&1; then
    cp "$INPUT_DIR"/*"$R1_SUFFIX" "$SCRATCHDIR/raw/"
    cp "$INPUT_DIR"/*"$R2_SUFFIX" "$SCRATCHDIR/raw/"
else
    echo "CRITICAL ERROR: No R1 files found in $INPUT_DIR"
    exit 1
fi

# --- Processing Loop (FastP -> BWA) ---
cd "$SCRATCHDIR/raw"
BAM_LIST=""

for R1_FILE in *"$R1_SUFFIX"; do
    R2_FILE="${R1_FILE/$R1_SUFFIX/$R2_SUFFIX}"
    if [ ! -f "$R2_FILE" ]; then continue; fi

    BASE_NAME=$(basename "$R1_FILE" "$R1_SUFFIX")

    TRIM_R1="$SCRATCHDIR/trimmed/${BASE_NAME}_trimmed_R1.fastq.gz"
    TRIM_R2="$SCRATCHDIR/trimmed/${BASE_NAME}_trimmed_R2.fastq.gz"

    fastp \
        --in1 "$R1_FILE" --in2 "$R2_FILE" \
        --out1 "$TRIM_R1" --out2 "$TRIM_R2" \
        --adapter_sequence "$P5_ADAPTER" --adapter_sequence_r2 "$P7_ADAPTER" \
        --cut_front --cut_tail \
        --length_required 25 --qualified_quality_phred 20 \
        --json "$SCRATCHDIR/qc_reports/${BASE_NAME}_fastp.json" \
        --html "$SCRATCHDIR/qc_reports/${BASE_NAME}_fastp.html" \
        --thread $THREADS_TOTAL

    TEMP_BAM="$SCRATCHDIR/bams/${BASE_NAME}.sorted.bam"
    
    bwa mem -t $THREADS_BWA -M \
        -R "@RG\tID:${BASE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:${SAMPLE_NAME}_Lib\tPU:${BASE_NAME}" \
        "$REF" "$TRIM_R1" "$TRIM_R2" | \
    samtools sort -@ $THREADS_SORT -m 2G -T "$SCRATCHDIR/tmp/${BASE_NAME}" -o "$TEMP_BAM" -

    BAM_LIST="$BAM_LIST $TEMP_BAM"
done

# --- Merge & Mark Duplicates ---
MERGED_BAM="$SCRATCHDIR/${SAMPLE_NAME}_merged.bam"
FINAL_BAM="$SCRATCHDIR/${SAMPLE_NAME}_final_markdup.bam"

samtools merge -@ $THREADS_TOTAL "$MERGED_BAM" $BAM_LIST
samtools index -@ $THREADS_TOTAL "$MERGED_BAM"

samtools collate -@ 4 -O -u "$MERGED_BAM" "$SCRATCHDIR/tmp/collate_temp" | \
samtools fixmate -@ 4 -m -u - - | \
samtools sort -@ 4 -m 2G -T "$SCRATCHDIR/tmp/sort_temp" -u - | \
samtools markdup -@ 4 - "$FINAL_BAM"

samtools index -@ $THREADS_TOTAL "$FINAL_BAM"

# --- Final Statistics & Cleanup ---
STATS_DIR="$SCRATCHDIR/final_stats"
mkdir -p "$STATS_DIR"

samtools stats -@ $THREADS_TOTAL "$FINAL_BAM" > "$STATS_DIR/${SAMPLE_NAME}_stats.txt"
samtools flagstat -@ $THREADS_TOTAL "$FINAL_BAM" > "$STATS_DIR/${SAMPLE_NAME}_flagstat.txt"
samtools coverage "$FINAL_BAM" > "$STATS_DIR/${SAMPLE_NAME}_coverage.txt"
samtools depth -Q 20 "$FINAL_BAM" > "$STATS_DIR/${SAMPLE_NAME}_depth.txt"

mkdir -p "$OUTPUT_DIR/qc_reports" "$OUTPUT_DIR/alignment_stats"

if [ -f "$REF.fai" ]; then cp "$REF.fai" "$GENOME_DIR/"; fi
if [ -f "${REF%.*}.dict" ]; then cp "${REF%.*}.dict" "$GENOME_DIR/"; fi
cp "$REF".bwt "$GENOME_DIR/" 2>/dev/null || true
cp "$REF".pac "$GENOME_DIR/" 2>/dev/null || true
cp "$REF".ann "$GENOME_DIR/" 2>/dev/null || true
cp "$REF".amb "$GENOME_DIR/" 2>/dev/null || true
cp "$REF".sa  "$GENOME_DIR/" 2>/dev/null || true

cp "$SCRATCHDIR/qc_reports/"* "$OUTPUT_DIR/qc_reports/"
cp "$FINAL_BAM" "$OUTPUT_DIR/"
cp "${FINAL_BAM}.bai" "$OUTPUT_DIR/"
cp "$STATS_DIR/"* "$OUTPUT_DIR/alignment_stats/"

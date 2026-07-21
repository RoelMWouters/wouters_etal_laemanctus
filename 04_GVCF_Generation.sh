#!/bin/bash
#PBS -N GVCF_Generation
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -J 1-50
#PBS -j oe

# ==============================================================================
# SCRIPT 4: GVCF GENERATION & COVERAGE METRICS (JOB ARRAY)
# ==============================================================================
# Description:
# Takes a BAM file and generates a Genomic VCF (GVCF) using GATK HaplotypeCaller. 
# Employs an internal scatter-gather approach (SplitIntervals) to parallelize 
# calling across regions to heavily reduce walltime. Also runs mosdepth for 
# rapid coverage calculation.
#
# Dependencies: gatk, samtools, mosdepth
# ==============================================================================

export TMPDIR=$SCRATCHDIR

BAM_INPUT_DIR="/path/to/cluster/project_dir/bam_output"
GVCF_OUTPUT_DIR="/path/to/cluster/project_dir/bam_output/gvcfs"
STATS_OUTPUT_DIR="/path/to/cluster/project_dir/bam_output/stats"

mkdir -p "$GVCF_OUTPUT_DIR"
mkdir -p "$STATS_OUTPUT_DIR"

SAMPLE_LIST="/path/to/cluster/project_dir/bam_output/samplelist.txt"
REF_GENOME="/path/to/cluster/project_dir/genome_dir/Genome.fasta"
MOSDEPTH_EXE="/path/to/cluster/tools/mosdepth"

JOBS_PARALLEL=16
MEM_JAVA="-Xmx6g"

# --- Select Sample & Stage Data ---
SAMPLE=$(sed -n "${PBS_ARRAY_INDEX}p" "$SAMPLE_LIST" | tr -d '\r')

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found at line $PBS_ARRAY_INDEX"
    exit 1
fi

ORIGINAL_BAM="$BAM_INPUT_DIR/${SAMPLE}.bam"
FINAL_GVCF="$GVCF_OUTPUT_DIR/${SAMPLE}.g.vcf.gz"

if [ ! -f "$ORIGINAL_BAM" ]; then
    echo "CRITICAL ERROR: Input BAM not found: $ORIGINAL_BAM"
    exit 1
fi

cp "$ORIGINAL_BAM" "$SCRATCHDIR/"

if [ -f "${ORIGINAL_BAM}.bai" ]; then
    cp "${ORIGINAL_BAM}.bai" "$SCRATCHDIR/"
elif [ -f "${ORIGINAL_BAM%.bam}.bai" ]; then
    cp "${ORIGINAL_BAM%.bam}.bai" "$SCRATCHDIR/"
else
    module add samtools
    samtools index "$SCRATCHDIR/${SAMPLE}.bam"
fi

INPUT_BAM="$SCRATCHDIR/${SAMPLE}.bam"

# --- Run Mosdepth Coverage ---
if [ ! -f "$STATS_OUTPUT_DIR/${SAMPLE}.per-base.bed.gz" ]; then
    "$MOSDEPTH_EXE" -t 4 --by 500 "$SCRATCHDIR/${SAMPLE}" "$INPUT_BAM"
    cp "$SCRATCHDIR/${SAMPLE}.per-base.bed.gz" "$STATS_OUTPUT_DIR/"
    cp "$SCRATCHDIR/${SAMPLE}.per-base.bed.gz.csi" "$STATS_OUTPUT_DIR/"
    cp "$SCRATCHDIR/${SAMPLE}.mosdepth.summary.txt" "$STATS_OUTPUT_DIR/"
fi

# --- Prepare Shards (GATK Native Scattering) ---
if [ -f "$FINAL_GVCF" ]; then
    echo "GVCF already exists. Cleaning up and exiting."
    clean_scratch
    exit 0
fi

module add gatk samtools

REF_BASE="${REF_GENOME%.*}"
cp ${REF_BASE}* "$SCRATCHDIR/"
REF_NAME=$(basename "$REF_GENOME")
REF_PATH="$SCRATCHDIR/$REF_NAME"

if [ ! -f "${REF_PATH}.fai" ]; then samtools faidx "$REF_PATH"; fi
if [ ! -f "${REF_PATH%.*}.dict" ]; then gatk CreateSequenceDictionary -R "$REF_PATH"; fi

awk '{print $1}' "${REF_PATH}.fai" > "$SCRATCHDIR/intervals.list"
ALL_INTERVAL_ARGS=""
while read -r CONTIG; do
    if [ -n "$CONTIG" ]; then
        ALL_INTERVAL_ARGS="$ALL_INTERVAL_ARGS -L $CONTIG"
    fi
done < "$SCRATCHDIR/intervals.list"

mkdir -p "$SCRATCHDIR/scatter_dir"

gatk SplitIntervals \
    -R "$REF_PATH" \
    $ALL_INTERVAL_ARGS \
    --scatter-count $JOBS_PARALLEL \
    -O "$SCRATCHDIR/scatter_dir"

CHUNKS=($(find "$SCRATCHDIR/scatter_dir" -type f))
NUM_CHUNKS=${#CHUNKS[@]}

if [ "$NUM_CHUNKS" -eq 0 ]; then
    echo "CRITICAL ERROR: No chunks created by SplitIntervals."
    exit 1
fi

# --- Parallel Haplotype Caller ---
run_shard() {
    local INTERVAL_FILE=$1
    local ID=$2
    local OUT_SHARD="$SCRATCHDIR/shards/${SAMPLE}_${ID}.g.vcf.gz"
    local SHARD_TMP="$SCRATCHDIR/tmp_${ID}"
    
    mkdir -p "$SCRATCHDIR/shards" "$SHARD_TMP"
    
    local RAW_COORD_ARGS=$(grep -v '^@' "$INTERVAL_FILE" | awk '{print "-L " $1 ":" $2 "-" $3}' | tr '\n' ' ')
    
    gatk --java-options "$MEM_JAVA" HaplotypeCaller \
        -R "$REF_PATH" \
        -I "$INPUT_BAM" \
        -O "$OUT_SHARD" \
        -ERC GVCF \
        ${RAW_COORD_ARGS} \
        --native-pair-hmm-threads 1 \
        --tmp-dir "$SHARD_TMP" \
        > "$SCRATCHDIR/shards/${SAMPLE}_${ID}.log" 2>&1
}

export -f run_shard
export SAMPLE INPUT_BAM REF_PATH MEM_JAVA SCRATCHDIR

i=0
for chunk in "${CHUNKS[@]}"; do
    run_shard "$chunk" "$i" &
    ((i++))
done
wait

# --- Merge Shards ---
find "$SCRATCHDIR/shards" -maxdepth 1 -name "${SAMPLE}_*.g.vcf.gz" | sort -V > "$SCRATCHDIR/merge.list"
NUM_SHARDS=$(cat "$SCRATCHDIR/merge.list" | wc -l)

if [ "$NUM_SHARDS" -eq 0 ]; then
    echo "CRITICAL ERROR: No shards were produced."
    exit 1
fi

gatk --java-options "-Xmx32g" MergeVcfs -I "$SCRATCHDIR/merge.list" -O "$FINAL_GVCF"
gatk IndexFeatureFile -I "$FINAL_GVCF"

clean_scratch

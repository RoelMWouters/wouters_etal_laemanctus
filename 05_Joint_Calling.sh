#!/bin/bash
#PBS -N Joint_Calling_Filtering
#PBS -l select=1:ncpus=16:mem=128gb:scratch_local=200gb
#PBS -l walltime=48:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 5: JOINT CALLING & PARALOG PURGING (FILTERING)
# ==============================================================================
# Description:
# Combines GVCF files from a population, performs Joint Genotyping using GATK, 
# and strictly filters the resulting VCF. Filtering removes multi-allelic sites, 
# non-SNPs, and purges paralogs based on abnormal depth thresholds (filtering out 
# sites falling outside 0.5x - 1.5x of the population average depth).
#
# Dependencies: gatk, bcftools, samtools, vcftools
# ==============================================================================

export TMPDIR=$SCRATCHDIR

GVCF_INPUT_DIR="/path/to/cluster/project_dir/bam_output/gvcfs"
OUTPUT_DIR="/path/to/cluster/project_dir/output_dir/clonality"
REF_GENOME="/path/to/cluster/project_dir/genome_dir/genome.fasta"

mkdir -p "$OUTPUT_DIR/vcfs" "$OUTPUT_DIR/stats" "$SCRATCHDIR/gvcfs_local" "$SCRATCHDIR/tmp"

SPECIES_ID="species"
JOBS_PARALLEL=16
MEM_JAVA="-Xmx7g"

# --- Setup & Stage Data ---
module add gatk bcftools samtools vcftools

REF_BASE="${REF_GENOME%.*}"
cp ${REF_BASE}* "$SCRATCHDIR/"
REF_NAME=$(basename "$REF_GENOME")
REF_PATH="$SCRATCHDIR/$REF_NAME"

if [ ! -f "${REF_PATH}.fai" ]; then samtools faidx "$REF_PATH"; fi
if [ ! -f "${REF_PATH%.*}.dict" ]; then gatk CreateSequenceDictionary -R "$REF_PATH"; fi

cp "$GVCF_INPUT_DIR"/*.g.vcf.gz "$SCRATCHDIR/gvcfs_local/"
cp "$GVCF_INPUT_DIR"/*.g.vcf.gz.tbi "$SCRATCHDIR/gvcfs_local/" 2>/dev/null || true

VARIANT_LIST="$SCRATCHDIR/variants.list"
ls "$SCRATCHDIR"/gvcfs_local/*.g.vcf.gz > "$VARIANT_LIST" 2>/dev/null

if [ ! -s "$VARIANT_LIST" ]; then
    echo "CRITICAL ERROR: No GVCF files found to merge."
    exit 1
fi

# --- Prepare Shards ---
awk '{print $1}' "${REF_PATH}.fai" > "$SCRATCHDIR/intervals.list"
ALL_INTERVAL_ARGS=""
while read -r CONTIG; do
    if [ -n "$CONTIG" ]; then ALL_INTERVAL_ARGS="$ALL_INTERVAL_ARGS -L $CONTIG"; fi
done < "$SCRATCHDIR/intervals.list"

mkdir -p "$SCRATCHDIR/scatter_dir"
gatk SplitIntervals -R "$REF_PATH" $ALL_INTERVAL_ARGS --scatter-count $JOBS_PARALLEL -O "$SCRATCHDIR/scatter_dir"

CHUNKS=($(find "$SCRATCHDIR/scatter_dir" -type f))

# --- Parallel Joint Calling ---
run_shard() {
    local INTERVAL_FILE=$1
    local ID=$2
    local OUT_COMBINED="$SCRATCHDIR/shards/combined_${ID}.g.vcf.gz"
    local OUT_RAW="$SCRATCHDIR/shards/raw_${ID}.vcf.gz"
    local SHARD_TMP="$SCRATCHDIR/tmp_${ID}"
    
    mkdir -p "$SCRATCHDIR/shards" "$SHARD_TMP"
    local RAW_COORD_ARGS=$(grep -v '^@' "$INTERVAL_FILE" | awk '{print "-L " $1 ":" $2 "-" $3}' | tr '\n' ' ')
    
    gatk --java-options "$MEM_JAVA" CombineGVCFs \
        -R "$REF_PATH" -V "$VARIANT_LIST" -O "$OUT_COMBINED" \
        ${RAW_COORD_ARGS} --tmp-dir "$SHARD_TMP" > "$SCRATCHDIR/shards/combine_${ID}.log" 2>&1
        
    gatk --java-options "$MEM_JAVA" GenotypeGVCFs \
        -R "$REF_PATH" -V "$OUT_COMBINED" -O "$OUT_RAW" \
        ${RAW_COORD_ARGS} --tmp-dir "$SHARD_TMP" > "$SCRATCHDIR/shards/genotype_${ID}.log" 2>&1
}

export -f run_shard
export REF_PATH VARIANT_LIST MEM_JAVA SCRATCHDIR

i=0
for chunk in "${CHUNKS[@]}"; do
    run_shard "$chunk" "$i" &
    pids[${i}]=$!
    ((i++))
done

FAIL=0
for pid in ${pids[*]}; do 
    wait $pid || { echo "CRITICAL ERROR!"; FAIL=1; }
done

if [ "$FAIL" -eq 1 ]; then exit 1; fi

# --- Merge Shards ---
RAW_VCF="$SCRATCHDIR/${SPECIES_ID}.raw.vcf.gz"
find "$SCRATCHDIR/shards" -maxdepth 1 -name "raw_*.vcf.gz" | sort -V > "$SCRATCHDIR/merge.list"

gatk --java-options "-Xmx32g" MergeVcfs -I "$SCRATCHDIR/merge.list" -O "$RAW_VCF"

# --- Strict Filtering (Paralog Purge) ---
SPLIT_VCF="$SCRATCHDIR/${SPECIES_ID}.split.vcf.gz"
SNP_VCF="$SCRATCHDIR/${SPECIES_ID}.snps.vcf.gz"
FINAL_VCF="$SCRATCHDIR/${SPECIES_ID}.strict_filtered.vcf.gz"

bcftools norm -m -any -f "$REF_PATH" "$RAW_VCF" -O z -o "$SPLIT_VCF"
bcftools view -v snps -m 2 -M 2 "$SPLIT_VCF" -O z -o "$SNP_VCF"

AVG_DP=$(bcftools query -f '[%DP\n]' "$SNP_VCF" | awk '
    BEGIN {sum=0; n=0} 
    {if ($1 != ".") {sum+=$1; n++}} 
    END {if (n>0) printf "%.0f\n", (sum/n); else print 15}
')

MIN_DP=$(awk "BEGIN {printf \"%.0f\", $AVG_DP * 0.5}")
MAX_DP=$(awk "BEGIN {printf \"%.0f\", $AVG_DP * 1.5}")

bcftools filter \
  -i "INFO/MQ > 40 && INFO/QD > 2.0 && INFO/FS < 60.0" \
  "$SNP_VCF" -O z -o "$SCRATCHDIR/mq_filtered.vcf.gz"

vcftools --gzvcf "$SCRATCHDIR/mq_filtered.vcf.gz" \
    --minDP $MIN_DP \
    --maxDP $MAX_DP \
    --min-alleles 2 \
    --max-alleles 2 \
    --maf 0.05 \
    --max-missing 1.0 \
    --recode --recode-INFO-all \
    --out "$SCRATCHDIR/${SPECIES_ID}_purged"

bgzip -c "$SCRATCHDIR/${SPECIES_ID}_purged.recode.vcf" > "$FINAL_VCF"
tabix -p vcf "$FINAL_VCF"

# --- Metrics & Save ---
bcftools gtcheck "$FINAL_VCF" > "$OUTPUT_DIR/stats/${SPECIES_ID}_gtcheck.txt"

cp "$FINAL_VCF" "$OUTPUT_DIR/vcfs/"
cp "$FINAL_VCF.tbi" "$OUTPUT_DIR/vcfs/"

clean_scratch

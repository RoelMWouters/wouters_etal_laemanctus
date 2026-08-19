#!/bin/bash
#PBS -N YAHS_Scaffolding
#PBS -l select=1:ncpus=32:mem=256gb:scratch_local=800gb
#PBS -l walltime=24:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 2A: GENOME SCAFFOLDING USING YAHS
# ==============================================================================
# Description:
# Scaffolds a draft assembly using Hi-C data. Maps reads with BWA-MEM, removes 
# duplicates, runs the YAHS scaffolder, evaluates the output with QUAST, and 
# generates a Juicer `.hic` map for manual curation in Juicebox.
#
# Dependencies: bwa, samtools, openjdk, mambaforge/conda (yahs, juicer_tools)
# ==============================================================================

# --- Configuration Block ---
PREFIX="species_uncollapsed"
PROJECT_DIR="/path/to/cluster/project_dir/yahs_uncollapsed_3D"
DRAFT_GENOME="/path/to/cluster/project_dir/hifiasm_assembly/final_results/Assembly_Output/primary.fasta"
HIC_READ_1="/path/to/cluster/project_dir/hic/trimmed_data/trimmed_R1.fq.gz"
HIC_READ_2="/path/to/cluster/project_dir/hic/trimmed_data/trimmed_R2.fq.gz"
RESTRICTION_MOTIFS="GATC,GANTC"

JUICER_JAR="/path/to/cluster/tools/juicer_tools.jar"
MAMBA_ENV="/path/to/cluster/envs/yahs_env"

FINAL_DIR="$PROJECT_DIR/final_results"
mkdir -p "$FINAL_DIR"

rescue_data() {
    echo "=== EMERGENCY EXIT: RESCUING DATA TO STORAGE ==="
    cp $SCRATCHDIR/*_scaffolds_final.fa "$FINAL_DIR/" 2>/dev/null || true
    cp $SCRATCHDIR/*.hic "$FINAL_DIR/" 2>/dev/null || true
    cp $SCRATCHDIR/*.agp "$FINAL_DIR/" 2>/dev/null || true
    cp $SCRATCHDIR/*.assembly "$FINAL_DIR/" 2>/dev/null || true
    cp $SCRATCHDIR/*.bin "$FINAL_DIR/" 2>/dev/null || true
    cp $SCRATCHDIR/*.fai "$FINAL_DIR/" 2>/dev/null || true
    cp $SCRATCHDIR/*_mapping_*.txt "$FINAL_DIR/" 2>/dev/null || true
    cp -r $SCRATCHDIR/QUAST_Output "$FINAL_DIR/" 2>/dev/null || true
}
trap rescue_data EXIT

cd "$SCRATCHDIR"
set -e
set -o pipefail

# --- Environment Setup ---
echo "=== INITIALIZING ENVIRONMENT ==="
module add bwa samtools openjdk mambaforge

source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$MAMBA_ENV"

export _JAVA_OPTIONS="-Djava.io.tmpdir=$SCRATCHDIR"

cp "$DRAFT_GENOME" "$SCRATCHDIR/draft_genome.fasta"
cp "$JUICER_JAR" "$SCRATCHDIR/juicer_tools.jar"

echo "=== INDEXING DRAFT GENOME ==="
samtools faidx draft_genome.fasta

# --- Step 1: Mapping & Deduplication ---
if [ -f "$FINAL_DIR/${PREFIX}_deduplicated.bam" ]; then
    echo "=== CHECKPOINT 1 FOUND: Skipping alignment. ==="
    cp "$FINAL_DIR/${PREFIX}_deduplicated.bam" "$SCRATCHDIR/"
else
    echo "=== CHECKPOINT 1 NOT FOUND: Running BWA-MEM and Deduplication ==="
    cp "$HIC_READ_1" "$SCRATCHDIR/read_1.fq.gz"
    cp "$HIC_READ_2" "$SCRATCHDIR/read_2.fq.gz"
    
    bwa index draft_genome.fasta

    bwa mem -5SP -t 32 draft_genome.fasta read_1.fq.gz read_2.fq.gz | \
        samtools view -bhS - | \
        samtools sort -n -@ 32 -m 2G -o namesorted.bam -

    samtools fixmate -m -@ 32 namesorted.bam fixmate.bam
    samtools sort -@ 32 -m 2G -o coordsorted.bam fixmate.bam
    samtools markdup -r -@ 32 coordsorted.bam ${PREFIX}_deduplicated.bam
    
    samtools flagstat -@ 32 ${PREFIX}_deduplicated.bam > ${PREFIX}_mapping_flagstat.txt
    samtools stats -@ 32 ${PREFIX}_deduplicated.bam > ${PREFIX}_mapping_stats.txt
    
    cp ${PREFIX}_deduplicated.bam "$FINAL_DIR/"
    cp ${PREFIX}_mapping_flagstat.txt "$FINAL_DIR/"
    cp ${PREFIX}_mapping_stats.txt "$FINAL_DIR/"
    
    rm namesorted.bam fixmate.bam coordsorted.bam read_1.fq.gz read_2.fq.gz
fi

# --- Step 2: YAHS Scaffolding ---
if [ -f "$FINAL_DIR/yahs.out.bin" ]; then
    echo "=== CHECKPOINT 2 FOUND: Skipping scaffolding. ==="
    cp "$FINAL_DIR/yahs.out_scaffolds_final.fa" "$SCRATCHDIR/"
    cp "$FINAL_DIR/yahs.out_scaffolds_final.agp" "$SCRATCHDIR/"
    cp "$FINAL_DIR/yahs.out.bin" "$SCRATCHDIR/"
else
    echo "=== CHECKPOINT 2 NOT FOUND: Running YAHS ==="
    yahs -e "$RESTRICTION_MOTIFS" draft_genome.fasta ${PREFIX}_deduplicated.bam
    cp "$SCRATCHDIR"/yahs.out* "$FINAL_DIR/"
fi

# --- Step 3: Juicer Pre (Heatmap Generation) ---
if [ -f "$FINAL_DIR/${PREFIX}_scaffolded.hic" ]; then
    echo "=== CHECKPOINT 3 FOUND: Skipping heatmap generation. ==="
    cp "$FINAL_DIR/${PREFIX}_scaffolded.hic" "$SCRATCHDIR/"
    cp "$FINAL_DIR/curatable_map.assembly" "$SCRATCHDIR/"
else
    echo "=== CHECKPOINT 3 NOT FOUND: Generating interactive heatmap files ==="
    samtools faidx yahs.out_scaffolds_final.fa
    juicer pre -a -o curatable_map yahs.out.bin yahs.out_scaffolds_final.agp draft_genome.fasta.fai

    # Dynamic scaling for large genomes
    TRUE_SIZE=$(awk '{sum+=$2} END {print sum}' yahs.out_scaffolds_final.fa.fai)
    if [ "$TRUE_SIZE" -gt 2147483647 ]; then
        echo "WARNING: Genome > 2.14 Gb detected. Scaling coordinates."
        SCALED_SIZE=$(( TRUE_SIZE / 2 ))
    else
        SCALED_SIZE=$TRUE_SIZE
    fi
    echo "assembly $SCALED_SIZE" > jbat_chrom.sizes

    java -Xmx150G -jar juicer_tools.jar pre -n curatable_map.txt ${PREFIX}_scaffolded.hic jbat_chrom.sizes

    cp ${PREFIX}_scaffolded.hic "$FINAL_DIR/"
    cp curatable_map.assembly "$FINAL_DIR/"
fi

# --- Step 4: QUAST Evaluation ---
if [ -d "$FINAL_DIR/QUAST_Output" ]; then
    echo "=== CHECKPOINT 4 FOUND: Skipping QUAST evaluation. ==="
else
    echo "=== CHECKPOINT 4 NOT FOUND: Running QUAST ==="
    quast.py draft_genome.fasta yahs.out_scaffolds_final.fa -o QUAST_Output -t 32 --large
    cp -r QUAST_Output "$FINAL_DIR/"
fi

echo "=== PIPELINE COMPLETE ==="
trap - EXIT

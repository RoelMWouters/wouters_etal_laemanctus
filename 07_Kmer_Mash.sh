#!/bin/bash
#PBS -N Kmer_Mash
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=50gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 7: FAST K-MER DISTANCE ESTIMATION (MASH)
# ==============================================================================
# Description:
# Computes a fast k-mer based distance matrix directly from raw fastq reads 
# to assess clustering/clonality prior to mapping. Uses Mash sketching.
#
# Dependencies: curl, tar, mash (downloads binary internally)
# ==============================================================================

export TMPDIR=$SCRATCHDIR

RAW_READS_DIR="/path/to/cluster/project_dir/raw_fastqs"
STATS_DIR="/path/to/cluster/project_dir/output_dir/clonality/stats"
SAMPLES=("S01" "S02" "S03" "S04" "S05")

mkdir -p "$STATS_DIR" "$SCRATCHDIR"
cd "$SCRATCHDIR"

# Download standalone Mash
curl -L -s https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar -o mash-Linux64.tar
tar -xf mash-Linux64.tar
MASH_BIN="$SCRATCHDIR/mash-Linux64-v2.3/mash"
chmod +x "$MASH_BIN"

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_FOLDER="${RAW_READS_DIR}/${SAMPLE}"
    
    shopt -s nullglob
    FASTQ_FILES=("$SAMPLE_FOLDER"/*.f*q*)
    shopt -u nullglob

    if [ ${#FASTQ_FILES[@]} -eq 0 ]; then continue; fi
    $MASH_BIN sketch -r -m 2 -k 31 -s 10000000 -o "$SCRATCHDIR/${SAMPLE}" "${FASTQ_FILES[@]}"
done

$MASH_BIN paste "$SCRATCHDIR/all_samples" "$SCRATCHDIR"/*.msh
$MASH_BIN triangle "$SCRATCHDIR/all_samples.msh" > "$STATS_DIR/kmer_distance_matrix.txt"

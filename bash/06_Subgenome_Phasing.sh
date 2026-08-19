#!/bin/bash
#PBS -N Subgenome_Phasing
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=50gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 6: SUBGENOME PHASING VIA MAPPING BIAS
# ==============================================================================
# Description:
# Segregates a draft genome assembly into distinct parental subgenomes (e.g., 
# for hybrids) by analyzing the read mapping depth/bias from the two parental 
# lineages. Calculates Reads Per Million (RPM) per contig.
#
# Dependencies: samtools, seqtk, awk
# ==============================================================================

export TMPDIR=$SCRATCHDIR

REF_GENOME="/path/to/cluster/project_dir/genome_dir/genome.fasta"
PARENT_A_BAM="/path/to/cluster/project_dir/bam_output/ParentA_final_markdup.bam"
PARENT_B_BAM="/path/to/cluster/project_dir/bam_output/ParentB_final_markdup.bam"
OUTPUT_DIR="/path/to/cluster/project_dir/genome_dir/haploid_genome"

BIAS_THRESHOLD=0.8
MIN_CONTIG_LEN=10000 

module add samtools seqtk

mkdir -p "$OUTPUT_DIR" "$SCRATCHDIR"
cd "$SCRATCHDIR" 

ln -s "$PARENT_A_BAM" "$SCRATCHDIR/parentA.bam"
cp "${PARENT_A_BAM}.bai" "$SCRATCHDIR/parentA.bam.bai"

ln -s "$PARENT_B_BAM" "$SCRATCHDIR/parentB.bam"
cp "${PARENT_B_BAM}.bai" "$SCRATCHDIR/parentB.bam.bai"

cp "$REF_GENOME" "$SCRATCHDIR/genome.fasta"

# Calculate lengths and mapped reads
samtools idxstats "$SCRATCHDIR/parentA.bam" | cut -f 1,2,3 | sort -k1,1 > "$SCRATCHDIR/stats_A.txt"
samtools idxstats "$SCRATCHDIR/parentB.bam" | cut -f 1,3   | sort -k1,1 > "$SCRATCHDIR/stats_B.txt"

join -1 1 -2 1 "$SCRATCHDIR/stats_A.txt" "$SCRATCHDIR/stats_B.txt" > "$SCRATCHDIR/combined_stats.txt"

TOTAL_A=$(awk '{sum+=$3} END {print sum}' "$SCRATCHDIR/combined_stats.txt")
TOTAL_B=$(awk '{sum+=$4} END {print sum}' "$SCRATCHDIR/combined_stats.txt")

# Phase contigs based on RPM bias
awk -v totA="$TOTAL_A" -v totB="$TOTAL_B" -v thresh="$BIAS_THRESHOLD" -v minLen="$MIN_CONTIG_LEN" '
BEGIN {
    OFS="\t"
    print "Contig", "Length", "Mapped_A", "Mapped_B", "RPM_A", "RPM_B", "Bias_A", "Assignment" > "final_phasing_report.tsv"
}
{
    contig=$1
    len=$2
    readsA=$3
    readsB=$4
    
    if (contig == "*" || len < minLen) next;
    
    rpmA = (readsA / totA) * 1000000
    rpmB = (readsB / totB) * 1000000
    totalRPM = rpmA + rpmB
    
    if (totalRPM == 0) {
        biasA = 0.5
        assign = "Unmapped"
    } else {
        biasA = rpmA / totalRPM
        
        if (biasA >= thresh) {
            assign = "Parent_A"
            print contig > "list_Parent_A.txt"
        } else if (biasA <= (1 - thresh)) {
            assign = "Parent_B"
            print contig > "list_Parent_B.txt"
        } else {
            assign = "Shared"
            print contig > "list_Shared.txt"
        }
    }
    print contig, len, readsA, readsB, rpmA, rpmB, biasA, assign > "final_phasing_report.tsv"
}' "$SCRATCHDIR/combined_stats.txt"

# Extract Subgenome FASTAs
cat "list_Parent_A.txt" "list_Shared.txt" > "$SCRATCHDIR/target_A.txt"
seqtk subseq "$SCRATCHDIR/genome.fasta" "$SCRATCHDIR/target_A.txt" > "$SCRATCHDIR/Subgenome_A.fasta"

cat "list_Parent_B.txt" "list_Shared.txt" > "$SCRATCHDIR/target_B.txt"
seqtk subseq "$SCRATCHDIR/genome.fasta" "$SCRATCHDIR/target_B.txt" > "$SCRATCHDIR/Subgenome_B.fasta"

cp "$SCRATCHDIR/final_phasing_report.tsv" "$OUTPUT_DIR/"
cp "$SCRATCHDIR/Subgenome_A.fasta" "$OUTPUT_DIR/"
cp "$SCRATCHDIR/Subgenome_B.fasta" "$OUTPUT_DIR/"
cp "$SCRATCHDIR"/list_*.txt "$OUTPUT_DIR/"

clean_scratch

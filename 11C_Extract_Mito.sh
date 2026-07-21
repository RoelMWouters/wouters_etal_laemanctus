#!/bin/bash
#PBS -N Extract_Mito_BAM
#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=20gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==============================================================================
# SCRIPT 11C: CONSENSUS EXTRACTION FROM BAM
# ==============================================================================
# Description:
# Parses BLAST coordinates from 11B, and extracts those target regions directly 
# from sample BAM files by creating a local consensus sequence. This effectively 
# recovers the mitochondrial markers for each sequenced individual.
#
# Dependencies: blast, bedtools, samtools, bcftools, seqtk
# ==============================================================================

export TMPDIR=$SCRATCHDIR

module add blast bedtools samtools bcftools seqtk

GENOME_DIR="/path/to/cluster/project_dir/genome_dir"
BAM_DIR="/path/to/cluster/project_dir/bam_output"
OUT_DIR="/path/to/cluster/project_dir/output_dir/mt_markers"
BAITS="${OUT_DIR}/markers.fasta"

cd "${OUT_DIR}"

makeblastdb -in "${GENOME_DIR}/genome.fasta" -dbtype nucl -out species_db
makeblastdb -in "${GENOME_DIR}/ref.fna" -dbtype nucl -out outgroup_db

blastn -task blastn -query "${BAITS}" -db species_db -outfmt "6 sseqid sstart send qseqid" -max_target_seqs 1 -max_hsps 1 > blast_species.txt
awk '{OFS="\t"; if ($2 < $3) print $1, $2-1, $3, $4; else print $1, $3-1, $2, $4}' blast_species.txt > coords_species.bed

for BAM in "${BAM_DIR}"/*.bam; do
    SAMPLE=$(basename "${BAM}" .bam)
    while read -r contig start end gene; do
        REGION="${contig}:${start}-${end}"
        samtools view -b "${BAM}" "${REGION}" | \
          bcftools mpileup -Ou -f "${GENOME_DIR}/genome.fasta" - | \
          bcftools call -c - | \
          vcfutils.pl vcf2fq | \
          seqtk seq -A > tmp_${SAMPLE}_${gene}.fasta
        sed -i "1s/.*/\>${SAMPLE}/" tmp_${SAMPLE}_${gene}.fasta
    done < coords_species.bed
done

grep "^>" "${BAITS}" | sed 's/>//' > gene_list.txt
while read -r gene; do
    cat tmp_*_${gene}.fasta >> compiled_${gene}.fasta
done < gene_list.txt

rm tmp_*.fasta blast_*.txt coords_*.bed gene_list.txt

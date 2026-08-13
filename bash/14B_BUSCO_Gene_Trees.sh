#!/bin/bash
#PBS -N busco_extract_align
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=100gb
#PBS -l walltime=60:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 14B (STREAMLINED): BUSCO EXTRACTION & ALIGNMENT PIPELINE
# ==============================================================================
# Description:
# Extracts BUSCO loci from a population VCF, filters out paralogous "ghost" 
# alleles using parental reference BLASTs, and aligns the true orthologs 
# with MAFFT. Applies a PI-site QC check and saves the clean alignments 
# for downstream gCF/Supermatrix analysis.
# ==============================================================================

export TMPDIR=$SCRATCHDIR

# --- Configuration Paths ---
BASE_DIR="/path/to/cluster/project_dir"
REF="${BASE_DIR}/genome_dir/target_genome.fasta"
VCF="${BASE_DIR}/output_dir/clonality/vcfs/population_final_filtered.vcf.gz"
BED="${BASE_DIR}/genome_dir/combined_buscos.bed"

OUTGROUP_GENOME="${BASE_DIR}/genome_dir/outgroup.fna"
PAT_GENOME="${BASE_DIR}/genome_dir/paternal_ref.fna" 
MAT_GENOME="${BASE_DIR}/genome_dir/haploid_genome/maternal_ref.fasta"

SAMPLES=("Sample1" "Sample2" "Sample3" "Sample4")

OUT_DIR="${BASE_DIR}/output_dir/busco_phylogeny"
# STREAMLINE: Changed output to alignments_clean and removed trees directory
mkdir -p "$OUT_DIR/fastas" "$OUT_DIR/alignments_clean"
mkdir -p "$OUT_DIR/intermediates/unaligned_with_refs"
mkdir -p "$SCRATCHDIR/unaligned_genes"

# Ensure scratch cleanup on exit
trap 'rm -rf "$SCRATCHDIR"/*' TERM EXIT

# ==============================================================================
# PHASE 1: EXTRACTION 
# ==============================================================================
echo -e "\n--- PHASE 1: Sequence Extraction ---"
module purge
module add bcftools bedtools samtools

for SAMPLE in "${SAMPLES[@]}"; do
    if [ -s "$OUT_DIR/fastas/${SAMPLE}_buscos.fa" ]; then
        echo " -> $SAMPLE FASTA exists. Skipping."
    else
        echo "Extracting alleles for $SAMPLE..."
        bcftools consensus -s "$SAMPLE" -H 1 "$VCF" < "$REF" > "$SCRATCHDIR/${SAMPLE}_temp.fa"
        bedtools getfasta -fi "$SCRATCHDIR/${SAMPLE}_temp.fa" -bed "$BED" -name -s -fo "$OUT_DIR/fastas/${SAMPLE}_buscos.fa"
        rm -f "$SCRATCHDIR/${SAMPLE}_temp.fa"
        sed -i "s/>/>${SAMPLE}_/g" "$OUT_DIR/fastas/${SAMPLE}_buscos.fa"
    fi
done

# ==============================================================================
# PHASE 2: GROUPING & HEADER CLEANING
# ==============================================================================
echo -e "\n--- PHASE 2: Grouping & Header Cleaning ---"
if ls "$OUT_DIR/intermediates/unaligned_with_refs/"*.fasta 1> /dev/null 2>&1; then
    echo " -> Saved files found. Copying to scratch..."
    cp "$OUT_DIR/intermediates/unaligned_with_refs/"*.fasta "$SCRATCHDIR/unaligned_genes/"
    SKIP_PHASE_3=true
else
    echo " -> Cleaning headers and grouping copies..."
    cat "$OUT_DIR/fastas"/*_buscos.fa > "$SCRATCHDIR/all_pool.fa"

    awk -v outdir="$SCRATCHDIR/unaligned_genes" '
      /^>/ { 
        split($0, a, "_"); 
        sample = a[1];
        gene = a[2];
        copy_part = a[3];
        split(copy_part, c_clean, ":");
        header = sample "_" gene "_" c_clean[1];
        filename = outdir "/" gene ".fasta";
        print header >> filename;
        next;
      }
      { print >> filename }
    ' "$SCRATCHDIR/all_pool.fa"
    SKIP_PHASE_3=false
fi

# ==============================================================================
# PHASE 2.5: GHOST BUSTER (MATERNAL/PATERNAL FILTERING)
# ==============================================================================
if [ "$SKIP_PHASE_3" = false ]; then
    echo -e "\n--- PHASE 2.5: Removing Ghost Alleles via Reference BLAST ---"
    module purge
    module add blast-plus

    makeblastdb -in "$MAT_GENOME" -dbtype nucl -out "$SCRATCHDIR/mat_db" > /dev/null
    makeblastdb -in "$PAT_GENOME" -dbtype nucl -out "$SCRATCHDIR/pat_db" > /dev/null

    for FILE in "$SCRATCHDIR/unaligned_genes"/*.fasta; do
        BASENAME=$(basename "$FILE" .fasta)
        
        # 1. Filter Maternal 
        grep -A 1 "MatSample_${BASENAME}" "$FILE" > "$SCRATCHDIR/mat_query.fa"
        blastn -query "$SCRATCHDIR/mat_query.fa" -db "$SCRATCHDIR/mat_db" -outfmt "6 qseqid bitscore" -max_target_seqs 1 > "$SCRATCHDIR/mat_hits.txt"
        
        if [ $(wc -l < "$SCRATCHDIR/mat_hits.txt") -eq 2 ]; then
            GHOST_MAT=$(sort -k2,2nr "$SCRATCHDIR/mat_hits.txt" | tail -n 1 | awk '{print $1}')
            awk -v ghost="$GHOST_MAT" 'BEGIN{RS=">";ORS=""} NR>1 {if($1!=ghost) print ">"$0}' "$FILE" > "$SCRATCHDIR/temp.fa" && mv "$SCRATCHDIR/temp.fa" "$FILE"
            echo "   -> Dropped Maternal Ghost: $GHOST_MAT"
        fi

        # 2. Filter Paternal
        grep -A 1 "PatSample_${BASENAME}" "$FILE" > "$SCRATCHDIR/pat_query.fa"
        blastn -query "$SCRATCHDIR/pat_query.fa" -db "$SCRATCHDIR/pat_db" -outfmt "6 qseqid bitscore" -max_target_seqs 1 > "$SCRATCHDIR/pat_hits.txt"
        
        if [ $(wc -l < "$SCRATCHDIR/pat_hits.txt") -eq 2 ]; then
            GHOST_PAT=$(sort -k2,2nr "$SCRATCHDIR/pat_hits.txt" | tail -n 1 | awk '{print $1}')
            awk -v ghost="$GHOST_PAT" 'BEGIN{RS=">";ORS=""} NR>1 {if($1!=ghost) print ">"$0}' "$FILE" > "$SCRATCHDIR/temp.fa" && mv "$SCRATCHDIR/temp.fa" "$FILE"
            echo "   -> Dropped Paternal Ghost: $GHOST_PAT"
        fi
    done
fi

# ==============================================================================
# PHASE 3: HARVESTING OUTGROUP
# ==============================================================================
if [ "$SKIP_PHASE_3" = false ]; then
    echo -e "\n--- PHASE 3: Fetching Outgroup Sequence ---"
    module purge
    module add blast-plus samtools

    cp "$OUTGROUP_GENOME" "$SCRATCHDIR/outgroup.fasta"
    makeblastdb -in "$SCRATCHDIR/outgroup.fasta" -dbtype nucl -out "$SCRATCHDIR/outgroup_db" > /dev/null
    samtools faidx "$SCRATCHDIR/outgroup.fasta"
    
    # STREAMLINE: Paternal reference harvesting removed entirely

    for FILE in "$SCRATCHDIR/unaligned_genes"/*.fasta; do
        BASENAME=$(basename "$FILE" .fasta)
        head -n 2 "$FILE" > "$SCRATCHDIR/query_${BASENAME}.fa"
        
        blastn -task blastn -query "$SCRATCHDIR/query_${BASENAME}.fa" -db "$SCRATCHDIR/outgroup_db" -outfmt "6 sseqid sstart send" -max_target_seqs 1 -max_hsps 1 > "$SCRATCHDIR/hit_out_${BASENAME}.txt"
        
        if [ -s "$SCRATCHDIR/hit_out_${BASENAME}.txt" ]; then
            read SSEQID SSTART SEND < "$SCRATCHDIR/hit_out_${BASENAME}.txt"
            if [ "$SSTART" -lt "$SEND" ]; then START=$SSTART; END=$SEND; REV=false; else START=$SEND; END=$SSTART; REV=true; fi
            samtools faidx "$SCRATCHDIR/outgroup.fasta" "${SSEQID}:${START}-${END}" | grep -v "^>" | tr -d '\n' > "$SCRATCHDIR/out_seq.txt"
            if [ "$REV" = true ]; then cat "$SCRATCHDIR/out_seq.txt" | tr 'ACGTacgt' 'TGCAtgca' | rev > "$SCRATCHDIR/out_seq_rev.txt"; mv "$SCRATCHDIR/out_seq_rev.txt" "$SCRATCHDIR/out_seq.txt"; fi
            if [ -s "$SCRATCHDIR/out_seq.txt" ]; then echo -e "\n>OUTGROUP" >> "$FILE"; cat "$SCRATCHDIR/out_seq.txt" >> "$FILE"; echo "" >> "$FILE"; fi
        fi
    done
fi

# ==============================================================================
# PHASE 4: ALIGNMENT & DRY-RUN QC
# ==============================================================================
echo -e "\n--- PHASE 4: Alignment & Quality Control ---"
module purge
module add mafft 
module add iqtree/3.0.1-mpi

for FILE in "$SCRATCHDIR/unaligned_genes"/*.fasta; do
    BASENAME=$(basename "$FILE" .fasta)
    ALIGNED_FILE="$OUT_DIR/alignments_clean/${BASENAME}_aligned.fasta"
    
    if [ -s "$ALIGNED_FILE" ]; then
        echo "   -> Skipping $BASENAME (Already aligned)"
        continue
    fi

    echo ">> Aligning: $BASENAME"
    mafft --auto "$FILE" > "$ALIGNED_FILE" 2>/dev/null
    
    # The Goldilocks QC Check (Dry run IQ-TREE to parse SNPs)
    iqtree3-mpi -s "$ALIGNED_FILE" -m JC -n 0 -pre "$SCRATCHDIR/qc_${BASENAME}" -quiet
    PI_SITES=$(grep "Number of parsimony informative sites:" "$SCRATCHDIR/qc_${BASENAME}.iqtree" | awk '{print $6}')
    
    if [ -z "$PI_SITES" ]; then PI_SITES=0; fi
    echo "   -> Parsimony Informative Sites: $PI_SITES"

    # STREAMLINE: Discard bad alignments, but skip actual tree generation
    if [ "$PI_SITES" -lt 10 ]; then
        echo "   -> SKIPPING: Too conserved (< 10 PI sites). Deleting alignment."
        rm -f "$ALIGNED_FILE"
        continue
    elif [ "$PI_SITES" -gt 150 ]; then
        echo "   -> SKIPPING: Too noisy/divergent (> 150 PI sites). Deleting alignment."
        rm -f "$ALIGNED_FILE"
        continue
    fi

    echo "   -> Alignment passed QC and saved to alignments_clean."
done

# ==============================================================================
# PHASE 5: SAVE INTERMEDIATES
# ==============================================================================
if [ "$SKIP_PHASE_3" = false ]; then
    cp -r "$SCRATCHDIR/unaligned_genes/"* "$OUT_DIR/intermediates/unaligned_with_refs/"
fi

echo "=================================================="
echo "PIPELINE COMPLETE - DATA READY FOR gCF SCRIPT"
echo "=================================================="

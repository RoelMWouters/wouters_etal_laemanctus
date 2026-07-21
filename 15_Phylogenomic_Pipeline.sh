#!/bin/bash
#PBS -N Phylogenomic_Tree_Pipeline
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=50gb
#PBS -l walltime=48:00:00
#PBS -M your.email@institution.edu
#PBS -m abe
#PBS -j oe

# ==============================================================================
# SCRIPT 15: UNIFIED PHYLOGENOMIC TREE INFERENCE PIPELINE
# ==============================================================================
# Description:
# This superscript executes a complete phylogenomic workflow from grouped markers:
#   PHASE 1: Harvests and appends an outgroup sequence via BLAST.
#   PHASE 2: Aligns sequences and runs QC (filtering by PI sites).
#   PHASE 3: Infers Maximum Likelihood gene trees using IQ-TREE.
#   PHASE 4: Cleans tree tip labels and generates an IQ-TREE consensus tree.
#   PHASE 5: Infers the final coalescent Species Tree using ASTRAL.
#
# Dependencies: blast, samtools, mafft, iqtree, java, unzip, wget
# ==============================================================================

export TMPDIR=$SCRATCHDIR

# --- Configuration Paths ---
WORK_DIR="/path/to/cluster/project_dir/output_dir/busco_phylogeny"
INPUT_FASTAS="$WORK_DIR/grouped_fastas"       # Folder with your pre-grouped gene sequences
OUTGROUP_GENOME="/path/to/cluster/project_dir/genome_dir/outgroup.fna"
QUERY_SAMPLE="F01_HapA"                       # The sample ID used as the "query dart" to find outgroups

ALIGN_DIR="$WORK_DIR/alignments"
TREE_DIR="$WORK_DIR/trees"
FINAL_DIR="$WORK_DIR/final_species_trees"

mkdir -p "$ALIGN_DIR" "$TREE_DIR" "$FINAL_DIR"
cd "$SCRATCHDIR" || exit 1

# ==============================================================================
# PHASE 1: OUTGROUP HARVESTER
# ==============================================================================
echo -e "\n=== PHASE 1: HARVESTING OUTGROUP SEQUENCES ==="
module purge
module add blast-plus samtools

cp "$OUTGROUP_GENOME" "$SCRATCHDIR/outgroup.fasta"

echo "Building Outgroup BLAST Database..."
makeblastdb -in "$SCRATCHDIR/outgroup.fasta" -dbtype nucl -out "$SCRATCHDIR/outgroup_db" >/dev/null
samtools faidx "$SCRATCHDIR/outgroup.fasta"

# Loop through all available FASTA gene groups
for MASTER_FASTA in "$INPUT_FASTAS"/*.fasta; do
    [ -e "$MASTER_FASTA" ] || continue
    GENE_NAME=$(basename "$MASTER_FASTA" .fasta)
    
    echo "Processing Outgroup for $GENE_NAME..."
    
    # Extract the query sequence for the specific sample to act as our BLAST dart
    awk -v target="$QUERY_SAMPLE" 'BEGIN{RS=">";ORS=""} $1 ~ target {print ">"$0}' "$MASTER_FASTA" > "$SCRATCHDIR/query.fa"
    
    if [ ! -s "$SCRATCHDIR/query.fa" ]; then
        echo "   -> WARNING: Query sample $QUERY_SAMPLE not found in $GENE_NAME. Skipping outgroup."
        cp "$MASTER_FASTA" "$SCRATCHDIR/${GENE_NAME}_ready.fasta"
        continue
    fi

    blastn -query "$SCRATCHDIR/query.fa" -db "$SCRATCHDIR/outgroup_db" \
        -outfmt "6 sseqid sstart send evalue" \
        -max_target_seqs 1 -max_hsps 1 > "$SCRATCHDIR/top_hit.txt"
        
    if [ -s "$SCRATCHDIR/top_hit.txt" ]; then
        read OUT_CONTIG OUT_START OUT_END EVALUE < "$SCRATCHDIR/top_hit.txt"
        
        # Handle reverse strand logic
        if [ "$OUT_START" -gt "$OUT_END" ]; then
            REAL_START=$OUT_END; REAL_END=$OUT_START; STRAND="minus"
        else
            REAL_START=$OUT_START; REAL_END=$OUT_END; STRAND="plus"
        fi
        
        samtools faidx "$SCRATCHDIR/outgroup.fasta" "${OUT_CONTIG}:${REAL_START}-${REAL_END}" > "$SCRATCHDIR/outgroup_raw.fa"
        
        cp "$MASTER_FASTA" "$SCRATCHDIR/${GENE_NAME}_ready.fasta"
        echo -e "\n>OUTGROUP" >> "$SCRATCHDIR/${GENE_NAME}_ready.fasta"
        
        if [ "$STRAND" == "minus" ]; then
            grep -v ">" "$SCRATCHDIR/outgroup_raw.fa" | tr -d '\n' | rev | tr 'ACGTacgt' 'TGCAtgca' >> "$SCRATCHDIR/${GENE_NAME}_ready.fasta"
        else
            grep -v ">" "$SCRATCHDIR/outgroup_raw.fa" >> "$SCRATCHDIR/${GENE_NAME}_ready.fasta"
        fi
        echo -e "\n   -> Outgroup successfully added."
    else
        echo "   -> WARNING: No outgroup hit found."
        cp "$MASTER_FASTA" "$SCRATCHDIR/${GENE_NAME}_ready.fasta"
    fi
done

# ==============================================================================
# PHASE 2 & 3: ALIGNMENT, QC, & GENE TREE INFERENCE
# ==============================================================================
echo -e "\n=== PHASE 2 & 3: ALIGNMENT & IQ-TREE INFERENCE ==="
module purge
module add mafft iqtree/3.0.1-mpi  # Adjust iqtree module name based on your cluster

for FILE in "$SCRATCHDIR"/*_ready.fasta; do
    [ -e "$FILE" ] || continue
    BASENAME=$(basename "$FILE" _ready.fasta)

    echo ">> Aligning: $BASENAME"
    mafft --auto "$FILE" > "$ALIGN_DIR/${BASENAME}.aln" 2>/dev/null
    
    # 1. Goldilocks QC Check (Dry run IQ-TREE to parse Parsimony Informative sites)
    iqtree3-mpi -s "$ALIGN_DIR/${BASENAME}.aln" -m JC -n 0 -pre "$SCRATCHDIR/qc_${BASENAME}" -quiet
    
    PI_SITES=$(grep "Number of parsimony informative sites:" "$SCRATCHDIR/qc_${BASENAME}.iqtree" | awk '{print $6}')
    PI_SITES=${PI_SITES:-0}

    echo "   -> Parsimony Informative Sites: $PI_SITES"

    # 2. Apply Thresholds (Between 10 and 150 PI sites)
    if [ "$PI_SITES" -lt 10 ]; then
        echo "   -> SKIPPING: Too conserved (< 10 PI sites)."
        continue
    elif [ "$PI_SITES" -gt 150 ]; then
        echo "   -> SKIPPING: Too noisy/divergent (> 150 PI sites)."
        continue
    fi

    # 3. Safety Check: Outgroup Routing
    if grep -q ">OUTGROUP" "$ALIGN_DIR/${BASENAME}.aln"; then
        ROOT_FLAG="-o OUTGROUP"
    else
        echo "   -> WARNING: No OUTGROUP found. Building unrooted tree."
        ROOT_FLAG=""
    fi

    # 4. Infer Maximum Likelihood Tree with ModelFinder (MFP)
    echo "   -> Passing QC. Building Tree..."
    iqtree3-mpi -s "$ALIGN_DIR/${BASENAME}.aln" -m MFP -B 1000 -T AUTO $ROOT_FLAG --prefix "$TREE_DIR/${BASENAME}" -quiet
done

# ==============================================================================
# PHASE 4: LABEL CLEANING & IQ-TREE CONSENSUS
# ==============================================================================
echo -e "\n=== PHASE 4: TIP CLEANING & CONSENSUS TREE ==="
module purge
module add iqtree  # Assuming standard iqtree2 is available for the consensus step

# Concatenate all generated treefiles and clean the BUSCO ID suffixes from the tips
# e.g., strips "C01_12345at6789" down to "C01" or "C01_Copy1" depending on your exact format.
cat "$TREE_DIR"/*.treefile | sed -E 's/_[0-9]+at[0-9]+//g' > "$FINAL_DIR/all_genetrees_clean.trees"

echo "Calculating Majority-Rule Consensus Tree..."
iqtree2 -t "$FINAL_DIR/all_genetrees_clean.trees" -con -pre "$FINAL_DIR/Species_Consensus" -quiet

# ==============================================================================
# PHASE 5: ASTRAL SPECIES TREE INFERENCE
# ==============================================================================
echo -e "\n=== PHASE 5: ASTRAL SPECIES TREE INFERENCE ==="
module purge
module add openjdk

cd "$SCRATCHDIR"
if [ ! -f "Astral/astral.5.7.8.jar" ]; then
    echo "Downloading official ASTRAL engine..."
    wget -q https://github.com/smirarab/ASTRAL/raw/master/Astral.5.7.8.zip
    unzip -q -o Astral.5.7.8.zip
fi

echo "Running ASTRAL Coalescent Model..."
java -jar Astral/astral.5.7.8.jar \
    -i "$FINAL_DIR/all_genetrees_clean.trees" \
    -o "$FINAL_DIR/ASTRAL_Species_Tree.tre" \
    2> "$FINAL_DIR/astral_run.log"

echo "=================================================="
echo "PHYLOGENOMIC PIPELINE FINISHED!"
echo "Master ASTRAL Tree: $FINAL_DIR/ASTRAL_Species_Tree.tre"
echo "Consensus Tree:     $FINAL_DIR/Species_Consensus.contree"
echo "=================================================="

#!/bin/bash
# Automatic Sample Detection and Complete Exome Analysis
# Automatically detects FASTQ pairs and runs full pipeline
# Usage: bash auto_detect_and_analyze.sh [data_directory] [threads]

set -e

# Configuration
DATA_DIR=${1:-~/NGS/data}
THREADS=${2:-16}
WORK_DIR=~/NGS

echo "╔════════════════════════════════════════════════════════════╗"
echo "║     AUTOMATIC EXOME ANALYSIS - SAMPLE AUTO-DETECTION      ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Data Directory: $DATA_DIR"
echo "  Threads: $THREADS"
echo "  Work Directory: $WORK_DIR"
echo ""

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "❌ Data directory not found: $DATA_DIR"
    echo ""
    echo "Please create directory and add FASTQ files:"
    echo "  mkdir -p $DATA_DIR"
    echo "  # Upload your *_R1.fastq.gz and *_R2.fastq.gz files"
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "STEP 1: Detecting FASTQ Samples"
echo "════════════════════════════════════════════════════════════"
echo ""

# Find all R1 files (various naming patterns)
cd $DATA_DIR
R1_FILES=$(find . -maxdepth 1 -type f \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" -o -name "*R1*.fastq.gz" \) | sort)

if [ -z "$R1_FILES" ]; then
    echo "❌ No FASTQ files found in $DATA_DIR"
    echo ""
    echo "Looking for files matching:"
    echo "  *_R1*.fastq.gz or *_R1*.fq.gz"
    echo "  *_1.fastq.gz or *_1.fq.gz"
    echo "  *R1*.fastq.gz"
    echo ""
    echo "Current files in directory:"
    ls -lh
    exit 1
fi

# Parse sample names and find pairs
declare -A SAMPLES
for r1 in $R1_FILES; do
    # Get base name without R1/R2 suffix
    base=$(basename "$r1")
    
    # Try different patterns to find R2
    r2=""
    if [[ $base == *"_R1"* ]]; then
        r2=$(echo "$base" | sed 's/_R1/_R2/g')
    elif [[ $base == *"_1."* ]]; then
        r2=$(echo "$base" | sed 's/_1\./_2\./g')
    elif [[ $base == *"R1"* ]]; then
        r2=$(echo "$base" | sed 's/R1/R2/g')
    fi
    
    # Check if R2 exists
    if [ -n "$r2" ] && [ -f "$r2" ]; then
        # Extract sample name (remove R1/R2 and extensions)
        sample_name=$(echo "$base" | sed -E 's/[._-]*(R1|_1|_r1)[._-]*.*//' | sed 's/\.fastq\.gz$//' | sed 's/\.fq\.gz$//')
        
        # Store the pair
        SAMPLES["$sample_name"]="$r1,$r2"
        echo "✅ Found pair: $sample_name"
        echo "   R1: $r1"
        echo "   R2: $r2"
        echo ""
    else
        echo "⚠️  Warning: No R2 pair found for $r1"
        echo ""
    fi
done

# Count samples
SAMPLE_COUNT=${#SAMPLES[@]}

if [ $SAMPLE_COUNT -eq 0 ]; then
    echo "❌ No valid FASTQ pairs found"
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "Summary: Found $SAMPLE_COUNT sample(s)"
echo "════════════════════════════════════════════════════════════"
echo ""

# Ask for confirmation
echo "Samples to be analyzed:"
for sample in "${!SAMPLES[@]}"; do
    echo "  • $sample"
done
echo ""
read -p "Continue with analysis? (y/n): " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Analysis cancelled."
    exit 0
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "STEP 2: Running Complete Analysis Pipeline"
echo "════════════════════════════════════════════════════════════"
echo ""

# Process each sample
counter=1
for sample in "${!SAMPLES[@]}"; do
    IFS=',' read -r r1 r2 <<< "${SAMPLES[$sample]}"
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Processing Sample $counter of $SAMPLE_COUNT: $sample"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "R1: $r1"
    echo "R2: $r2"
    echo ""
    
    # Get absolute paths
    R1_FULL="$DATA_DIR/$r1"
    R2_FULL="$DATA_DIR/$r2"
    
    # Check if already analyzed
    if [ -f "$WORK_DIR/results/$sample/annovar/separated_by_type/SNPs.txt" ]; then
        echo "⚠️  Sample $sample already analyzed. Skipping..."
        echo "   (Delete results/$sample/ to re-analyze)"
        echo ""
        counter=$((counter + 1))
        continue
    fi
    
    # Run complete pipeline
    if [ -f "$WORK_DIR/run_complete_analysis_smart.sh" ]; then
        # Use smart pipeline (with parallel ANNOVAR install)
        bash $WORK_DIR/run_complete_analysis_smart.sh "$R1_FULL" "$R2_FULL" "$sample" $THREADS
    elif [ -f "$WORK_DIR/run_complete_analysis.sh" ]; then
        # Use standard pipeline
        bash $WORK_DIR/run_complete_analysis.sh "$R1_FULL" "$R2_FULL" "$sample" $THREADS
    else
        echo "❌ Pipeline script not found!"
        exit 1
    fi
    
    # Annotate if not already done
    if [ ! -f "$WORK_DIR/results/$sample/annovar/annotated_${sample}_clinical_plus.hg19_multianno.txt" ]; then
        if [ -f "$WORK_DIR/annotate_clinical_plus.sh" ]; then
            echo ""
            echo "Annotating with clinical databases..."
            bash $WORK_DIR/annotate_clinical_plus.sh "$sample"
        fi
    fi
    
    # Separate by variant type
    if [ -f "$WORK_DIR/separate_by_variant_type.sh" ]; then
        echo ""
        echo "Separating variants by type..."
        bash $WORK_DIR/separate_by_variant_type.sh "$sample"
    fi
    
    echo ""
    echo "✅ Sample $sample complete!"
    echo ""
    
    counter=$((counter + 1))
done

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              ✅ ALL ANALYSES COMPLETE!                     ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Completed: $(date)"
echo ""

# Generate summary report
echo "════════════════════════════════════════════════════════════"
echo "SUMMARY REPORT"
echo "════════════════════════════════════════════════════════════"
echo ""

for sample in "${!SAMPLES[@]}"; do
    echo "Sample: $sample"
    echo "────────────────────────────────────────────────────────────"
    
    RESULT_DIR="$WORK_DIR/results/$sample"
    
    if [ -d "$RESULT_DIR" ]; then
        # Check for key files
        if [ -f "$RESULT_DIR/filtered/filtered_variants.vcf" ]; then
            VARIANTS=$(grep -v "^#" "$RESULT_DIR/filtered/filtered_variants.vcf" | wc -l)
            PASS_VARIANTS=$(grep -v "^#" "$RESULT_DIR/filtered/filtered_variants.vcf" | grep -w "PASS" | wc -l)
            echo "  Variants called: $VARIANTS"
            echo "  PASS variants: $PASS_VARIANTS"
        fi
        
        # Check annotation
        ANNOT_FILE=$(find "$RESULT_DIR/annovar" -name "*_clinical_plus.hg19_multianno.txt" -o -name "*_PASS.hg19_multianno.txt" 2>/dev/null | head -1)
        if [ -f "$ANNOT_FILE" ]; then
            ANNOT_SIZE=$(du -h "$ANNOT_FILE" | cut -f1)
            echo "  Annotated: Yes ($ANNOT_SIZE)"
            
            # Count pathogenic
            PATHOGENIC=$(grep -i "pathogenic" "$ANNOT_FILE" 2>/dev/null | wc -l)
            echo "  Pathogenic variants: $PATHOGENIC"
        fi
        
        # Check separation
        if [ -d "$RESULT_DIR/annovar/separated_by_type" ]; then
            echo "  Separated by type: Yes"
            if [ -f "$RESULT_DIR/annovar/separated_by_type/SNPs.txt" ]; then
                SNPS=$(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/SNPs.txt") - 1))
                echo "    SNPs: $SNPS"
            fi
            if [ -f "$RESULT_DIR/annovar/separated_by_type/Insertions.txt" ]; then
                INS=$(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/Insertions.txt") - 1))
                echo "    Insertions: $INS"
            fi
            if [ -f "$RESULT_DIR/annovar/separated_by_type/Deletions.txt" ]; then
                DEL=$(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/Deletions.txt") - 1))
                echo "    Deletions: $DEL"
            fi
        fi
        
        # Storage
        DIR_SIZE=$(du -sh "$RESULT_DIR" | cut -f1)
        echo "  Storage used: $DIR_SIZE"
    else
        echo "  Status: ❌ Failed or not found"
    fi
    
    echo ""
done

echo "════════════════════════════════════════════════════════════"
echo "Results Location: $WORK_DIR/results/"
echo ""
echo "Download files:"
echo "  • results/SAMPLE/annovar/*.hg19_multianno.txt (Excel)"
echo "  • results/SAMPLE/annovar/separated_by_type/*.txt (By type)"
echo ""
echo "Total storage used:"
du -sh $WORK_DIR/results
echo "════════════════════════════════════════════════════════════"


#!/bin/bash
# Smart NGS Analysis Pipeline with Parallel ANNOVAR Installation
# Usage: bash run_complete_analysis_smart.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> [threads]

set -e

if [ $# -lt 3 ]; then
    echo "Usage: bash run_complete_analysis_smart.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> [threads]"
    echo ""
    echo "Example:"
    echo "  bash run_complete_analysis_smart.sh data/R1.fastq.gz data/R2.fastq.gz patient_001 16"
    echo ""
    exit 1
fi

# Configuration
R1_FILE=$1
R2_FILE=$2
OUTPUT_NAME=$3
THREADS=${4:-8}

WORK_DIR=~/NGS
REFERENCE=$WORK_DIR/reference/hg19.fa
ANNOVAR_DIR=$WORK_DIR/tools/annovar
OUTPUT_DIR=$WORK_DIR/results/$OUTPUT_NAME
GATK=/opt/gatk-4.6.2.0/gatk

ANNOVAR_LOCK=~/NGS/.annovar_installing
ANNOVAR_SUCCESS=~/NGS/.annovar_installed

# Create output directories
mkdir -p $OUTPUT_DIR/{fastqc,trimmed,aligned,sorted,dedup,variants,filtered,annovar}

# Log file
LOG=$OUTPUT_DIR/pipeline.log
exec > >(tee -a $LOG) 2>&1

echo "╔════════════════════════════════════════════════════════════╗"
echo "║      NGS EXOME ANALYSIS - SMART PIPELINE WITH ANNOVAR      ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Sample Name: $OUTPUT_NAME"
echo "  R1 File:     $R1_FILE"
echo "  R2 File:     $R2_FILE"
echo "  Threads:     $THREADS"
echo ""

# Check ANNOVAR status
ANNOVAR_READY=false
ANNOVAR_INSTALLING=false

if [ -f "$ANNOVAR_SUCCESS" ] && [ -d "$ANNOVAR_DIR" ]; then
    echo "✅ ANNOVAR: Already installed"
    ANNOVAR_READY=true
elif [ -f "$ANNOVAR_LOCK" ]; then
    echo "⏳ ANNOVAR: Installation in progress (background)"
    ANNOVAR_INSTALLING=true
else
    echo "⚙️  ANNOVAR: Starting installation in background..."
    nohup bash install_annovar_background.sh > $OUTPUT_DIR/annovar_install.log 2>&1 &
    ANNOVAR_PID=$!
    echo "   Installation PID: $ANNOVAR_PID"
    echo "   Log: $OUTPUT_DIR/annovar_install.log"
    ANNOVAR_INSTALLING=true
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "Starting Main Pipeline"
echo "════════════════════════════════════════════════════════════"
echo ""

# Function to print step
step() {
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "STEP $1: $2"
    echo "Time: $(date)"
    echo "════════════════════════════════════════════════════════════"
}

# Check input files
if [ ! -f "$R1_FILE" ]; then
    echo "❌ ERROR: R1 file not found: $R1_FILE"
    exit 1
fi

if [ ! -f "$R2_FILE" ]; then
    echo "❌ ERROR: R2 file not found: $R2_FILE"
    exit 1
fi

# STEP 1: Quality Control
step 1 "Quality Control (FastQC)"
fastqc -t $THREADS -o $OUTPUT_DIR/fastqc $R1_FILE $R2_FILE
echo "✅ FastQC complete"

# STEP 2: Trimming
step 2 "Read Trimming (fastp)"
fastp \
    -i $R1_FILE \
    -I $R2_FILE \
    -o $OUTPUT_DIR/trimmed/r1_trimmed.fastq.gz \
    -O $OUTPUT_DIR/trimmed/r2_trimmed.fastq.gz \
    -h $OUTPUT_DIR/trimmed/fastp_report.html \
    -j $OUTPUT_DIR/trimmed/fastp_report.json \
    --thread $THREADS \
    --detect_adapter_for_pe \
    --length_required 50 \
    --qualified_quality_phred 20
echo "✅ Trimming complete"

# STEP 3: Alignment
step 3 "Read Alignment (BWA-MEM)"
bwa mem \
    -t $THREADS \
    -R "@RG\tID:${OUTPUT_NAME}\tSM:${OUTPUT_NAME}\tPL:ILLUMINA\tLB:lib_${OUTPUT_NAME}\tPU:unit1" \
    $REFERENCE \
    $OUTPUT_DIR/trimmed/r1_trimmed.fastq.gz \
    $OUTPUT_DIR/trimmed/r2_trimmed.fastq.gz \
    > $OUTPUT_DIR/aligned/aligned.sam
echo "✅ Alignment complete"

# STEP 4: SAM to BAM
step 4 "SAM to BAM Conversion"
samtools view -@ $THREADS -bS $OUTPUT_DIR/aligned/aligned.sam > $OUTPUT_DIR/aligned/aligned.bam
rm $OUTPUT_DIR/aligned/aligned.sam
echo "✅ Conversion complete"

# STEP 5: Sort BAM
step 5 "BAM Sorting"
samtools sort -@ $THREADS -o $OUTPUT_DIR/sorted/sorted.bam $OUTPUT_DIR/aligned/aligned.bam
rm $OUTPUT_DIR/aligned/aligned.bam
samtools index $OUTPUT_DIR/sorted/sorted.bam
echo "✅ Sorting complete"

# STEP 6: Mark Duplicates
step 6 "Mark Duplicates (GATK)"
$GATK MarkDuplicates \
    -I $OUTPUT_DIR/sorted/sorted.bam \
    -O $OUTPUT_DIR/dedup/dedup.bam \
    -M $OUTPUT_DIR/dedup/metrics.txt \
    --CREATE_INDEX true
echo "✅ Duplicates marked"

# STEP 7: Variant Calling
step 7 "Variant Calling (GATK HaplotypeCaller)"
$GATK HaplotypeCaller \
    -R $REFERENCE \
    -I $OUTPUT_DIR/dedup/dedup.bam \
    -O $OUTPUT_DIR/variants/raw_variants.vcf \
    --native-pair-hmm-threads $THREADS
echo "✅ Variant calling complete"

RAW_VARIANTS=$(grep -v "^#" $OUTPUT_DIR/variants/raw_variants.vcf | wc -l)
echo "   Total variants called: $RAW_VARIANTS"

# STEP 8: Variant Filtering
step 8 "Variant Filtering"
$GATK VariantFiltration \
    -R $REFERENCE \
    -V $OUTPUT_DIR/variants/raw_variants.vcf \
    -O $OUTPUT_DIR/filtered/filtered_variants.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"
echo "✅ Filtering complete"

FILTERED_VARIANTS=$(grep -v "^#" $OUTPUT_DIR/filtered/filtered_variants.vcf | grep "PASS" | wc -l)
echo "   Variants passing filters: $FILTERED_VARIANTS"

# Compress filtered VCF
bgzip -c $OUTPUT_DIR/filtered/filtered_variants.vcf > $OUTPUT_DIR/filtered/filtered_variants.vcf.gz
tabix -p vcf $OUTPUT_DIR/filtered/filtered_variants.vcf.gz

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║          MAIN PIPELINE COMPLETE!                           ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Pipeline time: $SECONDS seconds ($(($SECONDS / 60)) minutes)"
echo ""

# Now handle ANNOVAR annotation
step 9 "ANNOVAR Annotation"

# Check if ANNOVAR is ready
if [ -f "$ANNOVAR_SUCCESS" ] && [ -d "$ANNOVAR_DIR" ]; then
    echo "✅ ANNOVAR is ready! Starting annotation..."
    
    bash annotate_results.sh $OUTPUT_NAME
    
    echo ""
    echo "✅ Complete analysis finished with ANNOVAR annotation!"
    
elif [ -f "$ANNOVAR_LOCK" ]; then
    echo "⏳ ANNOVAR installation is still running in background"
    echo ""
    echo "You have 3 options:"
    echo ""
    echo "1. WAIT for ANNOVAR installation to complete (recommended)"
    echo "   Then annotation will run automatically"
    echo ""
    echo "2. RUN annotation later manually:"
    echo "   bash annotate_results.sh $OUTPUT_NAME"
    echo ""
    echo "3. SKIP annotation for now"
    echo ""
    read -p "Choose option (1/2/3): " choice
    
    case $choice in
        1)
            echo ""
            echo "Waiting for ANNOVAR installation..."
            echo "You can monitor progress in another terminal:"
            echo "  tail -f $OUTPUT_DIR/annovar_install.log"
            echo ""
            
            # Wait for installation (check every 30 seconds)
            while [ -f "$ANNOVAR_LOCK" ] && [ ! -f "$ANNOVAR_SUCCESS" ]; do
                echo -n "."
                sleep 30
            done
            
            echo ""
            if [ -f "$ANNOVAR_SUCCESS" ]; then
                echo "✅ ANNOVAR installation complete!"
                bash annotate_results.sh $OUTPUT_NAME
            else
                echo "⚠️  ANNOVAR installation failed or was cancelled"
                echo "Run annotation later: bash annotate_results.sh $OUTPUT_NAME"
            fi
            ;;
        2)
            echo ""
            echo "Skipping annotation for now."
            echo "Run later with: bash annotate_results.sh $OUTPUT_NAME"
            ;;
        3)
            echo ""
            echo "Skipping annotation."
            ;;
        *)
            echo ""
            echo "Invalid choice. Skipping annotation."
            echo "Run later with: bash annotate_results.sh $OUTPUT_NAME"
            ;;
    esac
else
    echo "⚠️  ANNOVAR not installed and installation not started"
    echo ""
    echo "To annotate results:"
    echo "1. Install ANNOVAR: bash setup_annovar.sh"
    echo "2. Run annotation: bash annotate_results.sh $OUTPUT_NAME"
fi

# Generate summary
echo ""
echo "════════════════════════════════════════════════════════════"
echo "ANALYSIS SUMMARY"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Sample:           $OUTPUT_NAME"
echo "Total time:       $(($SECONDS / 60)) minutes"
echo "Raw variants:     $RAW_VARIANTS"
echo "Filtered variants: $FILTERED_VARIANTS"
echo ""
echo "Results location: $OUTPUT_DIR"
echo ""
echo "Main files:"
echo "  • Filtered VCF: $OUTPUT_DIR/filtered/filtered_variants.vcf.gz"

if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt" ]; then
    echo "  • Annotated TXT: $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt"
    echo "  • Annotated VCF: $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf.gz"
fi

echo ""
echo "════════════════════════════════════════════════════════════"


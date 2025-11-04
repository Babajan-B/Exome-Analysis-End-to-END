#!/bin/bash
# NGS Exome Analysis Pipeline - Single Command Execution
# Usage: bash run_pipeline.sh <R1.fastq.gz> <R2.fastq.gz> <output_name> [threads]

set -e  # Exit on error

# Configuration
R1_FILE=$1
R2_FILE=$2
OUTPUT_NAME=${3:-"sample"}
THREADS=${4:-8}  # Default 8 threads for cloud instance

# Paths
WORK_DIR=~/NGS
REFERENCE=$WORK_DIR/reference/hg19.fa
SNPEFF_DIR=$WORK_DIR/tools/snpEff
OUTPUT_DIR=$WORK_DIR/results/$OUTPUT_NAME
GATK=/opt/gatk-4.6.2.0/gatk

# Create output directory
mkdir -p $OUTPUT_DIR/{fastqc,trimmed,aligned,sorted,dedup,bqsr,variants,filtered,annotated}

# Log file
LOG=$OUTPUT_DIR/pipeline.log
exec > >(tee -a $LOG) 2>&1

echo "=========================================="
echo "NGS Exome Analysis Pipeline"
echo "Started: $(date)"
echo "=========================================="
echo "Configuration:"
echo "  R1: $R1_FILE"
echo "  R2: $R2_FILE"
echo "  Output: $OUTPUT_NAME"
echo "  Threads: $THREADS"
echo "  Reference: $REFERENCE"
echo "=========================================="
echo ""

# Check input files
if [ ! -f "$R1_FILE" ]; then
    echo "ERROR: R1 file not found: $R1_FILE"
    exit 1
fi

if [ ! -f "$R2_FILE" ]; then
    echo "ERROR: R2 file not found: $R2_FILE"
    exit 1
fi

# Function to print step
step() {
    echo ""
    echo "=========================================="
    echo "STEP $1: $2"
    echo "Time: $(date)"
    echo "=========================================="
}

# STEP 1: Quality Check
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
    -R "@RG\tID:$OUTPUT_NAME\tSM:$OUTPUT_NAME\tPL:ILLUMINA\tLB:lib_$OUTPUT_NAME\tPU:unit1" \
    $REFERENCE \
    $OUTPUT_DIR/trimmed/r1_trimmed.fastq.gz \
    $OUTPUT_DIR/trimmed/r2_trimmed.fastq.gz \
    > $OUTPUT_DIR/aligned/aligned.sam
echo "✅ Alignment complete"

# STEP 4: SAM to BAM conversion
step 4 "SAM to BAM Conversion"
samtools view -@ $THREADS -bS $OUTPUT_DIR/aligned/aligned.sam > $OUTPUT_DIR/aligned/aligned.bam
rm $OUTPUT_DIR/aligned/aligned.sam  # Remove large SAM file
echo "✅ SAM to BAM conversion complete"

# STEP 5: Sort BAM
step 5 "BAM Sorting"
samtools sort -@ $THREADS -o $OUTPUT_DIR/sorted/sorted.bam $OUTPUT_DIR/aligned/aligned.bam
rm $OUTPUT_DIR/aligned/aligned.bam  # Remove unsorted BAM
echo "✅ Sorting complete"

# STEP 6: Mark Duplicates
step 6 "Mark Duplicates (GATK)"
$GATK MarkDuplicates \
    -I $OUTPUT_DIR/sorted/sorted.bam \
    -O $OUTPUT_DIR/dedup/dedup.bam \
    -M $OUTPUT_DIR/dedup/metrics.txt \
    --CREATE_INDEX true
rm $OUTPUT_DIR/sorted/sorted.bam  # Remove pre-dedup BAM
echo "✅ Duplicate marking complete"

# STEP 7: Base Quality Score Recalibration (optional - requires known sites)
# Skipping BQSR if no known sites available
echo ""
echo "⏭️  Skipping BQSR (no known sites configured)"

# STEP 8: Variant Calling
step 7 "Variant Calling (GATK HaplotypeCaller)"
$GATK HaplotypeCaller \
    -R $REFERENCE \
    -I $OUTPUT_DIR/dedup/dedup.bam \
    -O $OUTPUT_DIR/variants/raw_variants.vcf \
    --native-pair-hmm-threads $THREADS
echo "✅ Variant calling complete"

# STEP 9: Variant Filtering
step 8 "Variant Filtering"
$GATK VariantFiltration \
    -R $REFERENCE \
    -V $OUTPUT_DIR/variants/raw_variants.vcf \
    -O $OUTPUT_DIR/filtered/filtered_variants.vcf \
    --filter-name "QD_filter" --filter-expression "QD < 2.0" \
    --filter-name "FS_filter" --filter-expression "FS > 60.0" \
    --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
    --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
    --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0"
echo "✅ Variant filtering complete"

# STEP 10: Variant Annotation
step 9 "Variant Annotation (snpEff)"
cd $SNPEFF_DIR
java -Xmx8g -jar snpEff.jar \
    GRCh37.75 \
    $OUTPUT_DIR/filtered/filtered_variants.vcf \
    > $OUTPUT_DIR/annotated/annotated_variants.vcf
cp snpEff_summary.html $OUTPUT_DIR/annotated/
echo "✅ Annotation complete"

# STEP 11: Compress and Index VCF
step 10 "VCF Compression and Indexing"
bgzip -c $OUTPUT_DIR/annotated/annotated_variants.vcf > $OUTPUT_DIR/annotated/annotated_variants.vcf.gz
tabix -p vcf $OUTPUT_DIR/annotated/annotated_variants.vcf.gz
echo "✅ Compression complete"

# Generate summary
step 11 "Generating Summary"
cd $OUTPUT_DIR

cat > summary.txt << EOF
========================================
NGS Exome Analysis Summary
========================================
Analysis: $OUTPUT_NAME
Completed: $(date)
Reference: hg19
Threads: $THREADS

File Locations:
  ├── FastQC Reports:    $OUTPUT_DIR/fastqc/
  ├── Trimming Report:   $OUTPUT_DIR/trimmed/fastp_report.html
  ├── Final BAM:         $OUTPUT_DIR/dedup/dedup.bam
  ├── Raw Variants:      $OUTPUT_DIR/variants/raw_variants.vcf
  ├── Filtered VCF:      $OUTPUT_DIR/filtered/filtered_variants.vcf
  └── Annotated VCF:     $OUTPUT_DIR/annotated/annotated_variants.vcf.gz

Variant Statistics:
  Total variants: $(grep -v "^#" $OUTPUT_DIR/variants/raw_variants.vcf | wc -l)
  Filtered variants: $(grep -v "^#" $OUTPUT_DIR/filtered/filtered_variants.vcf | grep PASS | wc -l)
  
Disk Usage:
$(du -sh $OUTPUT_DIR)

========================================
EOF

cat summary.txt

echo ""
echo "=========================================="
echo "✅ PIPELINE COMPLETE!"
echo "=========================================="
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo "Total time: $SECONDS seconds"
echo ""
echo "Next steps:"
echo "  1. Review summary: cat $OUTPUT_DIR/summary.txt"
echo "  2. View FastQC: firefox $OUTPUT_DIR/fastqc/*.html"
echo "  3. Download VCF: $OUTPUT_DIR/annotated/annotated_variants.vcf.gz"
echo ""
echo "=========================================="


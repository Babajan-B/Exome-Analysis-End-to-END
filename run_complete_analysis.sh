#!/bin/bash
# Complete NGS Analysis Pipeline with ANNOVAR Annotation
# Single command to run: QC → Trim → Align → Variant Call → Annotate
# Usage: bash run_complete_analysis.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> [threads]

set -e  # Exit on error

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: bash run_complete_analysis.sh <R1.fastq.gz> <R2.fastq.gz> <sample_name> [threads]"
    echo ""
    echo "Example:"
    echo "  bash run_complete_analysis.sh data/R1.fastq.gz data/R2.fastq.gz patient_001 16"
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
SNPEFF_DIR=$WORK_DIR/tools/snpEff
ANNOVAR_DIR=$WORK_DIR/tools/annovar
OUTPUT_DIR=$WORK_DIR/results/$OUTPUT_NAME
GATK=/opt/gatk-4.6.2.0/gatk

# Create output directories
mkdir -p $OUTPUT_DIR/{fastqc,trimmed,aligned,sorted,dedup,bqsr,variants,filtered,annotated,annovar}

# Log file
LOG=$OUTPUT_DIR/pipeline.log
exec > >(tee -a $LOG) 2>&1

echo "╔════════════════════════════════════════════════════════════╗"
echo "║         NGS EXOME ANALYSIS - COMPLETE PIPELINE             ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Started: $(date)"
echo ""
echo "Configuration:"
echo "  Sample Name: $OUTPUT_NAME"
echo "  R1 File:     $R1_FILE"
echo "  R2 File:     $R2_FILE"
echo "  Threads:     $THREADS"
echo "  Reference:   $REFERENCE"
echo "  Output Dir:  $OUTPUT_DIR"
echo ""
echo "Pipeline Steps:"
echo "  1. Quality Control (FastQC)"
echo "  2. Read Trimming (fastp)"
echo "  3. Alignment (BWA-MEM)"
echo "  4. SAM to BAM Conversion"
echo "  5. BAM Sorting"
echo "  6. Mark Duplicates (GATK)"
echo "  7. Base Quality Recalibration (BQSR)"
echo "  8. Variant Calling (GATK HaplotypeCaller)"
echo "  9. Variant Filtering"
echo "  10. Variant Annotation (ANNOVAR)"
echo "  11. VCF Compression"
echo ""
echo "════════════════════════════════════════════════════════════"
echo ""

# Check input files
if [ ! -f "$R1_FILE" ]; then
    echo "❌ ERROR: R1 file not found: $R1_FILE"
    exit 1
fi

if [ ! -f "$R2_FILE" ]; then
    echo "❌ ERROR: R2 file not found: $R2_FILE"
    exit 1
fi

# Function to print step
step() {
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "STEP $1: $2"
    echo "Time: $(date)"
    echo "════════════════════════════════════════════════════════════"
}

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

# STEP 7: Base Quality Recalibration (skip if no known sites)
step 7 "Base Quality Score Recalibration (Skipped - optional)"
echo "⏭️  BQSR skipped (optional step)"
FINAL_BAM=$OUTPUT_DIR/dedup/dedup.bam

# STEP 8: Variant Calling
step 8 "Variant Calling (GATK HaplotypeCaller)"
$GATK HaplotypeCaller \
    -R $REFERENCE \
    -I $FINAL_BAM \
    -O $OUTPUT_DIR/variants/raw_variants.vcf \
    --native-pair-hmm-threads $THREADS
echo "✅ Variant calling complete"

RAW_VARIANTS=$(grep -v "^#" $OUTPUT_DIR/variants/raw_variants.vcf | wc -l)
echo "   Total variants called: $RAW_VARIANTS"

# STEP 9: Variant Filtering
step 9 "Variant Filtering"
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

# STEP 10: ANNOVAR Annotation
step 10 "Variant Annotation (ANNOVAR)"

if [ -d "$ANNOVAR_DIR" ] && [ -f "$ANNOVAR_DIR/table_annovar.pl" ]; then
    echo "Running ANNOVAR annotation with comprehensive databases..."
    
    perl $ANNOVAR_DIR/table_annovar.pl \
        $OUTPUT_DIR/filtered/filtered_variants.vcf \
        $ANNOVAR_DIR/humandb/ \
        -buildver hg19 \
        -out $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME} \
        -remove \
        -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70,icgc28 \
        -operation g,g,g,f,f,f,f,f,f \
        -nastring . \
        -vcfinput \
        -polish \
        2>&1 | grep -E "(NOTICE|done)" || true
    
    if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf" ]; then
        echo "✅ ANNOVAR annotation complete"
        
        ANNOTATED_VARIANTS=$(grep -v "^#" $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf | wc -l)
        echo "   Variants annotated: $ANNOTATED_VARIANTS"
        
        # Check for pathogenic variants
        if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt" ]; then
            PATHOGENIC=$(grep -i "pathogenic" $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt | wc -l)
            echo "   Pathogenic/Likely pathogenic variants: $PATHOGENIC"
        fi
    else
        echo "⚠️  ANNOVAR annotation failed or incomplete"
    fi
else
    echo "⚠️  ANNOVAR not found - skipping annotation"
    echo "   Install with: bash setup_annovar.sh"
fi

# STEP 11: VCF Compression
step 11 "VCF Compression and Indexing"

if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf" ]; then
    # Compress annotated VCF
    bgzip -c $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf > \
        $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf.gz
    tabix -p vcf $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf.gz
    echo "✅ Annotated VCF compressed"
fi

# Also compress filtered VCF
bgzip -c $OUTPUT_DIR/filtered/filtered_variants.vcf > $OUTPUT_DIR/filtered/filtered_variants.vcf.gz
tabix -p vcf $OUTPUT_DIR/filtered/filtered_variants.vcf.gz
echo "✅ Filtered VCF compressed"

# Generate summary
step 12 "Generating Analysis Summary"

cat > $OUTPUT_DIR/ANALYSIS_SUMMARY.txt << EOF
╔════════════════════════════════════════════════════════════╗
║           NGS EXOME ANALYSIS - SUMMARY REPORT              ║
╚════════════════════════════════════════════════════════════╝

Sample Name:      $OUTPUT_NAME
Analysis Date:    $(date)
Reference Genome: hg19
Processing Time:  $SECONDS seconds

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

VARIANT STATISTICS:

  Raw Variants Called:          $RAW_VARIANTS
  Variants Passing Filters:     $FILTERED_VARIANTS
  Filtering Rate:               $(awk "BEGIN {printf \"%.1f\", ($FILTERED_VARIANTS/$RAW_VARIANTS)*100}")%

EOF

if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt" ]; then
    PATHOGENIC=$(grep -i "pathogenic" $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt 2>/dev/null | wc -l)
    cat >> $OUTPUT_DIR/ANALYSIS_SUMMARY.txt << EOF
  Annotated Variants:           $ANNOTATED_VARIANTS
  Pathogenic Variants Found:    $PATHOGENIC

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ANNOTATION DATABASES USED:

  ✓ RefGene, UCSC, Ensembl      (Gene annotations)
  ✓ dbSNP 150                    (Variant IDs)
  ✓ gnomAD v2.1.1                (Population frequencies)
  ✓ ClinVar 2024                 (Clinical significance)
  ✓ dbNSFP 4.2a                  (Functional predictions)
  ✓ COSMIC 70                    (Cancer mutations)
  ✓ ICGC 28                      (Cancer somatic)

EOF
fi

cat >> $OUTPUT_DIR/ANALYSIS_SUMMARY.txt << EOF
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

OUTPUT FILES:

  Quality Control:
    • $OUTPUT_DIR/fastqc/

  Alignment Files:
    • $OUTPUT_DIR/dedup/dedup.bam
    • $OUTPUT_DIR/dedup/dedup.bai

  Variant Files:
    • $OUTPUT_DIR/variants/raw_variants.vcf
    • $OUTPUT_DIR/filtered/filtered_variants.vcf.gz

EOF

if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf" ]; then
    cat >> $OUTPUT_DIR/ANALYSIS_SUMMARY.txt << EOF
  Annotated Variants (ANNOVAR):
    • $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf.gz
    • $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt

EOF
fi

cat >> $OUTPUT_DIR/ANALYSIS_SUMMARY.txt << EOF
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

NEXT STEPS:

  1. Download annotated results:
     • VCF file: annotated_${OUTPUT_NAME}.hg19_multianno.vcf.gz
     • TXT file: annotated_${OUTPUT_NAME}.hg19_multianno.txt

  2. Open TXT file in Excel and filter for:
     • CLNSIG column contains "Pathogenic"
     • gnomAD_exome_ALL < 0.01 (rare variants)
     • SIFT_score < 0.05 (damaging predictions)

  3. Review pathogenic variants and validate with literature

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STORAGE USAGE:

$(du -sh $OUTPUT_DIR)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EOF

cat $OUTPUT_DIR/ANALYSIS_SUMMARY.txt

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              ✅ ANALYSIS COMPLETE!                         ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Total time: $SECONDS seconds ($(($SECONDS / 60)) minutes)"
echo ""
echo "Results location: $OUTPUT_DIR"
echo ""
echo "Key files to download:"
if [ -f "$OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt" ]; then
    echo "  1. $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.txt (Excel)"
    echo "  2. $OUTPUT_DIR/annovar/annotated_${OUTPUT_NAME}.hg19_multianno.vcf.gz"
else
    echo "  1. $OUTPUT_DIR/filtered/filtered_variants.vcf.gz"
fi
echo "  3. $OUTPUT_DIR/ANALYSIS_SUMMARY.txt"
echo "  4. $OUTPUT_DIR/fastqc/ (QC reports)"
echo ""
echo "════════════════════════════════════════════════════════════"


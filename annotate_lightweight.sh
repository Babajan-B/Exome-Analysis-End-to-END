#!/bin/bash
# Lightweight ANNOVAR annotation - Essential databases only
# Uses only 3 databases instead of 9 - much smaller files!
# Usage: bash annotate_lightweight.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash annotate_lightweight.sh <sample_name>"
    echo ""
    echo "This uses ONLY essential databases:"
    echo "  • refGene          - Gene annotations"
    echo "  • clinvar_20240917 - Clinical significance"
    echo "  • gnomad211_exome  - Population frequencies"
    echo ""
    echo "Result: 10-20MB files (much smaller!)"
    exit 1
fi

SAMPLE_NAME=$1
ANNOVAR_DIR=~/NGS/tools/annovar
OUTPUT_DIR=~/NGS/results/$SAMPLE_NAME
VCF_FILE=$OUTPUT_DIR/filtered/filtered_variants.vcf

# Check files
if [ ! -f "$VCF_FILE" ]; then
    echo "❌ Filtered VCF not found: $VCF_FILE"
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "Lightweight ANNOVAR Annotation (Essential Databases Only)"
echo "Sample: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""

# Create PASS-only VCF
PASS_VCF=$OUTPUT_DIR/filtered/filtered_PASS_only.vcf

echo "Step 1: Extracting PASS variants..."
grep "^#" $VCF_FILE > $PASS_VCF
grep -v "^#" $VCF_FILE | grep -w "PASS" >> $PASS_VCF

PASS_COUNT=$(grep -v "^#" $PASS_VCF | wc -l)
echo "  PASS variants: $PASS_COUNT"
echo ""

# Annotate with minimal databases
echo "Step 2: Annotating with essential databases only..."
echo "  Using: refGene, clinvar_20240917, gnomad211_exome"
echo ""

mkdir -p $OUTPUT_DIR/annovar

perl $ANNOVAR_DIR/table_annovar.pl \
    $PASS_VCF \
    $ANNOVAR_DIR/humandb/ \
    -buildver hg19 \
    -out $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight \
    -remove \
    -protocol refGene,clinvar_20240917,gnomad211_exome \
    -operation g,f,f \
    -nastring . \
    -vcfinput \
    -polish

if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.vcf" ]; then
    # Compress
    bgzip -f $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.vcf
    tabix -p vcf $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.vcf.gz
    
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "✅ Lightweight Annotation Complete!"
    echo "════════════════════════════════════════════════════════════"
    echo ""
    
    VCF_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.vcf.gz | cut -f1)
    TXT_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.txt | cut -f1)
    
    echo "Output files:"
    echo "  VCF: $VCF_SIZE (compressed)"
    echo "  TXT: $TXT_SIZE"
    echo "  Variants: $PASS_COUNT"
    echo ""
    
    # Count columns
    COLS=$(head -1 $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.txt | tr '\t' '\n' | wc -l)
    echo "  Columns: $COLS (vs ~150 with all databases)"
    echo ""
    
    # Count pathogenic
    if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.txt" ]; then
        PATHOGENIC=$(grep -i "pathogenic" $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.txt 2>/dev/null | wc -l)
        echo "  Pathogenic variants: $PATHOGENIC"
    fi
    
    echo ""
    echo "Key columns included:"
    echo "  From refGene:"
    echo "    • Gene name"
    echo "    • Function (exonic, intronic, etc.)"
    echo "    • Exonic function (nonsynonymous, frameshift, etc.)"
    echo "    • Amino acid change"
    echo ""
    echo "  From ClinVar:"
    echo "    • Clinical significance (Pathogenic, Benign, etc.)"
    echo "    • Review status"
    echo "    • Disease associations"
    echo ""
    echo "  From gnomAD:"
    echo "    • Allele frequency in general population"
    echo "    • Frequency by ethnicity"
    echo ""
    echo "Download: $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_lightweight.hg19_multianno.txt"
    echo ""
    echo "════════════════════════════════════════════════════════════"
else
    echo "❌ Annotation failed"
    exit 1
fi


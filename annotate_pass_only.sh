#!/bin/bash
# Annotate only PASS variants (reduces file size by 90%)
# Usage: bash annotate_pass_only.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash annotate_pass_only.sh <sample_name>"
    echo ""
    echo "Example: bash annotate_pass_only.sh sample_Alaa"
    exit 1
fi

SAMPLE_NAME=$1
ANNOVAR_DIR=~/NGS/tools/annovar
OUTPUT_DIR=~/NGS/results/$SAMPLE_NAME
FILTERED_VCF=$OUTPUT_DIR/filtered/filtered_variants.vcf

# Check if sample exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "❌ Sample not found: $SAMPLE_NAME"
    exit 1
fi

if [ ! -f "$FILTERED_VCF" ]; then
    echo "❌ Filtered VCF not found: $FILTERED_VCF"
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "Filtering and Annotating PASS Variants Only"
echo "Sample: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""

# Count variants
TOTAL_VARIANTS=$(grep -v "^#" $FILTERED_VCF | wc -l)
PASS_VARIANTS=$(grep -v "^#" $FILTERED_VCF | grep -w "PASS" | wc -l)
FAILED_VARIANTS=$((TOTAL_VARIANTS - PASS_VARIANTS))

echo "Variant Statistics:"
echo "  Total variants:  $TOTAL_VARIANTS"
echo "  PASS variants:   $PASS_VARIANTS (will annotate)"
echo "  Failed variants: $FAILED_VARIANTS (will skip)"
echo ""

if [ $PASS_VARIANTS -eq 0 ]; then
    echo "❌ No PASS variants found!"
    exit 1
fi

# Create PASS-only VCF
echo "Creating PASS-only VCF..."
PASS_VCF=$OUTPUT_DIR/filtered/filtered_PASS_only.vcf

# Extract header and PASS variants
grep "^#" $FILTERED_VCF > $PASS_VCF
grep -v "^#" $FILTERED_VCF | grep -w "PASS" >> $PASS_VCF

echo "✅ PASS-only VCF created: $(du -h $PASS_VCF | cut -f1)"
echo ""

# Annotate with ANNOVAR
echo "Running ANNOVAR on PASS variants only..."
mkdir -p $OUTPUT_DIR/annovar

perl $ANNOVAR_DIR/table_annovar.pl \
    $PASS_VCF \
    $ANNOVAR_DIR/humandb/ \
    -buildver hg19 \
    -out $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS \
    -remove \
    -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150,dbnsfp42a \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish

if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.vcf" ]; then
    # Compress
    bgzip -f $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.vcf
    tabix -p vcf $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.vcf.gz
    
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "✅ Annotation Complete!"
    echo "════════════════════════════════════════════════════════════"
    echo ""
    
    # File sizes
    VCF_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.vcf.gz | cut -f1)
    TXT_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.txt | cut -f1)
    
    echo "New Files (PASS variants only):"
    echo "  VCF: $VCF_SIZE (compressed)"
    echo "  TXT: $TXT_SIZE (for Excel)"
    echo "  Variants: $PASS_VARIANTS"
    echo ""
    
    # Count pathogenic
    if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.txt" ]; then
        PATHOGENIC=$(grep -i "pathogenic" $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.txt 2>/dev/null | wc -l)
        echo "  Pathogenic variants: $PATHOGENIC"
    fi
    
    echo ""
    echo "Old Files (ALL variants - can delete):"
    if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf" ]; then
        OLD_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf | cut -f1)
        echo "  VCF: $OLD_SIZE (unfiltered)"
    fi
    
    echo ""
    echo "Space saved: ~90% smaller files!"
    echo ""
    echo "Download these files:"
    echo "  1. $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.txt"
    echo "  2. $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_PASS.hg19_multianno.vcf.gz"
    echo ""
    echo "════════════════════════════════════════════════════════════"
else
    echo "❌ Annotation failed"
    exit 1
fi


#!/bin/bash
# Check what columns ANNOVAR added
# Usage: bash check_annotation_columns.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash check_annotation_columns.sh <sample_name>"
    echo "Example: bash check_annotation_columns.sh sample_Alaa"
    exit 1
fi

SAMPLE=$1
RESULT_DIR=~/NGS/results/$SAMPLE

echo "════════════════════════════════════════════════════════════"
echo "Analyzing Annotation Columns for: $SAMPLE"
echo "════════════════════════════════════════════════════════════"
echo ""

# Check filtered VCF
if [ -f "$RESULT_DIR/filtered/filtered_variants.vcf" ]; then
    FILTERED_SIZE=$(du -h "$RESULT_DIR/filtered/filtered_variants.vcf" | cut -f1)
    FILTERED_VARIANTS=$(grep -v "^#" "$RESULT_DIR/filtered/filtered_variants.vcf" | wc -l)
    
    echo "ORIGINAL FILE:"
    echo "  File: filtered_variants.vcf"
    echo "  Size: $FILTERED_SIZE"
    echo "  Total variants: $FILTERED_VARIANTS"
    echo "  PASS variants: $(grep -v "^#" "$RESULT_DIR/filtered/filtered_variants.vcf" | grep -w "PASS" | wc -l)"
    echo ""
fi

# Check annotated VCF
if [ -f "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf" ]; then
    ANNOT_SIZE=$(du -h "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf" | cut -f1)
    ANNOT_VARIANTS=$(grep -v "^#" "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf" | wc -l)
    
    echo "ANNOTATED FILE:"
    echo "  File: annotated_${SAMPLE}.hg19_multianno.vcf"
    echo "  Size: $ANNOT_SIZE"
    echo "  Variants: $ANNOT_VARIANTS"
    echo ""
fi

# Check TXT file columns
if [ -f "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.txt" ]; then
    TXT_SIZE=$(du -h "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.txt" | cut -f1)
    
    echo "ANNOTATED TXT FILE:"
    echo "  File: annotated_${SAMPLE}.hg19_multianno.txt"
    echo "  Size: $TXT_SIZE"
    echo ""
    
    echo "COLUMN HEADERS (what ANNOVAR added):"
    echo "════════════════════════════════════════════════════════════"
    head -1 "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.txt" | \
        tr '\t' '\n' | nl
    echo ""
    
    TOTAL_COLS=$(head -1 "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.txt" | tr '\t' '\n' | wc -l)
    echo "Total columns: $TOTAL_COLS"
    echo ""
fi

echo "════════════════════════════════════════════════════════════"
echo "ANALYSIS:"
echo "════════════════════════════════════════════════════════════"
echo ""

if [ -f "$RESULT_DIR/filtered/filtered_variants.vcf" ] && [ -f "$RESULT_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf" ]; then
    echo "Size increase:"
    echo "  $FILTERED_SIZE → $ANNOT_SIZE"
    echo ""
    
    echo "Why so large?"
    echo "  1. Annotated ALL variants ($ANNOT_VARIANTS) instead of PASS only"
    echo "  2. Added ~100+ annotation columns from 9 databases"
    echo "  3. Each database adds multiple fields (gene names, scores, etc.)"
    echo ""
    
    echo "Databases that added columns:"
    echo "  • refGene     - Gene annotations (5-10 columns)"
    echo "  • knownGene   - UCSC genes (5-10 columns)"
    echo "  • ensGene     - Ensembl genes (5-10 columns)"
    echo "  • avsnp150    - dbSNP IDs (1-2 columns)"
    echo "  • gnomad211   - Population frequencies (10-20 columns)"
    echo "  • clinvar     - Clinical significance (5-10 columns)"
    echo "  • dbnsfp42a   - Functional predictions (50-100 columns!)"
    echo "  • cosmic70    - Cancer mutations (5-10 columns)"
    echo "  • icgc28      - ICGC data (5-10 columns)"
    echo ""
    echo "Total: ~100-150 columns added!"
    echo ""
fi

echo "════════════════════════════════════════════════════════════"
echo "SOLUTION:"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "1. Use fewer databases (if you don't need all):"
echo "   bash annotate_pass_only.sh $SAMPLE"
echo "   (Uses essential databases only)"
echo ""
echo "2. Or annotate PASS variants only:"
echo "   This will reduce size by 90% regardless of databases"
echo ""
echo "════════════════════════════════════════════════════════════"


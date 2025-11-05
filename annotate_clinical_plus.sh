#!/bin/bash
# Clinical Plus ANNOVAR annotation - Essential + Predictions
# Includes: Genes, ClinVar, gnomAD, dbSNP, and Prediction scores
# Usage: bash annotate_clinical_plus.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash annotate_clinical_plus.sh <sample_name>"
    echo ""
    echo "Databases included (5):"
    echo "  1. refGene          - Gene annotations"
    echo "  2. clinvar_20240917 - Clinical significance"
    echo "  3. gnomad211_exome  - Population frequencies"
    echo "  4. avsnp150         - dbSNP IDs (rs numbers)"
    echo "  5. dbnsfp42a        - Prediction scores (SIFT, PolyPhen, CADD, etc.)"
    echo ""
    echo "Result: 30-50MB files (clinical + predictions)"
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
echo "Clinical Plus ANNOVAR Annotation"
echo "Sample: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Databases:"
echo "  ✅ refGene          - Gene annotations"
echo "  ✅ clinvar_20240917 - Clinical significance"
echo "  ✅ gnomad211_exome  - Population frequencies"
echo "  ✅ avsnp150         - dbSNP IDs"
echo "  ✅ dbnsfp42a        - Prediction scores"
echo ""

# Create PASS-only VCF
PASS_VCF=$OUTPUT_DIR/filtered/filtered_PASS_only.vcf

echo "Step 1: Extracting PASS variants..."
grep "^#" $VCF_FILE > $PASS_VCF
grep -v "^#" $VCF_FILE | grep -w "PASS" >> $PASS_VCF

PASS_COUNT=$(grep -v "^#" $PASS_VCF | wc -l)
echo "  PASS variants: $PASS_COUNT"
echo ""

# Annotate with clinical + prediction databases
echo "Step 2: Annotating with 5 databases..."
echo ""

mkdir -p $OUTPUT_DIR/annovar

perl $ANNOVAR_DIR/table_annovar.pl \
    $PASS_VCF \
    $ANNOVAR_DIR/humandb/ \
    -buildver hg19 \
    -out $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus \
    -remove \
    -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150,dbnsfp42a \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish

if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.vcf" ]; then
    # Compress
    bgzip -f $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.vcf
    tabix -p vcf $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.vcf.gz
    
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "✅ Clinical Plus Annotation Complete!"
    echo "════════════════════════════════════════════════════════════"
    echo ""
    
    VCF_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.vcf.gz | cut -f1)
    TXT_SIZE=$(du -h $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt | cut -f1)
    
    echo "Output files:"
    echo "  VCF: $VCF_SIZE (compressed)"
    echo "  TXT: $TXT_SIZE"
    echo "  Variants: $PASS_COUNT"
    echo ""
    
    # Count columns
    COLS=$(head -1 $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt | tr '\t' '\n' | wc -l)
    echo "  Columns: $COLS"
    echo ""
    
    # Show column headers
    echo "Column headers included:"
    echo "════════════════════════════════════════════════════════════"
    head -1 $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt | tr '\t' '\n' | nl | head -30
    
    if [ $COLS -gt 30 ]; then
        echo "... and $((COLS - 30)) more columns"
    fi
    
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo ""
    
    # Count pathogenic
    if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt" ]; then
        PATHOGENIC=$(grep -i "pathogenic" $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt 2>/dev/null | wc -l)
        echo "  Pathogenic variants found: $PATHOGENIC"
        echo ""
    fi
    
    echo "Key information included:"
    echo ""
    echo "From refGene:"
    echo "  • Gene name, Function, Amino acid changes"
    echo ""
    echo "From ClinVar:"
    echo "  • Clinical significance (Pathogenic/Benign)"
    echo "  • Disease associations"
    echo ""
    echo "From gnomAD:"
    echo "  • Population frequencies (overall + by ethnicity)"
    echo ""
    echo "From dbSNP (avsnp150):"
    echo "  • rs numbers (e.g., rs80357906)"
    echo "  • Variant IDs for literature search"
    echo ""
    echo "From dbNSFP (prediction scores):"
    echo "  • SIFT_score (< 0.05 = damaging)"
    echo "  • Polyphen2_HDIV_score (> 0.85 = damaging)"
    echo "  • CADD_phred (> 20 = deleterious)"
    echo "  • MetaSVM_pred (D = damaging, T = tolerated)"
    echo "  • REVEL_score (> 0.5 = pathogenic)"
    echo "  • And 40+ more prediction algorithms"
    echo ""
    echo "Download file:"
    echo "  $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt"
    echo ""
    echo "════════════════════════════════════════════════════════════"
else
    echo "❌ Annotation failed"
    exit 1
fi


#!/bin/bash
# Annotate existing pipeline results with ANNOVAR
# Usage: bash annotate_results.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash annotate_results.sh <sample_name>"
    echo ""
    echo "Example: bash annotate_results.sh patient_001"
    exit 1
fi

SAMPLE_NAME=$1
ANNOVAR_DIR=~/NGS/tools/annovar
OUTPUT_DIR=~/NGS/results/$SAMPLE_NAME
VCF_FILE=$OUTPUT_DIR/filtered/filtered_variants.vcf

# Check if sample exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "❌ Sample not found: $SAMPLE_NAME"
    echo "Available samples:"
    ls -1 ~/NGS/results/ 2>/dev/null | head -10
    exit 1
fi

# Check if VCF exists
if [ ! -f "$VCF_FILE" ]; then
    echo "❌ Filtered VCF not found: $VCF_FILE"
    echo "Pipeline may not be complete yet."
    exit 1
fi

# Check ANNOVAR installation
if [ ! -d "$ANNOVAR_DIR" ] || [ ! -f "$ANNOVAR_DIR/table_annovar.pl" ]; then
    echo "❌ ANNOVAR not installed"
    echo ""
    echo "Install with: bash setup_annovar.sh"
    echo "Or wait if background installation is running..."
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "ANNOVAR Annotation for: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""

mkdir -p $OUTPUT_DIR/annovar

# Run annotation
perl $ANNOVAR_DIR/table_annovar.pl \
    $VCF_FILE \
    $ANNOVAR_DIR/humandb/ \
    -buildver hg19 \
    -out $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME} \
    -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70,icgc28 \
    -operation g,g,g,f,f,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish

if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf" ]; then
    # Compress
    bgzip -f $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf
    tabix -p vcf $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf.gz
    
    echo ""
    echo "✅ Annotation complete!"
    echo ""
    echo "Results:"
    echo "  - $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf.gz"
    echo "  - $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.txt"
    
    # Count pathogenic
    if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.txt" ]; then
        PATHOGENIC=$(grep -i "pathogenic" $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.txt 2>/dev/null | wc -l)
        echo "  - Pathogenic variants found: $PATHOGENIC"
    fi
else
    echo "❌ Annotation failed"
    exit 1
fi


#!/bin/bash
# Complete workflow: Annotate + Separate by variant type
# Usage: bash run_complete_with_separation.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║   Complete Analysis: Annotation + Variant Type Separation ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Usage: bash run_complete_with_separation.sh <sample_name>"
    echo ""
    echo "This will:"
    echo "  1. Annotate variants with 5 databases"
    echo "  2. Separate into SNPs, Insertions, Deletions, MNPs"
    echo "  3. Generate summary statistics"
    echo ""
    echo "Example:"
    echo "  bash run_complete_with_separation.sh sample_Alaa"
    echo ""
    exit 1
fi

SAMPLE_NAME=$1

echo "╔════════════════════════════════════════════════════════════╗"
echo "║        Complete Annotation + Variant Type Separation      ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Sample: $SAMPLE_NAME"
echo "Started: $(date)"
echo ""

# Step 1: Annotate
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: ANNOVAR Annotation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

bash annotate_clinical_plus.sh $SAMPLE_NAME

if [ $? -ne 0 ]; then
    echo ""
    echo "❌ Annotation failed. Please check the error above."
    exit 1
fi

echo ""
echo "✅ Annotation complete!"
echo ""
sleep 2

# Step 2: Separate by type
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2: Separating by Variant Type"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

bash separate_by_variant_type.sh $SAMPLE_NAME

if [ $? -ne 0 ]; then
    echo ""
    echo "❌ Separation failed. Please check the error above."
    exit 1
fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                 ✅ COMPLETE SUCCESS!                       ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Completed: $(date)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Your Results:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📁 Location: ~/NGS/results/$SAMPLE_NAME/annovar/"
echo ""
echo "📊 Main files:"
echo "  • annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.txt"
echo "  • annotated_${SAMPLE_NAME}_clinical_plus.hg19_multianno.vcf.gz"
echo ""
echo "📊 Separated by type:"
echo "  • separated_by_type/SNPs.txt"
echo "  • separated_by_type/Insertions.txt"
echo "  • separated_by_type/Deletions.txt"
echo "  • separated_by_type/MNPs.txt"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""


#!/bin/bash
# Complete Analysis: Annotate with Zygosity + Separate + Zip
# Usage: bash run_full_analysis_and_zip.sh <sample_name> [vcf_path]

if [ $# -lt 1 ]; then
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  COMPLETE ANALYSIS: Annotate + Zygosity + Separate + Zip  ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Usage: bash run_full_analysis_and_zip.sh <sample_name> [vcf_path]"
    echo ""
    echo "Examples:"
    echo "  bash run_full_analysis_and_zip.sh sample_Alaa"
    echo "  bash run_full_analysis_and_zip.sh Bur ~/NGS/Bur/filtered_variants.vcf"
    echo ""
    exit 1
fi

SAMPLE_NAME=$1
VCF_PATH=$2

echo "╔════════════════════════════════════════════════════════════╗"
echo "║        COMPLETE ANALYSIS WORKFLOW WITH ZYGOSITY           ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Sample: $SAMPLE_NAME"
echo "Started: $(date)"
echo ""

# Step 1: Annotate with zygosity
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 1: ANNOVAR Annotation with Zygosity"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ -n "$VCF_PATH" ]; then
    bash annotate_with_zygosity.sh $SAMPLE_NAME $VCF_PATH
else
    bash annotate_with_zygosity.sh $SAMPLE_NAME
fi

if [ $? -ne 0 ]; then
    echo "❌ Annotation failed"
    exit 1
fi

echo ""
sleep 1

# Step 2: Advanced separation
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 2: Advanced Variant Separation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

bash separate_advanced.sh $SAMPLE_NAME

if [ $? -ne 0 ]; then
    echo "❌ Separation failed"
    exit 1
fi

echo ""
sleep 1

# Step 3: Create zip file
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "STEP 3: Creating Zip File (Small Files Only)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Determine output directory
if [ -d "~/NGS/results/$SAMPLE_NAME/annovar" ]; then
    ANNOT_DIR=~/NGS/results/$SAMPLE_NAME/annovar
elif [ -d "~/NGS/$SAMPLE_NAME" ]; then
    ANNOT_DIR=~/NGS/$SAMPLE_NAME
else
    # Try to find it
    ANNOT_DIR=$(find ~/NGS -type d -name "$SAMPLE_NAME" -o -path "*/results/$SAMPLE_NAME/annovar" | head -1)
fi

cd ~/NGS

echo "Zipping results for $SAMPLE_NAME..."
echo ""

zip -r -q "${SAMPLE_NAME}_COMPLETE.zip" \
    $(find . -path "*/$SAMPLE_NAME/annovar/*_with_zygosity.*.txt" 2>/dev/null) \
    $(find . -path "*/$SAMPLE_NAME/annovar/by_type/*.txt" 2>/dev/null) \
    $(find . -path "*/$SAMPLE_NAME/annovar/by_function/*.txt" 2>/dev/null) \
    $(find . -path "*/$SAMPLE_NAME/fastqc/*.html" 2>/dev/null) \
    $(find . -path "*/$SAMPLE_NAME/trimmed/fastp_report.html" 2>/dev/null) \
    2>/dev/null

# For custom paths like Bur
if [ ! -f "${SAMPLE_NAME}_COMPLETE.zip" ] || [ $(stat -f%z "${SAMPLE_NAME}_COMPLETE.zip" 2>/dev/null || stat -c%s "${SAMPLE_NAME}_COMPLETE.zip" 2>/dev/null || echo 0) -lt 1000 ]; then
    # Try alternative zipping for custom directories
    if [ -d "$SAMPLE_NAME" ]; then
        zip -r -q "${SAMPLE_NAME}_COMPLETE.zip" \
            "$SAMPLE_NAME/"*_with_zygosity*.txt \
            "$SAMPLE_NAME/by_type/" \
            "$SAMPLE_NAME/by_function/" \
            2>/dev/null
    fi
fi

if [ -f "${SAMPLE_NAME}_COMPLETE.zip" ]; then
    echo "✅ Zip file created!"
    echo ""
    ls -lh "${SAMPLE_NAME}_COMPLETE.zip"
    echo ""
    
    # Show contents
    echo "Zip contents:"
    unzip -l "${SAMPLE_NAME}_COMPLETE.zip" | tail -20
else
    echo "❌ Zip creation failed"
    exit 1
fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              ✅ COMPLETE WORKFLOW FINISHED!                ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Completed: $(date)"
echo ""
echo "Download file:"
echo "  ~/NGS/${SAMPLE_NAME}_COMPLETE.zip"
echo ""
echo "Contains:"
echo "  ✅ Main annotated file with zygosity"
echo "  ✅ Separated by type (SNPs, Insertions, Deletions)"
echo "  ✅ SNPs separated by function (Exonic vs Non-Exonic)"
echo "  ✅ Exonic SNPs categorized (Nonsynonymous, Synonymous, etc.)"
echo "  ✅ QC reports (HTML)"
echo ""
echo "File categories:"
echo "  • SNPs_Exonic.txt - Coding changes (⭐ Most important)"
echo "  • SNPs_Nonsynonymous.txt - Amino acid changes"
echo "  • SNPs_NonExonic.txt - Non-coding variants"
echo "  • Insertions_all.txt - Inserted bases"
echo "  • Deletions_all.txt - Deleted bases"
echo ""
echo "════════════════════════════════════════════════════════════"


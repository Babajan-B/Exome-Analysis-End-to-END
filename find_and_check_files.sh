#!/bin/bash
# Find annotation files and check their headers
# Usage: bash find_and_check_files.sh

echo "════════════════════════════════════════════════════════════"
echo "Searching for Annotation Files"
echo "════════════════════════════════════════════════════════════"
echo ""

# Search for all annotation files
echo "Searching for .hg19_multianno files..."
find ~/NGS/results -name "*.hg19_multianno.*" -type f 2>/dev/null

echo ""
echo "════════════════════════════════════════════════════════════"
echo "Checking sample_Alaa directory:"
echo "════════════════════════════════════════════════════════════"
ls -lh ~/NGS/results/sample_Alaa/annovar/ 2>/dev/null || echo "Directory not found"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "Checking sample_6 directory:"
echo "════════════════════════════════════════════════════════════"
ls -lh ~/NGS/results/sample_6/annovar/ 2>/dev/null || echo "Directory not found"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "If files found, check headers with:"
echo "════════════════════════════════════════════════════════════"
echo ""

# Try to find and show headers
TXT_FILE=$(find ~/NGS/results -name "annotated_sample_Alaa.hg19_multianno.txt" -type f 2>/dev/null | head -1)
VCF_FILE=$(find ~/NGS/results -name "annotated_sample_Alaa.hg19_multianno.vcf" -type f 2>/dev/null | head -1)

if [ -f "$TXT_FILE" ]; then
    echo "✅ Found TXT file: $TXT_FILE"
    echo ""
    echo "Column headers:"
    head -1 "$TXT_FILE" | tr '\t' '\n' | nl
elif [ -f "$VCF_FILE" ]; then
    echo "✅ Found VCF file: $VCF_FILE"
    echo ""
    echo "INFO fields (what ANNOVAR added):"
    grep "^##INFO" "$VCF_FILE" | head -30
else
    echo "❌ Files not found - they may have been deleted"
    echo ""
    echo "To recreate, first pull latest scripts:"
    echo "  cd ~/NGS && git pull"
    echo ""
    echo "Then re-annotate:"
    echo "  bash annotate_pass_only.sh sample_Alaa"
fi


#!/bin/bash
# Separate annotated variants by type (SNPs, Indels, etc.)
# Usage: bash separate_by_variant_type.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash separate_by_variant_type.sh <sample_name>"
    echo ""
    echo "This will separate variants into:"
    echo "  • SNPs (single nucleotide polymorphisms)"
    echo "  • Insertions"
    echo "  • Deletions"
    echo "  • Complex variants"
    echo ""
    echo "Example: bash separate_by_variant_type.sh sample_Alaa"
    exit 1
fi

SAMPLE_NAME=$1
OUTPUT_DIR=~/NGS/results/$SAMPLE_NAME/annovar

# Find the annotated file (try different naming patterns)
ANNOT_FILE=$(find $OUTPUT_DIR -name "*_clinical_plus.hg19_multianno.txt" -o -name "*_PASS.hg19_multianno.txt" -o -name "*_lightweight.hg19_multianno.txt" 2>/dev/null | head -1)

if [ ! -f "$ANNOT_FILE" ]; then
    echo "❌ Annotated file not found for sample: $SAMPLE_NAME"
    echo ""
    echo "Looking in: $OUTPUT_DIR"
    echo ""
    echo "Available files:"
    ls -lh $OUTPUT_DIR/*.txt 2>/dev/null || echo "  No .txt files found"
    echo ""
    echo "Run annotation first:"
    echo "  bash annotate_clinical_plus.sh $SAMPLE_NAME"
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "Separating Variants by Type"
echo "Sample: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Input file: $(basename $ANNOT_FILE)"
echo "Size: $(du -h $ANNOT_FILE | cut -f1)"
echo ""

# Create output directory for separated files
SEPARATED_DIR=$OUTPUT_DIR/separated_by_type
mkdir -p $SEPARATED_DIR

# Get header
HEADER=$(head -1 $ANNOT_FILE)

# Count total variants
TOTAL=$(tail -n +2 $ANNOT_FILE | wc -l)
echo "Total variants: $TOTAL"
echo ""

# Separate by type
echo "Step 1: Identifying variant types..."
echo ""

# SNPs: Ref and Alt are same length AND length=1
echo "  Extracting SNPs..."
echo "$HEADER" > $SEPARATED_DIR/SNPs.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4)==1 && length($5)==1' >> $SEPARATED_DIR/SNPs.txt
SNP_COUNT=$(($(wc -l < $SEPARATED_DIR/SNPs.txt) - 1))

# Insertions: Ref shorter than Alt
echo "  Extracting Insertions..."
echo "$HEADER" > $SEPARATED_DIR/Insertions.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4) < length($5)' >> $SEPARATED_DIR/Insertions.txt
INS_COUNT=$(($(wc -l < $SEPARATED_DIR/Insertions.txt) - 1))

# Deletions: Ref longer than Alt
echo "  Extracting Deletions..."
echo "$HEADER" > $SEPARATED_DIR/Deletions.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4) > length($5)' >> $SEPARATED_DIR/Deletions.txt
DEL_COUNT=$(($(wc -l < $SEPARATED_DIR/Deletions.txt) - 1))

# MNPs (Multiple Nucleotide Polymorphisms): Same length but not SNPs
echo "  Extracting MNPs (multi-nucleotide polymorphisms)..."
echo "$HEADER" > $SEPARATED_DIR/MNPs.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4)==length($5) && length($4)>1' >> $SEPARATED_DIR/MNPs.txt
MNP_COUNT=$(($(wc -l < $SEPARATED_DIR/MNPs.txt) - 1))

echo ""
echo "════════════════════════════════════════════════════════════"
echo "✅ Separation Complete!"
echo "════════════════════════════════════════════════════════════"
echo ""

# Summary
echo "Variant Type Summary:"
echo "────────────────────────────────────────────────────────────"
printf "%-20s %10s %8s\n" "Type" "Count" "Percent"
echo "────────────────────────────────────────────────────────────"
printf "%-20s %10d %7.1f%%\n" "SNPs" $SNP_COUNT $(awk "BEGIN {printf \"%.1f\", ($SNP_COUNT/$TOTAL)*100}")
printf "%-20s %10d %7.1f%%\n" "Insertions" $INS_COUNT $(awk "BEGIN {printf \"%.1f\", ($INS_COUNT/$TOTAL)*100}")
printf "%-20s %10d %7.1f%%\n" "Deletions" $DEL_COUNT $(awk "BEGIN {printf \"%.1f\", ($DEL_COUNT/$TOTAL)*100}")
printf "%-20s %10d %7.1f%%\n" "MNPs" $MNP_COUNT $(awk "BEGIN {printf \"%.1f\", ($MNP_COUNT/$TOTAL)*100}")
echo "────────────────────────────────────────────────────────────"
printf "%-20s %10d %7s\n" "Total" $TOTAL "100%"
echo "════════════════════════════════════════════════════════════"
echo ""

# File details
echo "Output files:"
echo ""
echo "  SNPs:"
echo "    File: $SEPARATED_DIR/SNPs.txt"
echo "    Size: $(du -h $SEPARATED_DIR/SNPs.txt | cut -f1)"
echo "    Variants: $SNP_COUNT"
echo ""
echo "  Insertions:"
echo "    File: $SEPARATED_DIR/Insertions.txt"
echo "    Size: $(du -h $SEPARATED_DIR/Insertions.txt | cut -f1)"
echo "    Variants: $INS_COUNT"
echo ""
echo "  Deletions:"
echo "    File: $SEPARATED_DIR/Deletions.txt"
echo "    Size: $(du -h $SEPARATED_DIR/Deletions.txt | cut -f1)"
echo "    Variants: $DEL_COUNT"
echo ""
echo "  MNPs:"
echo "    File: $SEPARATED_DIR/MNPs.txt"
echo "    Size: $(du -h $SEPARATED_DIR/MNPs.txt | cut -f1)"
echo "    Variants: $MNP_COUNT"
echo ""

# Pathogenic counts by type
echo "════════════════════════════════════════════════════════════"
echo "Pathogenic Variants by Type:"
echo "════════════════════════════════════════════════════════════"
echo ""

SNP_PATH=$(grep -i "pathogenic" $SEPARATED_DIR/SNPs.txt 2>/dev/null | wc -l)
INS_PATH=$(grep -i "pathogenic" $SEPARATED_DIR/Insertions.txt 2>/dev/null | wc -l)
DEL_PATH=$(grep -i "pathogenic" $SEPARATED_DIR/Deletions.txt 2>/dev/null | wc -l)
MNP_PATH=$(grep -i "pathogenic" $SEPARATED_DIR/MNPs.txt 2>/dev/null | wc -l)

printf "%-20s %10d\n" "SNPs" $SNP_PATH
printf "%-20s %10d\n" "Insertions" $INS_PATH
printf "%-20s %10d\n" "Deletions" $DEL_PATH
printf "%-20s %10d\n" "MNPs" $MNP_PATH
echo ""

echo "════════════════════════════════════════════════════════════"
echo "Next Steps:"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "1. Download separated files:"
echo "   $SEPARATED_DIR/"
echo ""
echo "2. Open in Excel for analysis:"
echo "   • SNPs.txt - Most common variant type"
echo "   • Insertions.txt - Added nucleotides"
echo "   • Deletions.txt - Removed nucleotides"
echo "   • MNPs.txt - Multiple nucleotide changes"
echo ""
echo "3. Filter each type separately:"
echo "   • By clinical significance"
echo "   • By prediction scores"
echo "   • By population frequency"
echo ""
echo "════════════════════════════════════════════════════════════"


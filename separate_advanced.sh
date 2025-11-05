#!/bin/bash
# Advanced Variant Separation with Functional Categorization
# Separates by type, then SNPs by ExonicFunc
# Usage: bash separate_advanced.sh <sample_name> [annot_file]

if [ $# -lt 1 ]; then
    echo "Usage: bash separate_advanced.sh <sample_name> [annotation_file]"
    echo "Example: bash separate_advanced.sh sample_Alaa"
    exit 1
fi

SAMPLE_NAME=$1

# Find annotation file
if [ $# -eq 2 ]; then
    ANNOT_FILE=$2
else
    # Try to find file with zygosity first
    ANNOT_FILE=$(find ~/NGS/results/$SAMPLE_NAME/annovar -name "*_with_zygosity.hg19_multianno.txt" 2>/dev/null | head -1)
    
    # If not found, try other patterns
    if [ ! -f "$ANNOT_FILE" ]; then
        ANNOT_FILE=$(find ~/NGS/results/$SAMPLE_NAME/annovar ~/NGS/$SAMPLE_NAME -name "*.hg19_multianno.txt" 2>/dev/null | head -1)
    fi
fi

if [ ! -f "$ANNOT_FILE" ]; then
    echo "❌ Annotation file not found for: $SAMPLE_NAME"
    exit 1
fi

OUTPUT_BASE=$(dirname "$ANNOT_FILE")

echo "════════════════════════════════════════════════════════════"
echo "Advanced Variant Separation"
echo "Sample: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Input: $(basename $ANNOT_FILE)"
echo "Size: $(du -h $ANNOT_FILE | cut -f1)"
echo ""

# Create directories
mkdir -p $OUTPUT_BASE/by_type
mkdir -p $OUTPUT_BASE/by_function

HEADER=$(head -1 $ANNOT_FILE)

# Find ExonicFunc column (usually column 9)
EXONIC_COL=$(head -1 $ANNOT_FILE | tr '\t' '\n' | grep -n "ExonicFunc.refGene" | cut -d: -f1)

if [ -z "$EXONIC_COL" ]; then
    echo "⚠️  ExonicFunc.refGene column not found, using column 9"
    EXONIC_COL=9
fi

echo "Step 1: Separating by variant type..."
echo ""

# Separate by type
echo "  Creating SNPs.txt..."
echo "$HEADER" > $OUTPUT_BASE/by_type/SNPs_all.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4)==1 && length($5)==1' >> $OUTPUT_BASE/by_type/SNPs_all.txt

echo "  Creating Insertions.txt..."
echo "$HEADER" > $OUTPUT_BASE/by_type/Insertions_all.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4) < length($5)' >> $OUTPUT_BASE/by_type/Insertions_all.txt

echo "  Creating Deletions.txt..."
echo "$HEADER" > $OUTPUT_BASE/by_type/Deletions_all.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4) > length($5)' >> $OUTPUT_BASE/by_type/Deletions_all.txt

echo "  Creating MNPs.txt..."
echo "$HEADER" > $OUTPUT_BASE/by_type/MNPs_all.txt
tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4)==length($5) && length($4)>1' >> $OUTPUT_BASE/by_type/MNPs_all.txt

# Count variants
SNP_COUNT=$(($(wc -l < $OUTPUT_BASE/by_type/SNPs_all.txt) - 1))
INS_COUNT=$(($(wc -l < $OUTPUT_BASE/by_type/Insertions_all.txt) - 1))
DEL_COUNT=$(($(wc -l < $OUTPUT_BASE/by_type/Deletions_all.txt) - 1))
MNP_COUNT=$(($(wc -l < $OUTPUT_BASE/by_type/MNPs_all.txt) - 1))

echo ""
echo "Step 2: Separating SNPs by ExonicFunc..."
echo ""

# Separate SNPs by ExonicFunc
echo "  Creating SNPs_Exonic.txt (has ExonicFunc value)..."
echo "$HEADER" > $OUTPUT_BASE/by_function/SNPs_Exonic.txt
tail -n +2 $OUTPUT_BASE/by_type/SNPs_all.txt | awk -F'\t' -v col=$EXONIC_COL '$col != "."' >> $OUTPUT_BASE/by_function/SNPs_Exonic.txt

echo "  Creating SNPs_NonExonic.txt (ExonicFunc = '.')..."
echo "$HEADER" > $OUTPUT_BASE/by_function/SNPs_NonExonic.txt
tail -n +2 $OUTPUT_BASE/by_type/SNPs_all.txt | awk -F'\t' -v col=$EXONIC_COL '$col == "."' >> $OUTPUT_BASE/by_function/SNPs_NonExonic.txt

SNP_EXONIC=$(($(wc -l < $OUTPUT_BASE/by_function/SNPs_Exonic.txt) - 1))
SNP_NONEXONIC=$(($(wc -l < $OUTPUT_BASE/by_function/SNPs_NonExonic.txt) - 1))

# Further categorize exonic SNPs
echo ""
echo "Step 3: Categorizing exonic SNPs by type..."
echo ""

echo "  Nonsynonymous SNVs..."
echo "$HEADER" > $OUTPUT_BASE/by_function/SNPs_Nonsynonymous.txt
tail -n +2 $OUTPUT_BASE/by_function/SNPs_Exonic.txt | grep "nonsynonymous SNV" >> $OUTPUT_BASE/by_function/SNPs_Nonsynonymous.txt

echo "  Synonymous SNVs..."
echo "$HEADER" > $OUTPUT_BASE/by_function/SNPs_Synonymous.txt
tail -n +2 $OUTPUT_BASE/by_function/SNPs_Exonic.txt | grep "synonymous SNV" >> $OUTPUT_BASE/by_function/SNPs_Synonymous.txt

echo "  Stopgain variants..."
echo "$HEADER" > $OUTPUT_BASE/by_function/SNPs_Stopgain.txt
tail -n +2 $OUTPUT_BASE/by_function/SNPs_Exonic.txt | grep "stopgain" >> $OUTPUT_BASE/by_function/SNPs_Stopgain.txt

NONSYN=$(($(wc -l < $OUTPUT_BASE/by_function/SNPs_Nonsynonymous.txt) - 1))
SYN=$(($(wc -l < $OUTPUT_BASE/by_function/SNPs_Synonymous.txt) - 1))
STOP=$(($(wc -l < $OUTPUT_BASE/by_function/SNPs_Stopgain.txt) - 1))

echo ""
echo "════════════════════════════════════════════════════════════"
echo "✅ Advanced Separation Complete!"
echo "════════════════════════════════════════════════════════════"
echo ""

# Summary table
echo "VARIANT TYPE SUMMARY:"
echo "────────────────────────────────────────────────────────────"
printf "%-25s %10s %8s\n" "Category" "Count" "Percent"
echo "────────────────────────────────────────────────────────────"

TOTAL=$((SNP_COUNT + INS_COUNT + DEL_COUNT + MNP_COUNT))

printf "%-25s %10d %7.1f%%\n" "Total Variants" $TOTAL 100.0
echo ""
printf "%-25s %10d %7.1f%%\n" "SNPs (all)" $SNP_COUNT $(awk "BEGIN {printf \"%.1f\", ($SNP_COUNT/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "  - Exonic SNPs" $SNP_EXONIC $(awk "BEGIN {printf \"%.1f\", ($SNP_EXONIC/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "    • Nonsynonymous" $NONSYN $(awk "BEGIN {printf \"%.1f\", ($NONSYN/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "    • Synonymous" $SYN $(awk "BEGIN {printf \"%.1f\", ($SYN/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "    • Stopgain" $STOP $(awk "BEGIN {printf \"%.1f\", ($STOP/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "  - Non-Exonic SNPs" $SNP_NONEXONIC $(awk "BEGIN {printf \"%.1f\", ($SNP_NONEXONIC/$TOTAL)*100}")
echo ""
printf "%-25s %10d %7.1f%%\n" "Insertions" $INS_COUNT $(awk "BEGIN {printf \"%.1f\", ($INS_COUNT/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "Deletions" $DEL_COUNT $(awk "BEGIN {printf \"%.1f\", ($DEL_COUNT/$TOTAL)*100}")
printf "%-25s %10d %7.1f%%\n" "MNPs" $MNP_COUNT $(awk "BEGIN {printf \"%.1f\", ($MNP_COUNT/$TOTAL)*100}")
echo "════════════════════════════════════════════════════════════"
echo ""

# Check zygosity distribution
if grep -q "Zygosity" "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" 2>/dev/null; then
    HET=$(grep -c "Heterozygous" "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" 2>/dev/null || echo 0)
    HOM=$(grep -c "Homozygous" "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" 2>/dev/null || echo 0)
    
    echo "ZYGOSITY DISTRIBUTION:"
    echo "────────────────────────────────────────────────────────────"
    printf "%-25s %10d %7.1f%%\n" "Heterozygous (0/1)" $HET $(awk "BEGIN {printf \"%.1f\", ($HET/$TOTAL)*100}")
    printf "%-25s %10d %7.1f%%\n" "Homozygous (1/1)" $HOM $(awk "BEGIN {printf \"%.1f\", ($HOM/$TOTAL)*100}")
    echo "════════════════════════════════════════════════════════════"
    echo ""
fi

echo "OUTPUT FILES:"
echo "────────────────────────────────────────────────────────────"
echo "Main file with zygosity:"
echo "  • annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt"
echo ""
echo "By variant type:"
echo "  • by_type/SNPs_all.txt ($(du -h $OUTPUT_BASE/by_type/SNPs_all.txt | cut -f1))"
echo "  • by_type/Insertions_all.txt ($(du -h $OUTPUT_BASE/by_type/Insertions_all.txt | cut -f1))"
echo "  • by_type/Deletions_all.txt ($(du -h $OUTPUT_BASE/by_type/Deletions_all.txt | cut -f1))"
echo ""
echo "By function (SNPs only):"
echo "  • by_function/SNPs_Exonic.txt ($(du -h $OUTPUT_BASE/by_function/SNPs_Exonic.txt | cut -f1))"
echo "  • by_function/SNPs_NonExonic.txt ($(du -h $OUTPUT_BASE/by_function/SNPs_NonExonic.txt | cut -f1))"
echo "  • by_function/SNPs_Nonsynonymous.txt ($(du -h $OUTPUT_BASE/by_function/SNPs_Nonsynonymous.txt | cut -f1))"
echo "  • by_function/SNPs_Synonymous.txt ($(du -h $OUTPUT_BASE/by_function/SNPs_Synonymous.txt | cut -f1))"
echo "  • by_function/SNPs_Stopgain.txt ($(du -h $OUTPUT_BASE/by_function/SNPs_Stopgain.txt | cut -f1))"
echo ""
echo "════════════════════════════════════════════════════════════"


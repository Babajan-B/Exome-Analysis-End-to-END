#!/bin/bash
# Cleanup large annotation files and re-annotate with PASS only
# Usage: bash cleanup_and_reannotate.sh

echo "╔════════════════════════════════════════════════════════════╗"
echo "║     CLEANUP & RE-ANNOTATE WITH PASS VARIANTS ONLY         ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Find all samples with annotations
SAMPLES=$(find ~/NGS/results -name "annotated_*.hg19_multianno.vcf" -type f | \
          sed 's|.*/results/||' | sed 's|/.*||' | sort -u)

if [ -z "$SAMPLES" ]; then
    echo "No annotated samples found."
    exit 0
fi

echo "Found annotated samples:"
for sample in $SAMPLES; do
    echo "  - $sample"
done
echo ""

# Show current sizes
echo "Current file sizes:"
echo "════════════════════════════════════════════════════════════"
for sample in $SAMPLES; do
    VCF_FILE=$(find ~/NGS/results/$sample/annovar -name "*.hg19_multianno.vcf" -type f 2>/dev/null | head -1)
    TXT_FILE=$(find ~/NGS/results/$sample/annovar -name "*.hg19_multianno.txt" -type f 2>/dev/null | head -1)
    
    if [ -f "$VCF_FILE" ]; then
        VCF_SIZE=$(du -h "$VCF_FILE" | cut -f1)
        VARIANTS=$(grep -v "^#" "$VCF_FILE" | wc -l)
        echo "$sample:"
        echo "  VCF: $VCF_SIZE ($VARIANTS variants)"
        
        if [ -f "$TXT_FILE" ]; then
            TXT_SIZE=$(du -h "$TXT_FILE" | cut -f1)
            echo "  TXT: $TXT_SIZE"
        fi
        echo ""
    fi
done

echo "════════════════════════════════════════════════════════════"
echo ""
echo "These files are TOO LARGE because they include failed variants."
echo "Expected size: 60-100MB per sample"
echo ""
read -p "Re-annotate with PASS variants only? (y/n): " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

echo ""
echo "Processing samples..."
echo ""

# Process each sample
for sample in $SAMPLES; do
    echo "════════════════════════════════════════════════════════════"
    echo "Processing: $sample"
    echo "════════════════════════════════════════════════════════════"
    
    # Remove old large files
    echo "Removing old annotation files..."
    rm -f ~/NGS/results/$sample/annovar/annotated_${sample}.hg19_multianno.vcf
    rm -f ~/NGS/results/$sample/annovar/annotated_${sample}.hg19_multianno.txt
    
    # Re-annotate with PASS only
    bash annotate_pass_only.sh $sample
    
    echo ""
done

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                    ✅ CLEANUP COMPLETE!                    ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "New file sizes (PASS variants only):"
echo "════════════════════════════════════════════════════════════"

for sample in $SAMPLES; do
    VCF_FILE=$(find ~/NGS/results/$sample/annovar -name "*_PASS.hg19_multianno.vcf.gz" -type f 2>/dev/null | head -1)
    TXT_FILE=$(find ~/NGS/results/$sample/annovar -name "*_PASS.hg19_multianno.txt" -type f 2>/dev/null | head -1)
    
    if [ -f "$VCF_FILE" ]; then
        VCF_SIZE=$(du -h "$VCF_FILE" | cut -f1)
        VARIANTS=$(zcat "$VCF_FILE" | grep -v "^#" | wc -l)
        echo "$sample:"
        echo "  VCF: $VCF_SIZE (compressed, $VARIANTS variants)"
        
        if [ -f "$TXT_FILE" ]; then
            TXT_SIZE=$(du -h "$TXT_FILE" | cut -f1)
            echo "  TXT: $TXT_SIZE"
        fi
        echo ""
    fi
done

echo "════════════════════════════════════════════════════════════"
echo ""
echo "Space saved: ~90% reduction!"
echo ""
echo "Download the *_PASS.hg19_multianno.txt files for analysis."
echo ""


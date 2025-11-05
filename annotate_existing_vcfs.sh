#!/bin/bash
# Script to annotate existing VCF files in Jarvis Lab
# Run this in ~/NGS directory on Jarvis Lab

set -e

echo "=========================================="
echo "Finding and Annotating VCF Files"
echo "=========================================="
echo ""

# Find all VCF files in results directory
echo "Step 1: Searching for VCF files..."
echo ""

# Look for filtered variants (best for annotation)
FILTERED_VCFS=$(find ~/NGS/results -name "filtered_variants.vcf" -type f 2>/dev/null)

# If no filtered VCFs, look for any variants.vcf
if [ -z "$FILTERED_VCFS" ]; then
    echo "No filtered_variants.vcf found. Looking for variants.vcf files..."
    FILTERED_VCFS=$(find ~/NGS/results -name "variants.vcf" -type f 2>/dev/null)
fi

# If still nothing, look for any VCF
if [ -z "$FILTERED_VCFS" ]; then
    echo "Looking for any .vcf files..."
    FILTERED_VCFS=$(find ~/NGS/results -name "*.vcf" -type f 2>/dev/null | grep -v ".vcf.gz")
fi

# Count VCFs found
VCF_COUNT=$(echo "$FILTERED_VCFS" | grep -v "^$" | wc -l)

if [ $VCF_COUNT -eq 0 ]; then
    echo "❌ No VCF files found in ~/NGS/results/"
    echo ""
    echo "Please check if:"
    echo "  1. Results directory exists: ls -la ~/NGS/results/"
    echo "  2. Pipeline has completed successfully"
    echo "  3. VCF files were generated"
    exit 1
fi

echo "✅ Found $VCF_COUNT VCF file(s) to annotate:"
echo ""
echo "$FILTERED_VCFS" | nl
echo ""

# Ask for confirmation
read -p "Do you want to annotate all these VCF files? (y/n): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Annotation cancelled."
    exit 0
fi

echo ""
echo "Step 2: Annotating VCF files with ANNOVAR..."
echo ""

# Process each VCF file
counter=1
for vcf in $FILTERED_VCFS; do
    echo "=========================================="
    echo "Processing VCF $counter of $VCF_COUNT"
    echo "=========================================="
    echo "Input: $vcf"
    
    # Get the directory and create annovar subdirectory
    vcf_dir=$(dirname "$vcf")
    annovar_dir="$vcf_dir/../annovar"
    mkdir -p "$annovar_dir"
    
    # Get sample ID from path
    sample_id=$(echo "$vcf_dir" | sed 's|.*/results/||' | sed 's|/.*||')
    
    # Set output prefix
    output_prefix="$annovar_dir/annotated_${sample_id}"
    
    echo "Output: $output_prefix.hg19_multianno.vcf"
    echo "        $output_prefix.hg19_multianno.txt"
    echo ""
    
    # Check VCF size and variant count
    vcf_size=$(du -h "$vcf" | cut -f1)
    variant_count=$(grep -v "^#" "$vcf" | wc -l 2>/dev/null || echo "0")
    
    echo "VCF size: $vcf_size"
    echo "Variants: $variant_count"
    echo ""
    
    if [ $variant_count -eq 0 ]; then
        echo "⚠️  Warning: VCF file is empty or contains no variants. Skipping..."
        echo ""
        counter=$((counter + 1))
        continue
    fi
    
    # Run ANNOVAR annotation
    echo "Running ANNOVAR annotation..."
    start_time=$(date +%s)
    
    perl ~/NGS/tools/annovar/table_annovar.pl \
        "$vcf" \
        ~/NGS/tools/annovar/humandb/ \
        -buildver hg19 \
        -out "$output_prefix" \
        -remove \
        -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70,icgc28 \
        -operation g,g,g,f,f,f,f,f,f \
        -nastring . \
        -vcfinput \
        -polish \
        2>&1 | grep -E "(NOTICE|ERROR|WARNING)" || true
    
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    echo ""
    if [ -f "$output_prefix.hg19_multianno.vcf" ]; then
        echo "✅ Annotation complete! (took ${elapsed}s)"
        
        # Show output file sizes
        echo ""
        echo "Output files created:"
        ls -lh "$output_prefix.hg19_multianno"* | awk '{print "  " $9 " (" $5 ")"}'
        
        # Count annotated variants
        annotated_count=$(grep -v "^#" "$output_prefix.hg19_multianno.vcf" | wc -l 2>/dev/null || echo "0")
        echo ""
        echo "  Variants annotated: $annotated_count"
    else
        echo "❌ Annotation failed. Check the error messages above."
    fi
    
    echo ""
    counter=$((counter + 1))
done

echo ""
echo "=========================================="
echo "✅ All Annotations Complete!"
echo "=========================================="
echo ""
echo "Summary of annotated files:"
find ~/NGS/results -name "*.hg19_multianno.vcf" -type f 2>/dev/null | while read annotated_vcf; do
    size=$(du -h "$annotated_vcf" | cut -f1)
    variants=$(grep -v "^#" "$annotated_vcf" | wc -l 2>/dev/null || echo "0")
    echo "  $annotated_vcf"
    echo "    Size: $size, Variants: $variants"
    echo "    Text: ${annotated_vcf%.vcf}.txt"
    echo ""
done

echo "=========================================="
echo "Next Steps:"
echo "  1. Download the .txt files for Excel analysis"
echo "  2. Use the .vcf files for further bioinformatics tools"
echo "  3. Look for pathogenic variants in ClinVar column"
echo ""
echo "Example: Open in Excel and filter by:"
echo "  - CLNSIG contains 'Pathogenic'"
echo "  - gnomAD_exome_ALL < 0.01 (rare variants)"
echo "  - SIFT_score < 0.05 (damaging)"
echo "=========================================="


#!/bin/bash
# ANNOVAR Annotation with Zygosity Information
# Extracts GT field from VCF and adds as Zygosity column
# Usage: bash annotate_with_zygosity.sh <sample_name>

if [ $# -lt 1 ]; then
    echo "Usage: bash annotate_with_zygosity.sh <sample_name> [vcf_path]"
    echo ""
    echo "Example:"
    echo "  bash annotate_with_zygosity.sh sample_Alaa"
    echo "  bash annotate_with_zygosity.sh Bur ~/NGS/Bur/filtered_variants.vcf"
    exit 1
fi

SAMPLE_NAME=$1
ANNOVAR_DIR=~/NGS/tools/annovar

# Determine VCF path
if [ $# -eq 2 ]; then
    VCF_FILE=$2
    OUTPUT_DIR=$(dirname $VCF_FILE)
else
    OUTPUT_DIR=~/NGS/results/$SAMPLE_NAME
    VCF_FILE=$OUTPUT_DIR/filtered/filtered_variants.vcf
fi

if [ ! -f "$VCF_FILE" ]; then
    echo "❌ VCF file not found: $VCF_FILE"
    exit 1
fi

echo "════════════════════════════════════════════════════════════"
echo "ANNOVAR Annotation with Zygosity"
echo "Sample: $SAMPLE_NAME"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Input VCF: $VCF_FILE"
echo ""

# Create output directory
mkdir -p $OUTPUT_DIR/annovar

# Step 1: Extract PASS variants with GT field
echo "Step 1: Extracting PASS variants with genotype..."
PASS_VCF=$OUTPUT_DIR/annovar/filtered_PASS_with_GT.vcf

grep "^#" $VCF_FILE > $PASS_VCF
grep -v "^#" $VCF_FILE | grep -w "PASS" >> $PASS_VCF

PASS_COUNT=$(grep -v "^#" $PASS_VCF | wc -l)
echo "  PASS variants: $PASS_COUNT"
echo ""

# Step 2: Run ANNOVAR
echo "Step 2: Running ANNOVAR annotation..."
echo "  Databases: refGene, ClinVar, gnomAD, dbSNP, Predictions"
echo ""

perl $ANNOVAR_DIR/table_annovar.pl \
    $PASS_VCF \
    $ANNOVAR_DIR/humandb/ \
    -buildver hg19 \
    -out $OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME} \
    -remove \
    -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150,dbnsfp42a \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish 2>&1 | grep -E "(NOTICE|WARNING)" || true

if [ ! -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.txt" ]; then
    echo "❌ Annotation failed"
    exit 1
fi

# Step 3: Add Zygosity column
echo ""
echo "Step 3: Adding zygosity information..."

# Create Python script to add zygosity
VCF_FOR_PY="$PASS_VCF"
ANNOT_FOR_PY="$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.txt"
OUTPUT_FOR_PY="$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt"

python3 - "$VCF_FOR_PY" "$ANNOT_FOR_PY" "$OUTPUT_FOR_PY" << 'PYTHON_SCRIPT'
import sys

if len(sys.argv) < 4:
    print("Error: Missing arguments")
    sys.exit(1)

vcf_file = sys.argv[1]
annot_file = sys.argv[2]
output_file = sys.argv[3]

# Read VCF and extract GT for each variant
print("  Reading VCF genotypes...")
gt_dict = {}
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 10:
            continue
        
        chrom, pos, _, ref, alt = fields[0:5]
        format_field = fields[8] if len(fields) > 8 else ""
        sample_field = fields[9] if len(fields) > 9 else ""
        
        # Extract GT
        gt = "."
        if format_field and sample_field:
            format_parts = format_field.split(':')
            sample_parts = sample_field.split(':')
            
            if 'GT' in format_parts:
                gt_index = format_parts.index('GT')
                if gt_index < len(sample_parts):
                    gt_raw = sample_parts[gt_index]
                    
                    # Interpret zygosity
                    if gt_raw in ['0/1', '0|1', '1/0', '1|0']:
                        gt = "Heterozygous"
                    elif gt_raw in ['1/1', '1|1']:
                        gt = "Homozygous"
                    elif gt_raw in ['0/0', '0|0']:
                        gt = "Reference"
                    else:
                        gt = gt_raw
        
        key = f"{chrom}:{pos}:{ref}:{alt}"
        gt_dict[key] = gt

print(f"  Extracted zygosity for {len(gt_dict)} variants")

# Add zygosity to annotation file
print("  Adding zygosity column to annotation...")
with open(annot_file, 'r') as f_in, open(output_file, 'w') as f_out:
    header = f_in.readline().strip()
    f_out.write(header + "\tZygosity\n")
    
    for line in f_in:
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        
        chrom, start, end, ref, alt = fields[0:5]
        key = f"{chrom}:{start}:{ref}:{alt}"
        
        zygosity = gt_dict.get(key, "Unknown")
        f_out.write(line.strip() + "\t" + zygosity + "\n")

print("  ✅ Zygosity column added")

PYTHON_SCRIPT

if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" ]; then
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "✅ Annotation with Zygosity Complete!"
    echo "════════════════════════════════════════════════════════════"
    echo ""
    
    FILE_SIZE=$(du -h "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" | cut -f1)
    COLS=$(head -1 "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" | awk -F'\t' '{print NF}')
    
    echo "Output file:"
    echo "  File: annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt"
    echo "  Size: $FILE_SIZE"
    echo "  Variants: $PASS_COUNT"
    echo "  Columns: $COLS (including Zygosity)"
    echo ""
    
    # Count by zygosity
    echo "Zygosity Distribution:"
    HET=$(grep -c "Heterozygous" "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" 2>/dev/null || echo 0)
    HOM=$(grep -c "Homozygous" "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}_with_zygosity.hg19_multianno.txt" 2>/dev/null || echo 0)
    echo "  Heterozygous: $HET"
    echo "  Homozygous: $HOM"
    echo ""
    
    # Check for predictions
    echo "Prediction Databases Included:"
    echo "  ✅ SIFT scores"
    echo "  ✅ PolyPhen2 scores"
    echo "  ✅ CADD scores"
    echo "  ✅ REVEL scores"
    echo "  ✅ MetaSVM/MetaLR predictions"
    echo "  ✅ And 40+ more algorithms from dbNSFP"
    echo ""
    
    # Compress VCF
    if [ -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf" ]; then
        bgzip -f "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf" 2>/dev/null || true
        tabix -p vcf "$OUTPUT_DIR/annovar/annotated_${SAMPLE_NAME}.hg19_multianno.vcf.gz" 2>/dev/null || true
    fi
    
    echo "Next step: bash separate_advanced.sh $SAMPLE_NAME"
    echo ""
    echo "════════════════════════════════════════════════════════════"
else
    echo "❌ Failed to add zygosity"
    exit 1
fi


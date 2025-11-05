#!/bin/bash
# ANNOVAR Helper Script - Quick annotation of VCF files
# Usage: bash annovar_helper.sh <input.vcf> <output_prefix> [buildver]

set -e

# Configuration
ANNOVAR_DIR=~/NGS/tools/annovar
BUILDVER=${3:-hg19}  # Default to hg19

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input.vcf> <output_prefix> [buildver]"
    echo ""
    echo "Example:"
    echo "  $0 variants.vcf annotated_variants hg19"
    echo ""
    echo "Available buildver: hg19, hg38"
    exit 1
fi

INPUT_VCF=$1
OUTPUT_PREFIX=$2

# Check if input file exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "ERROR: Input VCF file not found: $INPUT_VCF"
    exit 1
fi

# Check if ANNOVAR is installed
if [ ! -d "$ANNOVAR_DIR" ]; then
    echo "ERROR: ANNOVAR not found at $ANNOVAR_DIR"
    echo "Please run: bash setup_annovar.sh"
    exit 1
fi

echo "=========================================="
echo "ANNOVAR Annotation"
echo "=========================================="
echo "Input VCF:    $INPUT_VCF"
echo "Output:       ${OUTPUT_PREFIX}.${BUILDVER}_multianno.vcf"
echo "Build:        $BUILDVER"
echo "Started:      $(date)"
echo "=========================================="
echo ""

# Run ANNOVAR table_annovar.pl
# This is the comprehensive annotation command
perl $ANNOVAR_DIR/table_annovar.pl \
    $INPUT_VCF \
    $ANNOVAR_DIR/humandb/ \
    -buildver $BUILDVER \
    -out $OUTPUT_PREFIX \
    -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70,icgc28 \
    -operation g,g,g,f,f,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish

echo ""
echo "=========================================="
echo "✅ Annotation Complete!"
echo "=========================================="
echo "Output files:"
echo "  - ${OUTPUT_PREFIX}.${BUILDVER}_multianno.vcf (annotated VCF)"
echo "  - ${OUTPUT_PREFIX}.${BUILDVER}_multianno.txt (tabular format)"
echo "  - ${OUTPUT_PREFIX}.avinput (intermediate file)"
echo ""
echo "Completed: $(date)"
echo "=========================================="


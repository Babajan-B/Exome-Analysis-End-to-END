#!/bin/bash
# Post-Processing Script: Add snpEff, Advanced Separation, Zygosity, and Zip
# Runs after MASTER_PIPELINE.sh completes
# Usage: bash post_process_complete.sh

set -e

echo "╔════════════════════════════════════════════════════════════╗"
echo "║          POST-PROCESSING: snpEff + Advanced Analysis       ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

RESULTS_DIR=~/NGS/results
SNPEFF_DIR=~/NGS/tools/snpEff
SNPEFF_DB="GRCh37.75"

# Check if results exist
if [ ! -d "$RESULTS_DIR" ]; then
    echo "❌ No results directory found. Run MASTER_PIPELINE.sh first."
    exit 1
fi

# Find all samples
SAMPLES=($(ls -d $RESULTS_DIR/*/ 2>/dev/null | xargs -n 1 basename))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "❌ No samples found in $RESULTS_DIR"
    exit 1
fi

echo "Found ${#SAMPLES[@]} sample(s): ${SAMPLES[@]}"
echo ""
echo "Will process:"
echo "  1. Add snpEff annotation"
echo "  2. Add zygosity information"
echo "  3. Advanced functional separation"
echo "  4. Create final ZIP archive"
echo ""
read -p "Continue? (y/n): " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Processing cancelled."
    exit 0
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "PROCESSING ALL SAMPLES"
echo "════════════════════════════════════════════════════════════"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  Processing: $SAMPLE"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    SAMPLE_DIR=$RESULTS_DIR/$SAMPLE
    ANNOVAR_DIR=$SAMPLE_DIR/annovar
    FILTERED_VCF=$SAMPLE_DIR/filtered/filtered_variants.vcf
    PASS_VCF=$SAMPLE_DIR/filtered/filtered_PASS_only.vcf
    
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 1: snpEff Annotation
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "STEP 1: snpEff Annotation"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    mkdir -p $ANNOVAR_DIR/snpeff
    
    if [ -f "$PASS_VCF" ]; then
        echo "  Running snpEff on PASS variants..."
        java -Xmx8g -jar $SNPEFF_DIR/snpEff.jar \
            -v $SNPEFF_DB \
            -stats $ANNOVAR_DIR/snpeff/${SAMPLE}_snpEff_summary.html \
            -csvStats $ANNOVAR_DIR/snpeff/${SAMPLE}_snpEff_summary.csv \
            $PASS_VCF \
            > $ANNOVAR_DIR/snpeff/${SAMPLE}_snpEff_annotated.vcf
        
        echo "  ✅ snpEff annotation complete"
        echo "     Output: $ANNOVAR_DIR/snpeff/${SAMPLE}_snpEff_annotated.vcf"
        echo "     Report: $ANNOVAR_DIR/snpeff/${SAMPLE}_snpEff_summary.html"
    else
        echo "  ⚠️  PASS VCF not found, skipping snpEff"
    fi
    
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 2: Add Zygosity Information
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "STEP 2: Adding Zygosity Information"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    ANNOT_TXT=$ANNOVAR_DIR/annotated_${SAMPLE}.hg19_multianno.txt
    ANNOT_WITH_ZYG=$ANNOVAR_DIR/annotated_${SAMPLE}_with_zygosity.txt
    
    if [ -f "$ANNOT_TXT" ] && [ -f "$PASS_VCF" ]; then
        echo "  Extracting zygosity from VCF..."
        
        python3 << 'PYTHON_SCRIPT'
import sys

vcf_file = sys.argv[1]
annot_file = sys.argv[2]
output_file = sys.argv[3]

# Extract GT from VCF
gt_dict = {}
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 9:
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            format_field = parts[8].split(':')
            sample_field = parts[9].split(':')
            
            if 'GT' in format_field:
                gt_index = format_field.index('GT')
                if gt_index < len(sample_field):
                    genotype = sample_field[gt_index]
                    
                    if genotype in ['0/1', '1/0']:
                        zygosity = "Heterozygous"
                    elif genotype == '1/1':
                        zygosity = "Homozygous"
                    elif genotype == '0/0':
                        zygosity = "Reference"
                    else:
                        zygosity = "Unknown"
                    
                    key = f"{chrom}:{pos}:{ref}:{alt}"
                    gt_dict[key] = zygosity

# Add zygosity column
with open(annot_file, 'r') as f_in, open(output_file, 'w') as f_out:
    header = f_in.readline()
    f_out.write(header.strip() + "\tZygosity\n")
    
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            chrom, start, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            key = f"{chrom}:{start}:{ref}:{alt}"
            zygosity = gt_dict.get(key, "Unknown")
            f_out.write(line.strip() + "\t" + zygosity + "\n")

print("  ✅ Zygosity column added")
PYTHON_SCRIPT
        
        python3 -c "$(cat)" "$PASS_VCF" "$ANNOT_TXT" "$ANNOT_WITH_ZYG" << 'PYTHON_SCRIPT'
import sys

vcf_file = sys.argv[1]
annot_file = sys.argv[2]
output_file = sys.argv[3]

# Extract GT from VCF
gt_dict = {}
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 9:
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            format_field = parts[8].split(':')
            sample_field = parts[9].split(':')
            
            if 'GT' in format_field:
                gt_index = format_field.index('GT')
                if gt_index < len(sample_field):
                    genotype = sample_field[gt_index]
                    
                    if genotype in ['0/1', '1/0']:
                        zygosity = "Heterozygous"
                    elif genotype == '1/1':
                        zygosity = "Homozygous"
                    elif genotype == '0/0':
                        zygosity = "Reference"
                    else:
                        zygosity = "Unknown"
                    
                    key = f"{chrom}:{pos}:{ref}:{alt}"
                    gt_dict[key] = zygosity

# Add zygosity column
with open(annot_file, 'r') as f_in, open(output_file, 'w') as f_out:
    header = f_in.readline()
    f_out.write(header.strip() + "\tZygosity\n")
    
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            chrom, start, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            key = f"{chrom}:{start}:{ref}:{alt}"
            zygosity = gt_dict.get(key, "Unknown")
            f_out.write(line.strip() + "\t" + zygosity + "\n")

print("  ✅ Zygosity column added")
PYTHON_SCRIPT
        
    else
        echo "  ⚠️  Required files not found, skipping zygosity"
    fi
    
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # STEP 3: Advanced Functional Separation
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "STEP 3: Advanced Functional Separation"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    FUNC_DIR=$ANNOVAR_DIR/functional_classification
    mkdir -p $FUNC_DIR
    
    INPUT_FILE=$ANNOT_WITH_ZYG
    if [ ! -f "$INPUT_FILE" ]; then
        INPUT_FILE=$ANNOT_TXT
    fi
    
    if [ -f "$INPUT_FILE" ]; then
        echo "  Separating by functional effect..."
        
        HEADER=$(head -1 $INPUT_FILE)
        
        # SNPs - Exonic
        echo "$HEADER" > $FUNC_DIR/SNPs_Exonic.txt
        tail -n +2 $INPUT_FILE | awk -F'\t' 'length($4)==1 && length($5)==1 && $6=="exonic"' >> $FUNC_DIR/SNPs_Exonic.txt
        
        # SNPs - Non-Exonic
        echo "$HEADER" > $FUNC_DIR/SNPs_NonExonic.txt
        tail -n +2 $INPUT_FILE | awk -F'\t' 'length($4)==1 && length($5)==1 && $6!="exonic"' >> $FUNC_DIR/SNPs_NonExonic.txt
        
        # Exonic - Nonsynonymous
        echo "$HEADER" > $FUNC_DIR/Exonic_Nonsynonymous.txt
        tail -n +2 $INPUT_FILE | awk -F'\t' '$6=="exonic" && $9 ~ /nonsynonymous/' >> $FUNC_DIR/Exonic_Nonsynonymous.txt
        
        # Exonic - Synonymous
        echo "$HEADER" > $FUNC_DIR/Exonic_Synonymous.txt
        tail -n +2 $INPUT_FILE | awk -F'\t' '$6=="exonic" && $9 ~ /synonymous/' >> $FUNC_DIR/Exonic_Synonymous.txt
        
        # Exonic - Stopgain
        echo "$HEADER" > $FUNC_DIR/Exonic_Stopgain.txt
        tail -n +2 $INPUT_FILE | awk -F'\t' '$6=="exonic" && $9 ~ /stopgain/' >> $FUNC_DIR/Exonic_Stopgain.txt
        
        # Exonic - Frameshift
        echo "$HEADER" > $FUNC_DIR/Exonic_Frameshift.txt
        tail -n +2 $INPUT_FILE | awk -F'\t' '$6=="exonic" && $9 ~ /frameshift/' >> $FUNC_DIR/Exonic_Frameshift.txt
        
        # Count results
        SNP_EXONIC=$(($(wc -l < $FUNC_DIR/SNPs_Exonic.txt) - 1))
        SNP_NONEXONIC=$(($(wc -l < $FUNC_DIR/SNPs_NonExonic.txt) - 1))
        NONSYN=$(($(wc -l < $FUNC_DIR/Exonic_Nonsynonymous.txt) - 1))
        SYN=$(($(wc -l < $FUNC_DIR/Exonic_Synonymous.txt) - 1))
        STOP=$(($(wc -l < $FUNC_DIR/Exonic_Stopgain.txt) - 1))
        FRAME=$(($(wc -l < $FUNC_DIR/Exonic_Frameshift.txt) - 1))
        
        echo "  ✅ Functional classification complete:"
        echo "     SNPs Exonic: $SNP_EXONIC"
        echo "     SNPs Non-Exonic: $SNP_NONEXONIC"
        echo "     Nonsynonymous: $NONSYN"
        echo "     Synonymous: $SYN"
        echo "     Stopgain: $STOP"
        echo "     Frameshift: $FRAME"
    else
        echo "  ⚠️  Annotation file not found, skipping"
    fi
    
    echo ""
    echo "✅ Sample $SAMPLE: POST-PROCESSING COMPLETE!"
    echo ""
done

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 4: Create Final ZIP Archive
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║             CREATING FINAL ZIP ARCHIVE                     ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

cd ~/NGS

ZIP_NAME="NGS_Results_Complete_$(date +%Y%m%d_%H%M%S).zip"

echo "Creating archive: $ZIP_NAME"
echo ""

# Include essential files only
zip -r $ZIP_NAME \
    results/*/annovar/*.txt \
    results/*/annovar/snpeff/*.html \
    results/*/annovar/snpeff/*.csv \
    results/*/annovar/separated_by_type/*.txt \
    results/*/annovar/functional_classification/*.txt \
    results/*/fastqc/*.html \
    results/*/trimmed/fastp_report.html \
    MASTER_ANALYSIS_SUMMARY.txt \
    -x "*.bam" "*.sam" "*.fastq.gz" "*.vcf" "*.avinput" "*_dropped" "*_filtered" \
    2>/dev/null

ZIP_SIZE=$(du -sh $ZIP_NAME | cut -f1)

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              ✅ POST-PROCESSING COMPLETE!                  ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Archive created: $ZIP_NAME ($ZIP_SIZE)"
echo ""
echo "Contents:"
echo "  ✅ ANNOVAR annotation (.txt files)"
echo "  ✅ ANNOVAR with zygosity"
echo "  ✅ snpEff annotation + reports"
echo "  ✅ Variant type separation (SNPs, Insertions, Deletions)"
echo "  ✅ Functional classification (Exonic, Nonsynonymous, etc.)"
echo "  ✅ Quality control reports (FastQC, fastp)"
echo "  ✅ Summary report"
echo ""
echo "Download this file:"
echo "  ~/NGS/$ZIP_NAME"
echo ""
echo "════════════════════════════════════════════════════════════"
echo ""


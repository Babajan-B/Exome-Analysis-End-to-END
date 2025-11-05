#!/bin/bash
# MASTER PIPELINE - Complete End-to-End Exome Analysis
# Automatically detects samples and runs: QC → Alignment → Variant Calling → Annotation → Separation
# Usage: bash MASTER_PIPELINE.sh [data_directory] [threads]

set -e

# Configuration
DATA_DIR=${1:-~/NGS/data}
THREADS=${2:-16}
WORK_DIR=~/NGS
REFERENCE=$WORK_DIR/reference/hg19.fa
ANNOVAR_DIR=$WORK_DIR/tools/annovar
GATK=/opt/gatk-4.6.2.0/gatk

echo "╔════════════════════════════════════════════════════════════╗"
echo "║           MASTER EXOME ANALYSIS PIPELINE                  ║"
echo "║        Auto-Detection → Complete Analysis                 ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Pipeline Steps:"
echo "  1. Auto-detect FASTQ samples"
echo "  2. Quality Control (FastQC)"
echo "  3. Read Trimming (fastp)"
echo "  4. Alignment (BWA-MEM)"
echo "  5. BAM Processing (sort, mark duplicates)"
echo "  6. Variant Calling (GATK HaplotypeCaller)"
echo "  7. Variant Filtering"
echo "  8. ANNOVAR Annotation (5 databases)"
echo "  9. Variant Type Separation (SNPs, Indels, etc.)"
echo "  10. Summary Report Generation"
echo ""
echo "Started: $(date)"
echo ""

# Check data directory
if [ ! -d "$DATA_DIR" ]; then
    echo "❌ Data directory not found: $DATA_DIR"
    echo ""
    echo "Setup:"
    echo "  mkdir -p $DATA_DIR"
    echo "  # Upload your FASTQ files:"
    echo "  #   sample1_R1.fastq.gz, sample1_R2.fastq.gz"
    echo "  #   sample2_R1.fastq.gz, sample2_R2.fastq.gz"
    echo "  # Then run: bash MASTER_PIPELINE.sh"
    exit 1
fi

echo "Configuration:"
echo "  Data Directory: $DATA_DIR"
echo "  Threads: $THREADS"
echo "  Reference: $REFERENCE"
echo ""

# Function to detect FASTQ pairs
detect_samples() {
    local data_dir=$1
    declare -gA SAMPLE_PAIRS
    
    cd "$data_dir"
    
    # Find R1 files
    for r1_file in *R1*.fastq.gz *R1*.fq.gz *_1.fastq.gz *_1.fq.gz 2>/dev/null; do
        [ -f "$r1_file" ] || continue
        
        # Generate expected R2 filename
        r2_file=$(echo "$r1_file" | sed -e 's/R1/R2/g' -e 's/_1\./_2\./g')
        
        if [ -f "$r2_file" ]; then
            # Extract sample name
            sample_name=$(echo "$r1_file" | sed -E 's/[._-]*(R1|_1)[._-]*.*//' | sed -E 's/\.(fastq|fq)\.gz$//')
            
            SAMPLE_PAIRS["$sample_name"]="$data_dir/$r1_file,$data_dir/$r2_file"
        fi
    done
}

# Function to run complete analysis for one sample
analyze_sample() {
    local sample_name=$1
    local r1_path=$2
    local r2_path=$3
    local threads=$4
    
    local output_dir=$WORK_DIR/results/$sample_name
    
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  Analyzing: $sample_name"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Create output directories
    mkdir -p $output_dir/{fastqc,trimmed,aligned,sorted,dedup,variants,filtered,annovar}
    
    # Log file
    LOG=$output_dir/pipeline.log
    exec > >(tee -a $LOG) 2>&1
    
    # Step function
    step() {
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "STEP $1: $2"
        echo "Time: $(date)"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    }
    
    # 1. FastQC
    step 1 "Quality Control"
    fastqc -t $threads -o $output_dir/fastqc $r1_path $r2_path
    echo "✅ QC complete"
    
    # 2. Trimming
    step 2 "Read Trimming"
    fastp -i $r1_path -I $r2_path \
        -o $output_dir/trimmed/r1_trimmed.fastq.gz \
        -O $output_dir/trimmed/r2_trimmed.fastq.gz \
        -h $output_dir/trimmed/fastp_report.html \
        -j $output_dir/trimmed/fastp_report.json \
        --thread $threads \
        --detect_adapter_for_pe \
        --length_required 50 \
        --qualified_quality_phred 20
    echo "✅ Trimming complete"
    
    # 3. Alignment
    step 3 "Read Alignment"
    bwa mem -t $threads \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:lib_${sample_name}\tPU:unit1" \
        $REFERENCE \
        $output_dir/trimmed/r1_trimmed.fastq.gz \
        $output_dir/trimmed/r2_trimmed.fastq.gz \
        > $output_dir/aligned/aligned.sam
    echo "✅ Alignment complete"
    
    # 4. SAM to BAM
    step 4 "SAM to BAM"
    samtools view -@ $threads -bS $output_dir/aligned/aligned.sam > $output_dir/aligned/aligned.bam
    rm $output_dir/aligned/aligned.sam
    echo "✅ Conversion complete"
    
    # 5. Sort
    step 5 "BAM Sorting"
    samtools sort -@ $threads -o $output_dir/sorted/sorted.bam $output_dir/aligned/aligned.bam
    rm $output_dir/aligned/aligned.bam
    samtools index $output_dir/sorted/sorted.bam
    echo "✅ Sorting complete"
    
    # 6. Mark Duplicates
    step 6 "Mark Duplicates"
    $GATK MarkDuplicates \
        -I $output_dir/sorted/sorted.bam \
        -O $output_dir/dedup/dedup.bam \
        -M $output_dir/dedup/metrics.txt \
        --CREATE_INDEX true
    echo "✅ Duplicates marked"
    
    # 7. Variant Calling
    step 7 "Variant Calling"
    $GATK HaplotypeCaller \
        -R $REFERENCE \
        -I $output_dir/dedup/dedup.bam \
        -O $output_dir/variants/raw_variants.vcf \
        --native-pair-hmm-threads $threads
    
    RAW_COUNT=$(grep -v "^#" $output_dir/variants/raw_variants.vcf | wc -l)
    echo "✅ Called $RAW_COUNT variants"
    
    # 8. Filtering
    step 8 "Variant Filtering"
    $GATK VariantFiltration \
        -R $REFERENCE \
        -V $output_dir/variants/raw_variants.vcf \
        -O $output_dir/filtered/filtered_variants.vcf \
        --filter-expression "QD < 2.0" --filter-name "QD2" \
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \
        --filter-expression "FS > 60.0" --filter-name "FS60" \
        --filter-expression "SOR > 3.0" --filter-name "SOR3"
    
    PASS_COUNT=$(grep -v "^#" $output_dir/filtered/filtered_variants.vcf | grep -w "PASS" | wc -l)
    echo "✅ $PASS_COUNT variants passed filters"
    
    # 9. ANNOVAR Annotation
    step 9 "ANNOVAR Annotation"
    
    if [ -d "$ANNOVAR_DIR" ]; then
        # Create PASS-only VCF
        PASS_VCF=$output_dir/filtered/filtered_PASS_only.vcf
        grep "^#" $output_dir/filtered/filtered_variants.vcf > $PASS_VCF
        grep -v "^#" $output_dir/filtered/filtered_variants.vcf | grep -w "PASS" >> $PASS_VCF
        
        perl $ANNOVAR_DIR/table_annovar.pl \
            $PASS_VCF \
            $ANNOVAR_DIR/humandb/ \
            -buildver hg19 \
            -out $output_dir/annovar/annotated_${sample_name} \
            -remove \
            -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150,dbnsfp42a \
            -operation g,f,f,f,f \
            -nastring . \
            -vcfinput \
            -polish
        
        # Compress
        if [ -f "$output_dir/annovar/annotated_${sample_name}.hg19_multianno.vcf" ]; then
            bgzip -f $output_dir/annovar/annotated_${sample_name}.hg19_multianno.vcf
            tabix -p vcf $output_dir/annovar/annotated_${sample_name}.hg19_multianno.vcf.gz
            echo "✅ Annotation complete"
        fi
    else
        echo "⚠️  ANNOVAR not found - skipping annotation"
    fi
    
    # 10. Separate by variant type
    step 10 "Variant Type Separation"
    
    ANNOT_FILE=$output_dir/annovar/annotated_${sample_name}.hg19_multianno.txt
    if [ -f "$ANNOT_FILE" ]; then
        SEPARATED_DIR=$output_dir/annovar/separated_by_type
        mkdir -p $SEPARATED_DIR
        
        HEADER=$(head -1 $ANNOT_FILE)
        
        # SNPs
        echo "$HEADER" > $SEPARATED_DIR/SNPs.txt
        tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4)==1 && length($5)==1' >> $SEPARATED_DIR/SNPs.txt
        
        # Insertions
        echo "$HEADER" > $SEPARATED_DIR/Insertions.txt
        tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4) < length($5)' >> $SEPARATED_DIR/Insertions.txt
        
        # Deletions
        echo "$HEADER" > $SEPARATED_DIR/Deletions.txt
        tail -n +2 $ANNOT_FILE | awk -F'\t' 'length($4) > length($5)' >> $SEPARATED_DIR/Deletions.txt
        
        SNP_COUNT=$(($(wc -l < $SEPARATED_DIR/SNPs.txt) - 1))
        INS_COUNT=$(($(wc -l < $SEPARATED_DIR/Insertions.txt) - 1))
        DEL_COUNT=$(($(wc -l < $SEPARATED_DIR/Deletions.txt) - 1))
        
        echo "✅ Separated: SNPs=$SNP_COUNT, Insertions=$INS_COUNT, Deletions=$DEL_COUNT"
    else
        echo "⚠️  Annotation file not found - skipping separation"
    fi
    
    echo ""
    echo "✅ Sample $sample_name: COMPLETE!"
    echo ""
}

# Main execution
echo "════════════════════════════════════════════════════════════"
echo "STEP 1: Auto-Detecting Samples"
echo "════════════════════════════════════════════════════════════"
echo ""

detect_samples "$DATA_DIR"

if [ ${#SAMPLE_PAIRS[@]} -eq 0 ]; then
    echo "❌ No FASTQ pairs found in $DATA_DIR"
    echo ""
    echo "Expected file naming:"
    echo "  sample_R1.fastq.gz + sample_R2.fastq.gz"
    echo "  OR sample_1.fastq.gz + sample_2.fastq.gz"
    echo ""
    echo "Files found:"
    ls -lh $DATA_DIR/
    exit 1
fi

echo "✅ Detected ${#SAMPLE_PAIRS[@]} sample(s):"
echo ""
for sample in "${!SAMPLE_PAIRS[@]}"; do
    IFS=',' read -r r1 r2 <<< "${SAMPLE_PAIRS[$sample]}"
    echo "  Sample: $sample"
    echo "    R1: $(basename $r1)"
    echo "    R2: $(basename $r2)"
    echo ""
done

echo "════════════════════════════════════════════════════════════"
read -p "Start analysis for all samples? (y/n): " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Analysis cancelled."
    exit 0
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "STEP 2: Running Complete Pipeline"
echo "════════════════════════════════════════════════════════════"
echo ""

# Process each sample
counter=1
for sample in "${!SAMPLE_PAIRS[@]}"; do
    IFS=',' read -r r1 r2 <<< "${SAMPLE_PAIRS[$sample]}"
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  Sample $counter of ${#SAMPLE_PAIRS[@]}: $sample"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    analyze_sample "$sample" "$r1" "$r2" $THREADS
    
    counter=$((counter + 1))
done

# Generate master summary
SUMMARY_FILE=$WORK_DIR/MASTER_ANALYSIS_SUMMARY.txt

cat > $SUMMARY_FILE << EOF
╔════════════════════════════════════════════════════════════╗
║          MASTER EXOME ANALYSIS - SUMMARY REPORT           ║
╚════════════════════════════════════════════════════════════╝

Analysis Date: $(date)
Total Samples: ${#SAMPLE_PAIRS[@]}
Threads Used: $THREADS
Reference: hg19

════════════════════════════════════════════════════════════
SAMPLE RESULTS:
════════════════════════════════════════════════════════════

EOF

for sample in "${!SAMPLE_PAIRS[@]}"; do
    RESULT_DIR="$WORK_DIR/results/$sample"
    
    cat >> $SUMMARY_FILE << EOF
Sample: $sample
────────────────────────────────────────────────────────────
EOF
    
    if [ -d "$RESULT_DIR" ]; then
        # Variant counts
        if [ -f "$RESULT_DIR/filtered/filtered_variants.vcf" ]; then
            RAW=$(grep -v "^#" "$RESULT_DIR/variants/raw_variants.vcf" 2>/dev/null | wc -l)
            PASS=$(grep -v "^#" "$RESULT_DIR/filtered/filtered_variants.vcf" | grep -w "PASS" | wc -l)
            
            cat >> $SUMMARY_FILE << EOF
  Raw variants called:     $RAW
  PASS variants:           $PASS
EOF
        fi
        
        # Annotation
        ANNOT_FILE=$(find "$RESULT_DIR/annovar" -name "*.hg19_multianno.txt" 2>/dev/null | head -1)
        if [ -f "$ANNOT_FILE" ]; then
            PATHOGENIC=$(grep -i "pathogenic" "$ANNOT_FILE" 2>/dev/null | wc -l)
            cat >> $SUMMARY_FILE << EOF
  Annotated:               Yes
  Pathogenic variants:     $PATHOGENIC
EOF
        fi
        
        # Separated types
        if [ -d "$RESULT_DIR/annovar/separated_by_type" ]; then
            SNPS=$(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/SNPs.txt" 2>/dev/null || echo 1) - 1))
            INS=$(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/Insertions.txt" 2>/dev/null || echo 1) - 1))
            DEL=$(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/Deletions.txt" 2>/dev/null || echo 1) - 1))
            
            cat >> $SUMMARY_FILE << EOF
  
  Variant Types:
    SNPs:                  $SNPS
    Insertions:            $INS
    Deletions:             $DEL
EOF
        fi
        
        # Storage
        STORAGE=$(du -sh "$RESULT_DIR" 2>/dev/null | cut -f1)
        cat >> $SUMMARY_FILE << EOF
  
  Storage used:            $STORAGE

EOF
    fi
done

cat >> $SUMMARY_FILE << EOF
════════════════════════════════════════════════════════════
OUTPUT LOCATIONS:
════════════════════════════════════════════════════════════

For each sample, download:
  
  📊 Main Results:
    results/SAMPLE/annovar/annotated_SAMPLE.hg19_multianno.txt
  
  📁 Separated by Type:
    results/SAMPLE/annovar/separated_by_type/SNPs.txt
    results/SAMPLE/annovar/separated_by_type/Insertions.txt
    results/SAMPLE/annovar/separated_by_type/Deletions.txt
  
  📈 Quality Reports:
    results/SAMPLE/fastqc/*.html
    results/SAMPLE/trimmed/fastp_report.html

════════════════════════════════════════════════════════════
TOTAL STORAGE:
════════════════════════════════════════════════════════════

$(du -sh $WORK_DIR/results 2>/dev/null)

════════════════════════════════════════════════════════════
EOF

# Display summary
cat $SUMMARY_FILE

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              🎉 MASTER PIPELINE COMPLETE! 🎉              ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Summary saved to: $SUMMARY_FILE"
echo ""
echo "Next steps:"
echo "  1. Download annotated files for Excel analysis"
echo "  2. Review separated variant types (SNPs, Indels)"
echo "  3. Filter for pathogenic variants"
echo ""
echo "════════════════════════════════════════════════════════════"


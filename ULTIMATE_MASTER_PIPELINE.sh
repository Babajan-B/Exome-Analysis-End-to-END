#!/bin/bash
# ULTIMATE MASTER PIPELINE - Complete End-to-End Exome Analysis
# Includes: QC → Alignment → Variant Calling → ANNOVAR → snpEff → Advanced Separation → ZIP
# Usage: bash ULTIMATE_MASTER_PIPELINE.sh [data_directory] [threads]

set -eo pipefail

# Configuration
DATA_DIR=${1:-~/NGS/data}
THREADS=${2:-16}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR=${NGS_WORK_DIR:-$SCRIPT_DIR}
REFERENCE=$WORK_DIR/reference/hg19.fa
# Known-sites are BUILD-SPECIFIC — hg38 known-sites are NOT compatible with an
# hg19 reference (contig lengths differ). Resolve the dir from the reference name
# so BQSR only ever uses matching-build sites (and skips if that build has none).
REF_BUILD=$(basename "$REFERENCE" .fa)   # e.g. hg19 or hg38
KNOWN_SITES_DIR=$WORK_DIR/reference/known-sites/$REF_BUILD
KNOWN_DBSNP=$KNOWN_SITES_DIR/dbsnp.vcf.gz
KNOWN_MILLS=$KNOWN_SITES_DIR/mills.vcf.gz
KNOWN_INDELS=$KNOWN_SITES_DIR/1000G_indels.vcf.gz
ANNOVAR_DIR=$WORK_DIR/tools/annovar
SNPEFF_DIR=$WORK_DIR/tools/snpEff
SNPEFF_DB="GRCh37.75"
GATK=${GATK:-$WORK_DIR/gatk-4.6.2.0/gatk}
if [ ! -x "$GATK" ] && [ -x /opt/gatk-4.6.2.0/gatk ]; then
    GATK=/opt/gatk-4.6.2.0/gatk
fi

echo "╔════════════════════════════════════════════════════════════╗"
echo "║        ULTIMATE EXOME ANALYSIS PIPELINE                   ║"
echo "║    Complete: Pipeline → Annotation → ZIP Results          ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Pipeline Steps:"
echo "  PART 1: CORE PIPELINE"
echo "    1. Auto-detect FASTQ samples"
echo "    2. Quality Control (FastQC)"
echo "    3. Read Trimming (fastp)"
echo "    4. Alignment (BWA-MEM)"
echo "    5. BAM Processing (sort, mark duplicates)"
echo "    6. Variant Calling (GATK HaplotypeCaller)"
echo "    7. Variant Filtering"
echo "    8. ANNOVAR Annotation (5 databases)"
echo "    9. Variant Type Separation (SNPs, Indels)"
echo ""
echo "  PART 2: ADVANCED ANNOTATION"
echo "    10. snpEff Annotation"
echo "    11. Add Zygosity Information"
echo "    12. Functional Classification"
echo ""
echo "  PART 3: FINAL PACKAGING"
echo "    13. Compress VCF files"
echo "    14. Create ZIP Archive"
echo ""
echo "Started: $(date)"
echo ""

# Check data directory
if [ ! -d "$DATA_DIR" ]; then
    echo "❌ Data directory not found: $DATA_DIR"
    echo ""
    echo "Setup:"
    echo "  mkdir -p $DATA_DIR"
    echo "  # Upload your FASTQ files"
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
    SAMPLE_NAMES=()
    SAMPLE_R1S=()
    SAMPLE_R2S=()
    
    cd "$data_dir"
    
    # Enable nullglob to handle no matches gracefully
    shopt -s nullglob
    
    # Find R1 files
    for r1_file in *R1*.fastq.gz *R1*.fq.gz *_1.fastq.gz *_1.fq.gz; do
        [ -f "$r1_file" ] || continue
        
        # Generate expected R2 filename
        r2_file=$(echo "$r1_file" | sed -e 's/R1/R2/g' -e 's/_1\./_2\./g')
        
        if [ -f "$r2_file" ]; then
            # Extract sample name
            sample_name=$(echo "$r1_file" | sed -E 's/[._-]*(R1|_1)[._-]*.*//' | sed -E 's/\.(fastq|fq)\.gz$//')
            
            SAMPLE_NAMES+=("$sample_name")
            SAMPLE_R1S+=("$data_dir/$r1_file")
            SAMPLE_R2S+=("$data_dir/$r2_file")
        fi
    done
    
    # Disable nullglob
    shopt -u nullglob
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
    mkdir -p $output_dir/{fastqc,trimmed,aligned,sorted,dedup,bqsr,variants,filtered,annovar/snpeff,annovar/functional_classification}
    
    # State hygiene: purge stale halt and reasoning artifacts from previous runs
    rm -f "$output_dir/supervisor_reasoning.json" \
          "$output_dir/stage2_alignment_reasoning.json" \
          "$output_dir/halt_report.json" \
          "$output_dir"/*.applied
    
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
    
    # ═══════════════════════════════════════════════════════════════
    # PART 1: CORE PIPELINE
    # ═══════════════════════════════════════════════════════════════
    
    # 1. FastQC
    step 1 "Quality Control"
    if [ -f "$output_dir/.override_qc_gate" ] && [ -d "$output_dir/fastqc" ] && [ "$(ls -A $output_dir/fastqc 2>/dev/null)" ]; then
        echo "  ⚠️  [SUPERVISOR] Operator Override active — FastQC reports already exist. Fast-forwarding..."
    else
        fastqc -t $threads -o $output_dir/fastqc $r1_path $r2_path
    fi
    echo "✅ QC complete"
    
    # 2. Read Trimming & Remediation — Layer 2 Execution Worker & Supervisor Gate
    # Universal execution of fastp with autonomous Supervisor remediation loop.
    # Supervisor evaluates fastp_report.json against qc.json:
    #   - 0: Approved (Clinical or Research Grade Flag & Continue) -> advance
    #   - 42: Remediation required -> re-runs fastp with tighter parameters (retry <= 2)
    #   - 1: Rejection floor or hardware failure -> halts pipeline before Alignment
    step 2 "Read Trimming & Remediation"
    
    # ── Check for Operator Override Fast-Forward ──
    # If the user clicked Override & Force Run and trimmed reads already exist, resume directly without re-computing
    if [ -f "$output_dir/.override_qc_gate" ] && [ -f "$output_dir/trimmed/r1_trimmed.fastq.gz" ] && [ -f "$output_dir/trimmed/fastp_report.json" ]; then
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing trimmed FASTQs."
        echo "  Executing gate evaluation to apply override and resuming straight to Alignment..."
        set +e
        node "$SCRIPT_DIR/scripts/stage1_qc_gate.js" "$output_dir" "$sample_name"
        gate_status=$?
        set -e
        if [ "$gate_status" -eq 0 ]; then
            echo "✅ Stage 1 QC Authorized via Operator Override. Proceeding directly to Alignment."
        else
            echo "❌ [SUPERVISOR] Gate evaluation failed with code $gate_status"
            exit $gate_status
        fi
    else
        # Fresh / Normal run: ensure no lingering retry flags or stale reports from previous runs
        rm -f "$output_dir/trimmed/.retry_count" "$output_dir/trimmed/.remediation_flags" "$output_dir/trimmed/fastp_report.json"
        extra_qc_flags=""
        qc_attempt=0
        max_qc_attempts=3
        while [ $qc_attempt -lt $max_qc_attempts ]; do
            qc_attempt=$((qc_attempt + 1))
            # Clear fastp_report.json before run to prevent stale report if tool crashes
            rm -f "$output_dir/trimmed/fastp_report.json"

            if [ -n "$extra_qc_flags" ]; then
                echo "  [QC AGENT] Executing fastp with Supervisor trimming parameters (Attempt $qc_attempt): $extra_qc_flags"
                fastp_args="--detect_adapter_for_pe $extra_qc_flags"
            else
                echo "  [SUPERVISOR] Profiling raw read quality & adapter content (pass-through)..."
                fastp_args="--detect_adapter_for_pe --disable_quality_filtering"
            fi
            
            set +e
            fastp -i "$r1_path" -I "$r2_path" \
                -o "$output_dir/trimmed/r1_trimmed.fastq.gz" \
                -O "$output_dir/trimmed/r2_trimmed.fastq.gz" \
                --unpaired1 "$output_dir/trimmed/r1_unpaired.fastq.gz" \
                --unpaired2 "$output_dir/trimmed/r2_unpaired.fastq.gz" \
                $fastp_args \
                --thread $threads \
                -h "$output_dir/trimmed/fastp_report.html" \
                -j "$output_dir/trimmed/fastp_report.json"
            fastp_status=$?
            set -e

            if [ $fastp_status -ne 0 ]; then
                echo "❌ [SUPERVISOR] fastp failed with exit code $fastp_status (binary crash or IO failure)."
            fi
            
            # ── Supervisor Gate & Two-Tier Reviewer Evaluation (qc.json) ──
            set +e
            node "$SCRIPT_DIR/scripts/stage1_qc_gate.js" "$output_dir" "$sample_name"
            gate_status=$?
            set -e
            
            if [ "$gate_status" -eq 0 ]; then
                echo "✅ Stage 1 QC Approved by Supervisor. Proceeding to Alignment."
                break
            elif [ "$gate_status" -eq 42 ]; then
                if [ -f "$output_dir/trimmed/.remediation_flags" ]; then
                    extra_qc_flags=$(cat "$output_dir/trimmed/.remediation_flags")
                    echo "  [SUPERVISOR] Quality in 20-30 range. Executing trimming loop (Attempt $qc_attempt of $max_qc_attempts): $extra_qc_flags"
                fi
                continue
            elif [ "$gate_status" -eq 2 ]; then
                echo "❌ [SUPERVISOR] Stage 1 QC Tool Failure: fastp crashed or failed to emit valid reports. Execution halted."
                exit 2
            else
                echo "❌ [SUPERVISOR] Stage 1 QC Critical Biological Rejection (Phred < 20). Operator opinion required."
                exit 1
            fi
        done

        if [ "$gate_status" -ne 0 ]; then
            echo "❌ [SUPERVISOR] Maximum QC remediation attempts ($max_qc_attempts) exceeded without reaching quality threshold."
            exit 1
        fi
    fi
    echo "✅ Trimming complete"
    
    # 3. Alignment & Coordinate Sorting (Piped Streaming)
    step 3 "Read Alignment & Coordinate Sorting"
    mkdir -p $output_dir/aligned $output_dir/sorted $output_dir/dedup
    
    # Layer 1 Storage Custodian: Remove any leftover orphaned SAM/BAM files from previous interrupted runs
    rm -f $output_dir/aligned/aligned.sam $output_dir/aligned/aligned.bam
    
    # ── Check for Operator Override Fast-Forward ──
    if [ -f "$output_dir/.override_align_gate" ] && [ -s "$output_dir/sorted/sorted.bam" ] && [ -s "$output_dir/aligned/alignment_flagstat.txt" ]; then
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing sorted BAM."
        echo "  Executing gate evaluation to apply override and resuming straight to Deduplication..."
        set +e
        node "$SCRIPT_DIR/scripts/stage2_align_gate.js" "$output_dir" "$sample_name"
        ALIGN_EXIT=$?
        set -e
        if [ $ALIGN_EXIT -eq 0 ]; then
            echo "✅ Stage 2 Alignment Authorized via Operator Override. Proceeding directly to Deduplication."
        else
            echo "❌ [SUPERVISOR] Alignment gate evaluation failed with code $ALIGN_EXIT"
            exit $ALIGN_EXIT
        fi
    else
        # Read Stage 1 Supervisor Downstream Directives
        BWA_FLAGS="-M -Y"
        if [ -f "$output_dir/supervisor_reasoning.json" ]; then
            CUSTOM_FLAGS=$(node -e '
                try {
                    const d = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8"));
                    const flags = d?.downstreamDirectives?.bwaFlags;
                    if (Array.isArray(flags) && flags.length) console.log(flags.join(" "));
                    else if (typeof d?.downstreamDirectives?.bwaMemFlags === "string") console.log(d.downstreamDirectives.bwaMemFlags);
                } catch(e) {}
            ' "$output_dir/supervisor_reasoning.json" 2>/dev/null)
            [ -n "$CUSTOM_FLAGS" ] && BWA_FLAGS="$CUSTOM_FLAGS"
        fi
        echo "  [SUPERVISOR DIRECTIVES] Injecting BWA-MEM flags: $BWA_FLAGS"
        echo "  [STREAMING ALIGNMENT] Piping BWA-MEM directly to samtools sort (0 GB uncompressed SAM on disk)..."
        
        set +e
        bwa mem -t $threads $BWA_FLAGS \
            -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:lib_${sample_name}\tPU:unit1" \
            $REFERENCE \
            $output_dir/trimmed/r1_trimmed.fastq.gz \
            $output_dir/trimmed/r2_trimmed.fastq.gz \
            | samtools sort -@ $threads -o $output_dir/sorted/sorted.bam -
        bwa_pipe_status=( "${PIPESTATUS[@]}" )
        set -e
        
        if [ ${bwa_pipe_status[0]:-0} -ne 0 ] || [ ${bwa_pipe_status[1]:-0} -ne 0 ] || [ ! -s "$output_dir/sorted/sorted.bam" ]; then
            echo "❌ Alignment or sorting failed (BWA exit: ${bwa_pipe_status[0]:-0}, samtools sort exit: ${bwa_pipe_status[1]:-0}). Check disk space and memory."
            exit 1
        fi
        
        samtools index $output_dir/sorted/sorted.bam
        echo "✅ Alignment & coordinate sorting complete"
        
        # 4. Alignment Metrics
        step 4 "Alignment Metrics"
        echo "  [ALIGNMENT METRICS] Computing flagstat, stats, and idxstats..."
        samtools flagstat -@ $threads $output_dir/sorted/sorted.bam > $output_dir/aligned/alignment_flagstat.txt
        samtools stats -@ $threads $output_dir/sorted/sorted.bam > $output_dir/aligned/alignment_stats.txt
        samtools idxstats $output_dir/sorted/sorted.bam > $output_dir/aligned/alignment_idxstats.txt
        echo "✅ Alignment metrics generated"
        
        # Stage 2: Alignment Cognitive Supervisor Gate
        echo "  [SUPERVISOR] Evaluating Stage 2 Alignment metrics against clinical policies (qc.json)..."
        set +e
        node "$SCRIPT_DIR/scripts/stage2_align_gate.js" "$output_dir" "$sample_name"
        ALIGN_EXIT=$?
        set -e
        
        if [ $ALIGN_EXIT -ne 0 ]; then
            if [ $ALIGN_EXIT -eq 1 ]; then
                echo ""
                echo "🛑 [PIPELINE HALTED] Alignment Quality Gate failed rejection floor (Mapping < 90% or Pairing < 85%)."
                echo "   Human-in-the-Loop Operator Opinion Gate is required."
                echo "   Use the Web Dashboard to either:"
                echo "     1. [Abort & Re-sequence] (Recommended clinical action)"
                echo "     2. [Override & Force Run] (High-Risk Research Mode)"
                echo "   Or run: touch \"$output_dir/.override_align_gate\" and restart pipeline."
                exit 1
            elif [ $ALIGN_EXIT -eq 2 ]; then
                echo "❌ [TOOL CRASH] Alignment output missing or corrupted. Pipeline halted."
                exit 2
            else
                echo "❌ [ERROR] Unknown alignment gate error ($ALIGN_EXIT)."
                exit 1
            fi
        fi
    fi
    echo "✅ Stage 2 Alignment Approved by Supervisor. Proceeding to Deduplication."
    
    # 5. Duplicate Marking (formerly step 6)
    step 5 "Mark Duplicates"
    $GATK MarkDuplicates \
        -I $output_dir/sorted/sorted.bam \
        -O $output_dir/dedup/dedup.bam \
        -M $output_dir/dedup/metrics.txt \
        --CREATE_INDEX true
    
    # Ensure index exists as both dedup.bai and dedup.bam.bai for universal tool compatibility
    if [ -f "$output_dir/dedup/dedup.bai" ] && [ ! -f "$output_dir/dedup/dedup.bam.bai" ]; then
        cp "$output_dir/dedup/dedup.bai" "$output_dir/dedup/dedup.bam.bai"
    elif [ -f "$output_dir/dedup/dedup.bam.bai" ] && [ ! -f "$output_dir/dedup/dedup.bai" ]; then
        cp "$output_dir/dedup/dedup.bam.bai" "$output_dir/dedup/dedup.bai"
    fi
    echo "✅ Duplicates marked"
    
    # Layer 1 Storage Custodian: Safe Space Reclamation
    # Once dedup.bam is verified (>10 KB) and indexed, retire intermediate sorted.bam to save ~15-20 GB
    DEDUP_SIZE=$(wc -c < "$output_dir/dedup/dedup.bam" 2>/dev/null || echo 0)
    if [ "$DEDUP_SIZE" -gt 10240 ] && [ -f "$output_dir/dedup/dedup.bam.bai" -o -f "$output_dir/dedup/dedup.bai" ]; then
        RECLAIM_KB=$(du -sk "$output_dir/sorted/sorted.bam" 2>/dev/null | awk '{print $1}')
        echo "🧹 [L1 STORAGE CUSTODIAN] Reclaiming disk space: retiring intermediate sorted.bam (${RECLAIM_KB:-0} KB)..."
        rm -f $output_dir/sorted/sorted.bam $output_dir/sorted/sorted.bam.bai $output_dir/sorted/sorted.bai
        echo "   Active validated BAM for Genome Viewer & BQSR: dedup/dedup.bam"
    fi
    
    # 7. Base Quality Score Recalibration (BQSR)
    step 7 "Base Quality Score Recalibration"
    BQSR_INPUT=$output_dir/dedup/dedup.bam
    if [ -f "$KNOWN_DBSNP" ]; then
        KS_ARGS="--known-sites $KNOWN_DBSNP"
        [ -f "$KNOWN_MILLS" ] && KS_ARGS="$KS_ARGS --known-sites $KNOWN_MILLS"
        [ -f "$KNOWN_INDELS" ] && KS_ARGS="$KS_ARGS --known-sites $KNOWN_INDELS"
        $GATK BaseRecalibrator \
            -I $output_dir/dedup/dedup.bam \
            -R $REFERENCE \
            $KS_ARGS \
            -O $output_dir/bqsr/recal_data.table
        $GATK ApplyBQSR \
            -I $output_dir/dedup/dedup.bam \
            -R $REFERENCE \
            --bqsr-recal-file $output_dir/bqsr/recal_data.table \
            -O $output_dir/bqsr/recal.bam
        BQSR_INPUT=$output_dir/bqsr/recal.bam
        echo "✅ BQSR complete — recalibrated BAM ready for calling"
    else
        echo "⚠️  Known-sites (dbSNP/Mills) not found in reference/known-sites/ — SKIPPING BQSR; using dedup.bam."
        echo "    Install known-sites VCFs to enable Base Quality Score Recalibration."
    fi
    
    # 7. Variant Calling
    step 8 "Variant Calling"
    $GATK HaplotypeCaller \
        -R $REFERENCE \
        -I $BQSR_INPUT \
        -O $output_dir/variants/raw_variants.vcf \
        --native-pair-hmm-threads $threads
    
    RAW_COUNT=$(grep -v "^#" $output_dir/variants/raw_variants.vcf | wc -l)
    echo "✅ Called $RAW_COUNT variants"
    
    # 8. Filtering
    step 9 "Variant Filtering"
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
    
    # Create PASS-only VCF
    PASS_VCF=$output_dir/filtered/filtered_PASS_only.vcf
    grep "^#" $output_dir/filtered/filtered_variants.vcf > $PASS_VCF
    grep -v "^#" $output_dir/filtered/filtered_variants.vcf | grep -w "PASS" >> $PASS_VCF
    
    # 9. ANNOVAR Annotation
    step 10 "ANNOVAR Annotation"
    
    if [ -d "$ANNOVAR_DIR" ]; then
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
        
        echo "✅ ANNOVAR annotation complete"
    else
        echo "⚠️  ANNOVAR not found - skipping"
    fi
    
    # 10. Basic Variant Type Separation
    step 11 "Variant Type Separation"
    
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
    fi
    
    # ═══════════════════════════════════════════════════════════════
    # PART 2: ADVANCED ANNOTATION
    # ═══════════════════════════════════════════════════════════════
    
    # 11. snpEff Annotation
    step 12 "snpEff Annotation"
    
    if [ -f "$SNPEFF_DIR/snpEff.jar" ] && [ -f "$PASS_VCF" ]; then
        java -Xmx8g -jar $SNPEFF_DIR/snpEff.jar \
            -v $SNPEFF_DB \
            -stats $output_dir/annovar/snpeff/${sample_name}_snpEff_summary.html \
            -csvStats $output_dir/annovar/snpeff/${sample_name}_snpEff_summary.csv \
            $PASS_VCF \
            > $output_dir/annovar/snpeff/${sample_name}_snpEff_annotated.vcf
        
        echo "✅ snpEff annotation complete"

        # 11b. SnpSift dbSNP CAF frequency annotation
        if [ -f "$SNPEFF_DIR/SnpSift.jar" ] && [ -f "$KNOWN_DBSNP" ]; then
            echo "  Annotating global allele frequency (CAF/COMMON) from dbSNP via SnpSift..."
            java -jar $SNPEFF_DIR/SnpSift.jar annotate \
                -tabix \
                -info CAF,COMMON \
                $KNOWN_DBSNP \
                $output_dir/annovar/snpeff/${sample_name}_snpEff_annotated.vcf \
                > $output_dir/annovar/snpeff/${sample_name}_final_annotated.vcf
            
            if [ -s "$output_dir/annovar/snpeff/${sample_name}_final_annotated.vcf" ]; then
                mv $output_dir/annovar/snpeff/${sample_name}_final_annotated.vcf $output_dir/annovar/snpeff/${sample_name}_snpEff_annotated.vcf
                echo "✅ SnpSift global frequency annotation complete"
            else
                rm -f $output_dir/annovar/snpeff/${sample_name}_final_annotated.vcf
                echo "⚠️  SnpSift output empty, keeping snpEff VCF"
            fi
        fi
    else
        echo "⚠️  snpEff not found - skipping"
    fi
    
    # 12. Add Zygosity Information
    step 13 "Adding Zygosity Information"
    
    ANNOT_TXT=$output_dir/annovar/annotated_${sample_name}.hg19_multianno.txt
    ANNOT_WITH_ZYG=$output_dir/annovar/annotated_${sample_name}_with_zygosity.txt
    
    if [ -f "$ANNOT_TXT" ] && [ -f "$PASS_VCF" ]; then
        cat > /tmp/add_zygosity_${sample_name}.py << 'PYTHON_SCRIPT'
import sys

if len(sys.argv) != 4:
    print("Usage: script.py vcf_file annot_file output_file")
    sys.exit(1)

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

print("✅ Zygosity column added")
PYTHON_SCRIPT
        
        python3 /tmp/add_zygosity_${sample_name}.py "$PASS_VCF" "$ANNOT_TXT" "$ANNOT_WITH_ZYG"
        rm -f /tmp/add_zygosity_${sample_name}.py
    fi
    
    # 13. Advanced Functional Separation
    step 14 "Advanced Functional Separation"
    
    FUNC_DIR=$output_dir/annovar/functional_classification
    mkdir -p $FUNC_DIR
    
    INPUT_FILE=$ANNOT_WITH_ZYG
    if [ ! -f "$INPUT_FILE" ]; then
        INPUT_FILE=$ANNOT_TXT
    fi
    
    if [ -f "$INPUT_FILE" ]; then
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
        
        SNP_EXONIC=$(($(wc -l < $FUNC_DIR/SNPs_Exonic.txt) - 1))
        SNP_NONEXONIC=$(($(wc -l < $FUNC_DIR/SNPs_NonExonic.txt) - 1))
        NONSYN=$(($(wc -l < $FUNC_DIR/Exonic_Nonsynonymous.txt) - 1))
        SYN=$(($(wc -l < $FUNC_DIR/Exonic_Synonymous.txt) - 1))
        STOP=$(($(wc -l < $FUNC_DIR/Exonic_Stopgain.txt) - 1))
        FRAME=$(($(wc -l < $FUNC_DIR/Exonic_Frameshift.txt) - 1))
        
        echo "✅ Functional classification:"
        echo "   SNPs Exonic: $SNP_EXONIC | Non-Exonic: $SNP_NONEXONIC"
        echo "   Nonsynonymous: $NONSYN | Synonymous: $SYN"
        echo "   Stopgain: $STOP | Frameshift: $FRAME"
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

if [ ${#SAMPLE_NAMES[@]} -eq 0 ]; then
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

echo "✅ Detected ${#SAMPLE_NAMES[@]} sample(s):"
echo ""
for idx in "${!SAMPLE_NAMES[@]}"; do
    sample="${SAMPLE_NAMES[$idx]}"
    r1="${SAMPLE_R1S[$idx]}"
    r2="${SAMPLE_R2S[$idx]}"
    echo "  Sample: $sample"
    echo "    R1: $(basename $r1)"
    echo "    R2: $(basename $r2)"
    echo ""
done

echo "════════════════════════════════════════════════════════════"
read -p "Start complete analysis for all samples? (y/n): " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Analysis cancelled."
    exit 0
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "STEP 2: Running Complete Pipeline for Each Sample"
echo "════════════════════════════════════════════════════════════"
echo ""

# Process each sample
counter=1
for idx in "${!SAMPLE_NAMES[@]}"; do
    sample="${SAMPLE_NAMES[$idx]}"
    r1="${SAMPLE_R1S[$idx]}"
    r2="${SAMPLE_R2S[$idx]}"
    
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  Sample $counter of ${#SAMPLE_NAMES[@]}: $sample"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    analyze_sample "$sample" "$r1" "$r2" $THREADS
    
    counter=$((counter + 1))
done

# ═══════════════════════════════════════════════════════════════
# PART 3: FINAL PACKAGING
# ═══════════════════════════════════════════════════════════════

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║             FINAL STEP: CREATING ZIP ARCHIVE               ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

cd $WORK_DIR

# Get list of samples for ZIP processing
SAMPLES=($(ls -d results/*/ 2>/dev/null | xargs -n 1 basename))

ZIP_NAME="NGS_Results_Complete_$(date +%Y%m%d_%H%M%S).zip"

echo "Step 1: Compressing annotated VCF files..."
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_DIR=$WORK_DIR/results/$SAMPLE
    
    # Compress snpEff VCF
    if [ -f "$SAMPLE_DIR/annovar/snpeff/${SAMPLE}_snpEff_annotated.vcf" ]; then
        echo "  Compressing ${SAMPLE} snpEff VCF..."
        bgzip -f $SAMPLE_DIR/annovar/snpeff/${SAMPLE}_snpEff_annotated.vcf
        tabix -p vcf $SAMPLE_DIR/annovar/snpeff/${SAMPLE}_snpEff_annotated.vcf.gz
    fi
    
    # Compress ANNOVAR VCF if not already compressed
    if [ -f "$SAMPLE_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf" ]; then
        echo "  Compressing ${SAMPLE} ANNOVAR VCF..."
        bgzip -f $SAMPLE_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf
        tabix -p vcf $SAMPLE_DIR/annovar/annotated_${SAMPLE}.hg19_multianno.vcf.gz
    fi
done

echo ""
echo "Step 2: Creating ZIP archive..."
echo ""

# Generate master summary
SUMMARY_FILE=$WORK_DIR/MASTER_ANALYSIS_SUMMARY.txt

cat > $SUMMARY_FILE << EOF
╔════════════════════════════════════════════════════════════╗
║       ULTIMATE EXOME ANALYSIS - COMPLETE SUMMARY           ║
╚════════════════════════════════════════════════════════════╝

Analysis Date: $(date)
Total Samples: ${#SAMPLES[@]}
Threads Used: $THREADS
Reference: hg19

════════════════════════════════════════════════════════════
SAMPLE RESULTS:
════════════════════════════════════════════════════════════

EOF

for SAMPLE in "${SAMPLES[@]}"; do
    RESULT_DIR="$WORK_DIR/results/$SAMPLE"
    
    cat >> $SUMMARY_FILE << EOF
Sample: $SAMPLE
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
  Annotated (ANNOVAR):     Yes
  Pathogenic variants:     $PATHOGENIC
EOF
        fi
        
        # snpEff
        if [ -f "$RESULT_DIR/annovar/snpeff/${SAMPLE}_snpEff_summary.html" ]; then
            cat >> $SUMMARY_FILE << EOF
  Annotated (snpEff):      Yes
EOF
        fi
        
        # Separated types
        if [ -d "$RESULT_DIR/annovar/separated_by_type" ]; then
            SNPS=$([ -f "$RESULT_DIR/annovar/separated_by_type/SNPs.txt" ] && echo $(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/SNPs.txt") - 1)) || echo 0)
            INS=$([ -f "$RESULT_DIR/annovar/separated_by_type/Insertions.txt" ] && echo $(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/Insertions.txt") - 1)) || echo 0)
            DEL=$([ -f "$RESULT_DIR/annovar/separated_by_type/Deletions.txt" ] && echo $(($(wc -l < "$RESULT_DIR/annovar/separated_by_type/Deletions.txt") - 1)) || echo 0)
            
            cat >> $SUMMARY_FILE << EOF
  
  Variant Types:
    SNPs:                  $SNPS
    Insertions:            $INS
    Deletions:             $DEL
EOF
        fi
        
        # Functional classification
        if [ -d "$RESULT_DIR/annovar/functional_classification" ]; then
            NONSYN=$([ -f "$RESULT_DIR/annovar/functional_classification/Exonic_Nonsynonymous.txt" ] && echo $(($(wc -l < "$RESULT_DIR/annovar/functional_classification/Exonic_Nonsynonymous.txt") - 1)) || echo 0)
            SYN=$([ -f "$RESULT_DIR/annovar/functional_classification/Exonic_Synonymous.txt" ] && echo $(($(wc -l < "$RESULT_DIR/annovar/functional_classification/Exonic_Synonymous.txt") - 1)) || echo 0)
            STOP=$([ -f "$RESULT_DIR/annovar/functional_classification/Exonic_Stopgain.txt" ] && echo $(($(wc -l < "$RESULT_DIR/annovar/functional_classification/Exonic_Stopgain.txt") - 1)) || echo 0)
            
            cat >> $SUMMARY_FILE << EOF
  
  Functional Classification:
    Nonsynonymous:         $NONSYN
    Synonymous:            $SYN
    Stopgain:              $STOP
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
OUTPUT LOCATIONS (in ZIP):
════════════════════════════════════════════════════════════

For each sample:
  
  📊 Main Results:
    results/SAMPLE/annovar/annotated_SAMPLE.hg19_multianno.txt
    results/SAMPLE/annovar/annotated_SAMPLE_with_zygosity.txt
  
  📁 VCF Files (compressed):
    results/SAMPLE/annovar/annotated_SAMPLE.hg19_multianno.vcf.gz
    results/SAMPLE/annovar/snpeff/SAMPLE_snpEff_annotated.vcf.gz
  
  📁 Separated by Type:
    results/SAMPLE/annovar/separated_by_type/*.txt
  
  📁 Functional Classification:
    results/SAMPLE/annovar/functional_classification/*.txt
  
  📈 Quality Reports:
    results/SAMPLE/fastqc/*.html
    results/SAMPLE/trimmed/fastp_report.html
  
  📊 snpEff Reports:
    results/SAMPLE/annovar/snpeff/SAMPLE_snpEff_summary.html

════════════════════════════════════════════════════════════
TOTAL STORAGE:
════════════════════════════════════════════════════════════

$(du -sh $WORK_DIR/results 2>/dev/null)

════════════════════════════════════════════════════════════
EOF

# Create ZIP
zip -r $ZIP_NAME \
    results/*/annovar/*.txt \
    results/*/annovar/snpeff/*.html \
    results/*/annovar/snpeff/*.csv \
    results/*/annovar/snpeff/*.vcf.gz \
    results/*/annovar/snpeff/*.vcf.gz.tbi \
    results/*/annovar/*.hg19_multianno.vcf.gz \
    results/*/annovar/*.hg19_multianno.vcf.gz.tbi \
    results/*/annovar/separated_by_type/*.txt \
    results/*/annovar/functional_classification/*.txt \
    results/*/fastqc/*.html \
    results/*/trimmed/fastp_report.html \
    MASTER_ANALYSIS_SUMMARY.txt \
    -x "*.bam" "*.sam" "*.fastq.gz" "*.avinput" "*_dropped" "*_filtered" "raw_variants.vcf" "filtered_variants.vcf" "filtered_PASS_only.vcf" \
    2>/dev/null

ZIP_SIZE=$(du -sh "$ZIP_NAME" 2>/dev/null | cut -f1 || echo "N/A")
[ -z "$ZIP_SIZE" ] && ZIP_SIZE="N/A"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║        🎉 ULTIMATE PIPELINE COMPLETE! 🎉                  ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Archive created: $ZIP_NAME ($ZIP_SIZE)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "WHAT'S INCLUDED:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "✅ ANNOVAR annotation (.txt files)"
echo "✅ ANNOVAR with zygosity information"
echo "✅ ANNOVAR annotated VCF (compressed + indexed)"
echo "✅ snpEff annotated VCF (compressed + indexed)"
echo "✅ snpEff reports (HTML + CSV)"
echo "✅ Variant type separation (SNPs, Insertions, Deletions)"
echo "✅ Functional classification:"
echo "   - SNPs Exonic / Non-Exonic"
echo "   - Nonsynonymous / Synonymous"
echo "   - Stopgain / Frameshift"
echo "✅ Quality control reports (FastQC, fastp)"
echo "✅ Complete summary report"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "DOWNLOAD:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  ~/NGS/$ZIP_NAME"
echo ""
echo "VCF files are compressed and ready for IGV viewing!"
echo ""
echo "════════════════════════════════════════════════════════════"
echo ""

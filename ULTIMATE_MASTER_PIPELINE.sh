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
    
    # State hygiene: purge stale halt reports and dotfile override markers
    rm -f "$output_dir/halt_report.json" "$output_dir"/.*.applied "$output_dir"/*.applied 2>/dev/null || true
    
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
    
    # ── Cascading Checkpoint Evaluators ──
    # Checks if Stage 5 (Variant Calling & Quality Filtering) is already approved
    is_stage5_approved() {
        local has_vcf=0
        if [ -s "$output_dir/filtered/filtered_variants.vcf" ] || [ -s "$output_dir/filtered/filtered_PASS_only.vcf" ]; then
            has_vcf=1
        fi
        if [ $has_vcf -eq 1 ]; then
            local s5_file="$output_dir/qc/stage5_variant_reasoning.json"
            [ ! -f "$s5_file" ] && s5_file="$output_dir/stage5_variant_reasoning.json"
            if [ -f "$s5_file" ]; then
                local s5_tier
                s5_tier=$(grep '"tier":' "$s5_file" 2>/dev/null || echo "")
                if [[ "$s5_tier" == *"TIER_1"* || "$s5_tier" == *"RESEARCH"* || "$s5_tier" == *"OVERRID"* ]]; then
                    return 0
                fi
            fi
        fi
        return 1
    }

    # Checks if Stage 4 (BQSR) is already approved (either recal.bam or dedup.bam if bypassed)
    is_stage4_approved() {
        if is_stage5_approved; then
            return 0
        fi
        local has_recal=0
        if [ -s "$output_dir/bqsr/recal.bam" ] && [ -s "$output_dir/bqsr/recal_data.table" ]; then
            if [ -f "$output_dir/bqsr/recal.bam.bai" ] || [ -f "$output_dir/bqsr/recal.bai" ]; then
                has_recal=1
            fi
        elif [ -f "$output_dir/.bqsr_bypassed" ] && [ -s "$output_dir/dedup/dedup.bam" ]; then
            has_recal=1
        fi
        if [ $has_recal -eq 1 ] && [ -f "$output_dir/stage4_bqsr_reasoning.json" ]; then
            local s4_tier
            s4_tier=$(grep '"tier":' "$output_dir/stage4_bqsr_reasoning.json" 2>/dev/null || echo "")
            if [[ "$s4_tier" == *"TIER_1"* || "$s4_tier" == *"RESEARCH"* || "$s4_tier" == *"OVERRID"* ]]; then
                return 0
            fi
        fi
        return 1
    }

    # Checks if Stage 3 (Deduplication) is already approved
    is_stage3_approved() {
        if is_stage4_approved || is_stage5_approved; then
            return 0
        fi
        local has_dedup=0
        if [ -s "$output_dir/dedup/dedup.bam" ] && [ -s "$output_dir/dedup/metrics.txt" ]; then
            if [ -f "$output_dir/dedup/dedup.bam.bai" ] || [ -f "$output_dir/dedup/dedup.bai" ]; then
                has_dedup=1
            fi
        fi
        if [ $has_dedup -eq 1 ] && [ -f "$output_dir/stage3_dedup_reasoning.json" ]; then
            local s3_tier
            s3_tier=$(grep '"tier":' "$output_dir/stage3_dedup_reasoning.json" 2>/dev/null || echo "")
            if [[ "$s3_tier" == *"TIER_1"* || "$s3_tier" == *"RESEARCH"* || "$s3_tier" == *"OVERRID"* ]]; then
                return 0
            fi
        fi
        return 1
    }

    # Checks if Stage 2 (Alignment) is already approved
    is_stage2_approved() {
        if is_stage3_approved || is_stage4_approved || is_stage5_approved; then
            return 0
        fi
        local has_alignment_bam=0
        if [ -s "$output_dir/sorted/sorted.bam" ] || [ -s "$output_dir/dedup/dedup.bam" ] || [ -s "$output_dir/bqsr/recal.bam" ]; then
            has_alignment_bam=1
        fi
        if [ $has_alignment_bam -eq 1 ] && [ -s "$output_dir/aligned/alignment_flagstat.txt" ]; then
            if [ -f "$output_dir/stage2_alignment_reasoning.json" ]; then
                local s2_tier
                s2_tier=$(grep '"tier":' "$output_dir/stage2_alignment_reasoning.json" 2>/dev/null || echo "")
                if [[ "$s2_tier" == *"TIER_1"* || "$s2_tier" == *"RESEARCH"* || "$s2_tier" == *"OVERRID"* ]]; then
                    return 0
                fi
            fi
        fi
        return 1
    }

    # Checks if Stage 1 (QC & Trimming) is already approved
    is_stage1_approved() {
        if is_stage2_approved || is_stage3_approved || is_stage4_approved || is_stage5_approved; then
            return 0
        fi
        if [ -s "$output_dir/trimmed/r1_trimmed.fastq.gz" ] && [ -s "$output_dir/trimmed/fastp_report.json" ]; then
            if [ -f "$output_dir/supervisor_reasoning.json" ]; then
                local s1_tier
                s1_tier=$(grep '"tier":' "$output_dir/supervisor_reasoning.json" 2>/dev/null || echo "")
                if [[ "$s1_tier" == *"TIER_1"* || "$s1_tier" == *"RESEARCH"* || "$s1_tier" == *"OVERRID"* ]]; then
                    return 0
                fi
            fi
        fi
        return 1
    }

    # ── Checkpoint: Stage 1 (QC & Trimming) ──
    stage1_approved=0
    if is_stage1_approved; then
        stage1_approved=1
    elif [ -f "$output_dir/.override_qc_gate" ] && [ -s "$output_dir/trimmed/r1_trimmed.fastq.gz" ] && [ -s "$output_dir/trimmed/fastp_report.json" ]; then
        # Run the gate so it records TIER_3_OPERATOR_OVERRIDDEN, audits the decision, and archives the marker;
        # otherwise the reasoning stays HALT and a later restart re-halts here.
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing trimmed FASTQs. Applying override via Stage 1 gate..."
        set +e
        node "$SCRIPT_DIR/scripts/stage1_qc_gate.js" "$output_dir" "$sample_name"
        S1_OVERRIDE_EXIT=$?
        set -e
        if [ $S1_OVERRIDE_EXIT -ne 0 ]; then
            echo "❌ [SUPERVISOR] Stage 1 gate failed to apply override (code $S1_OVERRIDE_EXIT)"
            exit $S1_OVERRIDE_EXIT
        fi
        stage1_approved=1
    fi

    if [ $stage1_approved -eq 1 ]; then
        echo "⏩ [CHECKPOINT] Stage 1 (QC & Trimming) already completed and approved. Skipping to Alignment..."
    else
        step 1 "Quality Control"
        fastqc -t $threads -o $output_dir/fastqc $r1_path $r2_path
        echo "✅ QC complete"

        step 2 "Read Trimming & Remediation"
        rm -f "$output_dir/trimmed/.retry_count" "$output_dir/trimmed/.remediation_flags" "$output_dir/trimmed/fastp_report.json"
        extra_qc_flags=""
        qc_attempt=0
        max_qc_attempts=3
        while [ $qc_attempt -lt $max_qc_attempts ]; do
            qc_attempt=$((qc_attempt + 1))
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
        echo "✅ Trimming complete"
    fi
    
    # ── Checkpoint: Stage 2 (Alignment & Coordinate Sorting) ──
    stage2_approved=0
    if is_stage2_approved; then
        stage2_approved=1
    elif [ -f "$output_dir/.override_align_gate" ] && [ -s "$output_dir/sorted/sorted.bam" ] && [ -s "$output_dir/aligned/alignment_flagstat.txt" ]; then
        # Run the gate so it records TIER_3_OPERATOR_OVERRIDDEN, audits the decision, and archives the marker;
        # otherwise the reasoning stays HALT and a later restart re-aligns and re-halts.
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing sorted BAM. Applying override via Stage 2 gate..."
        set +e
        node "$SCRIPT_DIR/scripts/stage2_align_gate.js" "$output_dir" "$sample_name"
        S2_OVERRIDE_EXIT=$?
        set -e
        if [ $S2_OVERRIDE_EXIT -ne 0 ]; then
            echo "❌ [SUPERVISOR] Stage 2 gate failed to apply override (code $S2_OVERRIDE_EXIT)"
            exit $S2_OVERRIDE_EXIT
        fi
        stage2_approved=1
    fi

    if [ $stage2_approved -eq 1 ]; then
        echo "⏩ [CHECKPOINT] Stage 2 (Alignment & Coordinate Sorting) already completed and approved. Skipping to Deduplication..."
    else
        step 3 "Read Alignment & Coordinate Sorting"
        mkdir -p $output_dir/aligned $output_dir/sorted $output_dir/dedup
        rm -f $output_dir/aligned/aligned.sam $output_dir/aligned/aligned.bam

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

        step 4 "Alignment Metrics"
        echo "  [ALIGNMENT METRICS] Computing flagstat, stats, and idxstats..."
        samtools flagstat -@ $threads $output_dir/sorted/sorted.bam > $output_dir/aligned/alignment_flagstat.txt
        samtools stats -@ $threads $output_dir/sorted/sorted.bam > $output_dir/aligned/alignment_stats.txt
        samtools idxstats $output_dir/sorted/sorted.bam > $output_dir/aligned/alignment_idxstats.txt
        echo "✅ Alignment metrics generated"

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
        echo "✅ Stage 2 Alignment Approved by Supervisor. Proceeding to Deduplication."
    fi
    
    # ── Checkpoint: Stage 3 (Deduplication) ──
    stage3_approved=0
    if is_stage3_approved; then
        stage3_approved=1
    elif [ -f "$output_dir/.override_dedup_gate" ] && [ -s "$output_dir/dedup/dedup.bam" ]; then
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing deduplicated BAM."
        echo "  Executing gate evaluation to apply override and resuming straight to BQSR..."
        set +e
        node "$SCRIPT_DIR/scripts/stage3_dedup_gate.js" "$output_dir" "$sample_name"
        DEDUP_EXIT=$?
        set -e
        if [ $DEDUP_EXIT -eq 0 ]; then
            echo "✅ Stage 3 Deduplication Authorized via Operator Override. Proceeding directly to BQSR."
            stage3_approved=1
        else
            echo "❌ [SUPERVISOR] Deduplication gate evaluation failed with code $DEDUP_EXIT"
            exit $DEDUP_EXIT
        fi
    fi

    if [ $stage3_approved -eq 1 ]; then
        echo "⏩ [CHECKPOINT] Stage 3 (Deduplication) already completed and approved. Skipping MarkDuplicates to BQSR..."
    else
        step 5 "Mark Duplicates"
        
        # Dynamic Optical Pixel Distance Detection (Patterned vs Non-patterned Flowcells)
        DETECTED_OPTICAL_DIST=2500
        SAMPLE_R1_FASTQ="$output_dir/trimmed/r1_trimmed.fastq.gz"
        [ ! -f "$SAMPLE_R1_FASTQ" ] && SAMPLE_R1_FASTQ="$r1_path"

        if [ -f "$SAMPLE_R1_FASTQ" ]; then
            FIRST_READ_HEADER=$(gzip -dc "$SAMPLE_R1_FASTQ" 2>/dev/null | head -n 1 || true)
            # Standard Illumina 7-field header: @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>
            if [[ "$FIRST_READ_HEADER" =~ ^@([^:]+):[0-9]+:([^:]+):[0-9]+:[0-9]+:[0-9]+:[0-9]+ ]]; then
                INSTRUMENT_ID="${BASH_REMATCH[1]}"
                # Non-patterned flowcell instruments: MiSeq (M), MiniSeq (MN), NextSeq 500/550 (NB, NS), HiSeq 2000/2500 (D, SN, HWI-), GAIIx (HWUSI-)
                if [[ "$INSTRUMENT_ID" =~ ^(M[0-9]|MN[0-9]|NB[0-9]|NS[0-9]|D[0-9]|SN[0-9]|HWI-|HWUSI-) ]]; then
                    DETECTED_OPTICAL_DIST=100
                    echo "  ℹ️  [FLOWCELL DETECTION] Detected non-patterned flowcell (Instrument: $INSTRUMENT_ID). Setting optical distance to 100 pixels."
                else
                    echo "  ℹ️  [FLOWCELL DETECTION] Detected patterned flowcell (Instrument: $INSTRUMENT_ID). Setting optical distance to 2500 pixels."
                fi
            elif [ -n "$FIRST_READ_HEADER" ]; then
                DETECTED_OPTICAL_DIST=100
                echo "  ℹ️  [FLOWCELL DETECTION] Non-standard / SRA read headers detected. Optical distance set to 100 (optical duplicates unmeasurable from coordinates)."
            fi
        fi

        # Read Stage 2 Downstream Directives (extra picard flags)
        PICARD_EXTRA=""
        if [ -f "$output_dir/stage2_alignment_reasoning.json" ]; then
            CUSTOM_PICARD=$(node -e '
                try {
                    const d = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8"));
                    const flags = d?.downstreamDirectives?.picardFlags;
                    if (Array.isArray(flags) && flags.length) console.log(flags.join(" "));
                } catch(e) {}
            ' "$output_dir/stage2_alignment_reasoning.json" 2>/dev/null)
            [ -n "$CUSTOM_PICARD" ] && PICARD_EXTRA="$CUSTOM_PICARD"
        fi
        # --CREATE_INDEX and the optical distance are set here only; GATK rejects any argument given twice,
        # and the flowcell-detected distance must win over any directive value.
        CLEAN_PICARD_EXTRA=$(echo "$PICARD_EXTRA" | sed -E 's/--CREATE_INDEX[ =]+(true|false)//g; s/--OPTICAL_DUPLICATE_PIXEL_DISTANCE[ =]+[0-9]+//g')
        CLEAN_PICARD_EXTRA="$CLEAN_PICARD_EXTRA --OPTICAL_DUPLICATE_PIXEL_DISTANCE $DETECTED_OPTICAL_DIST"
        echo "  [SUPERVISOR DIRECTIVES] Injecting MarkDuplicates parameters: $CLEAN_PICARD_EXTRA"
        
        $GATK MarkDuplicates \
            -I $output_dir/sorted/sorted.bam \
            -O $output_dir/dedup/dedup.bam \
            -M $output_dir/dedup/metrics.txt \
            --CREATE_INDEX true \
            $CLEAN_PICARD_EXTRA
        
        # Ensure index exists as both dedup.bai and dedup.bam.bai for universal tool compatibility
        if [ -f "$output_dir/dedup/dedup.bai" ] && [ ! -f "$output_dir/dedup/dedup.bam.bai" ]; then
            cp "$output_dir/dedup/dedup.bai" "$output_dir/dedup/dedup.bam.bai"
        elif [ -f "$output_dir/dedup/dedup.bam.bai" ] && [ ! -f "$output_dir/dedup/dedup.bai" ]; then
            cp "$output_dir/dedup/dedup.bam.bai" "$output_dir/dedup/dedup.bai"
        fi
        echo "✅ Duplicates marked"
        
        # Stage 3: Deduplication Cognitive Supervisor Gate
        echo "  [SUPERVISOR] Evaluating Stage 3 Deduplication metrics against clinical policies (qc.json)..."
        set +e
        node "$SCRIPT_DIR/scripts/stage3_dedup_gate.js" "$output_dir" "$sample_name"
        DEDUP_EXIT=$?
        set -e
        
        if [ $DEDUP_EXIT -ne 0 ]; then
            if [ $DEDUP_EXIT -eq 1 ]; then
                echo ""
                echo "🛑 [PIPELINE HALTED] Deduplication Quality Gate failed rejection floor (Duplication >= 25% or Library < 10M)."
                echo "   Human-in-the-Loop Operator Opinion Gate is required."
                echo "   Use the Web Dashboard to either:"
                echo "     1. [Abort & Re-prepare Library] (Recommended clinical action)"
                echo "     2. [Override & Force Run] (High-Risk Research Mode)"
                echo "   Or run: touch \"$output_dir/.override_dedup_gate\" and restart pipeline."
                exit 1
            elif [ $DEDUP_EXIT -eq 2 ]; then
                echo "❌ [TOOL CRASH] Deduplication output missing or corrupted. Pipeline halted."
                exit 2
            else
                echo "❌ [ERROR] Unknown deduplication gate error ($DEDUP_EXIT)."
                exit 1
            fi
        fi
        echo "✅ Stage 3 Deduplication Approved by Supervisor. Proceeding to BQSR."
    fi
    
    # Layer 1 Storage Custodian: Safe Space Reclamation
    # Once Stage 3 gate approves (standard or override) and dedup.bam is verified (>10 KB, indexed), safely retire intermediate sorted.bam
    DEDUP_SIZE=$(wc -c < "$output_dir/dedup/dedup.bam" 2>/dev/null || echo 0)
    if [ "$DEDUP_SIZE" -gt 10240 ] && [ -f "$output_dir/dedup/dedup.bam.bai" -o -f "$output_dir/dedup/dedup.bai" ]; then
        if [ -f "$output_dir/sorted/sorted.bam" ]; then
            RECLAIM_KB=$(du -sk "$output_dir/sorted/sorted.bam" 2>/dev/null | awk '{print $1}')
            echo "🧹 [L1 STORAGE CUSTODIAN] Reclaiming disk space: retiring intermediate sorted.bam (${RECLAIM_KB:-0} KB)..."
            rm -f $output_dir/sorted/sorted.bam $output_dir/sorted/sorted.bam.bai $output_dir/sorted/sorted.bai
            echo "   Active validated BAM for Genome Viewer & BQSR: dedup/dedup.bam"
        fi
    fi
    
    # ── Checkpoint: Stage 4 (Base Quality Score Recalibration) ──
    stage4_approved=0
    if is_stage4_approved; then
        stage4_approved=1
    elif [ -f "$output_dir/.override_bqsr_gate" ]; then
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing BQSR output."
        echo "  Executing gate evaluation to apply override and resuming straight to Variant Calling..."
        set +e
        node "$SCRIPT_DIR/scripts/stage4_bqsr_gate.js" "$output_dir" "$sample_name" "$REFERENCE"
        BQSR_EXIT=$?
        set -e
        if [ $BQSR_EXIT -eq 0 ]; then
            echo "✅ Stage 4 BQSR Authorized via Operator Override. Proceeding directly to Variant Calling."
            stage4_approved=1
        else
            echo "❌ [SUPERVISOR] BQSR gate evaluation failed with code $BQSR_EXIT"
            exit $BQSR_EXIT
        fi
    fi

    if [ $stage4_approved -eq 1 ]; then
        echo "⏩ [CHECKPOINT] Stage 4 (Base Quality Score Recalibration) already completed and approved. Skipping to Variant Calling..."
        if [ -s "$output_dir/bqsr/recal.bam" ]; then
            BQSR_INPUT=$output_dir/bqsr/recal.bam
        else
            BQSR_INPUT=$output_dir/dedup/dedup.bam
        fi
    else
        step 7 "Base Quality Score Recalibration"
        mkdir -p "$output_dir/bqsr"
        rm -f "$output_dir/bqsr/post_recal_data.table"
        echo "$REF_BUILD" > "$output_dir/bqsr/ref_build.txt"
        KNOWN_SITES_FILE="$output_dir/bqsr/known_sites.txt"
        > "$KNOWN_SITES_FILE"
        BQSR_INPUT=$output_dir/dedup/dedup.bam

        if [ -f "$KNOWN_DBSNP" ]; then
            rm -f "$output_dir/.bqsr_bypassed"
            echo "$KNOWN_DBSNP" >> "$KNOWN_SITES_FILE"
            KS_ARGS="--known-sites $KNOWN_DBSNP"
            if [ -f "$KNOWN_MILLS" ]; then
                echo "$KNOWN_MILLS" >> "$KNOWN_SITES_FILE"
                KS_ARGS="$KS_ARGS --known-sites $KNOWN_MILLS"
            fi
            if [ -f "$KNOWN_INDELS" ]; then
                echo "$KNOWN_INDELS" >> "$KNOWN_SITES_FILE"
                KS_ARGS="$KS_ARGS --known-sites $KNOWN_INDELS"
            fi

            echo "  [EXECUTION AGENT] Running GATK BaseRecalibrator with reference-matched known-sites..."
            $GATK BaseRecalibrator \
                -I $output_dir/dedup/dedup.bam \
                -R $REFERENCE \
                $KS_ARGS \
                -O $output_dir/bqsr/recal_data.table

            echo "  [EXECUTION AGENT] Running GATK ApplyBQSR to write recalibrated BAM..."
            $GATK ApplyBQSR \
                -I $output_dir/dedup/dedup.bam \
                -R $REFERENCE \
                --bqsr-recal-file $output_dir/bqsr/recal_data.table \
                -O $output_dir/bqsr/recal.bam

            # Ensure index exists as both recal.bai and recal.bam.bai
            if [ ! -f "$output_dir/bqsr/recal.bam.bai" ] && [ ! -f "$output_dir/bqsr/recal.bai" ]; then
                echo "  [INDEXING] Indexing recalibrated BAM..."
                samtools index "$output_dir/bqsr/recal.bam"
            fi
            if [ -f "$output_dir/bqsr/recal.bai" ] && [ ! -f "$output_dir/bqsr/recal.bam.bai" ]; then
                cp "$output_dir/bqsr/recal.bai" "$output_dir/bqsr/recal.bam.bai"
            elif [ -f "$output_dir/bqsr/recal.bam.bai" ] && [ ! -f "$output_dir/bqsr/recal.bai" ]; then
                cp "$output_dir/bqsr/recal.bam.bai" "$output_dir/bqsr/recal.bai"
            fi

            BQSR_INPUT=$output_dir/bqsr/recal.bam
            echo "✅ Base quality recalibration applied"

            # Conditional two-pass (qc.json post_recal_residual_error): when pre-recal drift exceeds the clinical
            # threshold, re-run BaseRecalibrator on recal.bam to measure the residual miscalibration BQSR left behind.
            NEEDS_SECOND_PASS=$(node -e '
                try {
                    const fs = require("fs");
                    const e = require(process.argv[1]);
                    const q = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
                    const m = e.parseRecalTable(fs.readFileSync(process.argv[3], "utf8"));
                    const limit = q?.tiers?.tier_bqsr_recalibration_qc?.mean_quality_drift?.clinical_grade ?? 4.0;
                    console.log(m.meanQualityDrift > limit ? 1 : 0);
                } catch (err) { console.log(1); }
            ' "$SCRIPT_DIR/../../../lib/exome/bqsr-triage-engine.js" "$SCRIPT_DIR/../../../qc.json" "$output_dir/bqsr/recal_data.table" 2>/dev/null || echo 1)
            if [ "$NEEDS_SECOND_PASS" = "1" ]; then
                echo "  [EXECUTION AGENT] Pre-recal drift above clinical threshold — running second BaseRecalibrator pass on recal.bam to measure residual error..."
                $GATK BaseRecalibrator \
                    -I $output_dir/bqsr/recal.bam \
                    -R $REFERENCE \
                    $KS_ARGS \
                    -O $output_dir/bqsr/post_recal_data.table
            fi
        else
            echo "⚠️  Known-sites (dbSNP/Mills) not found in reference/known-sites/ — BYPASSING BQSR; using dedup.bam."
            echo "    Install known-sites VCFs to enable full Base Quality Score Recalibration."
            touch "$output_dir/.bqsr_bypassed"
            BQSR_INPUT=$output_dir/dedup/dedup.bam
        fi

        # Stage 4: BQSR Cognitive Supervisor Gate
        echo "  [SUPERVISOR] Evaluating Stage 4 BQSR metrics against clinical policies (qc.json)..."
        set +e
        node "$SCRIPT_DIR/scripts/stage4_bqsr_gate.js" "$output_dir" "$sample_name" "$REFERENCE"
        BQSR_EXIT=$?
        set -e

        if [ $BQSR_EXIT -ne 0 ]; then
            if [ $BQSR_EXIT -eq 1 ]; then
                echo ""
                echo "🛑 [PIPELINE HALTED] BQSR Quality Gate failed rejection floor (Empirical Q < 15 or Observations < 50M)."
                echo "   Human-in-the-Loop Operator Opinion Gate is required."
                echo "   Use the Web Dashboard to either:"
                echo "     1. [Abort Pipeline] (Recommended clinical action)"
                echo "     2. [Override & Force Run] (High-Risk Research Mode)"
                echo "   Or run: touch \"$output_dir/.override_bqsr_gate\" and restart pipeline."
                exit 1
            elif [ $BQSR_EXIT -eq 2 ]; then
                echo "❌ [TOOL CRASH] BQSR output missing or corrupted. Pipeline halted."
                exit 2
            else
                echo "❌ [ERROR] Unknown BQSR gate error ($BQSR_EXIT)."
                exit 1
            fi
        fi

        # Assert active BAM agreement
        EXPECTED_BQSR_INPUT="$output_dir/bqsr/recal.bam"
        [ -f "$output_dir/.bqsr_bypassed" ] && EXPECTED_BQSR_INPUT="$output_dir/dedup/dedup.bam"
        if [ "$BQSR_INPUT" != "$EXPECTED_BQSR_INPUT" ]; then
            echo "  ℹ️  [ALIGNMENT AGREEMENT] Synchronizing calling BAM: $EXPECTED_BQSR_INPUT"
            BQSR_INPUT="$EXPECTED_BQSR_INPUT"
        fi

        echo "✅ Stage 4 BQSR Approved by Supervisor. Proceeding to Variant Calling."
    fi

    # Layer 1 Storage Custodian: Safe Space Reclamation
    # Once Stage 4 gate approves (standard or override) and recal.bam is verified (>10 KB, indexed), safely retire intermediate dedup.bam
    RECAL_SIZE=$(wc -c < "$output_dir/bqsr/recal.bam" 2>/dev/null || echo 0)
    if [ "$RECAL_SIZE" -gt 10240 ] && [ -f "$output_dir/bqsr/recal.bam.bai" -o -f "$output_dir/bqsr/recal.bai" ]; then
        if [ -f "$output_dir/dedup/dedup.bam" ]; then
            RECLAIM_DEDUP_KB=$(du -sk "$output_dir/dedup/dedup.bam" 2>/dev/null | awk '{print $1}')
            echo "🧹 [L1 STORAGE CUSTODIAN] Reclaiming disk space: retiring intermediate dedup.bam (${RECLAIM_DEDUP_KB:-0} KB)..."
            rm -f $output_dir/dedup/dedup.bam $output_dir/dedup/dedup.bam.bai $output_dir/dedup/dedup.bai
            echo "   Active validated BAM for Genome Viewer & Calling: bqsr/recal.bam"
        fi
    fi
    
    # ── Checkpoint: Stage 5 (Variant Calling & Quality Filtering) ──
    stage5_approved=0
    if is_stage5_approved; then
        stage5_approved=1
    elif [ -f "$output_dir/.override_variant_gate" ] && [ -s "$output_dir/filtered/filtered_variants.vcf" -o -s "$output_dir/variants/raw_variants.vcf" ]; then
        echo "  ⚠️  [SUPERVISOR] Operator Override active with existing VCF. Applying override via Stage 5 gate..."
        set +e
        node "$SCRIPT_DIR/scripts/stage5_variant_gate.js" "$output_dir" "$sample_name" "$REFERENCE"
        S5_OVERRIDE_EXIT=$?
        set -e
        if [ $S5_OVERRIDE_EXIT -ne 0 ]; then
            echo "❌ [SUPERVISOR] Stage 5 gate failed to apply override (code $S5_OVERRIDE_EXIT)"
            exit $S5_OVERRIDE_EXIT
        fi
        stage5_approved=1
    fi

    if [ $stage5_approved -eq 1 ]; then
        echo "⏩ [CHECKPOINT] Stage 5 (Variant Calling & Quality Filtering) already completed and approved. Skipping to Variant Annotation..."
        PASS_VCF=$output_dir/filtered/filtered_PASS_only.vcf
        if [ ! -s "$PASS_VCF" ]; then
            PASS_VCF=$output_dir/filtered/filtered_variants.vcf
        fi
    else
        # 7. Variant Calling
        step 8 "Variant Calling"
        mkdir -p $output_dir/variants $output_dir/filtered

        HC_EXTRA_ARGS=""
        if [ -f "$output_dir/stage4_bqsr_reasoning.json" -o -f "$output_dir/qc/stage4_bqsr_reasoning.json" ]; then
            S4_FILE="$output_dir/stage4_bqsr_reasoning.json"
            [ ! -f "$S4_FILE" ] && S4_FILE="$output_dir/qc/stage4_bqsr_reasoning.json"
            HC_EXTRA_ARGS=$(node -e '
                try {
                    const d = JSON.parse(require("fs").readFileSync(process.argv[1], "utf8"));
                    const args = [];
                    if (d.downstreamDirectives?.minPruning) args.push("--min-pruning " + d.downstreamDirectives.minPruning);
                    if (d.downstreamDirectives?.pcrIndelModel) args.push("--pcr-indel-model " + d.downstreamDirectives.pcrIndelModel);
                    console.log(args.join(" "));
                } catch(e) {}
            ' "$S4_FILE" 2>/dev/null)
        fi

        echo "  [EXECUTION AGENT] Running GATK HaplotypeCaller on recalibrated BAM..."
        $GATK HaplotypeCaller \
            -R $REFERENCE \
            -I $BQSR_INPUT \
            -O $output_dir/variants/raw_variants.vcf \
            --native-pair-hmm-threads $threads \
            --stand-call-conf 30.0 \
            $HC_EXTRA_ARGS

        RAW_COUNT=$(grep -v "^#" $output_dir/variants/raw_variants.vcf | wc -l | tr -d ' ')
        echo "✅ Called $RAW_COUNT raw variants"

        # 8. Filtering
        step 9 "Variant Filtering"
        echo "  [EXECUTION AGENT] Running GATK VariantFiltration (Broad Best Practices SNV & Indel Hard Filters)..."
        $GATK VariantFiltration \
            -R $REFERENCE \
            -V $output_dir/variants/raw_variants.vcf \
            -O $output_dir/filtered/filtered_variants.vcf \
            --filter-expression "vc.isSNP() && QD < 2.0" --filter-name "LowQD_SNP" \
            --filter-expression "vc.isSNP() && QUAL < 30.0" --filter-name "LowQUAL_SNP" \
            --filter-expression "vc.isSNP() && FS > 60.0" --filter-name "StrandBiasFS_SNP" \
            --filter-expression "vc.isSNP() && SOR > 3.0" --filter-name "HighSOR_SNP" \
            --filter-expression "vc.isSNP() && MQ < 40.0" --filter-name "LowMQ_SNP" \
            --filter-expression "vc.isSNP() && MQRankSum < -12.5" --filter-name "LowMQRankSum_SNP" \
            --filter-expression "vc.isSNP() && ReadPosRankSum < -8.0" --filter-name "ReadPosBias_SNP" \
            --filter-expression "!vc.isSNP() && QD < 2.0" --filter-name "LowQD_INDEL" \
            --filter-expression "!vc.isSNP() && QUAL < 30.0" --filter-name "LowQUAL_INDEL" \
            --filter-expression "!vc.isSNP() && FS > 200.0" --filter-name "StrandBiasFS_INDEL" \
            --filter-expression "!vc.isSNP() && SOR > 10.0" --filter-name "HighSOR_INDEL" \
            --filter-expression "!vc.isSNP() && ReadPosRankSum < -20.0" --filter-name "ReadPosBias_INDEL"

        # Create PASS-only VCF
        PASS_VCF=$output_dir/filtered/filtered_PASS_only.vcf
        grep "^#" $output_dir/filtered/filtered_variants.vcf > $PASS_VCF
        grep -v "^#" $output_dir/filtered/filtered_variants.vcf | grep -w "PASS" >> $PASS_VCF

        PASS_COUNT=$(grep -v "^#" $PASS_VCF | wc -l | tr -d ' ')
        echo "✅ Filtered variants: $PASS_COUNT / $RAW_COUNT passed filters"

        # Stage 5 Cognitive Supervisor Gate
        echo "  [SUPERVISOR] Evaluating Stage 5 Variant Calling & Biology metrics against clinical policies (qc.json)..."
        set +e
        node "$SCRIPT_DIR/scripts/stage5_variant_gate.js" "$output_dir" "$sample_name" "$REFERENCE"
        S5_EXIT=$?
        set -e

        if [ $S5_EXIT -ne 0 ]; then
            if [ $S5_EXIT -eq 1 ]; then
                echo ""
                echo "🛑 [PIPELINE HALTED] Stage 5 Variant Callset Gate failed biological rejection floor."
                echo "   Human-in-the-Loop Operator Opinion Gate is required."
                echo "   Use the Web Dashboard to either:"
                echo "     1. [Abort Pipeline] (Recommended clinical action)"
                echo "     2. [Override & Force Run] (Research Mode with mandatory flagging)"
                echo "   Or run: touch \"$output_dir/.override_variant_gate\" and restart pipeline."
                exit 1
            elif [ $S5_EXIT -eq 2 ]; then
                echo "❌ [TOOL CRASH] Variant calling output missing or corrupted. Pipeline halted."
                exit 2
            else
                echo "❌ [ERROR] Unknown Stage 5 gate error ($S5_EXIT)."
                exit 1
            fi
        fi

        echo "✅ Stage 5 Variant Calling & Callset Biology Approved by Supervisor. Proceeding to Annotation."
    fi
    
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
            RAW=$([ -f "$RESULT_DIR/variants/raw_variants.vcf" ] && (grep -v "^#" "$RESULT_DIR/variants/raw_variants.vcf" 2>/dev/null || true) | wc -l | tr -d ' ' || echo 0)
            PASS=$((grep -v "^#" "$RESULT_DIR/filtered/filtered_variants.vcf" 2>/dev/null | grep -w "PASS" 2>/dev/null || true) | wc -l | tr -d ' ')
            
            cat >> $SUMMARY_FILE << EOF
  Raw variants called:     $RAW
  PASS variants:           $PASS
EOF
        fi
        
        # Annotation
        ANNOT_FILE=$(find "$RESULT_DIR/annovar" -name "*.hg19_multianno.txt" 2>/dev/null | head -1 || true)
        if [ -f "$ANNOT_FILE" ]; then
            PATHOGENIC=$((grep -i "pathogenic" "$ANNOT_FILE" 2>/dev/null || true) | wc -l | tr -d ' ')
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
    2>/dev/null || true

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

#!/usr/bin/env bash
# Stream, compact, and index gnomAD v2.1.1 exomes per-chromosome
# 1. Downloads on-the-fly and strips unnecessary INFO tags
#    Retains ONLY: AF, AC, AN, nhomalt, AF_popmax, faf95, faf95_afr, faf95_amr, faf95_eas, faf95_nfe, faf95_sas
# 2. Renames contigs (1 -> chr1) to align with UCSC hg19
# 3. Fail-safe resume logic with .done sentinel markers
# 4. Concatenates into karyotypic master compact VCF: gnomad.exomes.r2.1.1.compact.vcf.gz (~3.2 GB)

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NGS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes"
OUT_DIR="$NGS_DIR/reference/known-sites/hg19/gnomad_chroms"
MASTER_GNOMAD="$NGS_DIR/reference/known-sites/hg19/gnomad.exomes.r2.1.1.compact.vcf.gz"
MASTER_DONE="$NGS_DIR/reference/known-sites/hg19/.gnomad_compact.done"

mkdir -p "$OUT_DIR"

# Target chromosomes (can be overridden via CLI args, e.g. ./stream_gnomad_chroms.sh 8)
if [ $# -gt 0 ] && [ "$1" != "--all" ]; then
    TARGET_CHRS=("$@")
else
    TARGET_CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
fi

echo "============================================================"
echo "  [L1 PROVISIONING] Resumable Streaming & Compaction for gnomAD v2.1.1"
echo "  Target Chromosomes: ${TARGET_CHRS[*]}"
echo "============================================================"

# Fields to retain (ClinGen SVI 2020 calibrated)
KEEP_FIELDS="^INFO/AF,INFO/AC,INFO/AN,INFO/nhomalt,INFO/AF_popmax,INFO/faf95,INFO/faf95_afr,INFO/faf95_amr,INFO/faf95_eas,INFO/faf95_nfe,INFO/faf95_sas"

for CHR in "${TARGET_CHRS[@]}"; do
    # Remove leading 'chr' if user passed chr8
    C_NUM="${CHR#chr}"
    DONE_FLAG="$OUT_DIR/chr${C_NUM}.done"
    FINAL_FILE="$OUT_DIR/chr${C_NUM}.compact.vcf.gz"
    TMP_FILE="$OUT_DIR/chr${C_NUM}.tmp.vcf.gz"

    if [ -f "$DONE_FLAG" ] && [ -s "$FINAL_FILE" ] && [ -s "$FINAL_FILE.tbi" ]; then
        echo "  [L1 ADVISOR] Chromosome chr${C_NUM} already verified. Skipping."
        continue
    fi

    echo "  [L1 ADVISOR] Streaming and compacting gnomAD chr${C_NUM}..."
    rm -f "$TMP_FILE" "$TMP_FILE.tbi" "$DONE_FLAG"

    RENAME_MAP=$(mktemp -t rename_map_XXXXXX)
    echo -e "${C_NUM}\tchr${C_NUM}" > "$RENAME_MAP"

    # Stream on-the-fly: download -> strip tags -> rename contig -> compress
    set +e
    bcftools annotate \
        -x "$KEEP_FIELDS" \
        --rename-chrs "$RENAME_MAP" \
        "$GNOMAD_BASE/gnomad.exomes.r2.1.1.sites.${C_NUM}.vcf.bgz" \
        -Oz -o "$TMP_FILE" 2>/dev/null
    STREAM_EXIT=$?
    rm -f "$RENAME_MAP"
    set -e

    if [ $STREAM_EXIT -ne 0 ] || [ ! -s "$TMP_FILE" ]; then
        echo "❌ [ERROR] Streaming failed for chr${C_NUM} (exit $STREAM_EXIT)."
        rm -f "$TMP_FILE" "$TMP_FILE.tbi"
        exit 1
    fi

    if ! tabix -p vcf "$TMP_FILE"; then
        echo "❌ [ERROR] Tabix indexing failed for chr${C_NUM}."
        rm -f "$TMP_FILE" "$TMP_FILE.tbi"
        exit 1
    fi

    mv "$TMP_FILE" "$FINAL_FILE"
    mv "$TMP_FILE.tbi" "$FINAL_FILE.tbi"
    touch "$DONE_FLAG"
    SIZE_MB=$(du -m "$FINAL_FILE" | cut -f1)
    echo "  ✅ Chromosome chr${C_NUM} compacted (${SIZE_MB} MB) and indexed."
done

# If all 24 chromosomes were requested, concatenate in karyotypic order
if [ ${#TARGET_CHRS[@]} -eq 24 ]; then
    echo "============================================================"
    echo "  [L1 ADVISOR] Concatenating all 24 chromosomes in karyotypic order..."
    echo "============================================================"
    
    KARYOTYPIC_FILES=()
    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
        KARYOTYPIC_FILES+=("$OUT_DIR/chr${CHR}.compact.vcf.gz")
    done

    bcftools concat -a "${KARYOTYPIC_FILES[@]}" -Oz -o "$MASTER_GNOMAD"
    tabix -p vcf "$MASTER_GNOMAD"
    touch "$MASTER_DONE"
    
    TOTAL_SIZE=$(du -h "$MASTER_GNOMAD" | cut -f1)
    echo "✅ [SUCCESS] Master compact gnomAD v2.1.1 ready: $MASTER_GNOMAD ($TOTAL_SIZE)"
fi

#!/usr/bin/env bash
# Prepare dbSNP GRCh37 VCF:
# 1. Update CAF and TOPMED header definitions to Number=R
# 2. Split multi-allelic records with bcftools norm -m -any
#    (Automatically slices CAF and TOPMED to [ref_freq, alt_freq] per record)
# 3. Build tabix index (.tbi)
# 4. Atomically swap into dbsnp.vcf.gz
#
# Traceability:
# Prepared: 2026-09-24 18:04:00+03:00
# Output SHA256: 37068eed696f7a9030b11ff5468405749baf9ae9e8ea286d84d7f37b46e62599 (1.5 GB)
# Metadata: reference/known-sites/hg19/dbsnp_metadata.json

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NGS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
KNOWN_SITES_DIR="$NGS_DIR/reference/known-sites/hg19"
INPUT_DBSNP="$KNOWN_SITES_DIR/dbsnp.vcf.gz"
DONE_MARKER="$KNOWN_SITES_DIR/.dbsnp_split.done"

if [ -f "$DONE_MARKER" ]; then
    echo "✅ [DBSNP] Local dbSNP GRCh37 VCF already split and indexed with per-allele CAF."
    exit 0
fi

if [ ! -f "$INPUT_DBSNP" ]; then
    echo "❌ [ERROR] $INPUT_DBSNP not found!"
    exit 1
fi

echo "============================================================"
echo "  [L1 PROVISIONING] Splitting dbSNP Multi-Allelics & Slicing CAF"
echo "============================================================"

TMP_OUT="$KNOWN_SITES_DIR/dbsnp.split.tmp.vcf.gz"
trap 'rm -f "$TMP_OUT" "$TMP_OUT.tbi"' EXIT

echo "  Streaming, re-heading Number=R, and splitting multi-allelics..."
bcftools view "$INPUT_DBSNP" \
  | sed -e 's/ID=CAF,Number=\./ID=CAF,Number=R/' -e 's/ID=TOPMED,Number=\./ID=TOPMED,Number=R/' \
  | bcftools norm -m -any -Oz -o "$TMP_OUT"

echo "  Building tabix index..."
tabix -f -p vcf "$TMP_OUT"

echo "  Atomically replacing dbsnp.vcf.gz..."
mv "$TMP_OUT" "$INPUT_DBSNP"
mv "$TMP_OUT.tbi" "${INPUT_DBSNP}.tbi"
touch "$DONE_MARKER"

echo "✅ [SUCCESS] dbSNP successfully prepared with per-allele CAF at $INPUT_DBSNP"

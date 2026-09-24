#!/usr/bin/env bash
# Prepare NCBI ClinVar GRCh37 VCF:
# 1. Download clinvar.vcf.gz (~188 MB) + .tbi from NCBI
# 2. Exclude MT contig (UCSC chrM is NC_001807 vs ClinVar MT NC_012920)
# 3. Rename contigs Ensembl -> UCSC (1 -> chr1, X -> chrX, Y -> chrY)
# 4. Normalize with bcftools norm against hg19.fa
# 5. Build tabix index (.tbi)

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NGS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REF_FASTA="$NGS_DIR/reference/hg19.fa"
KNOWN_SITES_DIR="$NGS_DIR/reference/known-sites/hg19"
TARGET_CLINVAR="$KNOWN_SITES_DIR/clinvar.vcf.gz"
TARGET_TBI="$KNOWN_SITES_DIR/clinvar.vcf.gz.tbi"
DONE_MARKER="$KNOWN_SITES_DIR/.clinvar.done"

mkdir -p "$KNOWN_SITES_DIR"

if [ -f "$DONE_MARKER" ] && [ -s "$TARGET_CLINVAR" ] && [ -s "$TARGET_TBI" ]; then
    echo "✅ [CLINVAR] Local ClinVar GRCh37 VCF already prepared and indexed: $TARGET_CLINVAR"
    exit 0
fi

echo "============================================================"
echo "  [L1 PROVISIONING] Preparing NCBI ClinVar GRCh37 VCF"
echo "============================================================"

TMP_DIR="$(mktemp -d -t clinvar_prep_XXXXXX)"
trap 'rm -rf "$TMP_DIR"' EXIT

RAW_CLINVAR="$TMP_DIR/clinvar_raw.vcf.gz"
RAW_TBI="$TMP_DIR/clinvar_raw.vcf.gz.tbi"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"

echo "  Downloading NCBI ClinVar GRCh37 VCF & tabix index (~188 MB)..."
curl -fSL "$CLINVAR_URL" -o "$RAW_CLINVAR"
curl -fSL "$CLINVAR_URL.tbi" -o "$RAW_TBI"

echo "  Verifying MD5 checksum against NCBI..."
EXPECTED_MD5=$(curl -fsSL "${CLINVAR_URL}.md5" | awk '{print $1}')
if command -v md5 &>/dev/null; then
    ACTUAL_MD5=$(md5 -q "$RAW_CLINVAR")
else
    ACTUAL_MD5=$(md5sum "$RAW_CLINVAR" | awk '{print $1}')
fi

if [ -n "$EXPECTED_MD5" ] && [ "$EXPECTED_MD5" != "$ACTUAL_MD5" ]; then
    echo "❌ [ERROR] ClinVar MD5 checksum mismatch! Expected: $EXPECTED_MD5, got: $ACTUAL_MD5"
    exit 1
fi
echo "  ✅ ClinVar MD5 verified ($ACTUAL_MD5)"

echo "  Building chromosome renaming map (1 -> chr1, excluding MT)..."
RENAME_MAP="$TMP_DIR/rename_chrs.tsv"
AUTOSOMAL_SEX_CHRS=()
for i in {1..22} X Y; do
    echo -e "${i}\tchr${i}" >> "$RENAME_MAP"
    AUTOSOMAL_SEX_CHRS+=("$i")
done

echo "  Extracting autosomes and sex chromosomes (skipping MT) and renaming contigs..."
FILTERED_CLINVAR="$TMP_DIR/clinvar_renamed.vcf.gz"

# Query only 1..22 X Y, pipe through annotate --rename-chrs
REGIONS_LIST=$(IFS=,; echo "${AUTOSOMAL_SEX_CHRS[*]}")
bcftools view -r "$REGIONS_LIST" "$RAW_CLINVAR" -Ou \
  | bcftools annotate --rename-chrs "$RENAME_MAP" -Oz -o "$FILTERED_CLINVAR"

tabix -p vcf "$FILTERED_CLINVAR"

echo "  Normalizing ClinVar indels and splitting multi-allelics against hg19.fa..."
TMP_FINAL="$TMP_DIR/clinvar_normalized.vcf.gz"
bcftools norm -m -any -c w -f "$REF_FASTA" "$FILTERED_CLINVAR" -Oz -o "$TMP_FINAL"
tabix -p vcf "$TMP_FINAL"

echo "  Installing to $TARGET_CLINVAR..."
mv "$TMP_FINAL" "$TARGET_CLINVAR"
mv "$TMP_FINAL.tbi" "$TARGET_TBI"
touch "$DONE_MARKER"

TOTAL_RECORDS=$(bcftools view -H "$TARGET_CLINVAR" | wc -l | tr -d ' ')
echo "✅ [CLINVAR READY] Installed $TOTAL_RECORDS normalized ClinVar records into $TARGET_CLINVAR"

#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# EIL installer (Layer 1) — provisions the components a WES run requires.
# Args: one or more action ids:
#   index-hg38 | index-hg19 | known-sites-hg38 | known-sites-hg19 |
#   tool:<brew-formula> | playwright
# Writes plaintext progress to reference/install_status.txt (polled by the UI).
# Best-effort + non-destructive: never overwrites an already-built index.
# ─────────────────────────────────────────────────────────────────────────────
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="$SCRIPT_DIR/reference"
KS="$REF/known-sites"
GATK="${GATK:-$SCRIPT_DIR/gatk-4.6.2.0/gatk}"
STATUS="$REF/install_status.txt"
mkdir -p "$REF" "$KS"

set_status() { # running action step done error
  { echo "running=$1"; echo "action=$2"; echo "step=$3"; echo "done=$4"; echo "error=$5"; echo "updated=$(date -u +%FT%TZ)"; } > "$STATUS"
}
fail() { set_status 0 "$1" "$2" 1 "$3"; exit 1; }

# GATK resource-bundle URLs — public HTTP bucket (assembly38 files match hg38.fa).
BROAD="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
DBSNP_HG38="$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS_HG38="$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

download() { # url dest label
  local url="$1" dest="$2" label="$3"
  [ -f "$dest" ] && { echo "  $label already present"; return 0; }
  set_status 1 "$CUR" "Downloading $label" 0 ""
  curl -fL --retry 3 -o "$dest.part" "$url" || return 1
  mv "$dest.part" "$dest"
  # GATK needs the tabix index next to the VCF — fetch it (best effort).
  curl -fsL --retry 2 -o "$dest.tbi" "$url.tbi" 2>/dev/null || true
}

for action in "$@"; do
  CUR="$action"
  case "$action" in
    index-hg38|index-hg19)
      build="${action#index-}"
      fa="$REF/${build}.fa"
      set_status 1 "$action" "Preparing ${build} FASTA" 0 ""
      if [ ! -f "$fa" ] && [ -f "$fa.gz" ]; then
        set_status 1 "$action" "Decompressing ${build}.fa.gz (~3 GB)" 0 ""
        gunzip -kc "$fa.gz" > "$fa" || fail "$action" "decompress" "gunzip failed"
      fi
      [ -f "$fa" ] || fail "$action" "missing FASTA" "No ${build}.fa or ${build}.fa.gz in reference/"
      if [ ! -f "$fa.bwt" ]; then
        set_status 1 "$action" "Building BWA index (can take 30–60 min)" 0 ""
        bwa index "$fa" || fail "$action" "bwa index" "bwa index failed"
      fi
      if [ ! -f "$fa.fai" ]; then
        set_status 1 "$action" "samtools faidx" 0 ""
        samtools faidx "$fa" || fail "$action" "faidx" "samtools faidx failed"
      fi
      if [ ! -f "$REF/${build}.dict" ]; then
        set_status 1 "$action" "GATK CreateSequenceDictionary" 0 ""
        "$GATK" CreateSequenceDictionary -R "$fa" -O "$REF/${build}.dict" || fail "$action" "dict" "CreateSequenceDictionary failed"
      fi
      echo "✅ ${build} indexed"
      ;;
    known-sites-hg38)
      mkdir -p "$KS/hg38"
      download "$DBSNP_HG38" "$KS/hg38/dbsnp.vcf.gz" "dbSNP (hg38)" || fail "$action" "dbSNP" "dbSNP download failed"
      download "$MILLS_HG38" "$KS/hg38/mills.vcf.gz" "Mills/1000G indels (hg38)" || fail "$action" "Mills" "Mills download failed"
      echo "✅ hg38 known-sites installed"
      ;;
    known-sites-hg19)
      # UCSC hg19 (chr-prefixed) needs matching-contig known-sites; the b37 bundle
      # will not line up. Flag rather than install the wrong contigs.
      fail "$action" "manual" "hg19 (UCSC) needs chr-prefixed dbSNP/Mills — place them in reference/known-sites/ (dbsnp.vcf.gz, mills.vcf.gz)."
      ;;
    tool:*)
      formula="${action#tool:}"
      if command -v brew >/dev/null 2>&1; then
        set_status 1 "$action" "brew install $formula" 0 ""
        brew install "$formula" || fail "$action" "brew" "brew install $formula failed"
      else
        fail "$action" "no brew" "Homebrew not found — install $formula manually."
      fi
      ;;
    playwright)
      set_status 1 "$action" "pip install playwright + chromium" 0 ""
      python3 -m pip install playwright >/dev/null 2>&1 && python3 -m playwright install chromium || fail "$action" "playwright" "playwright provisioning failed"
      ;;
    *)
      echo "unknown action: $action" ;;
  esac
done

set_status 0 "" "All requested components installed" 1 ""
echo "✅ EIL install complete"

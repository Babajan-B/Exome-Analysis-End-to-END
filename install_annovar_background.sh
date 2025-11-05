#!/bin/bash
# Background ANNOVAR Installation Script
# This runs in parallel with the main pipeline

set -e

ANNOVAR_URL="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz"
ANNOVAR_DIR=~/NGS/tools/annovar
LOCK_FILE=~/NGS/.annovar_installing
SUCCESS_FILE=~/NGS/.annovar_installed

# Create lock file
touch $LOCK_FILE

echo "════════════════════════════════════════════════════════════"
echo "ANNOVAR Background Installation Started"
echo "Time: $(date)"
echo "════════════════════════════════════════════════════════════"
echo ""

# Check if already installed
if [ -f "$SUCCESS_FILE" ] && [ -d "$ANNOVAR_DIR" ]; then
    echo "✅ ANNOVAR already installed!"
    rm -f $LOCK_FILE
    exit 0
fi

# Download and install ANNOVAR
echo "Downloading ANNOVAR..."
mkdir -p ~/NGS/tools
cd ~/NGS/tools

if [ ! -d "annovar" ]; then
    wget -q --show-progress -O annovar.latest.tar.gz "$ANNOVAR_URL" || {
        echo "❌ Download failed"
        rm -f $LOCK_FILE
        exit 1
    }
    
    tar -xzf annovar.latest.tar.gz
    rm annovar.latest.tar.gz
    chmod +x annovar/*.pl
fi

# Download databases
echo ""
echo "Downloading ANNOVAR databases (this takes 20-30 minutes)..."
cd $ANNOVAR_DIR

databases=(
    "refGene"
    "knownGene"
    "ensGene"
    "clinvar_20240917"
    "avsnp150"
    "dbnsfp42a"
    "gnomad211_exome"
    "cosmic70"
    "icgc28"
)

for db in "${databases[@]}"; do
    echo "  - Downloading $db..."
    perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar $db humandb/ 2>&1 | grep -E "(NOTICE|done)" || true
done

# Mark as complete
touch $SUCCESS_FILE
rm -f $LOCK_FILE

echo ""
echo "════════════════════════════════════════════════════════════"
echo "✅ ANNOVAR Installation Complete!"
echo "Time: $(date)"
echo "════════════════════════════════════════════════════════════"


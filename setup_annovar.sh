#!/bin/bash
# ANNOVAR Installation Script for Jarvis Lab / Cloud Instance
# This script downloads and installs ANNOVAR with required databases

set -e  # Exit on any error

echo "=========================================="
echo "ANNOVAR Setup for NGS Pipeline"
echo "=========================================="
echo ""

# Configuration
ANNOVAR_URL="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz"
INSTALL_DIR=~/NGS/tools
ANNOVAR_DIR=$INSTALL_DIR/annovar

# Create tools directory if it doesn't exist
echo "[1/5] Creating installation directory..."
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR

# Download ANNOVAR
echo "[2/5] Downloading ANNOVAR..."
if [ -d "$ANNOVAR_DIR" ]; then
    echo "⚠️  ANNOVAR directory already exists. Backing up to annovar.bak..."
    mv $ANNOVAR_DIR ${ANNOVAR_DIR}.bak.$(date +%Y%m%d_%H%M%S)
fi

wget -O annovar.latest.tar.gz "$ANNOVAR_URL"
echo "✅ Download complete"

# Extract ANNOVAR
echo "[3/5] Extracting ANNOVAR..."
tar -xzf annovar.latest.tar.gz
rm annovar.latest.tar.gz
echo "✅ Extraction complete"

# Make scripts executable
echo "[4/5] Setting permissions..."
chmod +x $ANNOVAR_DIR/*.pl
echo "✅ Permissions set"

# Download databases (common ones for hg19)
echo "[5/5] Downloading ANNOVAR databases for hg19..."
echo "This will download several databases (~2-5 GB total)"
echo ""

cd $ANNOVAR_DIR

# Essential databases for hg19
databases=(
    "refGene"           # RefSeq genes
    "knownGene"         # UCSC known genes
    "ensGene"           # Ensembl genes
    "clinvar_20240917"  # ClinVar (adjust date as needed)
    "avsnp150"          # dbSNP version 150
    "dbnsfp42a"         # dbNSFP for functional predictions
    "gnomad211_exome"   # gnomAD exome frequencies
    "cosmic70"          # COSMIC database
    "icgc28"            # ICGC somatic mutations
)

echo "Databases to be downloaded:"
for db in "${databases[@]}"; do
    echo "  - $db"
done
echo ""

# Download each database
for db in "${databases[@]}"; do
    echo "----------------------------------------"
    echo "Downloading: $db"
    echo "----------------------------------------"
    
    # Try to download, but don't fail if one database is unavailable
    if perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar $db humandb/; then
        echo "✅ $db downloaded successfully"
    else
        echo "⚠️  Failed to download $db (might be unavailable or require different command)"
    fi
    echo ""
done

echo ""
echo "=========================================="
echo "✅ ANNOVAR Setup Complete!"
echo "=========================================="
echo ""
echo "Installation Summary:"
echo "  ✅ ANNOVAR installed at: $ANNOVAR_DIR"
echo "  ✅ Databases installed at: $ANNOVAR_DIR/humandb/"
echo ""
echo "Usage Examples:"
echo ""
echo "1. Annotate a VCF file:"
echo "   perl $ANNOVAR_DIR/table_annovar.pl input.vcf \\"
echo "       $ANNOVAR_DIR/humandb/ \\"
echo "       -buildver hg19 \\"
echo "       -out output_prefix \\"
echo "       -remove \\"
echo "       -protocol refGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a \\"
echo "       -operation g,f,f,f,f \\"
echo "       -nastring . \\"
echo "       -vcfinput"
echo ""
echo "2. Convert VCF to ANNOVAR input format:"
echo "   perl $ANNOVAR_DIR/convert2annovar.pl -format vcf4 input.vcf > input.avinput"
echo ""
echo "3. Annotate with gene-based annotation:"
echo "   perl $ANNOVAR_DIR/annotate_variation.pl -geneanno \\"
echo "       -buildver hg19 \\"
echo "       input.avinput \\"
echo "       $ANNOVAR_DIR/humandb/"
echo ""
echo "Storage used: $(du -sh $ANNOVAR_DIR | cut -f1)"
echo "=========================================="
echo ""
echo "Next Steps:"
echo "  1. Test ANNOVAR: bash test_annovar.sh (if available)"
echo "  2. Integrate with pipeline: see integrate_annovar.md"
echo "  3. Add to PATH (optional):"
echo "     echo 'export PATH=\$PATH:$ANNOVAR_DIR' >> ~/.bashrc"
echo "     source ~/.bashrc"
echo ""


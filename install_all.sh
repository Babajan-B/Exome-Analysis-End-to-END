#!/bin/bash
# Complete NGS Pipeline Installation Script
# Single command to install everything: tools, reference genome, and ANNOVAR
# For Jarvis Lab / Cloud Instance

set -e  # Exit on any error

echo "╔════════════════════════════════════════════════════════════╗"
echo "║     NGS EXOME ANALYSIS PIPELINE - COMPLETE INSTALLATION    ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "This will install:"
echo "  ✓ System dependencies"
echo "  ✓ Bioinformatics tools (FastQC, BWA, SAMtools, GATK, etc.)"
echo "  ✓ Reference genome (hg19) with indices"
echo "  ✓ snpEff for annotation"
echo "  ✓ ANNOVAR with databases"
echo ""
echo "Storage required: ~15 GB"
echo "Time required: 45-60 minutes"
echo ""

# Ask for confirmation
read -p "Continue with installation? (y/n): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Installation cancelled."
    exit 0
fi

echo ""
echo "Starting installation..."
echo ""

# Update system
echo "[1/9] Updating system packages..."
sudo apt-get update -qq

# Install system dependencies
echo "[2/9] Installing system dependencies..."
sudo apt-get install -y -qq \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    openjdk-17-jdk \
    openjdk-17-jre \
    python3-pip \
    python3-venv \
    parallel \
    perl

# Set Java 17 as default
sudo update-alternatives --set java /usr/lib/jvm/java-17-openjdk-amd64/bin/java 2>/dev/null || true

# Install bioinformatics tools
echo "[3/9] Installing bioinformatics tools..."
sudo apt-get install -y -qq \
    fastqc \
    bwa \
    samtools \
    bcftools \
    tabix

# Install fastp
echo "[4/9] Installing fastp..."
if ! command -v fastp &> /dev/null; then
    wget -q http://opengene.org/fastp/fastp
    chmod a+x ./fastp
    sudo mv ./fastp /usr/local/bin/
fi

# Install GATK
echo "[5/9] Installing GATK 4.6.2.0..."
if [ ! -d "/opt/gatk-4.6.2.0" ]; then
    cd /tmp
    wget -q https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
    unzip -q gatk-4.6.2.0.zip
    sudo mv gatk-4.6.2.0 /opt/
    sudo ln -sf /opt/gatk-4.6.2.0/gatk /usr/local/bin/gatk
fi

# Download reference genome (hg19)
echo "[6/9] Setting up reference genome (hg19)..."
cd ~/NGS
mkdir -p reference
cd reference

if [ ! -f hg19.fa ]; then
    echo "   Downloading hg19 (~3GB, this will take a few minutes)..."
    wget -q --show-progress http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip hg19.fa.gz
fi

# Index reference genome
echo "[7/9] Indexing reference genome (15-30 minutes)..."
if [ ! -f hg19.fa.bwt ]; then
    echo "   Creating BWA index..."
    bwa index hg19.fa
fi

if [ ! -f hg19.fa.fai ]; then
    echo "   Creating samtools index..."
    samtools faidx hg19.fa
fi

if [ ! -f hg19.dict ]; then
    echo "   Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
fi

# Download snpEff
echo "[8/9] Setting up snpEff..."
cd ~/NGS
mkdir -p tools
cd tools

if [ ! -f snpEff/snpEff.jar ]; then
    wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -q snpEff_latest_core.zip
    cd snpEff
    echo "   Downloading snpEff database..."
    java -jar snpEff.jar download -v GRCh37.75
    cd ..
fi

# Install ANNOVAR
echo "[9/9] Installing ANNOVAR with databases (~5GB, 20-30 min)..."
cd ~/NGS

# Download ANNOVAR
ANNOVAR_URL="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz"
ANNOVAR_DIR=~/NGS/tools/annovar

if [ ! -d "$ANNOVAR_DIR" ]; then
    mkdir -p ~/NGS/tools
    cd ~/NGS/tools
    
    echo "   Downloading ANNOVAR..."
    wget -q -O annovar.latest.tar.gz "$ANNOVAR_URL"
    tar -xzf annovar.latest.tar.gz
    rm annovar.latest.tar.gz
    chmod +x annovar/*.pl
    
    echo "   Downloading ANNOVAR databases..."
    cd annovar
    
    # Essential databases
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
        echo "      - $db"
        perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar $db humandb/ 2>&1 | grep -E "(NOTICE|done)" || true
    done
else
    echo "   ANNOVAR already installed"
fi

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║              ✅ INSTALLATION COMPLETE!                     ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "Installation Summary:"
echo "  ✅ System tools installed"
echo "  ✅ Bioinformatics tools: FastQC, fastp, BWA, SAMtools, GATK"
echo "  ✅ Reference genome (hg19) downloaded and indexed"
echo "  ✅ snpEff configured with GRCh37.75 database"
echo "  ✅ ANNOVAR installed with 9 databases"
echo ""
echo "Verification:"
echo "  - FastQC:  $(which fastqc)"
echo "  - BWA:     $(which bwa)"
echo "  - SAMtools: $(which samtools)"
echo "  - GATK:    $(which gatk)"
echo "  - ANNOVAR: $ANNOVAR_DIR/table_annovar.pl"
echo ""
echo "Storage used: $(du -sh ~/NGS | cut -f1)"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Next Steps:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. Upload your FASTQ files to ~/NGS/data/"
echo "   mkdir -p ~/NGS/data"
echo "   # Upload R1.fastq.gz and R2.fastq.gz"
echo ""
echo "2. Run complete analysis (single command):"
echo "   bash run_complete_analysis.sh \\"
echo "       data/R1.fastq.gz \\"
echo "       data/R2.fastq.gz \\"
echo "       sample_name \\"
echo "       16"
echo ""
echo "3. Results will be in ~/NGS/results/sample_name/"
echo "   - Annotated VCF: results/sample_name/annovar/*.hg19_multianno.vcf"
echo "   - Excel table:   results/sample_name/annovar/*.hg19_multianno.txt"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""


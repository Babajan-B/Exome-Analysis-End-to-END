#!/bin/bash
# Jarvis Lab / Cloud Instance Setup Script
# Single command to install all dependencies and tools

set -e  # Exit on any error

echo "=========================================="
echo "NGS Exome Analysis - Cloud Setup"
echo "=========================================="
echo ""

# Update system
echo "[1/8] Updating system packages..."
sudo apt-get update -qq

# Install system dependencies
echo "[2/8] Installing system dependencies..."
sudo apt-get install -y -qq \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    default-jre \
    default-jdk \
    python3-pip \
    python3-venv \
    parallel

# Install bioinformatics tools
echo "[3/8] Installing bioinformatics tools..."
sudo apt-get install -y -qq \
    fastqc \
    bwa \
    samtools \
    bcftools \
    tabix

# Install fastp
echo "[4/8] Installing fastp..."
wget -q http://opengene.org/fastp/fastp
chmod a+x ./fastp
sudo mv ./fastp /usr/local/bin/

# Install GATK
echo "[5/8] Installing GATK 4.6.2.0..."
cd /tmp
wget -q https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip -q gatk-4.6.2.0.zip
sudo mv gatk-4.6.2.0 /opt/
sudo ln -sf /opt/gatk-4.6.2.0/gatk /usr/local/bin/gatk

# Download reference genome (hg19)
echo "[6/8] Downloading hg19 reference genome (~3GB, this will take a few minutes)..."
cd ~/NGS
mkdir -p reference
cd reference

# Download hg19
if [ ! -f hg19.fa ]; then
    wget -q --show-progress http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip hg19.fa.gz
fi

# Index reference genome
echo "[7/8] Indexing reference genome (this will take 15-30 minutes)..."
if [ ! -f hg19.fa.bwt ]; then
    bwa index hg19.fa
fi

if [ ! -f hg19.fa.fai ]; then
    samtools faidx hg19.fa
fi

if [ ! -f hg19.dict ]; then
    gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
fi

# Download snpEff for annotation
echo "[8/8] Setting up snpEff..."
cd ~/NGS
mkdir -p tools
cd tools

if [ ! -f snpEff/snpEff.jar ]; then
    wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -q snpEff_latest_core.zip
    cd snpEff
    java -jar snpEff.jar download -v GRCh37.75
fi

echo ""
echo "=========================================="
echo "✅ Setup Complete!"
echo "=========================================="
echo ""
echo "Installation Summary:"
echo "  ✅ System tools installed"
echo "  ✅ Bioinformatics tools installed"
echo "  ✅ Reference genome downloaded and indexed"
echo "  ✅ snpEff configured"
echo ""
echo "Next steps:"
echo "  1. Upload your FASTQ files to ~/NGS/data/"
echo "  2. Run: bash run_pipeline.sh data/r1.fastq.gz data/r2.fastq.gz output_name"
echo ""
echo "Storage used: $(du -sh ~/NGS | cut -f1)"
echo "=========================================="


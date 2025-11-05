#!/bin/bash
# Download files from LimeWire share links in Jarvis Lab
# Usage: bash download_from_limewire.sh

set -e

echo "=========================================="
echo "LimeWire File Downloader for Jarvis Lab"
echo "=========================================="
echo ""

# Create data directory
mkdir -p ~/NGS/data
cd ~/NGS/data

# LimeWire share links
LINK1="https://limewire.com/d/8xMLO#oJz6URSR83"  # Alaa_R1.fastq.gz (2.64GB)
LINK2="https://limewire.com/d/449bp#FnYO9xuzW4"  # Alaa_R2.fastq.gz (2.56GB)

echo "Attempting to download files from LimeWire..."
echo ""

# Method 1: Try with wget (usually doesn't work for LimeWire, but worth trying)
echo "[Method 1] Trying direct wget download..."
echo ""

# Download R1
echo "Downloading Alaa_R1.fastq.gz (2.64GB)..."
if wget --progress=bar:force:noscroll \
    --user-agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36" \
    --referer="https://limewire.com/" \
    -O Alaa_R1.fastq.gz \
    "$LINK1" 2>&1 | grep -q "200 OK\|saved"; then
    echo "✅ Alaa_R1.fastq.gz downloaded successfully!"
else
    echo "❌ Direct wget failed (expected - LimeWire requires authentication)"
fi

# Download R2
echo ""
echo "Downloading Alaa_R2.fastq.gz (2.56GB)..."
if wget --progress=bar:force:noscroll \
    --user-agent="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36" \
    --referer="https://limewire.com/" \
    -O Alaa_R2.fastq.gz \
    "$LINK2" 2>&1 | grep -q "200 OK\|saved"; then
    echo "✅ Alaa_R2.fastq.gz downloaded successfully!"
else
    echo "❌ Direct wget failed (expected - LimeWire requires authentication)"
fi

echo ""
echo "=========================================="
echo "⚠️  Direct download likely failed"
echo "=========================================="
echo ""
echo "LimeWire share links require browser authentication."
echo "Please use one of these methods:"
echo ""
echo "📥 METHOD A: Browser Download (Recommended)"
echo "  1. In Jarvis Lab, open JupyterLab"
echo "  2. Click 'File' → 'New' → 'Terminal'"
echo "  3. Open a new browser tab and paste these links:"
echo "     $LINK1"
echo "     $LINK2"
echo "  4. Download the files from your browser"
echo "  5. In JupyterLab, navigate to ~/NGS/data/"
echo "  6. Click 'Upload' button and select the downloaded files"
echo ""
echo "📥 METHOD B: Download Locally & Upload via SCP"
echo "  1. Download files on your local machine from the links above"
echo "  2. From your local machine, run:"
echo "     scp Alaa_R1.fastq.gz ubuntu@YOUR_JARVIS_IP:~/NGS/data/"
echo "     scp Alaa_R2.fastq.gz ubuntu@YOUR_JARVIS_IP:~/NGS/data/"
echo ""
echo "📥 METHOD C: Use Python requests with authentication (if you have login)"
echo "  Run this Python script:"
echo "  python3 << 'PYEOF'"
echo "import requests"
echo "session = requests.Session()"
echo "# Add your LimeWire cookies/authentication here"
echo "# Then download using session.get(url, stream=True)"
echo "PYEOF"
echo ""
echo "✅ Once files are in ~/NGS/data/, run:"
echo "   cd ~/NGS && bash run_pipeline.sh data/Alaa_R1.fastq.gz data/Alaa_R2.fastq.gz Alaa_sample 32"
echo ""
echo "=========================================="

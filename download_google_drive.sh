#!/bin/bash
# Download samples from Google Drive
# Usage: bash download_google_drive.sh

echo "════════════════════════════════════════════════════════════"
echo "Downloading Samples from Google Drive"
echo "════════════════════════════════════════════════════════════"
echo ""

# Install gdown if needed
if ! command -v gdown &> /dev/null; then
    echo "Installing gdown..."
    pip install gdown --quiet
fi

# Create data directory
mkdir -p ~/NGS/data
cd ~/NGS/data

echo "Sample 1:"
echo "  Downloading from: 15Pt-p3Y04NrczmBveKSPwgxWs4SQU2WM"
gdown --folder "https://drive.google.com/drive/folders/15Pt-p3Y04NrczmBveKSPwgxWs4SQU2WM?usp=sharing" -O sample_1

echo ""
echo "Sample 2:"
echo "  Downloading from: 1WJUF390XGbdzhM4GbsaqzbFbl26DepYq"
gdown --folder "https://drive.google.com/drive/folders/1WJUF390XGbdzhM4GbsaqzbFbl26DepYq?usp=sharing" -O sample_2

echo ""
echo "Sample 3:"
echo "  Downloading from: 1w03p51DXNaNXd6ZaGvrcKdV61kGS6rZ1"
gdown --folder "https://drive.google.com/drive/folders/1w03p51DXNaNXd6ZaGvrcKdV61kGS6rZ1?usp=sharing" -O sample_3

echo ""
echo "✅ Download complete!"
echo ""
echo "Files downloaded to: ~/NGS/data/"
ls -lh ~/NGS/data/


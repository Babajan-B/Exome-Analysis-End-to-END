#!/bin/bash
# Cleanup and Download Samples from Google Drive
# Usage: bash cleanup_and_download.sh

echo "╔════════════════════════════════════════════════════════════╗"
echo "║        CLEANUP & DOWNLOAD NEW SAMPLES                      ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Step 1: Check current storage
echo "════════════════════════════════════════════════════════════"
echo "STEP 1: Checking Current Storage"
echo "════════════════════════════════════════════════════════════"
echo ""

TOTAL_BEFORE=$(du -sh ~/NGS | cut -f1)
echo "Current NGS storage: $TOTAL_BEFORE"
echo ""
echo "Breakdown:"
printf "%-20s %10s\n" "Component" "Size"
echo "────────────────────────────────────────────"
du -sh ~/NGS/tools ~/NGS/reference ~/NGS/results ~/NGS/data ~/NGS/uploads ~/NGS/Bur 2>/dev/null | \
    awk '{printf "%-20s %10s\n", $2, $1}'
echo "────────────────────────────────────────────"
echo ""

# Check disk space
DISK_AVAIL=$(df -h ~ | tail -1 | awk '{print $4}')
DISK_TOTAL=$(df -h ~ | tail -1 | awk '{print $2}')
DISK_USED=$(df -h ~ | tail -1 | awk '{print $3}')
echo "Disk space:"
echo "  Total:  $DISK_TOTAL"
echo "  Used:   $DISK_USED"
echo "  Free:   $DISK_AVAIL"
echo ""

read -p "Continue with cleanup? (y/n): " confirm1
if [[ ! "$confirm1" =~ ^[Yy]$ ]]; then
    echo "Cleanup cancelled."
    exit 0
fi

# Step 2: Cleanup
echo ""
echo "════════════════════════════════════════════════════════════"
echo "STEP 2: Cleaning Up Old Data"
echo "════════════════════════════════════════════════════════════"
echo ""

echo "Removing old results..."
rm -rf ~/NGS/results 2>/dev/null && echo "  ✅ Results removed" || echo "  ⚠️  No results folder"

echo "Removing old data..."
rm -rf ~/NGS/data 2>/dev/null && echo "  ✅ Data removed" || echo "  ⚠️  No data folder"

echo "Removing uploads..."
rm -rf ~/NGS/uploads 2>/dev/null && echo "  ✅ Uploads removed" || echo "  ⚠️  No uploads folder"

echo "Removing Bur folder..."
rm -rf ~/NGS/Bur 2>/dev/null && echo "  ✅ Bur removed" || echo "  ⚠️  No Bur folder"

echo "Removing zip files..."
rm -f ~/NGS/*.zip 2>/dev/null && echo "  ✅ Zip files removed" || echo "  ⚠️  No zip files"

echo "Removing temp files..."
rm -f ~/NGS/.annovar_installing ~/NGS/.annovar_installed 2>/dev/null
echo "  ✅ Temp files removed"

TOTAL_AFTER=$(du -sh ~/NGS | cut -f1)
echo ""
echo "Storage after cleanup: $TOTAL_AFTER"
echo ""

# Step 3: Download samples from Google Drive
echo "════════════════════════════════════════════════════════════"
echo "STEP 3: Downloading Samples from Google Drive"
echo "════════════════════════════════════════════════════════════"
echo ""

# Install gdown if not available
if ! command -v gdown &> /dev/null; then
    echo "Installing gdown for Google Drive download..."
    pip install gdown --quiet
fi

# Create data directory
mkdir -p ~/NGS/data
cd ~/NGS/data

echo "Downloading Sample 1..."
echo "  Link: https://drive.google.com/drive/folders/15Pt-p3Y04NrczmBveKSPwgxWs4SQU2WM?usp=sharing"
gdown --folder "https://drive.google.com/drive/folders/15Pt-p3Y04NrczmBveKSPwgxWs4SQU2WM?usp=sharing" -O sample_1 || \
    echo "  ⚠️  Download failed - try manual download"

echo ""
echo "Downloading Sample 2..."
echo "  Link: https://drive.google.com/drive/folders/1WJUF390XGbdzhM4GbsaqzbFbl26DepYq?usp=sharing"
gdown --folder "https://drive.google.com/drive/folders/1WJUF390XGbdzhM4GbsaqzbFbl26DepYq?usp=sharing" -O sample_2 || \
    echo "  ⚠️  Download failed - try manual download"

echo ""
echo "Downloading Sample 3..."
echo "  Link: https://drive.google.com/drive/folders/1w03p51DXNaNXd6ZaGvrcKdV61kGS6rZ1?usp=sharing"
gdown --folder "https://drive.google.com/drive/folders/1w03p51DXNaNXd6ZaGvrcKdV61kGS6rZ1?usp=sharing" -O sample_3 || \
    echo "  ⚠️  Download failed - try manual download"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "Checking downloaded files..."
echo "════════════════════════════════════════════════════════════"
echo ""

for dir in sample_1 sample_2 sample_3; do
    if [ -d "$dir" ]; then
        echo "✅ $dir:"
        ls -lh "$dir" | head -5
        echo ""
    fi
done

echo "════════════════════════════════════════════════════════════"
echo "✅ Setup Complete!"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Next step:"
echo "  bash MASTER_PIPELINE.sh ~/NGS/data 16"
echo ""
echo "Or organize files first:"
echo "  Check ~/NGS/data/ for downloaded files"
echo "  Ensure R1 and R2 FASTQ files are paired correctly"
echo ""


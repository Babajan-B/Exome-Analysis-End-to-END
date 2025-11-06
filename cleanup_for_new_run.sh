#!/bin/bash
# Cleanup Script - Remove old results and data for fresh start
# Keeps: tools, reference genome, scripts
# Removes: results, data, uploads, temporary files

echo "╔════════════════════════════════════════════════════════════╗"
echo "║            CLEANUP FOR NEW ANALYSIS                        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "This will DELETE:"
echo "  ❌ All results (~/NGS/results/)"
echo "  ❌ All uploaded data (~/NGS/data/)"
echo "  ❌ All uploads (~/NGS/uploads/)"
echo "  ❌ Temporary files (~/NGS/Bur/, *.zip)"
echo ""
echo "This will KEEP:"
echo "  ✅ Tools (ANNOVAR, GATK, etc.)"
echo "  ✅ Reference genome (hg19)"
echo "  ✅ All scripts"
echo ""

# Show current sizes
echo "Current storage:"
echo "────────────────────────────────────────────────────────────"
du -sh ~/NGS/results 2>/dev/null && echo "  Results folder" || echo "  Results: Not found"
du -sh ~/NGS/data 2>/dev/null && echo "  Data folder" || echo "  Data: Not found"
du -sh ~/NGS/uploads 2>/dev/null && echo "  Uploads folder" || echo "  Uploads: Not found"
du -sh ~/NGS/Bur 2>/dev/null && echo "  Bur folder" || echo "  Bur: Not found"
du -sh ~/NGS/*.zip 2>/dev/null && echo "  Zip files" || echo "  Zip files: None"
echo ""

TOTAL_BEFORE=$(du -sh ~/NGS | cut -f1)
echo "Total NGS storage: $TOTAL_BEFORE"
echo ""

read -p "⚠️  Continue with cleanup? (y/n): " confirm

if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Cleanup cancelled."
    exit 0
fi

echo ""
echo "Starting cleanup..."
echo ""

# Remove results
if [ -d ~/NGS/results ]; then
    echo "Removing results folder..."
    rm -rf ~/NGS/results
    echo "✅ Results removed"
fi

# Remove data
if [ -d ~/NGS/data ]; then
    echo "Removing data folder..."
    rm -rf ~/NGS/data
    echo "✅ Data removed"
fi

# Remove uploads
if [ -d ~/NGS/uploads ]; then
    echo "Removing uploads folder..."
    rm -rf ~/NGS/uploads
    echo "✅ Uploads removed"
fi

# Remove Bur folder
if [ -d ~/NGS/Bur ]; then
    echo "Removing Bur folder..."
    rm -rf ~/NGS/Bur
    echo "✅ Bur removed"
fi

# Remove zip files
if ls ~/NGS/*.zip 1> /dev/null 2>&1; then
    echo "Removing zip files..."
    rm -f ~/NGS/*.zip
    echo "✅ Zip files removed"
fi

# Remove temporary files
echo "Removing temporary files..."
rm -f ~/NGS/.annovar_installing ~/NGS/.annovar_installed
rm -rf ~/NGS/test_data 2>/dev/null
echo "✅ Temp files removed"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "✅ Cleanup Complete!"
echo "════════════════════════════════════════════════════════════"
echo ""

TOTAL_AFTER=$(du -sh ~/NGS | cut -f1)
echo "Storage before: $TOTAL_BEFORE"
echo "Storage after:  $TOTAL_AFTER"
echo ""

echo "Ready for new analysis!"
echo ""
echo "Next steps:"
echo "  1. Create data folder: mkdir -p ~/NGS/data"
echo "  2. Upload new FASTQ files to ~/NGS/data/"
echo "  3. Run pipeline: bash MASTER_PIPELINE.sh"
echo ""
echo "════════════════════════════════════════════════════════════"


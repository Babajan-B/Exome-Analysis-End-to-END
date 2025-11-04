#!/bin/bash

# Script to download and install IGV (Integrative Genomics Viewer)
# This script will download and set up IGV for visualizing BAM files

echo "Setting up IGV (Integrative Genomics Viewer)..."

# Create IGV directory
mkdir -p tools/igv
cd tools/igv

# Download IGV (current version)
echo "Downloading IGV..."
wget https://data.broadinstitute.org/igv/projects/downloads/2.15/IGV_2.15.2.zip

# Unzip IGV
echo "Extracting IGV..."
unzip IGV_2.15.2.zip
rm IGV_2.15.2.zip

# Move back to the root directory
cd ../..

# Create a convenience script to launch IGV
cat > igv.sh << 'EOF'
#!/bin/bash
cd $(dirname $0)
java -Xmx4g -jar tools/igv/IGV_2.15.2/igv.jar "$@"
EOF

# Make script executable
chmod +x igv.sh

echo "IGV setup complete."
echo ""
echo "You can now launch IGV with the following command:"
echo "./igv.sh"
echo ""
echo "To visualize a BAM file, use one of these options:"
echo "1. Launch IGV and use File → Load from File to open your BAM file"
echo "2. Launch IGV with a specific BAM file: ./igv.sh /path/to/your/file.bam"
echo ""
echo "Note: Make sure your BAM file is indexed (has a .bai file) for proper visualization" 
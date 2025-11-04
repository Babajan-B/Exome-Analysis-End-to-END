#!/bin/bash

# Exit on any error
set -e

echo "Setting up NGS Exome Analysis Web Application"
echo "============================================="

# Check if running on macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS system"
    
    # Check if Homebrew is installed
    if ! command -v brew &> /dev/null; then
        echo "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    fi
    
    # Install dependencies using Homebrew
    echo "Installing bioinformatics tools with Homebrew..."
    brew install fastqc
    brew install fastp
    brew install bwa
    brew install samtools
    
    # Create a Python virtual environment if it doesn't exist
    if [ ! -d "venv" ]; then
        echo "Creating Python virtual environment..."
        python -m venv venv
    fi
    
    # Activate the virtual environment and install Python packages
    echo "Installing Python dependencies..."
    source venv/bin/activate
    pip install -r requirements.txt
    
else
    echo "This script is designed for macOS systems."
    echo "For Linux installation, you'll need to manually install the required tools:"
    echo "- FastQC"
    echo "- fastp"
    echo "- BWA"
    echo "- Samtools"
    echo "- GATK"
    echo ""
    echo "And then set up the Python environment with:"
    echo "python -m venv venv"
    echo "source venv/bin/activate"
    echo "pip install -r requirements.txt"
fi

# Create necessary directories if they don't exist
echo "Creating necessary directories..."
mkdir -p uploads
mkdir -p results
mkdir -p reference

echo ""
echo "Setup completed. You can start the application by running:"
echo "source venv/bin/activate"
echo "python run.py" 
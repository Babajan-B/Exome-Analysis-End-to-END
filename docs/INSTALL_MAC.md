# Installation Guide for macOS

This guide will walk you through installing the NGS Exome Analysis Pipeline on macOS (10.15 Catalina or later).

## 📋 Prerequisites

- macOS 10.15 (Catalina) or later
- Administrator access
- At least 50GB free disk space
- Internet connection

## 🚀 Step-by-Step Installation

### Step 1: Install Xcode Command Line Tools

The Xcode Command Line Tools provide essential development tools for macOS:

```bash
# Install Xcode Command Line Tools
xcode-select --install
```

Click "Install" in the popup dialog that appears.

**Verify Installation:**
```bash
xcode-select -p
# Should output: /Library/Developer/CommandLineTools
```

### Step 2: Install Homebrew

Homebrew is a package manager for macOS that makes installing tools easy:

```bash
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

**For Apple Silicon (M1/M2/M3) Macs**, add Homebrew to your PATH:

```bash
# Add to ~/.zshrc
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zshrc
source ~/.zshrc
```

**For Intel Macs**, Homebrew should automatically be in your PATH.

**Verify Installation:**
```bash
brew --version
```

### Step 3: Install Python 3

```bash
# Install Python 3
brew install python@3.10

# Verify installation
python3 --version
pip3 --version
```

### Step 4: Install Java

Java is required for GATK and other tools:

```bash
# Install OpenJDK
brew install openjdk@11

# Link Java
sudo ln -sfn /opt/homebrew/opt/openjdk@11/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-11.jdk

# For Intel Macs, use:
# sudo ln -sfn /usr/local/opt/openjdk@11/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-11.jdk
```

**Verify Installation:**
```bash
java -version
```

### Step 5: Install Bioinformatics Tools

Install all required bioinformatics tools using Homebrew:

```bash
# Add bioinformatics tap
brew tap brewsci/bio

# Install core tools
brew install fastqc fastp bwa samtools bcftools

# Verify installations
fastqc --version
fastp --version
bwa
samtools --version
```

### Step 6: Install GATK

```bash
# Install GATK
brew install gatk

# Verify installation
gatk --version
```

**Alternative: Manual Installation**

If you prefer to install GATK manually:

```bash
# Download GATK
cd ~/Downloads
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip

# Extract
unzip gatk-4.6.2.0.zip

# Move to a permanent location
sudo mv gatk-4.6.2.0 /usr/local/gatk-4.6.2.0

# Add to PATH (add to ~/.zshrc)
echo 'export PATH="/usr/local/gatk-4.6.2.0:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### Step 7: Clone the Repository

```bash
# Navigate to where you want to install
cd ~/Documents  # or any directory you prefer

# Clone the repository
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git

# Navigate to the project directory
cd Exome-Analysis-End-to-END
```

### Step 8: Set Up Python Virtual Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install Python dependencies
pip install -r requirements.txt
```

### Step 9: Run the Automated Setup Script

The project includes a setup script that automates the remaining configuration:

```bash
# Make the setup script executable
chmod +x setup.sh

# Run the setup script
./setup.sh
```

The script will:
- Check for required tools
- Download the reference genome (hg19)
- Index the reference genome
- Set up snpEff
- Create necessary directories

### Step 10: Download Reference Genome (if not done by setup.sh)

If the setup script didn't download the reference genome:

```bash
# Create reference directory
mkdir -p reference
cd reference

# Download hg19 reference genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# Extract
gunzip hg19.fa.gz

# Index the reference for BWA
bwa index hg19.fa

# Create FASTA index
samtools faidx hg19.fa

# Create sequence dictionary for GATK
gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict

# Return to project directory
cd ..
```

**Alternative: Use hg38 (newer reference)**

```bash
# Download hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index hg38.fa
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
```

### Step 11: Download Known Variants for GATK (Optional but Recommended)

Known variant sites improve BQSR accuracy:

```bash
# Create known sites directory
mkdir -p reference/known_sites
cd reference/known_sites

# Download dbSNP (for hg19)
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi

# Rename for convenience
mv 00-All.vcf.gz dbsnp_138.hg19.vcf.gz
mv 00-All.vcf.gz.tbi dbsnp_138.hg19.vcf.gz.tbi

cd ../..
```

### Step 12: Configure the Application (Optional)

If tools are not in standard locations, edit `app/pipeline/pipeline.py`:

```python
# Tool paths (usually automatic with Homebrew)
FASTQC_PATH = "fastqc"
FASTP_PATH = "fastp"
BWA_PATH = "bwa"
SAMTOOLS_PATH = "samtools"
GATK_PATH = "gatk"  # or full path if installed manually
```

### Step 13: Run the Application

```bash
# Activate virtual environment (if not already activated)
source venv/bin/activate

# Start the Flask server
python run.py
```

You should see output like:
```
 * Running on http://0.0.0.0:5008/ (Press CTRL+C to quit)
 * Restarting with stat
 * Debugger is active!
```

### Step 14: Access the Web Interface

Open your web browser and navigate to:
```
http://localhost:5008
```

## 🎯 Quick Test Run

Test your installation with sample data:

```bash
# Download test data
mkdir -p test_data
cd test_data

# Small test FASTQ files (example URLs, adjust as needed)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR000/ERR000589/ERR000589_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR000/ERR000589/ERR000589_2.fastq.gz

cd ..
```

Then upload these files through the web interface.

## 🍎 macOS-Specific Notes

### Apple Silicon (M1/M2/M3) Considerations

1. **Rosetta 2**: Some tools may require Rosetta 2
   ```bash
   softwareupdate --install-rosetta --agree-to-license
   ```

2. **Architecture-specific paths**: Homebrew installs to different locations:
   - Apple Silicon: `/opt/homebrew`
   - Intel: `/usr/local`

3. **Check architecture**:
   ```bash
   uname -m
   # arm64 = Apple Silicon
   # x86_64 = Intel
   ```

### Performance Optimization

1. **Increase file descriptor limits**:
   ```bash
   # Add to ~/.zshrc
   ulimit -n 4096
   ```

2. **Use SSD storage**: Store data on SSD for faster I/O

3. **Memory considerations**: Close unnecessary applications during analysis

### Security Settings

macOS may block some tools. If you get security warnings:

1. Go to **System Preferences → Security & Privacy**
2. Click "Allow Anyway" for blocked applications
3. Or use: `sudo spctl --master-disable` (not recommended for security)

## ❓ Troubleshooting

### Issue: "command not found" errors

**Solution**: Check if Homebrew is in your PATH:
```bash
echo $PATH
# Should include /opt/homebrew/bin or /usr/local/bin

# Add to PATH if missing
echo 'export PATH="/opt/homebrew/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

### Issue: Port 5008 already in use

**Solution**: Change port in `run.py`:
```python
app.run(debug=True, host='0.0.0.0', port=5009)
```

Or kill the process using the port:
```bash
lsof -ti:5008 | xargs kill
```

### Issue: Java not found

**Solution**: Ensure Java is properly linked:
```bash
java -version
# If not working, check:
/usr/libexec/java_home -V

# Set JAVA_HOME
export JAVA_HOME=$(/usr/libexec/java_home)
echo 'export JAVA_HOME=$(/usr/libexec/java_home)' >> ~/.zshrc
```

### Issue: Permission denied

**Solution**: Don't use `sudo` with pip in virtual environment:
```bash
# Activate venv first
source venv/bin/activate
# Then install without sudo
pip install -r requirements.txt
```

### Issue: Out of memory

**Solution**: 
- Close other applications
- Reduce number of threads in pipeline
- Use smaller test datasets first

### Issue: Reference genome download is slow

**Solution**: Use a mirror or download using a download manager:
```bash
# Use aria2c for faster downloads
brew install aria2
aria2c -x 16 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```

## 🔄 Updating the Application

To update to the latest version:

```bash
cd ~/Documents/Exome-Analysis-End-to-END
git pull origin main

# Activate virtual environment
source venv/bin/activate

# Update dependencies
pip install -r requirements.txt --upgrade
```

## 🧹 Uninstalling

To completely remove the application:

```bash
# Remove the project directory
rm -rf ~/Documents/Exome-Analysis-End-to-END

# Optionally remove installed tools
brew uninstall fastqc fastp bwa samtools gatk
```

## 📚 Additional Resources

- [Homebrew Documentation](https://docs.brew.sh/)
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
- [Samtools Documentation](http://www.htslib.org/doc/samtools.html)

## 🆘 Getting Help

If you encounter issues:
1. Check the troubleshooting section above
2. Review logs in `results/[analysis_id]/pipeline.log`
3. Check terminal output for error messages
4. Open an issue on [GitHub](https://github.com/Babajan-B/Exome-Analysis-End-to-END/issues)

---

**macOS Installation Complete!** 🎉

You can now:
- Access the web interface at http://localhost:5008
- Upload FASTQ files for analysis
- Monitor pipeline progress in real-time
- Download results and reports

Enjoy analyzing your exome data!


# Installation Guide for Windows

This guide will walk you through installing the NGS Exome Analysis Pipeline on Windows 10 or Windows 11.

## 📋 Prerequisites

Before starting, ensure you have administrator access on your Windows machine.

## 🚀 Step-by-Step Installation

### Step 1: Install Python

1. **Download Python**:
   - Go to [https://www.python.org/downloads/](https://www.python.org/downloads/)
   - Download Python 3.10 or higher for Windows
   - **Important**: During installation, check "Add Python to PATH"

2. **Verify Installation**:
   ```cmd
   python --version
   pip --version
   ```

### Step 2: Install Java

1. **Download Java**:
   - Go to [https://www.oracle.com/java/technologies/downloads/](https://www.oracle.com/java/technologies/downloads/)
   - Download Java SE Development Kit (JDK) 11 or higher
   - Run the installer with default settings

2. **Verify Installation**:
   ```cmd
   java -version
   ```

### Step 3: Install Git

1. **Download Git**:
   - Go to [https://git-scm.com/download/win](https://git-scm.com/download/win)
   - Download and install Git for Windows
   - Use default settings during installation

2. **Verify Installation**:
   ```cmd
   git --version
   ```

### Step 4: Install WSL2 (Windows Subsystem for Linux) - Recommended

Many bioinformatics tools work better on Linux. We recommend using WSL2:

1. **Enable WSL2**:
   ```powershell
   # Run PowerShell as Administrator
   wsl --install
   ```

2. **Restart your computer** when prompted

3. **Set up Ubuntu**:
   - Open "Ubuntu" from Start menu
   - Create a username and password
   - Update packages:
   ```bash
   sudo apt update && sudo apt upgrade -y
   ```

4. **Install bioinformatics tools in WSL2**:
   ```bash
   # Install dependencies
   sudo apt install -y build-essential wget curl default-jre default-jdk \
                       python3-pip fastqc samtools bwa bcftools tabix unzip

   # Install fastp
   wget http://opengene.org/fastp/fastp
   chmod a+x ./fastp
   sudo mv ./fastp /usr/local/bin/
   ```

### Step 5: Clone the Repository

**Option A: Using WSL2 (Recommended)**

```bash
# In WSL2 Ubuntu terminal
cd ~
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git
cd Exome-Analysis-End-to-END
```

**Option B: Using Windows Command Prompt**

```cmd
# In Windows Command Prompt or PowerShell
cd %USERPROFILE%\Documents
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git
cd Exome-Analysis-End-to-END
```

### Step 6: Set Up Python Virtual Environment

**Using WSL2:**

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Install Python dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

**Using Windows:**

```cmd
# Create virtual environment
python -m venv venv

# Activate virtual environment
venv\Scripts\activate

# Install Python dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

### Step 7: Download Reference Genome

**Using WSL2:**

```bash
# Create reference directory
mkdir -p reference

# Download hg19 reference (or use hg38)
cd reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# Index the reference
bwa index hg19.fa
samtools faidx hg19.fa

# Create sequence dictionary
java -jar ../tools/picard.jar CreateSequenceDictionary \
    R=hg19.fa O=hg19.dict

cd ..
```

**Using Windows (download manually):**

1. Download hg19.fa.gz from [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)
2. Extract using 7-Zip or similar tool
3. Place in the `reference/` folder
4. Index using WSL2 or install tools on Windows

### Step 8: Download and Set Up GATK

```bash
# Download GATK (if not included)
cd ~
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip

# The gatk folder should be in your project directory
# Or update paths in pipeline.py to point to GATK location
```

### Step 9: Configure the Application

**Edit `app/pipeline/pipeline.py`** if needed to update tool paths:

```python
# Example paths for WSL2
FASTQC_PATH = "fastqc"
FASTP_PATH = "fastp"
BWA_PATH = "bwa"
SAMTOOLS_PATH = "samtools"
GATK_PATH = "/mnt/c/Users/YourName/Documents/Exome-Analysis-End-to-END/gatk-4.6.2.0/gatk"
```

### Step 10: Run the Application

**Using WSL2:**

```bash
# Activate virtual environment
source venv/bin/activate

# Start the server
python run.py
```

**Using Windows:**

```cmd
# Activate virtual environment
venv\Scripts\activate

# Start the server
python run.py
```

### Step 11: Access the Web Interface

Open your web browser and navigate to:
```
http://localhost:5008
```

## 🔧 Alternative: Using Docker (Easier Option)

If you have Docker Desktop installed, you can use a containerized version:

```bash
# Pull and run the Docker container (if available)
docker pull babajanb/ngs-exome-analysis:latest
docker run -p 5008:5008 -v ${PWD}/uploads:/app/uploads -v ${PWD}/results:/app/results babajanb/ngs-exome-analysis:latest
```

## 📝 Windows-Specific Notes

### File Paths
- Windows uses backslashes (`\`) in paths
- WSL2 can access Windows files at `/mnt/c/Users/...`
- Convert paths when needed: `wslpath 'C:\Users\...'`

### Performance Tips
- Use WSL2 for better performance with bioinformatics tools
- Store data files on WSL2 filesystem (not Windows) for faster I/O
- Allocate at least 8GB RAM to WSL2 in `.wslconfig`

### WSL2 Configuration

Create `%USERPROFILE%\.wslconfig`:

```ini
[wsl2]
memory=16GB
processors=4
swap=8GB
```

## ❓ Troubleshooting

### Issue: "Command not found" in Windows
- **Solution**: Make sure tools are in your PATH or use full paths in configuration

### Issue: WSL2 is slow
- **Solution**: Store files on WSL2 filesystem, not Windows filesystem
- Move project to `~/` in WSL2 for better performance

### Issue: Port 5008 already in use
- **Solution**: Change port in `run.py`:
  ```python
  app.run(debug=True, host='0.0.0.0', port=5009)
  ```

### Issue: Java not found
- **Solution**: Verify Java installation and PATH:
  ```cmd
  where java
  echo %JAVA_HOME%
  ```

### Issue: Out of memory
- **Solution**: Increase WSL2 memory allocation in `.wslconfig`

## 📚 Additional Resources

- [WSL2 Documentation](https://docs.microsoft.com/en-us/windows/wsl/)
- [Python Windows Installation](https://docs.python.org/3/using/windows.html)
- [GATK Documentation](https://gatk.broadinstitute.org/hc/en-us)
- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)

## 🆘 Getting Help

If you encounter issues:
1. Check the troubleshooting section above
2. Review logs in `results/[analysis_id]/pipeline.log`
3. Open an issue on [GitHub](https://github.com/Babajan-B/Exome-Analysis-End-to-END/issues)

---

**Windows Installation Complete!** You can now run analyses through the web interface at http://localhost:5008


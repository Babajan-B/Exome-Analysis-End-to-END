# NGS Exome Analysis Pipeline - Web Application

A comprehensive web-based application for processing and analyzing Next-Generation Sequencing (NGS) exome data. This tool provides an intuitive interface for uploading paired-end FASTQ files and running a complete exome analysis pipeline.

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue)](https://www.python.org/)
[![Flask](https://img.shields.io/badge/Flask-2.0.1-green)](https://flask.palletsprojects.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 🚀 Features

- **User-Friendly Web Interface** - Simple drag-and-drop file upload
- **Complete Analysis Pipeline**:
  - ✅ Quality control using FastQC
  - ✅ Read trimming and filtering using fastp
  - ✅ Alignment to reference genome (hg19/hg38) using BWA
  - ✅ SAM to BAM conversion and sorting
  - ✅ Duplicate marking using GATK MarkDuplicates
  - ✅ Base Quality Score Recalibration (BQSR)
  - ✅ Variant calling using GATK HaplotypeCaller
  - ✅ Variant annotation using snpEff
- **Real-time Progress Tracking** - Monitor your analysis as it runs
- **Interactive Visualization** - Built-in IGV.js browser for BAM file visualization
- **Downloadable Results** - VCF files, QC reports, and more
- **Multiple Input Options** - Support for FASTQ, BAM, and VCF files
- **Resume Capability** - Continue interrupted analyses

## 📋 System Requirements

### Minimum Requirements
- **OS**: Windows 10/11, macOS 10.15+, or Linux (Ubuntu 20.04+)
- **RAM**: 8GB (16GB+ recommended for whole exome)
- **Storage**: 50GB free disk space (more for reference genomes and analysis results)
- **CPU**: 4+ cores recommended
- **Internet**: Required for downloading reference genomes, GATK, and dependencies

### Software Prerequisites
- Python 3.7 or higher
- Java Runtime Environment (JRE) 8 or higher
- Git (for cloning the repository)

### Large Files NOT Included in Repository

The following files are **NOT** included in the Git repository due to their size:

1. **Reference Genome** (~3GB) - Download from UCSC Genome Browser
2. **GATK JAR files** (~200MB) - Download from Broad Institute
3. **snpEff Database** (~700MB) - Download using snpEff
4. **Picard Tools** (~15MB) - Download from Broad Institute

**All these files will be downloaded automatically when you run the `setup.sh` script**, or you can download them manually following the instructions in the installation guides.

## 🛠️ Installation

### Quick Start

Choose your operating system:

- **[Windows Installation Guide](docs/INSTALL_WINDOWS.md)** - Step-by-step instructions for Windows
- **[macOS Installation Guide](docs/INSTALL_MAC.md)** - Step-by-step instructions for macOS
- **[Linux Installation Guide](docs/INSTALL_LINUX.md)** - Step-by-step instructions for Linux

### Quick Install (macOS/Linux)

```bash
# Clone the repository
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git
cd Exome-Analysis-End-to-END

# Run the automated setup script
chmod +x setup.sh
./setup.sh

# Activate virtual environment
source venv/bin/activate

# Start the application
python run.py
```

### Quick Install (Windows)

```cmd
# Clone the repository
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git
cd Exome-Analysis-End-to-END

# Run the setup script
setup_windows.bat

# Activate virtual environment
venv\Scripts\activate

# Start the application
python run.py
```

## 📦 Downloading Required Large Files

**Important:** These files are NOT included in the repository. You must download them separately.

### Option 1: Automatic Download (Recommended)

Run the setup script which will download all required files automatically:

```bash
# macOS/Linux
./setup.sh

# Windows (in WSL2)
bash setup.sh
```

### Option 2: Manual Download

If you prefer to download files manually:

#### 1. Download GATK

```bash
# Download GATK 4.6.2.0
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip
# Move to project directory if needed
```

Or download from: [GATK Releases](https://github.com/broadinstitute/gatk/releases)

#### 2. Download Reference Genome

**Option A: hg19 (GRCh37)**
```bash
cd reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# Index the reference
bwa index hg19.fa
samtools faidx hg19.fa
gatk CreateSequenceDictionary -R hg19.fa -O hg19.dict
```

**Option B: hg38 (GRCh38 - newer)**
```bash
cd reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Index the reference
bwa index hg38.fa
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
```

#### 3. Download snpEff JAR files

```bash
# Download snpEff
cd tools/snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
mv snpEff/* .
rmdir snpEff
```

Or download from: [snpEff Download](http://pcingola.github.io/SnpEff/download/)

#### 4. Download snpEff Database

```bash
cd tools/snpEff
# For hg19
java -jar snpEff.jar download -v GRCh37.75

# For hg38
java -jar snpEff.jar download -v GRCh38.99
```

#### 5. Download Picard Tools (Optional)

```bash
cd tools
wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
```

### File Sizes

- **hg19 reference genome**: ~3GB (compressed), ~3.2GB (uncompressed)
- **hg38 reference genome**: ~3.1GB (compressed), ~3.3GB (uncompressed)  
- **GATK package**: ~200MB
- **snpEff core**: ~50MB
- **snpEff database (GRCh37.75)**: ~700MB
- **Picard tools**: ~15MB

**Total space required**: ~10-15GB including indexes and temporary files

## 🎯 Usage

### Starting the Server

1. **Activate the virtual environment** (if not already activated):
   ```bash
   # macOS/Linux
   source venv/bin/activate
   
   # Windows
   venv\Scripts\activate
   ```

2. **Start the Flask server**:
   ```bash
   python run.py
   ```

3. **Open your browser** and navigate to:
   ```
   http://localhost:5008
   ```

### Running an Analysis

1. **Upload FASTQ Files**:
   - Click on the "Upload FASTQ Files" section
   - Select your R1 (forward) and R2 (reverse) FASTQ files
   - Optionally skip QC or trimming steps
   - Click "Start Analysis"

2. **Monitor Progress**:
   - You'll be redirected to a status page
   - View real-time logs and progress updates
   - Each pipeline step is tracked individually

3. **Download Results**:
   - Once complete, download your results:
     - `variants.vcf` - Called variants
     - `variants.ann.vcf` - Annotated variants
     - `fastp_report.html` - Trimming report
     - `*_fastqc.html` - Quality control reports
     - `*.bam` and `*.bai` - Alignment files

### Alternative Input Options

**Upload BAM File**: If you already have aligned reads
**Upload VCF File**: For annotation only
**Direct File Path**: Use files already on the server
**Visualize BAM**: Use the built-in IGV.js browser

## 📁 Project Structure

```
NGS/
├── app/
│   ├── __init__.py           # Flask app initialization
│   ├── routes.py             # Web routes and handlers
│   ├── pipeline/
│   │   └── pipeline.py       # Main analysis pipeline
│   ├── static/               # CSS, JavaScript
│   └── templates/            # HTML templates
├── reference/                # Reference genome files
│   ├── hg19.fa              # Human reference genome
│   ├── hg19.fa.fai          # FASTA index
│   └── hg19.dict            # Sequence dictionary
├── tools/                    # Bioinformatics tools
│   ├── snpEff/              # Variant annotation
│   └── picard.jar           # Picard tools
├── gatk-4.6.2.0/            # GATK toolkit
├── uploads/                  # User uploaded files
├── results/                  # Analysis results
├── run.py                    # Application entry point
├── requirements.txt          # Python dependencies
└── setup.sh                  # Automated setup script
```

## 🔬 Pipeline Details

### Step-by-Step Process

1. **Quality Control (FastQC)**
   - Analyzes raw sequencing data quality
   - Generates detailed QC reports

2. **Trimming (fastp)**
   - Removes adapter sequences
   - Filters low-quality reads
   - Trims low-quality bases

3. **Alignment (BWA-MEM)**
   - Aligns reads to reference genome
   - Produces SAM file

4. **SAM Processing**
   - Converts SAM to BAM format
   - Sorts by coordinate
   - Creates index

5. **Duplicate Marking (GATK)**
   - Identifies PCR and optical duplicates
   - Marks duplicates without removing

6. **Base Quality Score Recalibration (BQSR)**
   - Corrects systematic errors in quality scores
   - Uses known variant sites

7. **Variant Calling (GATK HaplotypeCaller)**
   - Calls SNPs and small indels
   - Generates VCF file

8. **Variant Annotation (snpEff)**
   - Annotates variants with functional information
   - Predicts variant effects

## 🐛 Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| Command not found | Ensure tools are installed and in PATH |
| Out of memory | Reduce thread count or use more RAM |
| Port already in use | Change port in `run.py` or stop other process |
| Missing reference genome | Run setup script or download manually |
| Java not found | Install Java 8+ and add to PATH |

### Getting Help

- Check the logs in `results/[analysis_id]/pipeline.log`
- Review terminal output for error messages
- Ensure all dependencies are properly installed
- Open an [issue on GitHub](https://github.com/Babajan-B/Exome-Analysis-End-to-END/issues)

## 📊 Output Files

| File | Description |
|------|-------------|
| `variants.vcf` | Raw variant calls (VCF format) |
| `variants.ann.vcf` | Annotated variants with functional predictions |
| `sample.bam` | Aligned and processed reads |
| `sample.bam.bai` | BAM index file |
| `fastp_report.html` | Trimming and filtering report |
| `*_fastqc.html` | Quality control reports |
| `pipeline.log` | Complete pipeline execution log |
| `progress.json` | Pipeline progress tracking |

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

This application uses the following excellent open-source tools:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control
- [fastp](https://github.com/OpenGene/fastp) - Read preprocessing
- [BWA](https://github.com/lh3/bwa) - Read alignment
- [Samtools](http://www.htslib.org/) - SAM/BAM manipulation
- [GATK](https://gatk.broadinstitute.org/) - Variant calling
- [snpEff](http://pcingola.github.io/SnpEff/) - Variant annotation
- [IGV.js](https://github.com/igvteam/igv.js/) - Genome visualization
- [Flask](https://flask.palletsprojects.com/) - Web framework

## 📧 Contact

For questions or support, please open an issue on GitHub or contact the maintainer.

## 🌟 Citation

If you use this tool in your research, please cite:

```
NGS Exome Analysis Pipeline
https://github.com/Babajan-B/Exome-Analysis-End-to-END
```

---

**Made with ❤️ for the bioinformatics community**

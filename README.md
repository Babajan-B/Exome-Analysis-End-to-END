# NGS Exome Analysis Pipeline

A complete, automated pipeline for Next-Generation Sequencing (NGS) exome analysis. Process FASTQ files to annotated variants with a single command.

[![Python](https://img.shields.io/badge/Python-3.7%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## ⚡ Quick Start (3 Commands)

```bash
# 1. Clone the repository (run once)
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS
cd ~/NGS

# 2. Install every dependency (tools, reference genome, databases)
bash install_all.sh

# 3. Run the full pipeline on all FASTQ pairs in ~/NGS/data/
bash ULTIMATE_MASTER_PIPELINE.sh ~/NGS/data 16
```

After the first run you only need step 3.  
The pipeline creates a timestamped ZIP in `~/NGS/` with every report, annotation, and summary.

---

## 🎯 What It Does

### Complete End-to-End Pipeline

**Input:** Paired-end FASTQ files  
**Output:** Annotated variants, functional classifications, quality reports

**Automated Steps:**
1. Quality Control (FastQC)
2. Read Trimming (fastp)
3. Alignment (BWA-MEM to hg19)
4. BAM Processing (SAMtools, GATK)
5. Variant Calling (GATK HaplotypeCaller)
6. Variant Filtering (PASS only)
7. **ANNOVAR Annotation** (5 databases)
   - refGene
   - ClinVar
   - gnomAD
   - dbSNP (avsnp150)
   - Prediction scores (dbnsfp42a)
8. **snpEff Annotation**
9. **Zygosity Classification** (Heterozygous/Homozygous)
10. **Functional Separation**
    - By Type: SNPs, Insertions, Deletions
    - By Location: Exonic, Non-Exonic
    - By Effect: Nonsynonymous, Synonymous, Stopgain, Frameshift
11. **VCF Compression & Indexing** (for IGV)
12. **ZIP Archive** with all essential results

---

## 🚀 Features

✅ **One-Command Installation** - Everything installed automatically  
✅ **One-Command Analysis** - From FASTQ to results in one run  
✅ **Dual Annotation** - Both ANNOVAR and snpEff  
✅ **Smart Filtering** - PASS variants only, reduces file sizes by 90%  
✅ **Functional Classification** - Automatic separation by variant type and effect  
✅ **Zygosity Information** - Het/Hom status for each variant  
✅ **IGV-Ready** - Compressed, indexed VCFs for genome browser  
✅ **Excel-Ready** - Tab-delimited TXT files for easy analysis  
✅ **Production-Ready** - Used in clinical research labs  

---

## 📋 System Requirements

### Minimum
- **OS**: Linux (Ubuntu 20.04+, CentOS 7+, or similar)
- **RAM**: 16 GB minimum, 32 GB recommended
- **Storage**: 100 GB free space
- **CPU**: 8+ cores (16+ recommended)
- **Internet**: Required for initial setup

### Software (Auto-Installed)
- Python 3.7+
- Java 21 (for GATK, snpEff)
- Perl (for ANNOVAR)
- All bioinformatics tools installed by `install_all.sh`

---

## 🛠️ Installation

### Step 1: Clone Repository

```bash
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS
cd ~/NGS
```

### Step 2: Run Installation Script

```bash
bash install_all.sh
```

**This installs:**
- FastQC, fastp, BWA, SAMtools, GATK, snpEff
- Reference genome (hg19) with indices
- ANNOVAR with 9 clinical databases
- All dependencies

**Time:** 45-60 minutes  
**Storage:** ~15 GB

### Step 3: Verify Installation

```bash
bash CHECK_INSTALLATION.sh
```

Confirms all tools and databases are ready.

---

## 🎯 Usage

### Standard Workflow

   ```bash
# 1. Add FASTQ files to data directory
mkdir -p ~/NGS/data
# Copy your *_R1.fastq.gz and *_R2.fastq.gz files here

# 2. Run complete pipeline
cd ~/NGS
bash ULTIMATE_MASTER_PIPELINE.sh ~/NGS/data 16
```

**Parameters:**
- First argument: Directory containing FASTQ files
- Second argument: Number of CPU threads (default: 16)

**The pipeline will:**
- Auto-detect all sample pairs
- Process each sample sequentially
- Generate comprehensive results
- Create final ZIP archive

---

## 📦 Output Structure

### Final Deliverable

**One ZIP File:** `NGS_Results_Complete_[timestamp].zip` (~400-600 MB)

**Contains:**

```
results/
├── Sample1/
│   ├── annovar/
│   │   ├── annotated_Sample1.hg19_multianno.txt          # ANNOVAR annotation
│   │   ├── annotated_Sample1_with_zygosity.txt           # With Het/Hom status
│   │   ├── annotated_Sample1.hg19_multianno.vcf.gz       # ANNOVAR VCF
│   │   ├── separated_by_type/
│   │   │   ├── SNPs.txt
│   │   │   ├── Insertions.txt
│   │   │   └── Deletions.txt
│   │   ├── functional_classification/
│   │   │   ├── SNPs_Exonic.txt
│   │   │   ├── SNPs_NonExonic.txt
│   │   │   ├── Exonic_Nonsynonymous.txt
│   │   │   ├── Exonic_Synonymous.txt
│   │   │   ├── Exonic_Stopgain.txt
│   │   │   └── Exonic_Frameshift.txt
│   │   └── snpeff/
│   │       ├── Sample1_snpEff_annotated.vcf.gz          # snpEff VCF
│   │       ├── Sample1_snpEff_summary.html              # snpEff report
│   │       └── Sample1_snpEff_summary.csv
│   ├── fastqc/
│   │   ├── Sample1_R1_fastqc.html                       # QC reports
│   │   └── Sample1_R2_fastqc.html
│   └── trimmed/
│       └── fastp_report.html                            # Trimming stats
└── MASTER_ANALYSIS_SUMMARY.txt                          # Overall summary
```

**Excluded from ZIP** (to keep size small):
- ❌ BAM files (large, 15-20 GB each)
- ❌ FASTQ files (large, input data)
- ❌ Intermediate files

---

## 📊 Understanding the Results

### Priority Files for Analysis

1. **`annotated_[SAMPLE]_with_zygosity.txt`**
   - Main annotation file with zygosity
   - Open in Excel/LibreOffice
   - Filter by:
     - ClinVar significance (Pathogenic/Likely Pathogenic)
     - gnomAD frequency (rare: < 0.01)
     - Functional effect (nonsynonymous, stopgain)
     - Zygosity (Heterozygous/Homozygous)

2. **`functional_classification/Exonic_Nonsynonymous.txt`**
   - Coding variants that change amino acids
   - High priority for pathogenicity analysis

3. **`functional_classification/Exonic_Stopgain.txt`**
   - Variants creating premature stop codons
   - Often deleterious

4. **`snpeff/[SAMPLE]_snpEff_summary.html`**
   - Visual summary of variant effects
   - Statistics and charts

5. **VCF Files** (`.vcf.gz`)
   - Load into IGV (Integrative Genomics Viewer)
   - Visualize variants in genomic context

---

## 🔬 Pipeline Details

### Annotation Databases

**ANNOVAR (5 databases):**
- **refGene**: Gene-based annotation
- **ClinVar**: Clinical significance
- **gnomAD**: Population frequencies (exomes)
- **dbSNP (avsnp150)**: rsIDs and allele frequencies
- **dbnsfp42a**: Prediction scores (SIFT, PolyPhen, CADD, etc.)

**snpEff:**
- Functional consequences (HIGH/MODERATE/LOW impact)
- Protein change predictions
- Transcript-level annotations

### Functional Classifications

**By Variant Type:**
- SNPs (Single Nucleotide Polymorphisms)
- Insertions
- Deletions
- MNPs (Multiple Nucleotide Polymorphisms)

**By Location:**
- Exonic (in coding regions) vs Non-Exonic
- UTR, intronic, intergenic

**By Effect:**
- Nonsynonymous (amino acid change)
- Synonymous (silent mutation)
- Stopgain (premature termination)
- Frameshift (indel changing reading frame)

---

## ⏱️ Performance

**Per Sample (Typical Exome):**
- Raw FASTQ: 6-8 GB
- Analysis time: 2.5-3 hours (16 CPUs)
- Storage used: 20-25 GB
- Final ZIP: 150-200 MB

**3 Samples:**
- Total time: 7-9 hours
- Runs sequentially (one completes → next starts)

---

## 🐛 Troubleshooting

### Installation Issues

```bash
# Verify installation
bash CHECK_INSTALLATION.sh

# Check disk space
df -h ~

# Check Java version (needs 21+)
java -version
```

### Common Errors

| Error | Solution |
|-------|----------|
| `command not found` | Re-run `install_all.sh` or check PATH |
| `Out of memory` | Reduce thread count or add more RAM |
| `No FASTQ pairs found` | Check file naming: `*_R1.fastq.gz` + `*_R2.fastq.gz` |
| `snpEff class version error` | Install Java 21: `sudo apt install openjdk-21-jdk` |

### Getting Logs

```bash
# Pipeline logs for each sample
cat ~/NGS/results/[SAMPLE]/pipeline.log

# Summary of all samples
cat ~/NGS/MASTER_ANALYSIS_SUMMARY.txt
```

---

## 📚 Documentation

- **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** - Detailed project layout
- **[RUN_ULTIMATE_PIPELINE.txt](RUN_ULTIMATE_PIPELINE.txt)** - Quick reference guide
- **[LICENSE](LICENSE)** - MIT License

---

## 🤝 Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

For major changes, open an issue first.

---

## 🙏 Acknowledgments

This pipeline integrates these excellent tools:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control
- [fastp](https://github.com/OpenGene/fastp) - Read preprocessing
- [BWA](https://github.com/lh3/bwa) - Read alignment
- [SAMtools](http://www.htslib.org/) - BAM/SAM manipulation
- [GATK](https://gatk.broadinstitute.org/) - Variant calling
- [ANNOVAR](http://annovar.openbioinformatics.org/) - Variant annotation
- [snpEff](http://pcingola.github.io/SnpEff/) - Functional annotation

---

## 📝 License

MIT License - see [LICENSE](LICENSE) file for details.

---

## 📧 Support

- **Issues**: [GitHub Issues](https://github.com/Babajan-B/Exome-Analysis-End-to-END/issues)
- **Documentation**: Check the docs/ folder
- **Updates**: Watch the repository for new releases

---

## 🌟 Citation

If you use this pipeline in your research:

```
NGS Exome Analysis Pipeline
https://github.com/Babajan-B/Exome-Analysis-End-to-END
```

---

**Three scripts. One command. Complete analysis.**

🧬 Made for researchers, by researchers.

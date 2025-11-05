# Ultimate Quick Start - Complete Exome Analysis

## 🎯 **ONE Command Does Everything**

Put your FASTQ files in `~/NGS/data/` then run:

```bash
bash MASTER_PIPELINE.sh
```

**That's it!** 🎉

---

## 📋 **Complete Setup + Run**

### Step 1: Upload Your FASTQ Files

```bash
mkdir -p ~/NGS/data
# Upload your files to ~/NGS/data/:
#   sample1_R1.fastq.gz, sample1_R2.fastq.gz
#   sample2_R1.fastq.gz, sample2_R2.fastq.gz
```

### Step 2: Run Master Pipeline

```bash
cd ~/NGS
bash MASTER_PIPELINE.sh
```

**What happens:**
1. ✅ Auto-detects all FASTQ pairs
2. ✅ Shows you what it found
3. ✅ Asks for confirmation
4. ✅ Runs complete analysis for all samples:
   - Quality Control
   - Trimming
   - Alignment
   - Variant Calling
   - Filtering
   - **ANNOVAR Annotation** (5 databases)
   - **Variant Type Separation** (SNPs, Indels)
5. ✅ Generates summary report

---

## ⏱️ **Time Estimates**

| Samples | Instance | Time |
|---------|----------|------|
| 1 sample | 16 vCPU | 1.5 hours |
| 2 samples | 16 vCPU | 3 hours |
| 1 sample | 32 vCPU | 45 min |

---

## 📊 **What You Get**

For each sample:

```
results/SAMPLE_NAME/
├── fastqc/                           (QC reports)
├── trimmed/                          (Trimming reports)
├── dedup/dedup.bam                   (Aligned reads)
├── filtered/filtered_variants.vcf    (Filtered variants)
└── annovar/
    ├── annotated_SAMPLE.hg19_multianno.txt        (Excel file) ⭐
    ├── annotated_SAMPLE.hg19_multianno.vcf.gz     (VCF file)
    └── separated_by_type/
        ├── SNPs.txt                   (SNP variants) ⭐
        ├── Insertions.txt             (Insertion variants) ⭐
        ├── Deletions.txt              (Deletion variants) ⭐
        └── MNPs.txt                   (Complex variants)
```

---

## 🗂️ **File Naming Patterns Supported**

The script auto-detects these patterns:

✅ `sample_R1.fastq.gz` + `sample_R2.fastq.gz`  
✅ `sample_1.fastq.gz` + `sample_2.fastq.gz`  
✅ `sampleR1.fastq.gz` + `sampleR2.fastq.gz`  
✅ `sample_R1.fq.gz` + `sample_R2.fq.gz`  

---

## 💡 **Example Workflow**

```bash
# Upload files to data folder
mkdir -p ~/NGS/data
# (Upload patient1_R1.fastq.gz, patient1_R2.fastq.gz via JupyterLab)

# Run master pipeline
cd ~/NGS
bash MASTER_PIPELINE.sh

# Script detects: patient1
# Asks: "Start analysis? (y/n)"
# You type: y
# 
# Pipeline runs automatically (~1.5 hours)
# 
# Results in:
#   results/patient1/annovar/separated_by_type/SNPs.txt
#   results/patient1/annovar/separated_by_type/Insertions.txt
#   results/patient1/annovar/separated_by_type/Deletions.txt
```

---

## 📥 **Download Results**

After completion:

```bash
# Check results
ls -lh ~/NGS/results/*/annovar/separated_by_type/

# Download via JupyterLab:
# Navigate to results/SAMPLE/annovar/separated_by_type/
# Right-click → Download
```

---

## 🔧 **Advanced Options**

### Custom Threads

```bash
bash MASTER_PIPELINE.sh ~/NGS/data 32
```

### Custom Data Directory

```bash
bash MASTER_PIPELINE.sh /path/to/fastq/files 16
```

---

## ✅ **Summary**

**One script does:**
1. ✅ Auto-detects samples
2. ✅ Complete pipeline (FASTQ → VCF)
3. ✅ ANNOVAR annotation (5 databases)
4. ✅ Variant type separation
5. ✅ Summary report

**Just run:**
```bash
bash MASTER_PIPELINE.sh
```

**That's it!** 🚀

---

**Made with 🧬 for fully automated exome analysis**


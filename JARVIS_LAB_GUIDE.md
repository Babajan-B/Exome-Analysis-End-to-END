# Jarvis Lab Cloud Instance Setup Guide

Complete guide to run NGS Exome Analysis Pipeline on Jarvis Lab cloud instances with maximum speed.

---

## 🎯 **QUICK REFERENCE - Copy & Paste This!**

### **In Jarvis Lab:**
1. **Select Template**: **PyTorch** (only suitable option)
2. **Instance Size**: 16 vCPU, 32 GB RAM, 100 GB Storage
3. **Launch & Open Terminal**

### **Command 1 - Setup (Run once, ~40 min):**
```bash
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS && cd ~/NGS && bash cloud_setup.sh
```

### **Command 2 - Upload Your Files:**
```bash
mkdir -p ~/NGS/data
# Upload R1.fastq.gz and R2.fastq.gz to ~/NGS/data/ using JupyterLab upload
```

### **Command 3 - Run Pipeline (~1-2 hours):**
```bash
bash run_pipeline.sh data/R1.fastq.gz data/R2.fastq.gz my_sample 16
```

### **Download Results:**
- Navigate to: `~/NGS/results/my_sample/annotated/`
- Download: `annotated_variants.vcf.gz`

**That's it! Just 3 commands!** ⚡

---

## 🚀 Jarvis Lab Instance Selection

### ⚠️ **Important: No Blank Ubuntu Template Available**

Since Jarvis Lab doesn't offer a blank Ubuntu option, we'll use the **PyTorch template** as the base.

### **SELECT: PyTorch Template** ✅

**Why PyTorch?**
- ✅ Has Ubuntu 20.04/22.04 base
- ✅ Pre-installed build tools (gcc, make)
- ✅ Python already configured
- ✅ We'll install bioinformatics tools on top
- ⚠️ Ignore all the ML/GPU stuff (we don't need it)

---

## 💻 **Recommended Instance Configurations**

### **Option 1: Standard Analysis (RECOMMENDED)** ⭐
- **Template**: **PyTorch**
- **Instance Size**: 16 vCPU, 32 GB RAM
- **Storage**: 100-150 GB SSD
- **GPU**: **Select CPU-only if available** (or smallest GPU if forced)
- **Estimated Cost**: ~$0.80-1.20/hour
- **Estimated Time**: **1-1.5 hours** for full exome

### **Option 2: Fast Processing** ⚡
- **Template**: **PyTorch**
- **Instance Size**: 32 vCPU, 64 GB RAM
- **Storage**: 150 GB SSD
- **GPU**: CPU-only preferred
- **Estimated Cost**: ~$1.50-2.00/hour
- **Estimated Time**: **45-60 minutes** for full exome

### **Option 3: Budget-Friendly** 💰
- **Template**: **PyTorch**
- **Instance Size**: 8 vCPU, 16 GB RAM
- **Storage**: 100 GB SSD
- **GPU**: CPU-only
- **Estimated Cost**: ~$0.40-0.60/hour
- **Estimated Time**: **2-3 hours** for full exome

---

## 📋 Quick Start - 3 Steps Only!

### Step 1: Launch Jarvis Lab Instance with PyTorch Template

1. **Go to**: [Jarvis Labs](https://jarvislabs.ai/) and sign in
2. **Click**: "Create New Instance" or "Launch Instance"
3. **Select Template**: **PyTorch** ⭐
4. **Configure**:
   - **vCPUs**: 16 (recommended) or 32 (faster)
   - **RAM**: 32 GB minimum
   - **Storage**: 100-150 GB SSD
   - **GPU**: Select **CPU-only** if available, or **smallest/cheapest GPU** if forced
5. **Launch** the instance
6. **Connect**: Click "Open JupyterLab" or use SSH
7. **Open Terminal**: In JupyterLab, click "Terminal" icon

**Note**: We're using PyTorch template only for the Ubuntu base - we'll install bioinformatics tools ourselves.

---

### Step 2: Clone Repository & Setup (Single Command)

```bash
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS && \
cd ~/NGS && \
bash cloud_setup.sh
```

**This single command will**:
- ✅ Clone the repository
- ✅ Install all system dependencies
- ✅ Install bioinformatics tools (FastQC, BWA, SAMtools, GATK, etc.)
- ✅ Download reference genome (hg19, ~3GB)
- ✅ Index reference genome
- ✅ Set up snpEff for annotation
- ⏱️ **Time**: 30-45 minutes (mostly downloading reference genome)

---

### Step 3: Run Pipeline (Single Command)

**Upload your FASTQ files first**, then run:

```bash
cd ~/NGS && \
bash run_pipeline.sh \
    /path/to/your/R1.fastq.gz \
    /path/to/your/R2.fastq.gz \
    my_sample_name \
    16
```

**Parameters**:
- `R1.fastq.gz` - Path to your R1 (forward) reads
- `R2.fastq.gz` - Path to your R2 (reverse) reads
- `my_sample_name` - Name for your analysis
- `16` - Number of threads (match your instance vCPUs)

---

## 🎯 Complete Example Workflow

```bash
# 1. Clone and setup (run once per instance)
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS
cd ~/NGS
bash cloud_setup.sh

# 2. Create data directory and upload your files
mkdir -p ~/NGS/data
# Upload your FASTQ files to ~/NGS/data/

# 3. Run the pipeline
bash run_pipeline.sh \
    data/sample_R1.fastq.gz \
    data/sample_R2.fastq.gz \
    patient_001 \
    16

# 4. Results will be in ~/NGS/results/patient_001/
```

---

## ⚡ Performance Optimization Tips

### For Maximum Speed:

1. **Use Maximum Threads**:
   ```bash
   # Find number of CPUs
   nproc
   
   # Use all available cores
   bash run_pipeline.sh data/R1.fq.gz data/R2.fq.gz sample $(nproc)
   ```

2. **Use tmpfs for Temporary Files** (if you have lots of RAM):
   ```bash
   # Create RAM disk for temp files
   sudo mkdir -p /mnt/ramdisk
   sudo mount -t tmpfs -o size=32G tmpfs /mnt/ramdisk
   export TMPDIR=/mnt/ramdisk
   ```

3. **Optimize for NVMe/SSD Storage**:
   ```bash
   # The script automatically uses SSD if available
   # Store data on fastest disk
   ```

4. **Skip Optional Steps**:
   Edit `run_pipeline.sh` to comment out:
   - FastQC (saves 3-5 min)
   - BQSR if no known sites (saves 10-15 min)

---

## 📊 Expected Timeline by Instance Type

### 8 vCPU Instance (cpu.4xlarge)
| Step | Time | Total |
|------|------|-------|
| FastQC | 5 min | 5 min |
| Trimming | 10 min | 15 min |
| Alignment | 60 min | 75 min |
| SAM→BAM+Sort | 15 min | 90 min |
| Mark Duplicates | 10 min | 100 min |
| Variant Calling | 45 min | 145 min |
| Filtering | 2 min | 147 min |
| Annotation | 5 min | 152 min |
| **Total** | **~2.5 hours** | |

### 16 vCPU Instance (cpu.8xlarge) - **RECOMMENDED**
| Step | Time | Total |
|------|------|-------|
| FastQC | 3 min | 3 min |
| Trimming | 5 min | 8 min |
| Alignment | 30 min | 38 min |
| SAM→BAM+Sort | 8 min | 46 min |
| Mark Duplicates | 5 min | 51 min |
| Variant Calling | 25 min | 76 min |
| Filtering | 1 min | 77 min |
| Annotation | 3 min | 80 min |
| **Total** | **~1.5 hours** | |

### 32 vCPU Instance (cpu.16xlarge)
| Step | Time | Total |
|------|------|-------|
| FastQC | 2 min | 2 min |
| Trimming | 3 min | 5 min |
| Alignment | 15 min | 20 min |
| SAM→BAM+Sort | 5 min | 25 min |
| Mark Duplicates | 3 min | 28 min |
| Variant Calling | 15 min | 43 min |
| Filtering | 1 min | 44 min |
| Annotation | 2 min | 46 min |
| **Total** | **~45 min** | |

---

## 💰 Cost Estimation

**For typical 6GB exome (R1 + R2 = 6GB total)**:

| Instance | Time | Cost/hr | Total Cost |
|----------|------|---------|------------|
| 8 vCPU | 2.5 hrs | $0.40 | **$1.00** |
| 16 vCPU | 1.5 hrs | $0.80 | **$1.20** |
| 32 vCPU | 0.75 hrs | $1.60 | **$1.20** |

**Recommendation**: **16 vCPU** instance offers best price/performance ratio!

---

## 📥 Uploading Data to Jarvis Lab

### Option 1: Using SCP (from your local machine)
```bash
# Upload FASTQ files
scp sample_R1.fastq.gz ubuntu@your-instance-ip:~/NGS/data/
scp sample_R2.fastq.gz ubuntu@your-instance-ip:~/NGS/data/
```

### Option 2: Using wget (if files are on a server)
```bash
# In Jarvis Lab terminal
cd ~/NGS/data
wget https://your-server.com/sample_R1.fastq.gz
wget https://your-server.com/sample_R2.fastq.gz
```

### Option 3: Using JupyterLab Upload
1. Open JupyterLab interface
2. Navigate to `~/NGS/data/`
3. Use upload button to transfer files

---

## 📥 Downloading Results from Jarvis Lab

### Option 1: Using SCP (to your local machine)
```bash
# Download all results
scp -r ubuntu@your-instance-ip:~/NGS/results/my_sample_name ./local_results/

# Download just the final VCF
scp ubuntu@your-instance-ip:~/NGS/results/my_sample_name/annotated/annotated_variants.vcf.gz ./
```

### Option 2: Using JupyterLab Download
1. Navigate to `~/NGS/results/` in JupyterLab
2. Right-click on files/folders
3. Select "Download"

---

## 🔧 Troubleshooting

### Issue: Out of Memory
```bash
# Monitor memory usage
watch -n 5 free -h

# If running out of memory, reduce threads
bash run_pipeline.sh data/R1.fq.gz data/R2.fq.gz sample 4
```

### Issue: Out of Disk Space
```bash
# Check disk space
df -h

# Clean up intermediate files manually in the script
# Or add --delete-intermediates flag (custom)
```

### Issue: Pipeline Fails Mid-Run
```bash
# Check the log
tail -50 ~/NGS/results/your_sample/pipeline.log

# Resume from specific step by editing run_pipeline.sh
# Comment out completed steps and re-run
```

---

## 🎯 Production-Ready Configuration

For running multiple samples efficiently:

```bash
#!/bin/bash
# Process multiple samples in parallel

SAMPLES=(
    "patient001:data/patient001_R1.fq.gz:data/patient001_R2.fq.gz"
    "patient002:data/patient002_R1.fq.gz:data/patient002_R2.fq.gz"
    "patient003:data/patient003_R1.fq.gz:data/patient003_R2.fq.gz"
)

# Process in parallel (one per 16 cores)
for sample in "${SAMPLES[@]}"; do
    IFS=':' read -r name r1 r2 <<< "$sample"
    bash run_pipeline.sh "$r1" "$r2" "$name" 16 &
done

wait
echo "All samples processed!"
```

---

## 📊 Monitor Progress

### Real-time monitoring:
```bash
# Watch log file
tail -f ~/NGS/results/your_sample/pipeline.log

# Check CPU and memory usage
htop

# Check disk usage
watch -n 10 df -h

# Check specific process
ps aux | grep -E "(bwa|samtools|gatk)"
```

---

## 🔐 Security Notes

1. **Don't leave instances running** - Stop when done
2. **Download results immediately** - Instances are temporary
3. **Delete sensitive data** after download
4. **Use encrypted transfer** (SCP with SSH keys)

---

## 📦 What Gets Installed

### System Packages:
- build-essential, wget, curl, git, unzip
- Java JDK 11
- Python 3

### Bioinformatics Tools:
- **FastQC** (v0.11.9+) - Quality control
- **fastp** (latest) - Read preprocessing
- **BWA** (v0.7.17+) - Read alignment
- **SAMtools** (v1.15+) - SAM/BAM manipulation
- **GATK** (v4.6.2.0) - Variant calling
- **snpEff** (latest) - Variant annotation
- **BCFtools** (v1.15+) - VCF manipulation
- **Tabix** (v1.15+) - VCF indexing

### Reference Data:
- **hg19** reference genome (~3.2 GB)
- **BWA indices** (~3.2 GB)
- **SAMtools index** (~5 MB)
- **GATK dictionary** (~5 MB)
- **snpEff database** GRCh37.75 (~700 MB)

**Total storage**: ~10 GB for tools + 20-50 GB per sample analysis

---

## 🎓 Advanced Usage

### Run with Custom Parameters

Edit `run_pipeline.sh` to customize:

```bash
# Alignment parameters
bwa mem -t 32 -k 19 -M -R "@RG..." ...

# Variant calling parameters
gatk HaplotypeCaller \
    --native-pair-hmm-threads 32 \
    --max-alternate-alleles 3 \
    --min-base-quality-score 20 \
    ...

# Filtering parameters
gatk VariantFiltration \
    --filter-expression "QD < 5.0" \
    ...
```

---

## 📞 Support

If you encounter issues:
1. Check `~/NGS/results/your_sample/pipeline.log`
2. Verify all tools installed: `which fastqc bwa samtools gatk`
3. Check disk space: `df -h`
4. Open issue on GitHub

---

## ✅ Summary

**Two Commands to Success:**

```bash
# 1. Setup (once per instance)
git clone https://github.com/Babajan-B/Exome-Analysis-End-to-END.git ~/NGS && \
cd ~/NGS && bash cloud_setup.sh

# 2. Run analysis (per sample)
bash run_pipeline.sh data/R1.fq.gz data/R2.fq.gz sample_name 16
```

**That's it!** 🎉

---

**Made with ⚡ for high-performance cloud computing**


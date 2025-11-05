# Smart Workflow Guide - Parallel ANNOVAR Installation

## 🚀 Problem Solved

ANNOVAR installation takes 20-30 minutes. Instead of waiting, this smart workflow:
1. **Starts pipeline immediately** (doesn't wait for ANNOVAR)
2. **Installs ANNOVAR in background** during pipeline run  
3. **Auto-annotates when both finish**
4. **Gives options if ANNOVAR still installing**

## ⚡ Quick Start

### Option 1: Smart Pipeline (Recommended)

```bash
cd ~/NGS
bash run_complete_analysis_smart.sh \
    data/R1.fastq.gz \
    data/R2.fastq.gz \
    patient_001 \
    16
```

**What happens:**
- ✅ Pipeline starts immediately
- ✅ ANNOVAR installs in background (if not already installed)
- ✅ When pipeline finishes, checks ANNOVAR status:
  - If ready → Auto-annotates
  - If installing → Gives you options:
    1. Wait for installation (recommended)
    2. Annotate later manually
    3. Skip annotation

### Option 2: Install ANNOVAR First (Traditional)

```bash
# 1. Install everything including ANNOVAR
bash install_all.sh

# 2. Run analysis (ANNOVAR already ready)
bash run_complete_analysis.sh data/R1.fastq.gz data/R2.fastq.gz patient_001 16
```

## 📋 Three Main Scripts

### 1. `run_complete_analysis_smart.sh` ⭐ SMART
- Runs pipeline immediately
- Handles ANNOVAR installation in parallel
- Auto-detects ANNOVAR status
- Gives intelligent options

### 2. `run_complete_analysis.sh` (Standard)
- Assumes ANNOVAR is already installed
- Simpler, faster if ANNOVAR ready
- Good for repeat analyses

### 3. `install_all.sh` (Setup)
- Installs everything upfront
- Takes 45-60 minutes
- Good for first-time setup

## 🔄 Workflow Scenarios

### Scenario A: First Time (No ANNOVAR)

```bash
cd ~/NGS
bash run_complete_analysis_smart.sh data/R1.fq.gz data/R2.fq.gz sample1 16
```

**Timeline:**
```
Time 0:     Pipeline starts + ANNOVAR installation starts in background
Time 60:    Pipeline completes
Time 60:    Check ANNOVAR status
  - If ready: Auto-annotate (5 min)
  - If not: Choose to wait or annotate later
Time 90:    ANNOVAR installation completes
```

### Scenario B: ANNOVAR Already Installed

```bash
cd ~/NGS
bash run_complete_analysis_smart.sh data/R1.fq.gz data/R2.fq.gz sample2 16
```

**Timeline:**
```
Time 0:     Detects ANNOVAR installed
Time 0:     Pipeline starts immediately
Time 60:    Pipeline completes
Time 60:    Auto-annotates (5 min)
Time 65:    Complete with annotation
```

### Scenario C: ANNOVAR Installing (Another Analysis Running)

```bash
cd ~/NGS
bash run_complete_analysis_smart.sh data/R1.fq.gz data/R2.fq.gz sample3 16
```

**Timeline:**
```
Time 0:     Detects ANNOVAR installation in progress
Time 0:     Pipeline starts (doesn't start duplicate installation)
Time 60:    Pipeline completes
Time 60:    Check ANNOVAR - still installing
Time 60:    Prompts for action:
              1. Wait (recommended)
              2. Annotate later
              3. Skip
```

## 🛠️ Additional Tools

### Annotate Existing Results

If you chose to annotate later:

```bash
bash annotate_results.sh sample_name
```

This will:
- Check if ANNOVAR is installed
- Find the sample's filtered VCF
- Run ANNOVAR annotation
- Compress and index results

### Background ANNOVAR Installation

Manual trigger (usually automatic):

```bash
nohup bash install_annovar_background.sh > annovar_install.log 2>&1 &
```

### Check ANNOVAR Status

```bash
# Check if installed
ls -la ~/NGS/.annovar_installed

# Check if installing
ls -la ~/NGS/.annovar_installing

# Check ANNOVAR directory
ls -la ~/NGS/tools/annovar/
```

## 📊 Timing Comparison

### Traditional (Sequential)
```
Install ANNOVAR:  30 min
Run Pipeline:     60 min
Annotate:         5 min
────────────────────────
TOTAL:            95 min
```

### Smart (Parallel)
```
Run Pipeline:     60 min (ANNOVAR installing in background)
Wait for ANNOVAR: 0 min (already done!)
Annotate:         5 min
────────────────────────
TOTAL:            65 min
```

**Time Saved: 30 minutes!** ⚡

## 🔍 Status Files

The smart system uses lock files:

- `~/NGS/.annovar_installed` - ANNOVAR ready to use
- `~/NGS/.annovar_installing` - Installation in progress
- `~/NGS/results/SAMPLE/annovar_install.log` - Installation log

## 💡 Tips

1. **First run**: Use smart pipeline, grab coffee while it installs ☕
2. **Subsequent runs**: Use either smart or standard pipeline
3. **Multiple samples**: Smart pipeline handles concurrent ANNOVAR installs
4. **Check logs**: Monitor `annovar_install.log` in another terminal
5. **Annotate later**: Always possible with `annotate_results.sh`

## ⚠️ Important Notes

- Only ONE ANNOVAR installation runs at a time (uses lock file)
- Pipeline NEVER waits for ANNOVAR (runs immediately)
- Safe to run multiple analyses simultaneously
- ANNOVAR installation ~5GB, takes 20-30 min
- Annotation per sample ~5-10 min

## 🎯 Recommended Workflow

### For Production (Multiple Samples)

```bash
# First sample (installs ANNOVAR in background)
bash run_complete_analysis_smart.sh data/sample1_R1.fq.gz data/sample1_R2.fq.gz sample1 16

# While sample1 running, start sample2 in another terminal
bash run_complete_analysis_smart.sh data/sample2_R1.fq.gz data/sample2_R2.fq.gz sample2 16

# Both pipelines run immediately
# ANNOVAR installation happens once in background
# Both auto-annotate when ready
```

## 📖 Quick Reference

```bash
# Smart run (handles everything)
bash run_complete_analysis_smart.sh R1.fq.gz R2.fq.gz sample 16

# Annotate later
bash annotate_results.sh sample

# Check ANNOVAR status
ls -la ~/NGS/.annovar_*

# Manual ANNOVAR install
bash setup_annovar.sh

# Standard run (ANNOVAR must be ready)
bash run_complete_analysis.sh R1.fq.gz R2.fq.gz sample 16
```

---

**Made with ⚡ for parallel processing efficiency**


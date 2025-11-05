# ANNOVAR Installation Summary for Jarvis Lab

## 📦 What Was Created

The following files have been added to your NGS pipeline for ANNOVAR integration:

### 🔧 Installation Scripts
1. **`setup_annovar.sh`** - Main ANNOVAR installation script
2. **`cloud_setup_with_annovar.sh`** - Extended cloud setup including ANNOVAR
3. **`test_annovar.sh`** - Verify ANNOVAR installation

### 🚀 Helper Scripts
4. **`annovar_helper.sh`** - Quick annotation wrapper script

### 📚 Documentation
5. **`ANNOVAR_QUICKSTART.md`** - Quick reference guide
6. **`docs/ANNOVAR_INTEGRATION.md`** - Complete integration guide
7. **`docs/ANNOTATION_TOOLS_COMPARISON.md`** - Comparison of annotation tools
8. **`ANNOVAR_INSTALLATION_SUMMARY.md`** - This file

---

## 🎯 Quick Start Commands for Jarvis Lab

### Step 1: Install ANNOVAR (One-Time Setup)

```bash
cd ~/NGS
bash setup_annovar.sh
```

**What this does:**
- Downloads ANNOVAR from the official source
- Installs to `~/NGS/tools/annovar/`
- Downloads essential databases for hg19:
  - refGene (RefSeq genes)
  - knownGene (UCSC genes)
  - ensGene (Ensembl genes)
  - avsnp150 (dbSNP)
  - gnomad211_exome (population frequencies)
  - clinvar_20240917 (pathogenic variants)
  - dbnsfp42a (functional predictions)
  - cosmic70 (cancer mutations)
  - icgc28 (cancer somatic)

**Time**: 20-30 minutes  
**Storage**: ~5 GB

---

### Step 2: Test Installation

```bash
bash test_annovar.sh
```

**Expected output:**
```
✅ ANNOVAR directory found
✅ Scripts are functional
✅ Databases are present
✅ Test annotation successful
```

---

### Step 3: Use ANNOVAR

#### Option A: Annotate Existing VCF

```bash
bash annovar_helper.sh input.vcf output_prefix
```

#### Option B: Annotate Pipeline Results

```bash
# After running: bash run_pipeline.sh data/R1.fastq.gz data/R2.fastq.gz sample_name 16
bash annovar_helper.sh \
    results/sample_name/filtered/filtered_variants.vcf \
    results/sample_name/annovar/annotated
```

---

## 📋 Complete Workflow Example

### Full Pipeline with ANNOVAR Annotation

```bash
# 1. Navigate to project directory
cd ~/NGS

# 2. Create data directory and upload FASTQ files
mkdir -p data
# Upload your R1.fastq.gz and R2.fastq.gz to data/ folder

# 3. Run the main NGS pipeline
bash run_pipeline.sh \
    data/R1.fastq.gz \
    data/R2.fastq.gz \
    patient_001 \
    16

# 4. Wait for pipeline to complete (~1-2 hours on 16 vCPU)

# 5. Annotate variants with ANNOVAR (~5-10 minutes)
bash annovar_helper.sh \
    results/patient_001/filtered/filtered_variants.vcf \
    results/patient_001/annovar/patient_001_annotated

# 6. Results will be in:
# - results/patient_001/annovar/patient_001_annotated.hg19_multianno.vcf
# - results/patient_001/annovar/patient_001_annotated.hg19_multianno.txt
```

---

## 🔗 File Locations After Installation

```
~/NGS/
├── tools/
│   └── annovar/                           # ANNOVAR installation
│       ├── annotate_variation.pl          # Main annotation script
│       ├── table_annovar.pl               # Table annotation script
│       ├── convert2annovar.pl             # Format converter
│       └── humandb/                       # Database directory
│           ├── hg19_refGene.txt           # RefSeq genes
│           ├── hg19_avsnp150.txt          # dbSNP
│           ├── hg19_gnomad211_exome.txt   # gnomAD
│           ├── hg19_clinvar_20240917.txt  # ClinVar
│           └── ... (more databases)
│
├── setup_annovar.sh                       # Installation script
├── annovar_helper.sh                      # Quick annotation script
├── test_annovar.sh                        # Test script
├── ANNOVAR_QUICKSTART.md                  # Quick reference
└── docs/
    ├── ANNOVAR_INTEGRATION.md             # Complete guide
    └── ANNOTATION_TOOLS_COMPARISON.md     # Tool comparison
```

---

## 💡 Common Usage Patterns

### 1. Quick Clinical Annotation

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out annotated \
    -protocol refGene,clinvar_20240917,gnomad211_exome \
    -operation g,f,f \
    -vcfinput
```

### 2. Comprehensive Research Annotation

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out annotated_full \
    -protocol refGene,knownGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70 \
    -operation g,g,f,f,f,f,f \
    -vcfinput \
    -polish
```

### 3. Cancer-Specific Annotation

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    somatic_variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out cancer_annotated \
    -protocol refGene,cosmic70,icgc28,clinvar_20240917 \
    -operation g,f,f,f \
    -vcfinput
```

---

## 📊 Expected Output Files

After running ANNOVAR, you'll get:

1. **`.hg19_multianno.vcf`** - Annotated VCF file
   - Same format as input VCF
   - Additional INFO fields with annotations
   - Can be used with IGV, GATK, etc.

2. **`.hg19_multianno.txt`** - Tab-delimited table
   - One variant per line
   - Easy to open in Excel
   - Easy to filter and sort

3. **`.avinput`** - Intermediate file (usually removed with `-remove` flag)

---

## 🎓 Understanding the Output

### Key Columns in the `.txt` File

| Column | Description | Example |
|--------|-------------|---------|
| `Chr` | Chromosome | chr1 |
| `Start` | Position | 12345 |
| `Ref` | Reference allele | A |
| `Alt` | Alternate allele | G |
| `Func.refGene` | Functional region | exonic |
| `Gene.refGene` | Gene name | BRCA1 |
| `ExonicFunc.refGene` | Variant type | nonsynonymous SNV |
| `AAChange.refGene` | Protein change | BRCA1:c.5266dupC:p.Q1756fs |
| `avsnp150` | dbSNP ID | rs80357906 |
| `gnomAD_exome_ALL` | Population frequency | 0.0001 |
| `CLNSIG` | Clinical significance | Pathogenic |
| `SIFT_score` | Deleteriousness | 0.001 (damaging) |
| `Polyphen2_HDIV_score` | Pathogenicity | 0.999 (damaging) |

### Interpreting Clinical Significance

**Pathogenic Variants:**
- `CLNSIG` = "Pathogenic" or "Likely_pathogenic"
- Low population frequency (< 0.01)
- Damaging predictions (SIFT < 0.05, PolyPhen > 0.85)

**Benign Variants:**
- `CLNSIG` = "Benign" or "Likely_benign"
- High population frequency (> 0.05)
- Tolerated predictions

**Variants of Uncertain Significance (VUS):**
- No ClinVar entry or `CLNSIG` = "Uncertain_significance"
- Rare but conflicting predictions

---

## 🔧 Troubleshooting

### Problem 1: "ANNOVAR not found"

**Solution:**
```bash
# Check installation
ls -la ~/NGS/tools/annovar/

# If not found, install
cd ~/NGS
bash setup_annovar.sh
```

### Problem 2: "Cannot open database"

**Solution:**
```bash
# Check databases
ls -lh ~/NGS/tools/annovar/humandb/

# Download missing database
cd ~/NGS/tools/annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
```

### Problem 3: "Out of memory"

**Solution:**
```bash
# Split large VCF by chromosome
bcftools view -r chr1 input.vcf > chr1.vcf
bash annovar_helper.sh chr1.vcf chr1_annotated

# Process each chromosome separately
```

### Problem 4: Download fails during setup

**Solution:**
```bash
# Try downloading databases one at a time
cd ~/NGS/tools/annovar

# Essential databases only
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20240917 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome humandb/
```

---

## 📥 Adding More Databases

### For hg19

```bash
cd ~/NGS/tools/annovar

# Additional population databases
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/

# Additional clinical databases
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar omim humandb/

# Conservation scores
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar phastConsElements100way humandb/
```

### For hg38 (if using hg38 reference)

```bash
cd ~/NGS/tools/annovar

# Core databases for hg38
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
```

---

## 🎯 Next Steps

1. **Test the installation**: `bash test_annovar.sh`
2. **Read the quick start**: `cat ANNOVAR_QUICKSTART.md`
3. **Run your first annotation**: `bash annovar_helper.sh test.vcf test_output`
4. **Integrate with pipeline**: Follow examples in `docs/ANNOVAR_INTEGRATION.md`

---

## 📚 Documentation Files

| File | Purpose |
|------|---------|
| `ANNOVAR_QUICKSTART.md` | Quick reference for common commands |
| `docs/ANNOVAR_INTEGRATION.md` | Complete integration guide with examples |
| `docs/ANNOTATION_TOOLS_COMPARISON.md` | Compare ANNOVAR vs snpEff vs others |
| `ANNOVAR_INSTALLATION_SUMMARY.md` | This file - installation overview |

---

## 🌐 Download Link Used

The installation script downloads ANNOVAR from:
```
http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
```

This is the personal download link you provided. If this link expires or changes, update the `ANNOVAR_URL` variable in `setup_annovar.sh`.

---

## 💻 System Requirements

| Component | Requirement |
|-----------|-------------|
| **Disk Space** | ~5 GB for ANNOVAR + databases |
| **Memory** | 2-4 GB for typical annotation |
| **Perl** | Version 5.10+ (usually pre-installed) |
| **Time** | 20-30 min setup, 5-10 min per annotation |

---

## ✅ Verification Checklist

After installation, verify:

- [ ] ANNOVAR directory exists: `ls ~/NGS/tools/annovar/`
- [ ] Scripts are executable: `ls -l ~/NGS/tools/annovar/*.pl`
- [ ] Databases downloaded: `ls ~/NGS/tools/annovar/humandb/`
- [ ] Test script passes: `bash test_annovar.sh`
- [ ] Can annotate test file: `bash annovar_helper.sh test.vcf test_output`

---

## 🆘 Getting Help

1. **Test Installation**: Run `bash test_annovar.sh`
2. **Check Documentation**: See `docs/ANNOVAR_INTEGRATION.md`
3. **ANNOVAR Website**: http://annovar.openbioinformatics.org/
4. **ANNOVAR Forum**: https://groups.google.com/g/annovar
5. **Project Issues**: Create an issue on your GitHub repository

---

**Installation Date**: $(date)  
**ANNOVAR Version**: Latest (downloaded from provided link)  
**Genome Build**: hg19 (GRCh37)

---

**Made with 🧬 for genomic variant annotation in Jarvis Lab**


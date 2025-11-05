# 🧬 ANNOVAR for NGS Pipeline - Complete Guide

**Quick Links**: [Installation](#-installation) | [Quick Start](#-quick-start) | [Usage Examples](#-usage-examples) | [Troubleshooting](#-troubleshooting)

---

## 📋 Overview

ANNOVAR (ANNOtate VARiation) is a powerful tool for functional annotation of genetic variants detected from diverse genomes. This integration allows you to annotate variants with:

- **Gene Information**: RefSeq, UCSC, Ensembl genes
- **Population Frequencies**: gnomAD, ExAC, 1000 Genomes
- **Clinical Databases**: ClinVar, COSMIC, ICGC
- **Functional Predictions**: SIFT, PolyPhen2, CADD scores
- **Conservation Scores**: PhyloP, GERP++

---

## 🚀 Installation

### For Jarvis Lab (Recommended)

```bash
cd ~/NGS
bash setup_annovar.sh
```

**Time**: 20-30 minutes  
**Storage**: ~5 GB  
**Downloads**: ANNOVAR + 9 essential databases for hg19

### Verify Installation

```bash
bash test_annovar.sh
```

Expected output: ✅ All tests pass

---

## ⚡ Quick Start

### 1. Annotate VCF File (Simplest)

```bash
bash annovar_helper.sh input.vcf output_prefix
```

### 2. Annotate Pipeline Results

```bash
bash annovar_helper.sh \
    results/sample_name/filtered/filtered_variants.vcf \
    results/sample_name/annovar/annotated
```

### 3. View Results

**VCF Format** (for tools):
```bash
less output_prefix.hg19_multianno.vcf
```

**Table Format** (for Excel):
```bash
# Download and open in Excel/Google Sheets
open output_prefix.hg19_multianno.txt
```

---

## 📚 Usage Examples

### Example 1: Clinical Annotation (Recommended)

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out clinical_annotated \
    -remove \
    -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150,dbnsfp42a \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
```

**Annotations added:**
- Gene names and transcript effects
- ClinVar pathogenicity
- Population frequencies (gnomAD)
- dbSNP IDs
- SIFT & PolyPhen predictions

### Example 2: Cancer Annotation

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    somatic_variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out cancer_annotated \
    -remove \
    -protocol refGene,cosmic70,icgc28,clinvar_20240917,gnomad211_exome \
    -operation g,f,f,f,f \
    -vcfinput
```

**Annotations added:**
- COSMIC cancer mutations
- ICGC somatic variants
- ClinVar cancer associations
- Population frequencies

### Example 3: Basic Gene Annotation (Fast)

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out basic_annotated \
    -remove \
    -protocol refGene \
    -operation g \
    -vcfinput
```

**Annotations added:**
- Gene names
- Exonic/intronic/intergenic
- Amino acid changes

### Example 4: Research (Comprehensive)

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    variants.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out research_annotated \
    -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70,icgc28 \
    -operation g,g,g,f,f,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
```

**Annotations added:**
- Multiple gene databases
- All clinical databases
- All population databases
- All prediction scores

---

## 📊 Understanding Output

### Output Files

| File | Description | Use For |
|------|-------------|---------|
| `.hg19_multianno.vcf` | Annotated VCF | IGV, GATK, downstream tools |
| `.hg19_multianno.txt` | Tab-delimited table | Excel, filtering, analysis |
| `.avinput` | Intermediate format | Debugging (removed with `-remove`) |

### Key Columns in `.txt` File

| Column | Meaning | Example |
|--------|---------|---------|
| `Func.refGene` | Functional location | exonic, intronic, UTR5 |
| `Gene.refGene` | Gene name | BRCA1 |
| `ExonicFunc.refGene` | Variant effect | nonsynonymous SNV |
| `AAChange.refGene` | Protein change | BRCA1:c.5266dupC:p.Q1756fs |
| `avsnp150` | dbSNP ID | rs80357906 |
| `gnomAD_exome_ALL` | Allele frequency | 0.0001 (0.01%) |
| `CLNSIG` | Clinical significance | Pathogenic |
| `SIFT_score` | Deleteriousness | 0.001 (damaging if < 0.05) |
| `Polyphen2_HDIV_score` | Pathogenicity | 0.999 (damaging if > 0.85) |

### Filtering for Important Variants

**Pathogenic Variants:**
```bash
# In Excel: Filter for
- CLNSIG contains "Pathogenic"
- gnomAD_exome_ALL < 0.01 (rare)
- SIFT_score < 0.05 (damaging)
```

**Novel Variants:**
```bash
# Filter for
- avsnp150 is empty (not in dbSNP)
- gnomAD_exome_ALL < 0.0001 (very rare)
```

---

## 🗄️ Available Databases

### Installed by Default (hg19)

| Database | Description | Size |
|----------|-------------|------|
| `refGene` | RefSeq genes | ~50 MB |
| `knownGene` | UCSC genes | ~60 MB |
| `ensGene` | Ensembl genes | ~50 MB |
| `avsnp150` | dbSNP 150 | ~1 GB |
| `gnomad211_exome` | gnomAD v2.1.1 | ~500 MB |
| `clinvar_20240917` | ClinVar | ~100 MB |
| `dbnsfp42a` | Functional predictions | ~3 GB |
| `cosmic70` | COSMIC mutations | ~200 MB |
| `icgc28` | ICGC somatic | ~150 MB |

### Download Additional Databases

```bash
cd ~/NGS/tools/annovar

# Download a database
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar DATABASE_NAME humandb/

# Examples:
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
```

---

## 🔄 Complete Workflow

### Full Pipeline + ANNOVAR

```bash
# Step 1: Setup (once per instance)
cd ~/NGS
bash setup_annovar.sh

# Step 2: Run NGS pipeline
bash run_pipeline.sh \
    data/R1.fastq.gz \
    data/R2.fastq.gz \
    patient_001 \
    16

# Step 3: Annotate with ANNOVAR
bash annovar_helper.sh \
    results/patient_001/filtered/filtered_variants.vcf \
    results/patient_001/annovar/patient_001

# Step 4: Download results
# - results/patient_001/annovar/patient_001.hg19_multianno.vcf
# - results/patient_001/annovar/patient_001.hg19_multianno.txt
```

**Total Time**: 
- Pipeline: 1-2 hours (16 vCPU)
- ANNOVAR: 5-10 minutes
- **Total: ~1.5-2.5 hours**

---

## 🛠️ Troubleshooting

### Issue 1: "ANNOVAR not found"

```bash
# Check installation
ls -la ~/NGS/tools/annovar/

# Install if missing
cd ~/NGS
bash setup_annovar.sh
```

### Issue 2: "Cannot open database file"

```bash
# Check databases
ls -lh ~/NGS/tools/annovar/humandb/

# Download missing database
cd ~/NGS/tools/annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
```

### Issue 3: "Out of memory"

```bash
# Split VCF by chromosome
bcftools view -r chr1 input.vcf > chr1.vcf
bash annovar_helper.sh chr1.vcf chr1_out

# Or reduce databases
# Use only: -protocol refGene,clinvar_20240917
```

### Issue 4: Slow annotation

```bash
# Use fewer databases
-protocol refGene,avsnp150,clinvar_20240917

# Or parallelize by chromosome
for chr in {1..22} X Y; do
    bcftools view -r chr${chr} input.vcf > chr${chr}.vcf
    bash annovar_helper.sh chr${chr}.vcf chr${chr}_out &
done
wait
```

### Issue 5: Download fails during setup

```bash
# Download databases individually
cd ~/NGS/tools/annovar

# Essential only
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20240917 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome humandb/
```

---

## 💡 Tips & Tricks

### Speed Up Annotation

1. **Use fewer databases**:
   ```bash
   -protocol refGene,clinvar_20240917  # Instead of 9 databases
   ```

2. **Skip removal of intermediate files** during debugging:
   ```bash
   # Remove the -remove flag
   ```

3. **Parallelize by chromosome**:
   ```bash
   for chr in {1..22}; do
       bash annovar_helper.sh chr${chr}.vcf chr${chr}_out &
   done
   ```

### Save Space

```bash
# Compress output VCF
bgzip output.hg19_multianno.vcf
tabix -p vcf output.hg19_multianno.vcf.gz

# Remove intermediate files
rm *.avinput
```

### Excel Tips

```bash
# The .txt file opens directly in Excel
# Use filters on columns:
# - CLNSIG for pathogenicity
# - gnomAD_exome_ALL for frequency
# - SIFT_score for deleteriousness
```

---

## 📈 Performance

### Typical Runtimes (16 vCPU)

| Variants | Time | Memory |
|----------|------|--------|
| 1,000 | 30 sec | 500 MB |
| 10,000 | 2 min | 1 GB |
| 100,000 | 10 min | 2 GB |
| 1,000,000 | 90 min | 4 GB |

### Storage Requirements

| Component | Size |
|-----------|------|
| ANNOVAR program | 50 MB |
| hg19 databases (9) | 4-5 GB |
| Per analysis | 50-100 MB |

---

## 📚 Documentation Files

| File | Purpose |
|------|---------|
| **`ANNOVAR_README.md`** | **This file - Complete guide** |
| `ANNOVAR_QUICKSTART.md` | Quick reference card |
| `ANNOVAR_INSTALLATION_SUMMARY.md` | Installation details |
| `docs/ANNOVAR_INTEGRATION.md` | Detailed integration guide |
| `docs/ANNOTATION_TOOLS_COMPARISON.md` | Compare with other tools |

---

## 🎓 Learning Resources

### Official Documentation
- **Website**: http://annovar.openbioinformatics.org/
- **Manual**: http://annovar.openbioinformatics.org/en/latest/
- **Database List**: http://annovar.openbioinformatics.org/en/latest/user-guide/download/

### Community
- **Forum**: https://groups.google.com/g/annovar
- **Citations**: 10,000+ publications use ANNOVAR

### Video Tutorials
- Search YouTube for "ANNOVAR tutorial"
- Recommended: OpenHelix ANNOVAR tutorials

---

## 🔬 Example: Finding Disease-Causing Variants

### Step-by-Step Analysis

```bash
# 1. Annotate VCF
bash annovar_helper.sh patient.vcf annotated

# 2. Open in Excel
open annotated.hg19_multianno.txt

# 3. Apply filters (in Excel):
Filter 1: Func.refGene = "exonic"
Filter 2: ExonicFunc.refGene = "nonsynonymous SNV" OR "frameshift" OR "stopgain"
Filter 3: gnomAD_exome_ALL < 0.01 (rare variants)
Filter 4: CLNSIG contains "Pathogenic" OR SIFT_score < 0.05

# 4. Review remaining variants:
# - Check OMIM for gene-disease associations
# - Review ClinVar entries
# - Consider inheritance pattern
# - Validate with Sanger sequencing
```

---

## 🆘 Getting Help

1. **Test Installation**: `bash test_annovar.sh`
2. **Quick Reference**: `cat ANNOVAR_QUICKSTART.md`
3. **Detailed Guide**: `cat docs/ANNOVAR_INTEGRATION.md`
4. **ANNOVAR Forum**: https://groups.google.com/g/annovar
5. **GitHub Issues**: Create issue in your repository

---

## 📝 Citation

If you use ANNOVAR in your research, please cite:

> Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research, 38:e164, 2010.

---

## ✅ Quick Checklist

Before running ANNOVAR:
- [ ] ANNOVAR installed: `ls ~/NGS/tools/annovar/`
- [ ] Databases downloaded: `ls ~/NGS/tools/annovar/humandb/`
- [ ] Test passed: `bash test_annovar.sh`
- [ ] VCF file ready
- [ ] Enough disk space (~100 MB per analysis)

---

## 🎯 Common Use Cases

### Use Case 1: Rare Disease Diagnosis
```bash
bash annovar_helper.sh patient.vcf patient_annotated
# Focus on: CLNSIG, OMIM, rare variants
```

### Use Case 2: Cancer Genomics
```bash
perl table_annovar.pl tumor.vcf humandb/ \
    -buildver hg19 -out tumor \
    -protocol refGene,cosmic70,icgc28 \
    -operation g,f,f -vcfinput
```

### Use Case 3: Population Study
```bash
perl table_annovar.pl cohort.vcf humandb/ \
    -buildver hg19 -out cohort \
    -protocol refGene,gnomad211_exome,1000g2015aug \
    -operation g,f,f -vcfinput
```

---

## 🔗 Integration with Other Tools

### With IGV (Visualization)
```bash
# Open annotated VCF in IGV
java -jar igv.jar annotated.hg19_multianno.vcf
```

### With bcftools (Filtering)
```bash
# Filter annotated VCF
bcftools view -i 'INFO/CLNSIG="Pathogenic"' annotated.vcf
```

### With R/Python (Analysis)
```R
# Read in R
library(data.table)
variants <- fread("annotated.hg19_multianno.txt")
pathogenic <- variants[CLNSIG %like% "Pathogenic"]
```

---

**Last Updated**: $(date)  
**ANNOVAR Version**: Latest  
**Genome Build**: hg19 (GRCh37)

---

**Made with 🧬 for genomic variant annotation**

**For questions**: Check documentation or create an issue


# ANNOVAR Annotation Options - File Size Comparison

## 📊 The Problem

**Original filtered VCF:** 85 MB  
**After ANNOVAR with 9 databases:** 1.5 GB (18x bigger!)  

**Why?** ANNOVAR added ~150 columns from 9 databases with detailed annotations.

---

## 🎯 Three Annotation Options

| Option | Databases | Columns | File Size | Use When |
|--------|-----------|---------|-----------|----------|
| **Lightweight** ⚡ | 3 | ~20 | **10-20 MB** | Clinical diagnosis |
| **Standard** 📊 | 5 | ~50 | **30-50 MB** | General research |
| **Full** 🔬 | 9 | ~150 | **300+ MB** | Comprehensive research |

---

## ⚡ Option 1: Lightweight (RECOMMENDED)

**Command:**
```bash
bash annotate_lightweight.sh sample_name
```

**Databases:**
- ✅ `refGene` - Gene annotations (essential)
- ✅ `clinvar_20240917` - Clinical significance (essential)
- ✅ `gnomad211_exome` - Population frequencies (essential)

**Output:**
- **VCF:** 5-10 MB (compressed)
- **TXT:** 10-20 MB
- **Columns:** ~20

**Best for:**
- Clinical diagnosis
- Finding pathogenic variants
- Quick analysis
- Easy Excel handling

**Key columns you get:**
```
1. Gene name
2. Function (exonic, intronic)
3. Exonic function (nonsynonymous, frameshift)
4. Amino acid change
5. Clinical significance (ClinVar)
6. Population frequency (gnomAD)
7. Disease associations
```

---

## 📊 Option 2: Standard

**Command:**
```bash
bash annotate_pass_only.sh sample_name
```

**Databases:**
- ✅ `refGene` - RefSeq genes
- ✅ `clinvar_20240917` - Clinical significance
- ✅ `gnomad211_exome` - Population frequencies
- ✅ `avsnp150` - dbSNP IDs
- ✅ `dbnsfp42a` - Functional predictions (SIFT, PolyPhen)

**Output:**
- **VCF:** 15-25 MB (compressed)
- **TXT:** 30-50 MB
- **Columns:** ~50

**Best for:**
- Research studies
- Functional impact analysis
- Prioritizing variants
- Publications

**Additional columns:**
```
+ dbSNP IDs (rs numbers)
+ SIFT scores (deleteriousness)
+ PolyPhen scores (pathogenicity)
+ CADD scores
+ Conservation scores
+ Multiple prediction algorithms
```

---

## 🔬 Option 3: Full (Current - TOO LARGE)

**Command:**
```bash
# This is what created your 1.5GB files
# NOT RECOMMENDED unless you need everything
```

**Databases:**
- All 9 databases including COSMIC, ICGC, multiple gene databases

**Output:**
- **VCF:** 1.5 GB (uncompressed) or 60-100 MB (compressed)
- **TXT:** 300-500 MB
- **Columns:** ~150

**Problems:**
- ❌ Too large to download easily
- ❌ Slow to open in Excel
- ❌ Contains redundant information
- ❌ Most columns unused

**Only use if:**
- Cancer research (need COSMIC, ICGC)
- Comparing multiple gene databases
- Very specialized research

---

## 📋 What Makes Files So Large?

### Database Contributions:

**refGene (Essential)** - ~5 columns
- Gene name
- Function
- Exonic function
- Amino acid change
- Transcript ID

**clinvar_20240917 (Essential)** - ~5 columns
- Clinical significance
- Review status
- Disease name
- Allele ID
- Last reviewed

**gnomad211_exome (Essential)** - ~10 columns
- Total allele frequency
- African frequency
- East Asian frequency
- European frequency
- South Asian frequency
- Latino frequency
- Other populations
- Homozygote counts

**avsnp150** - ~2 columns
- dbSNP ID (rs number)
- Validation status

**dbnsfp42a** - **50-100 columns!** 🚨
- SIFT score and prediction
- PolyPhen2 HDIV score and prediction
- PolyPhen2 HVAR score and prediction
- LRT score and prediction
- MutationTaster score and prediction
- MutationAssessor score and prediction
- FATHMM score
- PROVEAN score
- MetaSVM score
- MetaLR score
- CADD scores (raw and phred)
- DANN score
- PhyloP scores (multiple)
- PhastCons scores (multiple)
- GERP++ scores
- SiPhy scores
- ...and 70+ more!

**knownGene, ensGene** - ~10 columns each
- Duplicate gene annotations (often redundant)

**cosmic70, icgc28** - ~10 columns each
- Cancer-specific data
- Usually empty for germline variants

---

## 💡 Recommendation

### For Clinical Work:
```bash
bash annotate_lightweight.sh sample_name
```
**Result:** 10-20 MB - Perfect for Excel, has everything you need!

### For Research:
```bash
bash annotate_pass_only.sh sample_name
```
**Result:** 30-50 MB - Includes functional predictions

### Need Everything?
First try lightweight, then add more if needed:
```bash
# Start lightweight
bash annotate_lightweight.sh sample_name

# If you need functional predictions, upgrade to standard
bash annotate_pass_only.sh sample_name
```

---

## 🔍 Check What You Have

**Run this to see columns:**
```bash
bash check_annotation_columns.sh sample_name
```

**Sample output:**
```
COLUMN HEADERS:
1. Chr
2. Start
3. End
4. Ref
5. Alt
6. Func.refGene
7. Gene.refGene
8. GeneDetail.refGene
9. ExonicFunc.refGene
10. AAChange.refGene
...150 total columns
```

---

## 📊 File Size Comparison

**Your current situation:**

| Sample | Filtered VCF | Full Annotation | Lightweight | Savings |
|--------|--------------|-----------------|-------------|---------|
| sample_Alaa | 85 MB | 1.5 GB | 15 MB | **99%** |
| sample_6 | 85 MB | 1.6 GB | 16 MB | **99%** |

---

## 🚀 Quick Fix Commands

**For Jarvis Lab - Run these:**

```bash
cd ~/NGS

# Option 1: Lightweight (fastest, smallest)
bash annotate_lightweight.sh sample_Alaa
bash annotate_lightweight.sh sample_6

# Option 2: Standard (more info)
bash annotate_pass_only.sh sample_Alaa
bash annotate_pass_only.sh sample_6

# Check what was added
bash check_annotation_columns.sh sample_Alaa
```

---

## 📥 What to Download

**Lightweight:**
- `annotated_sample_Alaa_lightweight.hg19_multianno.txt` (15 MB)

**Standard:**
- `annotated_sample_Alaa_PASS.hg19_multianno.txt` (35 MB)

**Full (not recommended):**
- `annotated_sample_Alaa.hg19_multianno.txt` (314 MB)

---

## ✅ Summary

**The problem:** Used all 9 databases = too many columns = huge files

**The solution:** Use only essential databases = fewer columns = manageable files

**Best practice:**
1. Start with **lightweight** (3 databases, 15MB)
2. Upgrade to **standard** if you need predictions (5 databases, 35MB)
3. Never use full unless specifically needed (9 databases, 300MB+)

---

**Made with 📊 for optimal file sizes**


# Final Annotation Configuration

## ✅ **Your Requested Databases Included**

Based on your requirements, the annotation now includes:

### **5 Essential Databases:**

1. ✅ **refGene** - Gene annotations
2. ✅ **clinvar_20240917** - Clinical significance  
3. ✅ **gnomad211_exome** - Population frequencies
4. ✅ **avsnp150** - **dbSNP (rs numbers)** ← YOU REQUESTED
5. ✅ **dbnsfp42a** - **Prediction scores** ← YOU REQUESTED

---

## 📊 **What Each Database Provides**

### 1. **refGene** (Gene-based)
```
Columns added: ~5
Information:
- Gene name (BRCA1, TP53, etc.)
- Function (exonic, intronic, UTR, splicing)
- Exonic function (nonsynonymous, frameshift, stopgain)
- Amino acid change (p.V600E, p.R273H, etc.)
- Transcript IDs
```

### 2. **clinvar_20240917** (Clinical)
```
Columns added: ~5
Information:
- Clinical significance (Pathogenic, Benign, VUS)
- Disease associations
- Review status (stars)
- ClinVar IDs
- Last reviewed date
```

### 3. **gnomad211_exome** (Population)
```
Columns added: ~10
Information:
- Overall allele frequency
- African (AFR) frequency
- Latino (AMR) frequency
- Ashkenazi Jewish (ASJ) frequency
- East Asian (EAS) frequency
- Finnish (FIN) frequency
- European (NFE) frequency
- South Asian (SAS) frequency
- Homozygote counts
```

### 4. **avsnp150** (dbSNP) - YOUR REQUEST ⭐
```
Columns added: ~2
Information:
- rs numbers (e.g., rs80357906, rs1800497)
- Variant IDs for literature search
- PubMed reference linking

Why important:
✅ Search variants in literature (PubMed, Google Scholar)
✅ Compare with published studies
✅ Standardized variant nomenclature
✅ Cross-reference with other databases
```

### 5. **dbnsfp42a** (Predictions) - YOUR REQUEST ⭐
```
Columns added: ~40-50
Information:
- SIFT_score (< 0.05 = damaging)
- SIFT_pred (D = deleterious, T = tolerated)
- Polyphen2_HDIV_score (> 0.85 = probably damaging)
- Polyphen2_HDIV_pred (D = damaging, P = possibly damaging, B = benign)
- CADD_phred (> 20 = deleterious)
- REVEL_score (> 0.5 = likely pathogenic)
- MetaSVM_pred (D = damaging, T = tolerated)
- MetaLR_pred (D = damaging, T = tolerated)
- DANN_score
- FATHMM_pred
- MutationTaster_pred
- MutationAssessor_score
- PROVEAN_pred
- VEST4_score
- M-CAP_pred
- PrimateAI_pred
- DEOGEN2_pred
- BayesDel_addAF_pred
- And 30+ more algorithms!

Conservation scores:
- phyloP100way_vertebrate
- phastCons100way_vertebrate
- GERP++_RS
- SiPhy_29way_logOdds

Why important:
✅ Predict variant deleteriousness
✅ Multiple independent algorithms
✅ Consensus predictions more reliable
✅ Required for variant prioritization
✅ Needed for publications
```

---

## 📏 **File Size Comparison**

| Configuration | Databases | Columns | File Size | Your Case |
|---------------|-----------|---------|-----------|-----------|
| **Lightweight** | 3 | ~20 | 10-15 MB | Too minimal |
| **Clinical Plus** ⭐ | 5 | ~50 | **30-50 MB** | **✅ Perfect!** |
| **Full (old)** | 9 | ~150 | 300 MB+ | ❌ Too large |

---

## 🎯 **What You Get Now**

With these 5 databases, you can:

### Clinical Interpretation:
- ✅ Identify disease-causing genes
- ✅ Check ClinVar pathogenicity
- ✅ Assess variant rarity (gnomAD)
- ✅ **Search literature using rs numbers** (dbSNP)
- ✅ **Predict functional impact** (dbnsfp42a)

### Research & Publication:
- ✅ Multiple prediction algorithms for consensus
- ✅ Conservation scores for evolutionary importance
- ✅ Standardized variant IDs (rs numbers)
- ✅ Comprehensive functional annotations
- ✅ Ready for supplementary tables

### Variant Prioritization:
- ✅ Filter by prediction scores
- ✅ Combine multiple algorithms
- ✅ Check conservation across species
- ✅ Literature evidence via dbSNP
- ✅ Population frequency filtering

---

## 🚀 **Commands to Use**

### In Jarvis Lab:

```bash
cd ~/NGS
git pull

# Option 1: Use the new clinical_plus script (recommended)
bash annotate_clinical_plus.sh sample_Alaa
bash annotate_clinical_plus.sh sample_6

# Option 2: Use the updated standard script (same 5 databases)
bash annotate_pass_only.sh sample_Alaa
bash annotate_pass_only.sh sample_6

# Both commands now use the same 5 databases you requested!
```

---

## 📊 **Example Columns You'll See**

```
Column 1:  Chr
Column 2:  Start
Column 3:  End
Column 4:  Ref
Column 5:  Alt
Column 6:  Func.refGene
Column 7:  Gene.refGene
Column 8:  GeneDetail.refGene
Column 9:  ExonicFunc.refGene
Column 10: AAChange.refGene
Column 11: CLNSIG (ClinVar significance)
Column 12: CLNDN (ClinVar disease name)
Column 13: CLNREVSTAT (Review status)
Column 14: gnomAD_exome_ALL (Overall frequency)
Column 15: gnomAD_exome_AFR (African frequency)
Column 16: gnomAD_exome_EAS (East Asian frequency)
Column 17: gnomAD_exome_NFE (European frequency)
Column 18: avsnp150 (rs number) ⭐
Column 19: SIFT_score ⭐
Column 20: SIFT_pred ⭐
Column 21: Polyphen2_HDIV_score ⭐
Column 22: Polyphen2_HDIV_pred ⭐
Column 23: CADD_phred ⭐
Column 24: REVEL_score ⭐
Column 25: MetaSVM_pred ⭐
... and 30 more prediction columns!
```

---

## 💡 **Filtering Examples**

### In Excel - Find Pathogenic Variants:

**Filter 1: Clinical significance**
```
CLNSIG contains "Pathogenic"
```

**Filter 2: Rare variants**
```
gnomAD_exome_ALL < 0.01
```

**Filter 3: Predicted damaging**
```
SIFT_score < 0.05
AND Polyphen2_HDIV_score > 0.85
AND CADD_phred > 20
```

**Filter 4: Has dbSNP ID**
```
avsnp150 is not empty
```

### Command Line - Quick Counts:

```bash
# Count pathogenic variants
grep -i "pathogenic" annotated_sample_Alaa_clinical_plus.hg19_multianno.txt | wc -l

# Find variants with rs numbers
grep -v "^\." annotated_sample_Alaa_clinical_plus.hg19_multianno.txt | awk -F'\t' '$18 != "."' | wc -l

# Find damaging predictions (SIFT)
awk -F'\t' '$19 < 0.05' annotated_sample_Alaa_clinical_plus.hg19_multianno.txt | wc -l
```

---

## 🎓 **Key Prediction Score Thresholds**

| Score | Damaging Threshold | Interpretation |
|-------|-------------------|----------------|
| **SIFT** | < 0.05 | Deleterious |
| **PolyPhen2** | > 0.85 | Probably damaging |
| **CADD** | > 20 | Deleterious (> 30 = highly deleterious) |
| **REVEL** | > 0.5 | Likely pathogenic (> 0.75 = highly likely) |
| **MetaSVM** | D | Damaging |
| **GERP++** | > 4 | Conserved (> 6 = highly conserved) |

**Best practice:** Use multiple scores for consensus!

---

## ✅ **Summary**

### What Changed:
- ❌ Removed: knownGene, ensGene (redundant gene annotations)
- ❌ Removed: COSMIC, ICGC (cancer databases - mostly empty for germline)
- ✅ Kept: refGene (essential genes)
- ✅ Kept: ClinVar (essential clinical)
- ✅ Kept: gnomAD (essential frequencies)
- ✅ Added: **dbSNP** (your request - literature search)
- ✅ Added: **Predictions** (your request - functional impact)

### Result:
- **File size:** 30-50 MB (98% smaller than 1.5 GB!)
- **Information:** All essential + predictions + dbSNP
- **Quality:** Optimal for clinical diagnosis and research
- **Usability:** Excel-friendly, fast downloads

---

**Your annotation now includes everything you need!** ⭐

**Run:** `bash annotate_clinical_plus.sh sample_Alaa` in Jarvis Lab


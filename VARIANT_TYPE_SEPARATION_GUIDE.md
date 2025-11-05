# Variant Type Separation Guide

## 🎯 Overview

After annotation, variants can be separated by type for focused analysis:
- **SNPs** - Single Nucleotide Polymorphisms
- **Insertions** - Added nucleotides
- **Deletions** - Removed nucleotides
- **MNPs** - Multiple Nucleotide Polymorphisms

---

## 🚀 Quick Start

### Option 1: Complete Workflow (Annotate + Separate)

```bash
cd ~/NGS
git pull
bash run_complete_with_separation.sh sample_Alaa
bash run_complete_with_separation.sh sample_6
```

**This automatically:**
1. ✅ Annotates with 5 databases (30-50MB)
2. ✅ Separates by variant type
3. ✅ Generates statistics
4. ✅ Creates summary reports

**Time:** ~10-15 minutes per sample

---

### Option 2: Separate After Annotation

If you already have annotated files:

```bash
bash separate_by_variant_type.sh sample_Alaa
bash separate_by_variant_type.sh sample_6
```

---

## 📊 Variant Types Explained

### 1. **SNPs (Single Nucleotide Polymorphisms)**
```
Example:
Ref: A
Alt: G

Most common variant type (~70-80% of all variants)
Single base change: A→G, C→T, etc.
```

**Clinical significance:**
- Most common genetic variation
- Can cause amino acid changes (nonsynonymous)
- May be silent (synonymous)
- Easy to genotype and validate

### 2. **Insertions**
```
Example:
Ref: A
Alt: AGTC

Added nucleotides
Length(Alt) > Length(Ref)
```

**Clinical significance:**
- Can cause frameshifts if not multiple of 3
- May introduce premature stop codons
- Can disrupt protein function
- Common in repetitive regions

### 3. **Deletions**
```
Example:
Ref: AGTC
Alt: A

Removed nucleotides
Length(Ref) > Length(Alt)
```

**Clinical significance:**
- Often more pathogenic than SNPs
- Can cause frameshifts
- May remove functional domains
- Associated with many genetic diseases

### 4. **MNPs (Multiple Nucleotide Polymorphisms)**
```
Example:
Ref: AGT
Alt: TCA

Multiple bases changed
Same length but >1 nucleotide different
```

**Clinical significance:**
- Less common (~1-2% of variants)
- Can affect multiple amino acids
- May impact protein structure significantly

---

## 📁 Output Structure

After running, you'll get:

```
results/sample_Alaa/annovar/
├── annotated_sample_Alaa_clinical_plus.hg19_multianno.txt (35 MB)
├── annotated_sample_Alaa_clinical_plus.hg19_multianno.vcf.gz (20 MB)
└── separated_by_type/
    ├── SNPs.txt           (~25 MB, ~70% of variants)
    ├── Insertions.txt     (~4 MB, ~15% of variants)
    ├── Deletions.txt      (~5 MB, ~10% of variants)
    └── MNPs.txt           (~1 MB, ~5% of variants)
```

---

## 📊 Expected Distribution

**Typical exome (PASS variants):**

| Type | Count | Percentage | Clinical Impact |
|------|-------|------------|----------------|
| **SNPs** | ~35,000-42,000 | 70-80% | Variable |
| **Insertions** | ~3,000-9,000 | 6-18% | Often pathogenic |
| **Deletions** | ~3,000-9,000 | 6-18% | Often pathogenic |
| **MNPs** | ~500-2,000 | 1-4% | Rare |
| **Total** | ~50,000 | 100% | |

---

## 🔍 Analysis by Type

### Why Separate?

Different variant types have different characteristics:

**SNPs:**
- Most common, well-studied
- Large databases (dbSNP, ClinVar)
- Prediction tools optimized for SNPs
- Easier to validate

**Indels (Insertions + Deletions):**
- More likely to be pathogenic (if frameshift)
- Less common in population
- Harder to sequence/call
- Higher false positive rate
- Need extra validation

**MNPs:**
- Rare, less well-studied
- May represent sequencing errors
- Or true di/tri-nucleotide changes
- Require careful review

---

## 💡 Analysis Strategies

### For SNPs (SNPs.txt)

**Filter for:**
```
1. Nonsynonymous changes
   ExonicFunc.refGene = "nonsynonymous SNV"

2. Predicted damaging
   SIFT_score < 0.05
   Polyphen2_HDIV_score > 0.85

3. Rare variants
   gnomAD_exome_ALL < 0.01

4. Known pathogenic
   CLNSIG contains "Pathogenic"
```

### For Insertions/Deletions (Indels)

**Check for:**
```
1. Frameshift variants
   ExonicFunc.refGene = "frameshift insertion" or "frameshift deletion"

2. Very rare (often pathogenic if rare)
   gnomAD_exome_ALL < 0.001

3. In known disease genes
   Cross-reference with gene panels

4. ClinVar pathogenic indels
   CLNSIG contains "Pathogenic"
```

### For MNPs

**Review carefully:**
```
1. Manual inspection recommended
2. Check sequencing quality
3. May be sequencing artifacts
4. Or true complex variants
5. Literature search with rs numbers
```

---

## 📈 Statistics You'll Get

The separation script provides:

### Summary Table:
```
Type                    Count    Percent
────────────────────────────────────────
SNPs                   38,547     77.1%
Insertions              5,234     10.5%
Deletions               5,892     11.8%
MNPs                      391      0.6%
────────────────────────────────────────
Total                  50,064    100%
```

### Pathogenic Counts by Type:
```
Type                    Pathogenic Count
────────────────────────────────────────
SNPs                              42
Insertions                        18
Deletions                         23
MNPs                               2
```

---

## 🎓 Clinical Interpretation

### General Rules:

**For SNPs:**
- Check ClinVar first
- Use multiple prediction algorithms
- Consider conservation scores
- Look at population frequency

**For Insertions/Deletions:**
- Frameshifts are often pathogenic
- In-frame indels: depends on location
- Check if removes functional domain
- Very rare indels: suspicious for pathogenicity

**For MNPs:**
- Often false positives
- Validate with Sanger sequencing
- Check read depth and quality
- Search literature

---

## 📥 Download and Analyze

### In Excel:

**Open each type separately:**

1. **SNPs.txt**
   - Filter: SIFT_score < 0.05
   - Filter: gnomAD_exome_ALL < 0.01
   - Sort by CADD_phred descending

2. **Insertions.txt + Deletions.txt**
   - Filter: ExonicFunc contains "frameshift"
   - Filter: gnomAD_exome_ALL < 0.001
   - Check ClinVar column

3. **MNPs.txt**
   - Review all (usually small number)
   - Check read quality
   - Validate interesting ones

---

## 🔬 Advanced: Indel Size Distribution

You can further analyze indel sizes:

```bash
# For insertions - count by size
awk -F'\t' 'NR>1 {print length($5)-length($4)}' separated_by_type/Insertions.txt | sort -n | uniq -c

# For deletions - count by size
awk -F'\t' 'NR>1 {print length($4)-length($5)}' separated_by_type/Deletions.txt | sort -n | uniq -c
```

**Interpretation:**
- Size = 1 or 2: Often in homopolymer regions (sequencing errors)
- Size = 3, 6, 9: In-frame (less likely pathogenic)
- Size not multiple of 3: Frameshift (more likely pathogenic)

---

## 🎯 Prioritization Strategy

### Phase 1: Quick wins (SNPs)
1. Open SNPs.txt
2. Filter: CLNSIG = "Pathogenic"
3. Review these first (known pathogenic)

### Phase 2: Frameshift indels
1. Open Insertions.txt + Deletions.txt
2. Filter: ExonicFunc contains "frameshift"
3. Filter: gnomAD < 0.001
4. Check gene relevance

### Phase 3: Predicted damaging SNPs
1. Back to SNPs.txt
2. Filter: SIFT < 0.05 AND PolyPhen > 0.85
3. Filter: gnomAD < 0.01
4. Sort by CADD score

### Phase 4: Everything else
1. Review MNPs manually
2. Check in-frame indels
3. Look at splice variants

---

## 📞 Quick Commands

```bash
# Complete workflow
bash run_complete_with_separation.sh sample_Alaa

# Just separate (if already annotated)
bash separate_by_variant_type.sh sample_Alaa

# Check separation results
ls -lh ~/NGS/results/sample_Alaa/annovar/separated_by_type/

# Count pathogenic by type
grep -i "pathogenic" ~/NGS/results/sample_Alaa/annovar/separated_by_type/SNPs.txt | wc -l
grep -i "pathogenic" ~/NGS/results/sample_Alaa/annovar/separated_by_type/Insertions.txt | wc -l
grep -i "pathogenic" ~/NGS/results/sample_Alaa/annovar/separated_by_type/Deletions.txt | wc -l
```

---

## ✅ Summary

**Benefits of separation:**
1. ✅ Focused analysis by variant type
2. ✅ Different filtering strategies
3. ✅ Better interpretation
4. ✅ Easier to spot patterns
5. ✅ Type-specific validation

**Best workflow:**
```bash
cd ~/NGS && git pull
bash run_complete_with_separation.sh sample_Alaa
# Download all files in separated_by_type/ folder
# Analyze each type separately in Excel
```

---

**Made with 🧬 for comprehensive variant analysis**


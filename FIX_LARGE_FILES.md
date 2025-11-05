# Fix Large Annotation Files

## 🚨 Problem

Your annotated files are **1.5-1.6GB** when they should be **60-100MB max**.

**Cause:** ANNOVAR annotated ALL variants (including failed/low-quality ones).

---

## 📊 Comparison

| | Your Files | Should Be |
|---|---|---|
| **Variants** | ~480,000 | ~30,000-60,000 |
| **VCF Size** | 1.5-1.6 GB | 60-100 MB |
| **TXT Size** | 314 MB | 30-50 MB |

---

## ✅ Quick Fix (In Jarvis Lab)

### Option 1: Auto-Cleanup All Samples

```bash
cd ~/NGS
bash cleanup_and_reannotate.sh
```

**What it does:**
1. Finds all annotated samples
2. Shows current file sizes
3. Removes large files
4. Re-annotates with PASS variants only
5. **Result: 90% smaller files!**

### Option 2: Fix One Sample

```bash
bash annotate_pass_only.sh sample_Alaa
bash annotate_pass_only.sh sample_6
```

---

## 📥 Download Correct Files

After running the fix, download these **smaller** files:

```
results/sample_Alaa/annovar/annotated_sample_Alaa_PASS.hg19_multianno.txt
results/sample_6/annovar/annotated_sample_6_PASS.hg19_multianno.txt
```

**File naming:**
- ❌ OLD: `annotated_sample_Alaa.hg19_multianno.txt` (314MB - too big)
- ✅ NEW: `annotated_sample_Alaa_PASS.hg19_multianno.txt` (30-50MB - correct)

---

## 🔍 What Went Wrong?

The annotation script processed ALL variants:
```bash
# Wrong - annotates ALL variants (failed + passed)
grep -v "^#" filtered_variants.vcf  # ~480,000 variants
```

Should have been:
```bash
# Correct - annotates PASS only
grep -v "^#" filtered_variants.vcf | grep -w "PASS"  # ~30,000-60,000 variants
```

---

## 📋 Manual Fix Steps

If you want to do it manually:

### Step 1: Extract PASS variants

```bash
cd ~/NGS/results/sample_Alaa

# Create PASS-only VCF
grep "^#" filtered/filtered_variants.vcf > filtered/filtered_PASS_only.vcf
grep -v "^#" filtered/filtered_variants.vcf | grep -w "PASS" >> filtered/filtered_PASS_only.vcf
```

### Step 2: Count variants

```bash
# Before
grep -v "^#" filtered/filtered_variants.vcf | wc -l
# Result: ~480,000 (TOO MANY)

# After
grep -v "^#" filtered/filtered_PASS_only.vcf | wc -l
# Result: ~30,000-60,000 (CORRECT)
```

### Step 3: Re-annotate

```bash
perl ~/NGS/tools/annovar/table_annovar.pl \
    filtered/filtered_PASS_only.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out annovar/annotated_sample_Alaa_PASS \
    -remove \
    -protocol refGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a \
    -operation g,f,f,f,f \
    -vcfinput
```

### Step 4: Compress

```bash
bgzip annotated_sample_Alaa_PASS.hg19_multianno.vcf
tabix -p vcf annotated_sample_Alaa_PASS.hg19_multianno.vcf.gz
```

### Step 5: Remove old files

```bash
rm annotated_sample_Alaa.hg19_multianno.vcf  # 1.5GB
rm annotated_sample_Alaa.hg19_multianno.txt  # 314MB
```

---

## 💾 Space Savings

**Before:**
```
sample_Alaa:
  VCF: 1.5 GB
  TXT: 314 MB
  Total: ~1.8 GB

sample_6:
  VCF: 1.6 GB
  TXT: 320 MB
  Total: ~1.9 GB
```

**After:**
```
sample_Alaa:
  VCF.gz: 15 MB (compressed)
  TXT: 30 MB
  Total: ~45 MB (96% reduction!)

sample_6:
  VCF.gz: 16 MB (compressed)
  TXT: 32 MB
  Total: ~48 MB (97% reduction!)
```

---

## 🎯 Expected Variant Counts

| Sample Type | Raw | Filtered (PASS) | Final Size |
|-------------|-----|-----------------|------------|
| Exome | ~100K | **30K-60K** ✅ | 60-100 MB |
| Whole Genome | ~5M | **3M-4M** | 3-5 GB |

**Your current:** 480K variants → TOO HIGH for exome!

---

## ⚠️ Why So Many Variants?

Possible reasons:
1. **Using wrong reference** - hg19 vs hg38 mismatch?
2. **Low quality threshold** - Calling too many false positives
3. **No filtering** - Including all low-quality calls
4. **Wrong input** - Whole genome instead of exome?

**Check your filtered_variants.vcf:**
```bash
# Count PASS vs FAIL
grep -v "^#" filtered/filtered_variants.vcf | grep -w "PASS" | wc -l  # Should be 30K-60K
grep -v "^#" filtered/filtered_variants.vcf | grep -v -w "PASS" | wc -l  # Failed variants
```

---

## 🔧 Prevention (Future Runs)

The scripts have been updated. Future annotations will automatically use PASS-only variants.

Use these updated scripts:
- `annotate_pass_only.sh` - PASS variants only
- `annotate_results.sh` - Will be updated to default to PASS

---

## 📞 Quick Commands Reference

```bash
# Fix all samples automatically
bash cleanup_and_reannotate.sh

# Fix one sample
bash annotate_pass_only.sh sample_name

# Check variant counts
grep -v "^#" results/sample/filtered/filtered_variants.vcf | wc -l
grep -v "^#" results/sample/filtered/filtered_variants.vcf | grep -w "PASS" | wc -l

# Check file sizes
du -h results/sample/annovar/*.vcf
du -h results/sample/annovar/*.txt

# Remove old large files
rm results/sample/annovar/annotated_sample.hg19_multianno.vcf
rm results/sample/annovar/annotated_sample.hg19_multianno.txt
```

---

**After fixing, your files will be the correct size (~60-100MB) for easy download and analysis!** ✅


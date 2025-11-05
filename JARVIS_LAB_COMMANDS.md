# Commands to Run in Jarvis Lab

## 🎯 You Are Here: ANNOVAR Successfully Installed! ✅

Now let's annotate your VCF files.

---

## 📋 Step 1: Find Your VCF Files

Run these commands **in Jarvis Lab terminal**:

```bash
# Navigate to NGS directory
cd ~/NGS

# Find all VCF files
echo "=== Looking for VCF files ==="
find results -name "*.vcf" -type f 2>/dev/null | grep -v ".vcf.gz"
```

---

## 🚀 Step 2: Quick Automated Annotation (Recommended)

This script will find all VCF files and annotate them automatically:

```bash
cd ~/NGS
bash annotate_existing_vcfs.sh
```

**What this does:**
- ✅ Automatically finds all VCF files in `results/`
- ✅ Shows you what will be annotated
- ✅ Asks for confirmation
- ✅ Annotates each VCF with ANNOVAR
- ✅ Shows progress and results

**Time:** ~5-10 minutes per VCF file

---

## 💻 Step 3: Manual Annotation (If You Prefer)

If you want to annotate specific files manually:

### Find Your Sample IDs First

```bash
ls -la ~/NGS/results/
```

You should see directories like:
- `c284f797-37ff-4bd7-9849-669b936d9dad`
- `26a84d8c-eaf4-4d17-abde-12518537fb6f`

### Annotate First Sample

```bash
# Set your sample ID
SAMPLE_ID="c284f797-37ff-4bd7-9849-669b936d9dad"  # Change this

# Create output directory
mkdir -p ~/NGS/results/$SAMPLE_ID/annovar

# Run ANNOVAR
bash annovar_helper.sh \
    ~/NGS/results/$SAMPLE_ID/filtered/filtered_variants.vcf \
    ~/NGS/results/$SAMPLE_ID/annovar/annotated_${SAMPLE_ID}
```

### Annotate Second Sample

```bash
# Set your second sample ID
SAMPLE_ID="26a84d8c-eaf4-4d17-abde-12518537fb6f"  # Change this

# Create output directory
mkdir -p ~/NGS/results/$SAMPLE_ID/annovar

# Run ANNOVAR
bash annovar_helper.sh \
    ~/NGS/results/$SAMPLE_ID/filtered/filtered_variants.vcf \
    ~/NGS/results/$SAMPLE_ID/annovar/annotated_${SAMPLE_ID}
```

---

## 📊 Step 4: Check Results

After annotation completes:

```bash
# List all annotated files
find ~/NGS/results -name "*.hg19_multianno.*" -type f

# Check size and variant counts
for vcf in $(find ~/NGS/results -name "*.hg19_multianno.vcf" -type f); do
    echo "=== $vcf ==="
    echo "Size: $(du -h "$vcf" | cut -f1)"
    echo "Variants: $(grep -v "^#" "$vcf" | wc -l)"
    echo ""
done
```

---

## 📥 Step 5: Download Results

### From JupyterLab Interface

1. Navigate to `~/NGS/results/SAMPLE_ID/annovar/`
2. Right-click on files
3. Select "Download"

### Using SCP (from your local machine)

```bash
# Download annotated VCF
scp user@jarvis-lab:~/NGS/results/SAMPLE_ID/annovar/*.hg19_multianno.vcf ./

# Download TXT file (for Excel)
scp user@jarvis-lab:~/NGS/results/SAMPLE_ID/annovar/*.hg19_multianno.txt ./
```

---

## 🔍 Step 6: Analyze Results

### Open TXT File in Excel

The `.hg19_multianno.txt` file can be opened directly in Excel.

**Filter for pathogenic variants:**
1. Open file in Excel
2. Apply AutoFilter (Data → Filter)
3. Filter columns:
   - `CLNSIG` contains "Pathogenic"
   - `gnomAD_exome_ALL` < 0.01 (rare)
   - `ExonicFunc.refGene` = "nonsynonymous SNV" or "frameshift"

### Quick Command-Line Check

```bash
# Count total variants
grep -v "^#" annotated.hg19_multianno.vcf | wc -l

# Find pathogenic variants
grep "Pathogenic" annotated.hg19_multianno.txt | wc -l

# List genes with pathogenic variants
grep "Pathogenic" annotated.hg19_multianno.txt | cut -f7 | sort -u
```

---

## 📋 Expected Output Files

For each sample, you'll get:

```
results/SAMPLE_ID/annovar/
├── annotated_SAMPLE_ID.hg19_multianno.vcf  (VCF format - for tools)
└── annotated_SAMPLE_ID.hg19_multianno.txt  (Table format - for Excel)
```

---

## 🎓 Understanding Output Columns

**Key columns in the .txt file:**

| Column | Meaning | Example |
|--------|---------|---------|
| `Func.refGene` | Location | exonic, intronic |
| `Gene.refGene` | Gene name | BRCA1, TP53 |
| `ExonicFunc.refGene` | Variant type | nonsynonymous SNV |
| `AAChange.refGene` | Protein change | p.V600E |
| `avsnp150` | dbSNP ID | rs80357906 |
| `gnomAD_exome_ALL` | Frequency | 0.0001 (0.01%) |
| `CLNSIG` | Clinical sig | Pathogenic |
| `SIFT_score` | Deleteriousness | 0.001 (damaging) |
| `Polyphen2_HDIV_score` | Pathogenicity | 0.999 (damaging) |

---

## 🔧 Troubleshooting

### If VCF File Not Found

```bash
# Check exact location
find ~/NGS/results -name "*.vcf" -type f

# If in compressed format
find ~/NGS/results -name "*.vcf.gz" -type f

# Decompress if needed
gunzip ~/NGS/results/SAMPLE_ID/variants.vcf.gz
```

### If Annotation Fails

```bash
# Test ANNOVAR installation
bash test_annovar.sh

# Check database availability
ls -lh ~/NGS/tools/annovar/humandb/ | grep hg19

# Try with single database first
perl ~/NGS/tools/annovar/table_annovar.pl \
    input.vcf \
    ~/NGS/tools/annovar/humandb/ \
    -buildver hg19 \
    -out test \
    -protocol refGene \
    -operation g \
    -vcfinput
```

---

## ⚡ Quick Reference

```bash
# 1. Navigate to NGS directory
cd ~/NGS

# 2. Run automated annotation
bash annotate_existing_vcfs.sh

# 3. Check results
find results -name "*.hg19_multianno.*"

# 4. Download files (from your local machine)
scp user@jarvis:~/NGS/results/SAMPLE/annovar/*.txt ./
```

---

## 📞 Need Help?

- **Test installation**: `bash test_annovar.sh`
- **Check documentation**: `cat ANNOVAR_QUICKSTART.md`
- **Detailed guide**: `cat ANNOVAR_README.md`

---

**All commands above should be run in your Jarvis Lab terminal!** 🚀


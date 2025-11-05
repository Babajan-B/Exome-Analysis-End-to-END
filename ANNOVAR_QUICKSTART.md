# ANNOVAR Quick Start Guide for Jarvis Lab

## 🚀 Installation (One Command)

```bash
cd ~/NGS && bash setup_annovar.sh
```

**Time**: 20-30 minutes  
**Storage**: ~5 GB

---

## ⚡ Quick Usage

### 1. Annotate VCF File (Easiest)

```bash
bash annovar_helper.sh input.vcf output_prefix
```

### 2. Annotate After Pipeline Run

```bash
# After running the main pipeline
bash annovar_helper.sh \
    results/SAMPLE_NAME/filtered/filtered_variants.vcf \
    results/SAMPLE_NAME/annovar/annotated
```

---

## 📋 Common Commands

### Basic Annotation (RefSeq genes only)

```bash
cd ~/NGS/tools/annovar
perl table_annovar.pl \
    input.vcf \
    humandb/ \
    -buildver hg19 \
    -out output \
    -protocol refGene \
    -operation g \
    -vcfinput
```

### Clinical Annotation (Recommended)

```bash
perl table_annovar.pl \
    input.vcf \
    humandb/ \
    -buildver hg19 \
    -out output \
    -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150 \
    -operation g,f,f,f \
    -vcfinput
```

### Full Annotation (Research)

```bash
perl table_annovar.pl \
    input.vcf \
    humandb/ \
    -buildver hg19 \
    -out output \
    -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a,cosmic70 \
    -operation g,g,g,f,f,f,f,f \
    -vcfinput \
    -polish
```

---

## 🧬 Understanding Output

### Output Files

- **`output.hg19_multianno.vcf`** - Annotated VCF file
- **`output.hg19_multianno.txt`** - Tab-delimited table (easier to read in Excel)

### Key Columns in `.txt` File

| Column | Description |
|--------|-------------|
| `Func.refGene` | Function (exonic, intronic, UTR, etc.) |
| `Gene.refGene` | Gene name |
| `ExonicFunc.refGene` | Type (nonsynonymous, frameshift, etc.) |
| `AAChange.refGene` | Protein change (e.g., p.V600E) |
| `avsnp150` | dbSNP ID (rs number) |
| `gnomAD_exome_ALL` | Population frequency |
| `CLNSIG` | ClinVar significance |
| `SIFT_score` | Deleteriousness (< 0.05 = damaging) |
| `Polyphen2_HDIV_score` | Pathogenicity (> 0.85 = damaging) |

---

## 🎯 Find Clinically Relevant Variants

### Filter for Pathogenic Variants

```bash
# In Excel or command line:
grep "Pathogenic\|Likely_pathogenic" output.hg19_multianno.txt > pathogenic_variants.txt
```

### Filter for Rare Damaging Variants

```bash
# Frequency < 1%, predicted damaging
awk -F'\t' '$NF < 0.01 && $13 == "nonsynonymous SNV"' output.hg19_multianno.txt
```

---

## 📥 Download Additional Databases

```bash
cd ~/NGS/tools/annovar

# Download database
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar DATABASE_NAME humandb/

# Examples:
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
```

---

## 🔧 Troubleshooting

### Test Installation

```bash
bash test_annovar.sh
```

### Check Databases

```bash
ls -lh ~/NGS/tools/annovar/humandb/
```

### Common Issues

**Issue**: "Cannot open database"  
**Solution**: Download the database
```bash
cd ~/NGS/tools/annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
```

**Issue**: "Out of memory"  
**Solution**: Split VCF by chromosome
```bash
bcftools view -r chr1 input.vcf > chr1.vcf
bash annovar_helper.sh chr1.vcf chr1_annotated
```

---

## 📚 Available Databases

### Gene Annotations (operation: g)
- `refGene` - RefSeq genes ⭐ **Recommended**
- `knownGene` - UCSC genes
- `ensGene` - Ensembl genes

### Variant Frequencies (operation: f)
- `gnomad211_exome` - gnomAD v2.1.1 ⭐ **Recommended**
- `gnomad30_genome` - gnomAD v3.0
- `avsnp150` - dbSNP 150 ⭐ **Recommended**
- `exac03` - ExAC
- `1000g2015aug` - 1000 Genomes

### Clinical Databases (operation: f)
- `clinvar_20240917` - ClinVar ⭐ **Recommended**
- `cosmic70` - COSMIC mutations
- `icgc28` - ICGC cancer

### Functional Predictions (operation: f)
- `dbnsfp42a` - SIFT, PolyPhen, etc. ⭐ **Recommended**

---

## 💡 Tips

### Speed Up Annotation

Use fewer databases:
```bash
-protocol refGene,clinvar_20240917
```

### Compress Output

```bash
bgzip output.hg19_multianno.vcf
tabix -p vcf output.hg19_multianno.vcf.gz
```

### View in Excel

Open `output.hg19_multianno.txt` directly in Excel for easy filtering

---

## 🔗 Resources

- **Documentation**: [ANNOVAR_INTEGRATION.md](docs/ANNOVAR_INTEGRATION.md)
- **Official Site**: http://annovar.openbioinformatics.org/
- **Database List**: http://annovar.openbioinformatics.org/en/latest/user-guide/download/

---

## ⏱️ Typical Runtimes

| VCF Size | Runtime | Memory |
|----------|---------|--------|
| 1,000 variants | 30 sec | 500 MB |
| 10,000 variants | 2 min | 1 GB |
| 100,000 variants | 10 min | 2 GB |
| 1M variants | 1 hour | 4 GB |

---

**Quick Help**: For complete documentation, see [ANNOVAR_INTEGRATION.md](docs/ANNOVAR_INTEGRATION.md)


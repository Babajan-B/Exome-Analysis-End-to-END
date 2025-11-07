# ANNOVAR Integration Guide

Complete guide for using ANNOVAR with the NGS Exome Analysis Pipeline.

## 📋 Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Available Databases](#available-databases)
4. [Usage Examples](#usage-examples)
5. [Integration with Pipeline](#integration-with-pipeline)
6. [Troubleshooting](#troubleshooting)

---

## 🔧 Installation

### Option 1: Bundled installer (recommended)

```bash
cd ~/NGS
bash install_all.sh
```

This downloads ANNOVAR to `~/NGS/tools/annovar/` and populates the default hg19 databases used by the pipeline (refGene, ClinVar, gnomAD, dbSNP, dbNSFP, COSMIC, ICGC, etc.).

**Storage Required**: ~2–5 GB for databases  
**Time Required**: 15–30 minutes (depends on network speed)

### Option 2: Manual installation

```bash
cd ~/NGS/tools
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
mkdir -p annovar && tar -xzf annovar.latest.tar.gz --strip-components=1 -C annovar
rm annovar.latest.tar.gz
cd annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
# repeat annotate_variation to download any additional databases you need
```

---

## 🚀 Quick Start

### Basic command-line usage

```bash
cd ~/NGS/tools/annovar

perl table_annovar.pl     /path/to/variants.vcf     humandb/     -buildver hg19     -out annotated     -remove     -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150,dbnsfp42a     -operation g,f,f,f,f     -nastring .     -vcfinput     -polish
```

The bundled pipeline runs a similar command automatically on every sample after variant filtering.

---

## 📚 Available Databases

### Gene-Based Annotations (operation: g)

| Database | Description | Size |
|----------|-------------|------|
| `refGene` | RefSeq genes | ~50 MB |
| `knownGene` | UCSC known genes | ~60 MB |
| `ensGene` | Ensembl genes | ~50 MB |
| `wgEncodeGencodeBasicV19` | GENCODE genes | ~40 MB |

### Filter-Based Annotations (operation: f)

| Database | Description | Size |
|----------|-------------|------|
| `avsnp150` | dbSNP 150 | ~1 GB |
| `gnomad211_exome` | gnomAD v2.1.1 exome | ~500 MB |
| `gnomad30_genome` | gnomAD v3.0 genome | ~800 MB |
| `clinvar_20240917` | ClinVar pathogenic | ~100 MB |
| `dbnsfp42a` | Functional predictions | ~3 GB |
| `cosmic70` | COSMIC mutations | ~200 MB |
| `icgc28` | ICGC somatic | ~150 MB |
| `exac03` | ExAC database | ~400 MB |
| `esp6500siv2_all` | NHLBI ESP | ~50 MB |
| `1000g2015aug_all` | 1000 Genomes | ~200 MB |

### Region-Based Annotations (operation: r)

| Database | Description | Size |
|----------|-------------|------|
| `cytoBand` | Chromosome bands | ~1 MB |
| `genomicSuperDups` | Segmental duplications | ~5 MB |
| `dgvMerged` | Database of Genomic Variants | ~10 MB |
| `gwas` | GWAS catalog | ~20 MB |

---

## 💡 Usage Examples

### Example 1: Basic Gene Annotation

```bash
cd ~/NGS/tools/annovar

# Annotate with RefSeq genes only
perl table_annovar.pl \
    variants.vcf \
    humandb/ \
    -buildver hg19 \
    -out annotated_basic \
    -remove \
    -protocol refGene \
    -operation g \
    -vcfinput
```

### Example 2: Clinical Annotation

```bash
# Focus on clinically relevant databases
perl table_annovar.pl \
    variants.vcf \
    humandb/ \
    -buildver hg19 \
    -out annotated_clinical \
    -remove \
    -protocol refGene,clinvar_20240917,gnomad211_exome,avsnp150 \
    -operation g,f,f,f \
    -nastring . \
    -vcfinput
```

### Example 3: Comprehensive Research Annotation

```bash
# Full annotation with all databases
perl table_annovar.pl \
    variants.vcf \
    humandb/ \
    -buildver hg19 \
    -out annotated_full \
    -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,gnomad211_exome,gnomad30_genome,clinvar_20240917,dbnsfp42a,cosmic70,icgc28,exac03,esp6500siv2_all,1000g2015aug_all \
    -operation g,g,g,f,f,f,f,f,f,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
```

### Example 4: Cancer-Specific Annotation

```bash
# Annotate with cancer databases
perl table_annovar.pl \
    somatic_variants.vcf \
    humandb/ \
    -buildver hg19 \
    -out annotated_cancer \
    -remove \
    -protocol refGene,cosmic70,icgc28,clinvar_20240917,gnomad211_exome \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput
```

### Example 5: Convert VCF to ANNOVAR Format

```bash
# If you need avinput format
perl convert2annovar.pl \
    -format vcf4 \
    variants.vcf \
    > variants.avinput

# Then annotate
perl annotate_variation.pl \
    -geneanno \
    -buildver hg19 \
    variants.avinput \
    humandb/
```

---

## 🔗 Integration with Pipeline

### Option A: Standalone Annotation (Recommended for Now)

After running the main pipeline:

```bash
# 1. Run the main NGS pipeline
bash run_pipeline.sh data/R1.fastq.gz data/R2.fastq.gz sample_name 16

# 2. Wait for completion, then annotate with ANNOVAR
bash annovar_helper.sh \
    results/sample_name/filtered/filtered_variants.vcf \
    results/sample_name/annovar/annotated
```

### Option B: Manual Integration Point

In the pipeline, ANNOVAR would run after variant filtering:

```bash
# In run_pipeline.sh, add after Step 9 (Variant Filtering):

# STEP 10: Variant Annotation with ANNOVAR
step 10 "Variant Annotation (ANNOVAR)"
ANNOVAR_DIR=$WORK_DIR/tools/annovar
perl $ANNOVAR_DIR/table_annovar.pl \
    $OUTPUT_DIR/filtered/filtered_variants.vcf \
    $ANNOVAR_DIR/humandb/ \
    -buildver hg19 \
    -out $OUTPUT_DIR/annotated/annovar \
    -remove \
    -protocol refGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
echo "✅ ANNOVAR annotation complete"
```

### Option C: Add to Python Pipeline

To integrate into `app/pipeline/pipeline.py`, add this function:

```python
def annotate_with_annovar(input_vcf, output_prefix, buildver='hg19'):
    """Annotate variants using ANNOVAR."""
    annovar_dir = os.path.join(os.getcwd(), 'tools', 'annovar')
    
    cmd = f"""
    perl {annovar_dir}/table_annovar.pl \
        {input_vcf} \
        {annovar_dir}/humandb/ \
        -buildver {buildver} \
        -out {output_prefix} \
        -remove \
        -protocol refGene,avsnp150,gnomad211_exome,clinvar_20240917,dbnsfp42a \
        -operation g,f,f,f,f \
        -nastring . \
        -vcfinput \
        -polish
    """
    
    logger.info(f"Running ANNOVAR: {cmd}")
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if process.returncode != 0:
        logger.error(f"ANNOVAR failed: {process.stderr}")
        return False
    
    logger.info(f"ANNOVAR completed successfully")
    return True
```

---

## 🔍 Understanding ANNOVAR Output

### Output Files

1. **`.hg19_multianno.vcf`** - Annotated VCF with all INFO fields added
2. **`.hg19_multianno.txt`** - Tab-delimited table (easier to read)
3. **`.avinput`** - Intermediate ANNOVAR input format

### Key Annotation Fields (in VCF INFO)

- `Func.refGene` - Functional region (exonic, intronic, UTR, etc.)
- `Gene.refGene` - Gene name
- `ExonicFunc.refGene` - Type of exonic mutation (nonsynonymous, synonymous, etc.)
- `AAChange.refGene` - Amino acid change
- `avsnp150` - dbSNP rs ID
- `gnomAD_exome_ALL` - Population allele frequency
- `CLNSIG` - ClinVar clinical significance
- `SIFT_score` - SIFT prediction score
- `Polyphen2_HDIV_score` - PolyPhen-2 score

### Interpreting Results

**High-priority variants:**
- `CLNSIG=Pathogenic` or `Likely_pathogenic`
- `ExonicFunc=nonsynonymous_SNV` with low gnomAD frequency (<0.01)
- `SIFT_score<0.05` (damaging) and `Polyphen2>0.85` (probably damaging)

---

## 📥 Downloading Additional Databases

### For hg19

```bash
cd ~/NGS/tools/annovar

# Download a specific database
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <database_name> humandb/

# Examples:
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
```

### For hg38

```bash
# First time setup for hg38
cd ~/NGS/tools/annovar

# Download hg38 databases
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome humandb/
```

---

## 🔧 Troubleshooting

### Issue 1: Database Download Fails

```bash
# Try alternative download method
cd ~/NGS/tools/annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/

# If still fails, check ANNOVAR website for manual download links
```

### Issue 2: Perl Module Missing

```bash
# Install required Perl modules
sudo apt-get install -y libdbi-perl libdbd-mysql-perl
```

### Issue 3: "Cannot open database" Error

```bash
# Check if databases are downloaded
ls -lh ~/NGS/tools/annovar/humandb/

# Re-download missing database
cd ~/NGS/tools/annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <missing_db> humandb/
```

### Issue 4: Out of Memory

```bash
# For very large VCF files, split and process in chunks
bcftools view -r chr1 input.vcf > chr1.vcf
# Annotate each chromosome separately
bash annovar_helper.sh chr1.vcf chr1_annotated
```

---

## 📊 Performance Tips

### Speed Optimization

1. **Use fewer databases** for faster annotation:
   ```bash
   -protocol refGene,avsnp150,clinvar_20240917
   ```

2. **Parallelize by chromosome**:
   ```bash
   for chr in {1..22} X Y; do
       bcftools view -r chr$chr input.vcf | \
       bash annovar_helper.sh - output_chr${chr} &
   done
   wait
   ```

3. **Skip intermediate files**:
   ```bash
   -remove  # Removes avinput after completion
   ```

### Storage Optimization

```bash
# Compress output files
bgzip annotated.hg19_multianno.vcf
tabix -p vcf annotated.hg19_multianno.vcf.gz

# Remove intermediate files
rm *.avinput *.refGene *.log
```

---

## 📚 Additional Resources

- **ANNOVAR Website**: http://annovar.openbioinformatics.org/
- **Documentation**: http://annovar.openbioinformatics.org/en/latest/
- **Database List**: http://annovar.openbioinformatics.org/en/latest/user-guide/download/
- **Forum**: https://groups.google.com/g/annovar

---

## 🎯 Quick Reference

```bash
# Install ANNOVAR
bash setup_annovar.sh

# Annotate VCF (quick)
bash annovar_helper.sh input.vcf output

# Annotate VCF (custom)
perl ~/NGS/tools/annovar/table_annovar.pl \
    input.vcf humandb/ \
    -buildver hg19 \
    -out output \
    -protocol refGene,avsnp150 \
    -operation g,f \
    -vcfinput

# Download additional database
cd ~/NGS/tools/annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <db_name> humandb/

# List available databases
ls -lh ~/NGS/tools/annovar/humandb/
```

---

**Made with 🧬 for genomic variant annotation**


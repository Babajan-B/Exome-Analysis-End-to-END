# Methods

This document summarizes the computational workflow implemented by the **NGS Exome Analysis Pipeline** and highlights the major tools, databases, and decision points used for each phase. The goal is to provide a clear methodological record that can be referenced in manuscripts, SOPs, or collaborative work.

---

## Overview Flowchart

```mermaid
flowchart TD
    A[Raw FASTQ files] --> B[Quality Control<br/>(FastQC)]
    B --> C[Read Trimming<br/>(fastp)]
    C --> D[Alignment to hg19<br/>(BWA-MEM)]
    D --> E[SAM→BAM + Sorting<br/>(SAMtools)]
    E --> F[Mark Duplicates<br/>(GATK MarkDuplicates)]
    F --> G[Variant Calling<br/>(GATK HaplotypeCaller)]
    G --> H[Variant Filtering<br/>(GATK VariantFiltration)]
    H --> I[PASS-only VCF Extraction]
    I --> J[ANNOVAR Annotation
            • refGene
            • ClinVar
            • gnomAD
            • dbSNP (avsnp150)
            • dbNSFP (predictions)]
    I --> K[snpEff Annotation
            • Functional Effects
            • Impact Classification]
    J --> L[Zygosity Extraction
            (VCF GT field)]
    L --> M[Functional Classification
            • SNP vs Indel
            • Exonic vs Non-exonic
            • Nonsyn/ Syn/ Stopgain/ Frameshift]
    K --> N[snpEff Reports
            (HTML + CSV)]
    M --> O[Result Packaging
            – Annotated TXT/VCF
            – QC Reports (FastQC/fastp)
            – Summary ZIP archive]
    N --> O
```

---

## Detailed Methodology

### 1. Input Acquisition
- **Data type:** Paired-end FASTQ files (gzipped)
- **Naming pattern:** `Sample_R1.fastq.gz` and `Sample_R2.fastq.gz` (or `_1` / `_2` equivalents)
- **Storage:** Files uploaded or copied into `~/NGS/data/`

### 2. Quality Assessment (FastQC)
- **Tool:** [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **Purpose:** Evaluate raw read quality, GC distribution, adapter contamination, overrepresented sequences.
- **Output:** HTML reports stored under `results/<sample>/fastqc/`

### 3. Read Trimming and Filtering (fastp)
- **Tool:** [fastp](https://github.com/OpenGene/fastp)
- **Parameters:**
  - Adapter auto-detection
  - Minimum read length: 50 bp
  - Quality filtering threshold: Q20
  - Thread count: matches pipeline configuration
- **Output:** Trimmed FASTQ files + HTML/JSON summary under `results/<sample>/trimmed/`

### 4. Alignment (BWA-MEM)
- **Tool:** [BWA-MEM](https://github.com/lh3/bwa)
- **Reference genome:** hg19 (GRCh37) located in `~/NGS/reference/`
- **Read group tags:** ID, SM, PL, LB, PU automatically populated per sample
- **Output:** Aligned SAM file converted to BAM in subsequent steps

### 5. Post-alignment Processing (SAMtools + GATK)
1. **SAM→BAM conversion & sorting:** `samtools view` and `samtools sort`
2. **Indexing:** `samtools index`
3. **Duplicate marking:** `gatk MarkDuplicates` (creates deduplicated BAM + metrics)

### 6. Variant Calling (GATK HaplotypeCaller)
- **Tool:** [GATK](https://gatk.broadinstitute.org) v4.6.2.0
- **Mode:** Germline variant discovery in single-sample BAM
- **Output:** `raw_variants.vcf`

### 7. Variant Filtering (GATK VariantFiltration)
- **Filters applied:**
  - QD < 2.0 (quality by depth)
  - QUAL < 30.0
  - MQ < 40.0 (mapping quality)
  - FS > 60.0 (Fisher strand bias)
  - SOR > 3.0 (strand odds ratio)
- **Result:** `filtered_variants.vcf`

### 8. PASS-only Variant Extraction
- Extract header lines plus variants labeled `PASS` into `filtered_PASS_only.vcf`
- **Purpose:** Reduce annotation workload and output size while focusing on high-confidence calls

### 9. Annotation – ANNOVAR
- **Tool:** [ANNOVAR](http://annovar.openbioinformatics.org)
- **Script:** `table_annovar.pl`
- **Databases (hg19 build):**
  - `refGene` – Gene-based annotation
  - `clinvar_20240917` – Clinical significance
  - `gnomad211_exome` – Population frequencies
  - `avsnp150` – dbSNP identifiers
  - `dbnsfp42a` – Functional predictions (SIFT, PolyPhen, CADD, etc.)
- **Flags:** `-vcfinput`, `-polish`, `-remove`, `-nastring .`
- **Outputs:**
  - `annotated_<sample>.hg19_multianno.txt`
  - `annotated_<sample>.hg19_multianno.vcf`

### 10. Annotation – snpEff
- **Tool:** [snpEff](http://pcingola.github.io/SnpEff/), GRCh37.75 database
- **Command:** `snpEff.jar -stats ... filtered_PASS_only.vcf`
- **Outputs:**
  - Annotated VCF `*_snpEff_annotated.vcf`
  - HTML summary `*_snpEff_summary.html`
  - CSV statistics `*_snpEff_summary.csv`

### 11. Zygosity Extraction
- Python helper parses genotype (`GT`) fields from PASS VCF and appends a **Zygosity** column to the multianno TXT file:
  - `0/1` or `1/0` → *Heterozygous*
  - `1/1` → *Homozygous*
  - `0/0` → *Reference*
  - Other → *Unknown*

### 12. Functional Classification
- **Scripts:** AWK-based filters split the annotated table into:
  - SNPs vs Insertions vs Deletions
  - Exonic vs Non-exonic
  - Exonic subsets: Nonsynonymous, Synonymous, Stopgain, Frameshift
- **Output location:** `results/<sample>/annovar/functional_classification/`

### 13. Packaging and Summaries
- Compress and tabix-index the annotated VCFs (`bgzip`, `tabix`)
- Generate `MASTER_ANALYSIS_SUMMARY.txt` with per-sample statistics
- Bundle essential outputs into `NGS_Results_Complete_<timestamp>.zip` including:
  - ANNOVAR TXT/VCF (with zygosity)
  - snpEff annotated VCF + reports
  - Functional classification tables
  - FastQC & fastp HTML reports
  - Summary file

### 14. Automation & Reproducibility
- Entire workflow orchestrated by `ULTIMATE_MASTER_PIPELINE.sh`
- Installation script (`install_all.sh`) provisions:
  - Bioinformatics tools (FastQC, fastp, BWA, SAMtools, GATK, snpEff)
  - Reference genome and indices
  - ANNOVAR + default databases
- Verification script (`CHECK_INSTALLATION.sh`) confirms availability of executables, databases, and reference files.

---

## Key Software Versions
| Component | Version | Source |
|-----------|---------|--------|
| FastQC | 0.11.9 | Debian package / install_all.sh |
| fastp | 1.0.1 | Pre-built binary from OpenGene |
| BWA | 0.7.17 | Debian package |
| SAMtools | 1.13 | Debian package |
| GATK | 4.6.2.0 | Broad Institute release |
| ANNOVAR | latest (as downloaded) | OpenBioinformatics |
| snpEff | 5.3a (GRCh37.75 DB) | snpEff project |
| Python | 3.8+ | System supplied |
| Java | OpenJDK 21 | `openjdk-21-jdk` |

---

## Reproducibility Notes
- **Thread count:** Default to 16 but configurable when launching the pipeline.
- **Reference genome:** hg19; swap indices and annotation databases to use hg38 if required.
- **Temporary files:** Large intermediate BAM/FASTQ files remain outside the final ZIP to save space; re-run pipeline to regenerate if needed.
- **Logs:** Each sample has a cumulative log at `results/<sample>/pipeline.log`.
- **Version control:** All scripts maintained in the GitHub repository [Babajan-B/Exome-Analysis-End-to-END](https://github.com/Babajan-B/Exome-Analysis-End-to-END).

---

*Prepared automatically to document the pipeline methodology. Update this file if new tools, databases, or filtration criteria are introduced.*

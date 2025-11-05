# Complete Analysis Plan - ANNOVAR with Zygosity & Advanced Separation

## 🎯 **Your Requirements:**

1. ✅ Run ANNOVAR with prediction databases (dbnsfp42a)
2. ✅ Add zygosity information (GT field from VCF)
3. ✅ Separate by variant type (SNPs, Indels)
4. ✅ Within SNPs, separate by ExonicFunc.refGene:
   - **Exonic SNPs** (has value like "nonsynonymous SNV", "synonymous SNV", etc.)
   - **Non-Exonic SNPs** (has "." - intronic, intergenic, UTR, etc.)

---

## 📋 **Analysis Workflow:**

```
Input: filtered_variants.vcf
   ↓
Step 1: Extract PASS + Add Zygosity
   ↓
Step 2: ANNOVAR Annotation (5 databases with predictions)
   ↓
Step 3: Separate by Variant Type
   ├── SNPs
   ├── Insertions
   ├── Deletions
   └── MNPs
   ↓
Step 4: Within SNPs, Separate by ExonicFunc
   ├── SNPs_Exonic (has ExonicFunc value)
   └── SNPs_NonExonic (ExonicFunc = ".")
   ↓
Output: Multiple categorized files ready for analysis
```

---

## 🔬 **Understanding Zygosity:**

### **What is Zygosity?**

| Genotype (GT) | Zygosity | Meaning |
|---------------|----------|---------|
| 0/1 or 0\|1 | **Heterozygous** | One copy of variant |
| 1/1 or 1\|1 | **Homozygous** | Both copies have variant |
| 0/0 or 0\|0 | Reference | No variant (should not appear) |
| ./. | Unknown | No call |

### **Why Important?**

- **Heterozygous**: Carrier state, dominant diseases
- **Homozygous**: Recessive diseases, stronger effect
- **Clinical**: Essential for inheritance pattern analysis

---

## 📊 **ExonicFunc Categories:**

### **Exonic (Functional Impact):**
- `nonsynonymous SNV` - Changes amino acid ⭐ Most important
- `synonymous SNV` - Silent mutation
- `stopgain` - Creates stop codon (truncating)
- `stoploss` - Removes stop codon
- `frameshift insertion/deletion` - Shifts reading frame
- `nonframeshift insertion/deletion` - In-frame change

### **Non-Exonic (Usually Less Impact):**
- `.` (blank) - Can be:
  - Intronic (within gene but not coding)
  - Intergenic (between genes)
  - UTR (untranslated region)
  - Upstream/Downstream (near gene)
  - ncRNA (non-coding RNA)
  - Splicing (splice sites)

---

## 🎯 **File Organization Plan:**

```
results/SAMPLE/annovar/
├── annotated_with_zygosity.hg19_multianno.txt (full file)
│
├── by_type/
│   ├── SNPs_all.txt
│   ├── Insertions_all.txt
│   ├── Deletions_all.txt
│   └── MNPs_all.txt
│
└── by_function/
    ├── SNPs_Exonic.txt          (ExonicFunc has value)
    │   ├── nonsynonymous.txt    (amino acid changes)
    │   ├── synonymous.txt       (silent)
    │   ├── stopgain.txt         (truncating)
    │   └── frameshift.txt       (frameshifts)
    │
    └── SNPs_NonExonic.txt       (ExonicFunc = ".")
        ├── intronic.txt
        ├── UTR.txt
        └── intergenic.txt
```

---

## 📝 **Script Breakdown:**

### **Script 1: annotate_with_zygosity.sh**
- Extracts PASS variants
- Adds zygosity column from GT field
- Runs ANNOVAR with 5 databases
- Output: Full annotated file WITH zygosity

### **Script 2: separate_by_type_advanced.sh**
- Separates into SNPs, Insertions, Deletions, MNPs
- Each file includes zygosity info

### **Script 3: separate_snps_by_function.sh**
- Takes SNPs file
- Separates by ExonicFunc.refGene:
  - Exonic (has value)
  - Non-Exonic (has ".")
- Further categorizes exonic by type

### **Script 4: master_analysis_complete.sh**
- Runs all above in sequence
- Generates comprehensive report
- Creates organized folder structure

---

## 🎓 **Clinical Interpretation Guide:**

### **Priority 1: Exonic SNPs**
```
File: SNPs_Exonic.txt
Focus on: nonsynonymous SNV
Filter by:
  • CLNSIG = "Pathogenic"
  • SIFT_score < 0.05
  • Polyphen2 > 0.85
  • Zygosity = Homozygous (for recessive)
```

### **Priority 2: Indels**
```
Files: Insertions_all.txt, Deletions_all.txt
Focus on: frameshift variants
Filter by:
  • ExonicFunc = "frameshift"
  • gnomAD < 0.001 (very rare)
  • Zygosity
```

### **Priority 3: Non-Exonic SNPs**
```
File: SNPs_NonExonic.txt
Focus on: splicing variants
Look for: Func.refGene contains "splicing"
```

---

## 📦 **Final Zip Strategy:**

```
SAMPLE_COMPLETE.zip contains:
├── Full_Annotated_with_Zygosity.txt
├── by_type/
│   ├── SNPs_all.txt
│   ├── Insertions_all.txt
│   └── Deletions_all.txt
└── by_function/
    ├── SNPs_Exonic.txt (⭐ Most important)
    └── SNPs_NonExonic.txt
```

---

## ⏱️ **Estimated Sizes:**

| File | Approx Size |
|------|-------------|
| Full annotated | 30-50 MB |
| SNPs (all) | 20-35 MB |
| SNPs Exonic | 8-15 MB |
| SNPs Non-Exonic | 10-20 MB |
| Insertions | 3-5 MB |
| Deletions | 3-5 MB |
| **Total ZIP** | **50-80 MB** |

---

## 🚀 **Next Steps:**

I'll create 4 scripts:
1. `annotate_with_zygosity.sh` - Add zygosity + annotate
2. `separate_by_type_advanced.sh` - Separate with zygosity
3. `separate_snps_by_function.sh` - SNPs by ExonicFunc
4. `run_complete_categorization.sh` - Master script for all

**Shall I proceed to create these scripts?** ✅


#!/usr/bin/env python3
"""
Generate comprehensive tab-delimited annotation tables and functional classifications
from normalized snpEff & SnpSift (ClinVar + dbSNP) annotated VCF files.
Produces:
- annovar/annotated_<sample>.txt
- annovar/annotated_<sample>_with_zygosity.txt
- annovar/separated_by_type/{SNPs.txt, Insertions.txt, Deletions.txt}
- annovar/functional_classification/{Exonic_Nonsynonymous.txt, Exonic_Synonymous.txt, Exonic_Stopgain.txt, Exonic_Frameshift.txt, Splicing.txt, Requires_Sanger_Validation.txt}
"""

import sys
import os
import gzip

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def extract_genotype_info(parts):
    if len(parts) <= 9:
        return "Unknown", 0, 0, 1.0, "LOW_CONFIDENCE_GENOTYPE"
    format_keys = parts[8].split(":")
    sample_vals = parts[9].split(":")
    fmt_dict = dict(zip(format_keys, sample_vals))
    
    gt = fmt_dict.get("GT", "./.")
    dp = int(fmt_dict.get("DP", 0)) if fmt_dict.get("DP", "0").isdigit() else 0
    gq = int(fmt_dict.get("GQ", 0)) if fmt_dict.get("GQ", "0").isdigit() else 0
    ad = fmt_dict.get("AD", "0,0").split(",")
    ad0 = int(ad[0]) if len(ad) > 0 and ad[0].isdigit() else 0
    ad1 = int(ad[1]) if len(ad) > 1 and ad[1].isdigit() else 0
    
    total_informative = ad0 + ad1
    ab = (ad1 / total_informative) if total_informative > 0 else 0.0
    
    if gt in ["0/1", "1/0"]:
        zygosity = "Heterozygous"
    elif gt == "1/1":
        zygosity = "Homozygous"
    elif gt in ["1", "1/."]:
        zygosity = "Hemizygous"
    elif gt == "0/0":
        zygosity = "Reference"
    else:
        zygosity = "Unknown"
        
    # Genotype QC evaluation (qc.json):
    # DP >= 10, GQ >= 20, Het AB in [0.20, 0.80], Hom AB >= 0.90
    qc_pass = (dp >= 10 and gq >= 20)
    if zygosity == "Heterozygous":
        qc_pass = qc_pass and (0.20 <= ab <= 0.80)
    elif zygosity == "Homozygous":
        qc_pass = qc_pass and (ab >= 0.90)
        
    confidence = "HIGH_CONFIDENCE" if qc_pass else "LOW_CONFIDENCE_GENOTYPE"
    return zygosity, dp, gq, ab, confidence

def main():
    if len(sys.argv) < 3:
        print("Usage: generate_annotation_table.py <input_vcf> <output_dir> [sample_name]")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_dir = sys.argv[2]
    sample_name = sys.argv[3] if len(sys.argv) > 3 else "Sample"

    annovar_dir = os.path.join(output_dir, "annovar")
    sep_dir = os.path.join(annovar_dir, "separated_by_type")
    func_dir = os.path.join(annovar_dir, "functional_classification")

    os.makedirs(annovar_dir, exist_ok=True)
    os.makedirs(sep_dir, exist_ok=True)
    os.makedirs(func_dir, exist_ok=True)

    out_annot_txt = os.path.join(annovar_dir, f"annotated_{sample_name}.txt")
    out_annot_zyg = os.path.join(annovar_dir, f"annotated_{sample_name}_with_zygosity.txt")

    headers = [
        "Chr", "Start", "End", "Ref", "Alt", 
        "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene",
        "avsnp150", "gnomAD_exome_ALL", "Impact", "Consequence", "Transcript", "Zygosity",
        "ClinVar_Significance", "ClinVar_RevStat", "ClinVar_Disease", "Genotype_Confidence"
    ]

    header_str = "\t".join(headers) + "\n"

    # Category file handles
    f_snps = open(os.path.join(sep_dir, "SNPs.txt"), "w")
    f_ins = open(os.path.join(sep_dir, "Insertions.txt"), "w")
    f_del = open(os.path.join(sep_dir, "Deletions.txt"), "w")

    f_nonsyn = open(os.path.join(func_dir, "Exonic_Nonsynonymous.txt"), "w")
    f_syn = open(os.path.join(func_dir, "Exonic_Synonymous.txt"), "w")
    f_stop = open(os.path.join(func_dir, "Exonic_Stopgain.txt"), "w")
    f_frameshift = open(os.path.join(func_dir, "Exonic_Frameshift.txt"), "w")
    f_splicing = open(os.path.join(func_dir, "Splicing.txt"), "w")
    f_sanger = open(os.path.join(func_dir, "Requires_Sanger_Validation.txt"), "w")

    for f in [f_snps, f_ins, f_del, f_nonsyn, f_syn, f_stop, f_frameshift, f_splicing, f_sanger]:
        f.write(header_str)

    rows_zyg = []
    rows_no_zyg = []

    snp_count = 0
    ins_count = 0
    del_count = 0
    nonsyn_count = 0
    syn_count = 0
    stop_count = 0
    fs_count = 0
    sanger_count = 0
    clinvar_match_count = 0

    with open_vcf(input_vcf) as vf:
        for line in vf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0]
            start = parts[1]
            ref = parts[3]
            alt = parts[4]
            info = parts[7]

            # Calculate end position for deletion or SNV
            end = str(int(start) + len(ref) - 1) if len(ref) > 1 else start

            # Extract genotype and confidence
            zygosity, dp, gq, ab, confidence = extract_genotype_info(parts)

            # rsID
            avsnp150 = parts[2] if parts[2] != "." else "."

            # Population frequency from CAF or gnomAD
            gnomad = "."
            clnsig = "."
            clnrevstat = "."
            clndn = "."

            info_items = info.split(";")
            for item in info_items:
                if item.startswith("CAF="):
                    caf_vals = item[4:].split(",")
                    if len(caf_vals) > 1 and caf_vals[1] != ".":
                        gnomad = caf_vals[1]
                elif item.startswith("CLNSIG="):
                    clnsig = item[7:]
                elif item.startswith("CLNREVSTAT="):
                    clnrevstat = item[11:]
                elif item.startswith("CLNDN="):
                    clndn = item[6:]

            if clnsig != ".":
                clinvar_match_count += 1

            # Default consequence fields
            func_region = "intergenic"
            gene_name = "."
            exonic_func = "."
            aa_change = "."
            impact = "."
            consequence = "."
            feature_id = "."

            for item in info_items:
                if item.startswith("ANN="):
                    ann_entries = item[4:].split(",")
                    if ann_entries:
                        first_ann = ann_entries[0].split("|")
                        if len(first_ann) > 10:
                            consequence = first_ann[1] or "."
                            impact = first_ann[2] or "."
                            gene_name = first_ann[3] or "."
                            feature_id = first_ann[6] or "."
                            aa_change = first_ann[10] or "."

                            # Map consequence to region & exonic function
                            cq_lower = consequence.lower()
                            if "missense" in cq_lower:
                                func_region = "exonic"
                                exonic_func = "nonsynonymous SNV"
                            elif "synonymous" in cq_lower:
                                func_region = "exonic"
                                exonic_func = "synonymous SNV"
                            elif "stop_gained" in cq_lower:
                                func_region = "exonic"
                                exonic_func = "stopgain"
                            elif "frameshift" in cq_lower:
                                func_region = "exonic"
                                exonic_func = "frameshift insertion" if len(ref) < len(alt) else "frameshift deletion"
                            elif "splice" in cq_lower:
                                func_region = "splicing"
                            elif "intron" in cq_lower:
                                func_region = "intronic"
                            elif "utr" in cq_lower:
                                func_region = "UTR"
                            elif "upstream" in cq_lower or "downstream" in cq_lower:
                                func_region = "intergenic"
                    break

            row_data = [
                chrom, start, end, ref, alt,
                func_region, gene_name, ".", exonic_func, aa_change,
                avsnp150, gnomad, impact, consequence, feature_id, zygosity,
                clnsig, clnrevstat, clndn, confidence
            ]
            row_str = "\t".join(row_data) + "\n"
            rows_zyg.append(row_str)
            rows_no_zyg.append("\t".join(row_data[:-1]) + "\n")

            # Separation by variant type
            if len(ref) == 1 and len(alt) == 1:
                f_snps.write(row_str)
                snp_count += 1
            elif len(ref) < len(alt):
                f_ins.write(row_str)
                ins_count += 1
            elif len(ref) > len(alt):
                f_del.write(row_str)
                del_count += 1

            # Functional separation
            ex_lower = exonic_func.lower()
            if "nonsynonymous" in ex_lower:
                f_nonsyn.write(row_str)
                nonsyn_count += 1
            elif "synonymous" in ex_lower:
                f_syn.write(row_str)
                syn_count += 1
            elif "stopgain" in ex_lower:
                f_stop.write(row_str)
                stop_count += 1
            elif "frameshift" in ex_lower:
                f_frameshift.write(row_str)
                fs_count += 1

            if func_region == "splicing":
                f_splicing.write(row_str)

            # Check if this variant requires Stage 9 Sanger validation
            # (ClinVar P/LP or HIGH impact, but LOW_CONFIDENCE_GENOTYPE)
            is_plp = ("pathogenic" in clnsig.lower() or "likely_pathogenic" in clnsig.lower())
            is_high = (impact == "HIGH" or "stopgain" in ex_lower or "frameshift" in ex_lower)
            if (is_plp or is_high) and confidence == "LOW_CONFIDENCE_GENOTYPE":
                f_sanger.write(row_str)
                sanger_count += 1

    # Close handles
    for f in [f_snps, f_ins, f_del, f_nonsyn, f_syn, f_stop, f_frameshift, f_splicing, f_sanger]:
        f.close()

    # Write main tables
    with open(out_annot_zyg, "w") as f:
        f.write(header_str)
        f.writelines(rows_zyg)

    with open(out_annot_txt, "w") as f:
        f.write("\t".join(headers[:-1]) + "\n")
        f.writelines(rows_no_zyg)

    print(f"✅ Generated {len(rows_zyg)} annotated variants in {annovar_dir}")
    print(f"   - ClinVar Matches: {clinvar_match_count}")
    print(f"   - SNPs: {snp_count}, Insertions: {ins_count}, Deletions: {del_count}")
    print(f"   - Exonic Nonsynonymous: {nonsyn_count}, Synonymous: {syn_count}, Stopgain: {stop_count}, Frameshift: {fs_count}")
    print(f"   - Requiring Stage 9 Sanger Validation: {sanger_count}")

if __name__ == "__main__":
    main()

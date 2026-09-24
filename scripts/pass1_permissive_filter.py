#!/usr/bin/env python3
"""
Pass 1: Local Whole-Callset Permissive Filter & Genotype Confidence Tagging

Implements the Permissive Funnel (Version 2.5.0):
1. Precedence Rule 1 (Clinical Exemption):
   - ClinVar Pathogenic / Likely Pathogenic (P/LP) is retained UNCONDITIONALLY,
     even if common (e.g. HFE C282Y AF ~ 6%), regardless of consequence.
2. Precedence Rule 2 (Coding & Splice Permissive Net):
   - Variants with HIGH or MODERATE impact, or splice-region consequences,
     are retained IF population frequency is rare (AF < 0.01) OR unknown (AF_UNKNOWN).
3. Precedence Rule 3 (Polymorphism & Benign Filtration):
   - Proven common polymorphisms (population AF >= 0.01) without ClinVar P/LP are dropped.
   - Benign, synonymous, and non-splice intronic/intergenic variants are dropped.
4. Genotype Confidence Tagging:
   - DP >= 10, GQ >= 20, Het AB in [0.20, 0.80], Hom AB >= 0.90.
   - Low confidence calls are tagged LOW_CONFIDENCE_GENOTYPE.
   - Low confidence P/LP or High-impact calls are routed for Stage 9 Sanger validation.
"""

import sys
import os
import gzip
import json
import re

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def is_clinvar_plp(clnsig):
    if not clnsig or clnsig == ".":
        return False
    cln_lower = clnsig.lower()
    # Conflicting interpretations alone do not qualify as pure P/LP
    if "conflict" in cln_lower:
        return False
    return "pathogenic" in cln_lower

def extract_population_af(info_str):
    """
    Extracts true population allele frequency from dbSNP CAF or gnomAD.
    CRITICAL: Never reads GATK sample callset AF!
    """
    pop_af = None
    af_source = "NONE"

    for item in info_str.split(";"):
        if item.startswith("CAF="):
            val = item.split("=")[1].strip("[]")
            parts = val.split(",")
            # In dbSNP CAF: parts[0] = REF freq, parts[1] = ALT1 freq
            if len(parts) > 1 and parts[1] != ".":
                try:
                    pop_af = float(parts[1])
                    af_source = "dbSNP_CAF"
                    break
                except (ValueError, IndexError):
                    pass
        elif item.startswith("gnomAD_AF=") or item.startswith("gnomAD_exome_ALL="):
            try:
                val = item.split("=")[1]
                if val != ".":
                    pop_af = float(val)
                    af_source = "gnomAD"
                    break
            except (ValueError, IndexError):
                pass
        elif item.startswith("AF_popmax="):
            try:
                val = item.split("=")[1]
                if val != ".":
                    pop_af = float(val)
                    af_source = "gnomAD_popmax"
                    break
            except (ValueError, IndexError):
                pass

    return pop_af, af_source

def parse_snpeff_annotations(info_str):
    """
    Parses snpEff ANN field across all transcripts.
    Returns:
      has_high_mod (bool),
      is_splice (bool),
      max_impact (str),
      primary_gene (str),
      primary_transcript (str),
      primary_consequence (str),
      hgvs_c (str),
      hgvs_p (str)
    """
    ann_str = ""
    for item in info_str.split(";"):
        if item.startswith("ANN=") or item.startswith("CSQ="):
            ann_str = item.split("=")[1]
            break

    if not ann_str:
        return False, False, "MODIFIER", "Intergenic", ".", "unknown", ".", "."

    has_high_mod = False
    is_splice = False
    max_impact = "MODIFIER"
    primary_gene = "Intergenic"
    primary_transcript = "."
    primary_consequence = "unknown"
    hgvs_c = "."
    hgvs_p = "."

    impact_hierarchy = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
    current_max_score = 0

    for entry in ann_str.split(","):
        parts = entry.split("|")
        if len(parts) > 10:
            conseq = parts[1] if len(parts) > 1 else ""
            impact = parts[2] if len(parts) > 2 else "MODIFIER"
            gene = parts[3] if len(parts) > 3 else "Intergenic"
            feature_id = parts[6] if len(parts) > 6 else "."
            biotype = parts[7] if len(parts) > 7 else ""
            c_syntax = parts[9] if len(parts) > 9 else "."
            p_syntax = parts[10] if len(parts) > 10 else "."

            if impact in ["HIGH", "MODERATE"]:
                has_high_mod = True
            if "splice" in conseq.lower():
                is_splice = True

            score = impact_hierarchy.get(impact, 0)
            # Give slight preference to protein coding transcripts and RefSeq NM_
            if biotype == "protein_coding":
                score += 0.5
            if feature_id.startswith("NM_"):
                score += 0.2

            if score > current_max_score:
                current_max_score = score
                max_impact = impact
                primary_gene = gene if gene else "Intergenic"
                primary_transcript = feature_id
                primary_consequence = conseq.replace("&", " / ").replace("_", " ")
                hgvs_c = c_syntax if c_syntax else "."
                hgvs_p = p_syntax if p_syntax else "."

    return has_high_mod, is_splice, max_impact, primary_gene, primary_transcript, primary_consequence, hgvs_c, hgvs_p

def evaluate_genotype_qc(format_keys, sample_vals):
    fmt = dict(zip(format_keys, sample_vals))
    gt = fmt.get("GT", "./.")
    dp = int(fmt.get("DP", 0)) if fmt.get("DP", "0").isdigit() else 0
    gq = int(fmt.get("GQ", 0)) if fmt.get("GQ", "0").isdigit() else 0
    
    ad = fmt.get("AD", "0,0").split(",")
    ad0 = int(ad[0]) if len(ad) > 0 and ad[0].isdigit() else 0
    ad1 = int(ad[1]) if len(ad) > 1 and ad[1].isdigit() else 0
    tot = ad0 + ad1
    ab = (ad1 / tot) if tot > 0 else 0.0

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

    qc_pass = (dp >= 10 and gq >= 20)
    if zygosity == "Heterozygous":
        qc_pass = qc_pass and (0.20 <= ab <= 0.80)
    elif zygosity == "Homozygous":
        qc_pass = qc_pass and (ab >= 0.90)

    confidence = "HIGH_CONFIDENCE" if qc_pass else "LOW_CONFIDENCE_GENOTYPE"
    return zygosity, dp, gq, ad0, ad1, ab, confidence

def main():
    if len(sys.argv) < 3:
        print("Usage: pass1_permissive_filter.py <input_vcf> <output_dir> [sample_name]")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_dir = sys.argv[2]
    sample_name = sys.argv[3] if len(sys.argv) > 3 else "Sample"

    # Destination directories
    annot_dir = os.path.join(output_dir, "annotation")
    annovar_dir = os.path.join(output_dir, "annovar")
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(annot_dir, exist_ok=True)
    os.makedirs(annovar_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)

    out_vcf_path = os.path.join(annot_dir, f"{sample_name}_pass1_shortlist.vcf")
    out_tsv_path = os.path.join(annot_dir, f"{sample_name}_pass1_shortlist.tsv")
    out_sanger_path = os.path.join(annot_dir, "Requires_Sanger_Validation.txt")
    out_json_path = os.path.join(qc_dir, "pass1_filter_summary.json")

    # Metrics counters
    total_variants = 0
    retained_clinvar_plp = 0
    retained_high_moderate_splice = 0
    dropped_proven_common = 0
    dropped_benign_modifier = 0
    low_confidence_count = 0
    sanger_required_count = 0

    candidate_genes = set()
    retained_records = []
    sanger_records = []

    tsv_headers = [
        "Chr", "Pos", "Ref", "Alt", "Gene", "Transcript", "Impact", "Consequence",
        "HGVS_c", "HGVS_p", "Zygosity", "DP", "GQ", "AD", "AB", "Genotype_Confidence",
        "ClinVar_Significance", "ClinVar_RevStat", "ClinVar_Disease",
        "Pop_AF", "AF_Source", "Retention_Reason", "Requires_Sanger"
    ]

    header_lines = []
    
    with open_vcf(input_vcf) as f_in, \
         open(out_vcf_path, "w", encoding="utf-8") as f_vcf, \
         open(out_tsv_path, "w", encoding="utf-8") as f_tsv, \
         open(out_sanger_path, "w", encoding="utf-8") as f_sanger:

        f_tsv.write("\t".join(tsv_headers) + "\n")
        f_sanger.write("\t".join(tsv_headers) + "\n")

        for line in f_in:
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                # Inject Pass 1 INFO headers before #CHROM
                f_vcf.write('##INFO=<ID=PASS1_STATUS,Number=1,Type=String,Description="Pass 1 Permissive Shortlist retention status">\n')
                f_vcf.write('##INFO=<ID=PASS1_REASON,Number=1,Type=String,Description="Pass 1 retention or exemption rationale">\n')
                f_vcf.write('##INFO=<ID=GENOTYPE_CONF,Number=1,Type=String,Description="Genotype QC confidence rating (HIGH_CONFIDENCE or LOW_CONFIDENCE_GENOTYPE)">\n')
                f_vcf.write('##INFO=<ID=SANGER_REQUIRED,Number=1,Type=String,Description="Flag indicating variant requires Stage 9 Sanger validation">\n')
                for hl in header_lines:
                    f_vcf.write(hl)
                f_vcf.write(line)
                continue

            total_variants += 1
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, rsid, ref, alt, qual, filt, info = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7]

            # 1. ClinVar evaluation
            clnsig = "."
            clnrevstat = "."
            clndn = "."
            for item in info.split(";"):
                if item.startswith("CLNSIG="):
                    clnsig = item.split("=")[1]
                elif item.startswith("CLNREVSTAT="):
                    clnrevstat = item.split("=")[1]
                elif item.startswith("CLNDN="):
                    clndn = item.split("=")[1]

            clinvar_plp = is_clinvar_plp(clnsig)

            # 2. Population frequency evaluation
            pop_af, af_source = extract_population_af(info)
            is_proven_common = (pop_af is not None and pop_af >= 0.01)

            # 3. snpEff consequence & impact
            has_high_mod, is_splice, max_impact, gene, transcript, consequence, hgvs_c, hgvs_p = parse_snpeff_annotations(info)

            # 4. Genotype QC evaluation
            format_keys = cols[8].split(":") if len(cols) > 8 else []
            sample_vals = cols[9].split(":") if len(cols) > 9 else []
            zygosity, dp, gq, ad0, ad1, ab, confidence = evaluate_genotype_qc(format_keys, sample_vals)

            if confidence == "LOW_CONFIDENCE_GENOTYPE":
                low_confidence_count += 1

            # 5. Permissive Funnel Decision
            retained = False
            retention_reason = "DROPPED"
            pass1_status = "DROPPED"

            if clinvar_plp:
                retained = True
                retention_reason = "CLINVAR_PATHOGENIC_EXEMPTION"
                pass1_status = "RETAINED_CLINVAR_PLP"
                retained_clinvar_plp += 1
            elif (has_high_mod or is_splice) and not is_proven_common:
                retained = True
                retention_reason = "HIGH_OR_MODERATE_IMPACT_RARE_OR_UNKNOWN"
                pass1_status = "RETAINED_HIGH_MODERATE"
                retained_high_moderate_splice += 1
            elif is_proven_common:
                dropped_proven_common += 1
                retention_reason = "PROVEN_COMMON_POLYMORPHISM"
            else:
                dropped_benign_modifier += 1
                retention_reason = "BENIGN_OR_MODIFIER_IMPACT"

            if not retained:
                continue

            candidate_genes.add(gene)

            # Check Sanger confirmation routing
            sanger_req = False
            if confidence == "LOW_CONFIDENCE_GENOTYPE" and (clinvar_plp or max_impact == "HIGH"):
                sanger_req = True
                sanger_required_count += 1

            sanger_str = "true" if sanger_req else "false"

            # Append INFO tags
            pass1_info = f"PASS1_STATUS={pass1_status};PASS1_REASON={retention_reason};GENOTYPE_CONF={confidence};SANGER_REQUIRED={sanger_str}"
            cols[7] = f"{info};{pass1_info}"
            f_vcf.write("\t".join(cols) + "\n")

            # TSV row
            pop_af_str = f"{pop_af:.6f}" if pop_af is not None else "AF_UNKNOWN"
            tsv_row = [
                chrom, pos, ref, alt, gene, transcript, max_impact, consequence,
                hgvs_c, hgvs_p, zygosity, str(dp), str(gq), f"{ad0},{ad1}", f"{ab:.2f}", confidence,
                clnsig, clnrevstat, clndn,
                pop_af_str, af_source, retention_reason, sanger_str
            ]
            row_line = "\t".join(tsv_row) + "\n"
            f_tsv.write(row_line)

            if sanger_req:
                f_sanger.write(row_line)

    retained_total = retained_clinvar_plp + retained_high_moderate_splice

    summary_payload = {
        "timestamp": "2026-09-24T15:45:00Z",
        "sampleName": sample_name,
        "pass": "Pass 1 Permissive Filter",
        "totalInputVariants": total_variants,
        "retainedCandidatesCount": retained_total,
        "droppedVariantsCount": dropped_proven_common + dropped_benign_modifier,
        "retentionBreakdown": {
            "clinvarPlpExemptionRetained": retained_clinvar_plp,
            "highModerateSpliceRetained": retained_high_moderate_splice,
            "provenCommonDropped": dropped_proven_common,
            "benignModifierDropped": dropped_benign_modifier
        },
        "genotypeQc": {
            "lowConfidenceTagged": low_confidence_count,
            "sangerConfirmationRequired": sanger_required_count
        },
        "candidateGenesCount": len(candidate_genes),
        "topCandidateGenes": sorted(list(candidate_genes))[:25],
        "artifacts": {
            "shortlistVcf": os.path.relpath(out_vcf_path, output_dir),
            "shortlistTsv": os.path.relpath(out_tsv_path, output_dir),
            "sangerValidationTsv": os.path.relpath(out_sanger_path, output_dir)
        }
    }

    with open(out_json_path, "w", encoding="utf-8") as f_json:
        json.dump(summary_payload, f_json, indent=2)

    # Legacy compatibility symlink / copy to annovar folder
    legacy_shortlist_vcf = os.path.join(annovar_dir, f"{sample_name}_pass1_shortlist.vcf")
    legacy_shortlist_tsv = os.path.join(annovar_dir, f"{sample_name}_pass1_shortlist.tsv")
    try:
        if not os.path.exists(legacy_shortlist_vcf):
            os.symlink(out_vcf_path, legacy_shortlist_vcf)
        if not os.path.exists(legacy_shortlist_tsv):
            os.symlink(out_tsv_path, legacy_shortlist_tsv)
    except OSError:
        pass

    print(f"✅ [PASS 1 COMPLETE] Processed {total_variants} variants for {sample_name}:")
    print(f"   • Retained Candidate Shortlist: {retained_total} variants across {len(candidate_genes)} genes")
    print(f"     - ClinVar P/LP Exemptions: {retained_clinvar_plp}")
    print(f"     - Rare / Novel High/Moderate/Splice: {retained_high_moderate_splice}")
    print(f"   • Filtered Out: {dropped_proven_common + dropped_benign_modifier} variants")
    print(f"     - Proven Common Polymorphisms (AF >= 1%): {dropped_proven_common}")
    print(f"     - Benign / Synonymous / Modifier: {dropped_benign_modifier}")
    print(f"   • Genotype QC & Sanger Routing:")
    print(f"     - Low Confidence Genotypes Tagged: {low_confidence_count}")
    print(f"     - Routed for Stage 9 Sanger Confirmation: {sanger_required_count}")
    print(f"   • Shortlist Artifacts: {out_vcf_path}, {out_tsv_path}")

if __name__ == "__main__":
    main()

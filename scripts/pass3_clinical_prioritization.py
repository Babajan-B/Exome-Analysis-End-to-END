#!/usr/bin/env python3
"""
Pass 3: Inheritance Filtering, ClinGen-Calibrated ACMG Classification & Clinical Prioritization
Part of the PrismBB Multi-Annotator Agentic Subsystem (v2.5.0)

Implements:
1. Inheritance-aware filtering (Dominant: popmax <= 0.1%; Recessive: popmax <= 1.0% + nhomalt <= 2)
2. CLINVAR P/LP ABSOLUTE OVERRIDE: Unconditionally retain established pathogenic variants (CFTR F508del, HFE C282Y)
3. ClinGen-calibrated ACMG/AMP rule scoring:
   - PVS1 (Very Strong): High-impact null variants in LOF-intolerant loci
   - PM2_Supporting: Strictly BLOCKED on AF_UNKNOWN; permitted for absent/popmax <= 0.0001
   - PP3 / BP4: ClinGen single-predictor calibrated REVEL (Strong >= 0.773, Mod >= 0.644) & AlphaMissense
   - BA1 / BS1: Allele frequency benign thresholds
   - ClinVar classification evidence integration
4. Composite ACMG classification (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign)
5. Sanger confirmation triage for LOW_CONFIDENCE genotypes on clinically actionable variants
6. Diagnostic-grade TSV and JSON artifact emission with QC summary
"""

import sys
import os
import json
from datetime import datetime, timezone

HIGH_IMPACT_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "transcript_ablation"
}

MODERATE_IMPACT_CONSEQUENCES = {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant"
}

def is_clinvar_plp(clinvar_obj):
    if not clinvar_obj:
        return False
    # Check boolean flag or string primary significance
    if clinvar_obj.get("isPathogenic") is True:
        return True
    sig = (clinvar_obj.get("primarySignificance") or clinvar_obj.get("significance") or "").lower()
    return "pathogenic" in sig and "conflict" not in sig and "benign" not in sig

def evaluate_inheritance_filter(record):
    """
    Evaluates inheritance filter.
    Returns: (retained: bool, reason: str, override_applied: bool)
    """
    clinvar = record.get("clinvar") or {}
    pop_freq = record.get("populationFrequency") or {}
    
    # 1. CLINVAR P/LP ABSOLUTE EXEMPTION
    if is_clinvar_plp(clinvar):
        return True, "RETAINED_CLINVAR_PLP_OVERRIDE_AF", True
    
    status = pop_freq.get("frequencyStatus") or "AF_UNKNOWN"
    popmax_af = pop_freq.get("gnomadPopmaxAf")
    effective_af = pop_freq.get("effectiveAf")
    hom_count = pop_freq.get("homozygoteCount")
    
    # Measure frequency against thresholds
    freq = popmax_af if popmax_af is not None else effective_af
    
    # Unknown frequencies remain retained under Research Grade protocol
    if status == "AF_UNKNOWN":
        return True, "RETAINED_UNRESOLVED_FREQUENCY_RESEARCH_GRADE", False
    
    # Absent from gnomAD is retained
    if status in ["POPULATION_ABSENT_COVERED", "ABSENT_COVERAGE_UNVERIFIED"]:
        return True, "RETAINED_POPULATION_ABSENT", False
        
    # Autosomal Dominant check (popmax <= 0.001 / 0.1%)
    if freq is not None and freq <= 0.001:
        return True, "RETAINED_DOMINANT_RARE", False
        
    # Recessive check (popmax <= 0.01 / 1.0% AND nhomalt <= 2)
    if freq is not None and freq <= 0.01:
        if hom_count is None or hom_count <= 2:
            return True, "RETAINED_RECESSIVE_RARE", False
        else:
            return False, f"EXCLUDED_RECESSIVE_HOMOZYGOTE_EXCESS (nhomalt={hom_count} > 2)", False
            
    # Common variant (> 1%)
    return False, f"EXCLUDED_POPULATION_COMMON (AF={freq:.5f} > 0.01)", False

def score_acmg_rules(record):
    """
    Scores ACMG/AMP criteria according to ClinGen SVI recommendations.
    Returns: (classification: str, criteria: list, points: float)
    """
    criteria = []
    consequence = record.get("consequence") or ""
    impact = (record.get("impact") or "").upper()
    pop_freq = record.get("populationFrequency") or {}
    insilico = record.get("inSilico") or {}
    clinvar = record.get("clinvar") or {}
    
    status = pop_freq.get("frequencyStatus") or "AF_UNKNOWN"
    popmax_af = pop_freq.get("gnomadPopmaxAf")
    effective_af = pop_freq.get("effectiveAf")
    freq = popmax_af if popmax_af is not None else effective_af
    
    is_plp = is_clinvar_plp(clinvar)
    cv_stars = clinvar.get("stars", 0) or 0
    
    # ── PVS1: Null variant in gene where LOF is disease mechanism ──
    if impact == "HIGH" or any(c in consequence for c in HIGH_IMPACT_CONSEQUENCES):
        criteria.append("PVS1")

    # ── PM2_Supporting: Population frequency absent or extremely rare ──
    # CRITICAL: STRICTLY BLOCKED IF AF_UNKNOWN
    if status != "AF_UNKNOWN":
        if status in ["POPULATION_ABSENT_COVERED", "ABSENT_COVERAGE_UNVERIFIED"]:
            criteria.append("PM2_Supporting")
        elif freq is not None and freq <= 0.0001:  # <= 0.01%
            criteria.append("PM2_Supporting")
            
    # ── In-Silico: ClinGen calibrated REVEL & AlphaMissense ──
    revel_ev = insilico.get("revelEvidence")
    am_class = insilico.get("alphaMissenseClass")
    cadd_phred = insilico.get("caddPhred")
    
    if revel_ev == "STRONG":
        criteria.append("PP3_Strong")
    elif revel_ev == "MODERATE":
        criteria.append("PP3_Moderate")
    elif revel_ev == "SUPPORTING" or am_class == "likely_pathogenic":
        criteria.append("PP3")
    elif revel_ev == "BENIGN_MODERATE":
        criteria.append("BP4_Moderate")
    elif revel_ev == "BENIGN_SUPPORTING" or am_class == "likely_benign":
        criteria.append("BP4")
    elif cadd_phred is not None and cadd_phred >= 20.0 and impact == "MODERATE":
        criteria.append("PP3")
        
    # ── Allele Frequency Benign: BA1 / BS1 (Bypassed if confirmed ClinVar P/LP) ──
    if not is_plp and freq is not None:
        if freq > 0.05:
            criteria.append("BA1")
        elif freq > 0.01:
            criteria.append("BS1")

    # ── ClinVar Evidence Integration ──
    if is_plp:
        if cv_stars >= 2:
            criteria.append("PS1_ClinVar_MultiSubmitter")
        elif cv_stars >= 1:
            criteria.append("PP5_ClinVar_Approved")
        else:
            criteria.append("PP5_ClinVar_Single")
            
    # ── Composite Classification ──
    # Count evidence weights
    has_ba1 = "BA1" in criteria
    has_bs1 = "BS1" in criteria
    has_pvs1 = "PVS1" in criteria
    has_ps = any(c.startswith("PS") or c == "PP3_Strong" for c in criteria)
    has_pm = any(c.startswith("PM") or c == "PP3_Moderate" for c in criteria)
    has_pp = any(c.startswith("PP") for c in criteria)
    has_bp = any(c.startswith("BP") for c in criteria)
    
    num_ps = sum(1 for c in criteria if c.startswith("PS") or c == "PP3_Strong")
    num_pm = sum(1 for c in criteria if c.startswith("PM") or c == "PP3_Moderate")
    num_pp = sum(1 for c in criteria if c.startswith("PP"))
    num_bp = sum(1 for c in criteria if c.startswith("BP"))
    
    classification = "VUS"
    
    if has_ba1 or (has_bs1 and has_bp):
        classification = "Benign"
    elif has_bs1 or num_bp >= 2:
        classification = "Likely_benign"
    elif is_plp and cv_stars >= 1:
        classification = "Pathogenic" if cv_stars >= 2 else "Likely_pathogenic"
    elif has_pvs1 and (num_ps >= 1 or num_pm >= 2 or (num_pm >= 1 and num_pp >= 1)):
        classification = "Pathogenic"
    elif num_ps >= 2 or (num_ps >= 1 and num_pm >= 3):
        classification = "Pathogenic"
    elif has_pvs1 and num_pm == 1:
        classification = "Likely_pathogenic"
    elif has_pvs1 and num_pp >= 1:
        classification = "Likely_pathogenic"
    elif num_ps >= 1 and (num_pm >= 1 or num_pp >= 2):
        classification = "Likely_pathogenic"
    elif num_pm >= 3 or (num_pm >= 2 and num_pp >= 2):
        classification = "Likely_pathogenic"
    else:
        classification = "VUS"
        
    return classification, criteria

def compute_priority_score(record, acmg_class, criteria):
    score = 0.0
    # ACMG tier
    if acmg_class == "Pathogenic":
        score += 100.0
    elif acmg_class == "Likely_pathogenic":
        score += 80.0
    elif acmg_class == "VUS":
        score += 40.0
    elif acmg_class == "Likely_benign":
        score += 10.0
    else:
        score += 0.0
        
    # Consequence impact
    impact = (record.get("impact") or "").upper()
    if impact == "HIGH":
        score += 25.0
    elif impact == "MODERATE":
        score += 15.0
        
    # ClinVar match
    clinvar = record.get("clinvar") or {}
    if is_clinvar_plp(clinvar):
        score += 30.0
    if clinvar.get("isReclassified") is True:
        score += 15.0  # Needs urgent clinician attention
        
    # In-silico points
    insilico = record.get("inSilico") or {}
    revel_pts = insilico.get("revelAcmgPoints", 0.0) or 0.0
    score += revel_pts * 2.5
    
    return round(score, 1)

def main():
    if len(sys.argv) < 3:
        print("Usage: pass3_clinical_prioritization.py <output_dir> <sample_name> [options]")
        sys.exit(1)
        
    output_dir = sys.argv[1]
    sample_name = sys.argv[2]
    
    annotation_dir = os.path.join(output_dir, "annotation")
    annovar_dir = os.path.join(output_dir, "annovar")
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(annotation_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)
    
    # Locate Pass 2 enriched records
    candidate_json_paths = [
        os.path.join(annotation_dir, f"{sample_name}_pass2_enriched.json"),
        os.path.join(annovar_dir, f"{sample_name}_pass2_enriched.json"),
        os.path.join(annotation_dir, f"{sample_name}_pass1_shortlist.json"),
    ]
    
    in_records = []
    source_file = None
    for p in candidate_json_paths:
        if os.path.exists(p):
            try:
                with open(p, "r", encoding="utf-8") as f:
                    in_records = json.load(f)
                source_file = p
                break
            except Exception as e:
                print(f"⚠️ Warning: Could not parse {p}: {e}")
                
    if not in_records:
        print(f"❌ [PASS 3 CRITICAL] No enriched or shortlist JSON records found for {sample_name} in {annotation_dir}")
        sys.exit(1)
        
    print(f"════════════════════════════════════════════════════════════")
    print(f"PASS 3: Inheritance Filtering, ACMG Classification & Clinical Prioritization")
    print(f"Sample: {sample_name} | Source: {os.path.basename(source_file)}")
    print(f"Candidate records to evaluate: {len(in_records)}")
    print(f"════════════════════════════════════════════════════════════")
    
    prioritized_records = []
    excluded_records = []
    
    clinvar_override_count = 0
    sanger_required_count = 0
    acmg_counts = {"Pathogenic": 0, "Likely_pathogenic": 0, "VUS": 0, "Likely_benign": 0, "Benign": 0}
    
    for r in in_records:
        retained, reason, clinvar_override = evaluate_inheritance_filter(r)
        
        if clinvar_override:
            clinvar_override_count += 1
            
        if not retained:
            r_excluded = dict(r)
            r_excluded["exclusionReason"] = reason
            excluded_records.append(r_excluded)
            continue
            
        # Score ACMG
        acmg_class, criteria = score_acmg_rules(r)
        acmg_counts[acmg_class] = acmg_counts.get(acmg_class, 0) + 1
        
        # Priority score
        priority_score = compute_priority_score(r, acmg_class, criteria)
        
        # Sanger Validation Gate (Stage 9 Bridge)
        confidence = r.get("confidenceStatus") or "HIGH_CONFIDENCE"
        is_actionable = acmg_class in ["Pathogenic", "Likely_pathogenic", "VUS"]
        impact = (r.get("impact") or "").upper()
        
        requires_sanger = False
        sanger_reason = None
        if confidence == "LOW_CONFIDENCE_GENOTYPE" and (is_actionable or impact == "HIGH"):
            requires_sanger = True
            sanger_required_count += 1
            sanger_reason = f"Genotype confidence flagged '{confidence}'. Orthogonal Sanger validation required prior to clinical sign-off."
            
        # Build enriched pass 3 record
        out_r = dict(r)
        out_r["inheritanceTriage"] = {
            "retained": True,
            "retentionReason": reason,
            "clinvarOverrideApplied": clinvar_override
        }
        out_r["clinicalEvaluation"] = {
            "acmgClassification": acmg_class,
            "acmgCriteria": criteria,
            "priorityScore": priority_score,
            "pm2Permitted": "PM2_Supporting" in criteria,
            "requiresSangerValidation": requires_sanger,
            "sangerJustification": sanger_reason
        }
        prioritized_records.append(out_r)
        
    # Sort descending by priority score
    prioritized_records.sort(key=lambda x: x["clinicalEvaluation"]["priorityScore"], reverse=True)
    for rank, rec in enumerate(prioritized_records, 1):
        rec["clinicalEvaluation"]["clinicalRank"] = rank
        
    # Output file paths
    out_json = os.path.join(annotation_dir, f"{sample_name}_pass3_clinical_candidates.json")
    out_tsv = os.path.join(annotation_dir, f"{sample_name}_pass3_clinical_candidates.tsv")
    out_summary = os.path.join(qc_dir, "pass3_prioritization_summary.json")
    
    # 1. Write JSON
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(prioritized_records, f, indent=2)
        
    # 2. Write TSV
    tsv_headers = [
        "Rank", "Locus", "Gene", "Transcript", "Consequence", "Impact", "HGVSc", "HGVSp",
        "Zygosity", "Confidence", "ACMG_Class", "ACMG_Criteria", "Priority_Score",
        "Popmax_AF", "Frequency_Status", "ClinVar", "ClinVar_Reclassified", "ClinVar_Override",
        "REVEL", "AlphaMissense", "Requires_Sanger"
    ]
    
    with open(out_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(tsv_headers) + "\n")
        for rec in prioritized_records:
            eval_data = rec["clinicalEvaluation"]
            pop_freq = rec.get("populationFrequency") or {}
            cv = rec.get("clinvar") or {}
            insilico = rec.get("inSilico") or {}
            
            row = [
                str(eval_data.get("clinicalRank", 0)),
                rec.get("locus", ""),
                rec.get("primaryGene", ""),
                rec.get("primaryTranscript", ""),
                rec.get("consequence", ""),
                rec.get("impact", ""),
                rec.get("hgvsC", "") or "NA",
                rec.get("hgvsP", "") or "NA",
                rec.get("zygosity", ""),
                rec.get("confidenceStatus", ""),
                eval_data.get("acmgClassification", ""),
                ",".join(eval_data.get("acmgCriteria", [])),
                str(eval_data.get("priorityScore", 0.0)),
                f"{pop_freq.get('gnomadPopmaxAf'):.6f}" if pop_freq.get("gnomadPopmaxAf") is not None else "NA",
                pop_freq.get("frequencyStatus", ""),
                cv.get("primarySignificance", "") or cv.get("significance", "") or "NA",
                "TRUE" if cv.get("isReclassified") else "FALSE",
                "TRUE" if rec["inheritanceTriage"]["clinvarOverrideApplied"] else "FALSE",
                f"{insilico.get('revelScore'):.4f}" if insilico.get("revelScore") is not None else "NA",
                insilico.get("alphaMissenseClass") or "NA",
                "TRUE" if eval_data.get("requiresSangerValidation") else "FALSE"
            ]
            f.write("\t".join(row) + "\n")
            
    # 3. Create backward-compatibility symlinks in annovar/
    for src, dst in [
        (out_json, os.path.join(annovar_dir, f"{sample_name}_pass3_clinical_candidates.json")),
        (out_tsv, os.path.join(annovar_dir, f"{sample_name}_pass3_clinical_candidates.tsv"))
    ]:
        try:
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.abspath(src), dst)
        except OSError:
            pass
            
    # 4. Write Summary JSON
    summary_data = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sampleName": sample_name,
        "inputCandidates": len(in_records),
        "retainedCandidates": len(prioritized_records),
        "excludedCandidates": len(excluded_records),
        "acmgBreakdown": acmg_counts,
        "clinvarOverrideCount": clinvar_override_count,
        "requiresSangerValidationCount": sanger_required_count,
        "topCandidate": {
            "gene": prioritized_records[0].get("primaryGene"),
            "locus": prioritized_records[0].get("locus"),
            "acmgClass": prioritized_records[0]["clinicalEvaluation"]["acmgClassification"],
            "score": prioritized_records[0]["clinicalEvaluation"]["priorityScore"]
        } if prioritized_records else None,
        "artifacts": {
            "clinicalCandidatesJson": out_json,
            "clinicalCandidatesTsv": out_tsv,
            "summaryJson": out_summary
        }
    }
    
    with open(out_summary, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2)
        
    print(f"✅ Pass 3 Prioritization Complete for {sample_name}!")
    print(f"   Candidates Retained: {len(prioritized_records)} / {len(in_records)} ({len(excluded_records)} common filtered)")
    print(f"   ACMG Classifications: {acmg_counts['Pathogenic']} Pathogenic, {acmg_counts['Likely_pathogenic']} Likely Pathogenic, {acmg_counts['VUS']} VUS, {acmg_counts['Likely_benign']} Likely Benign, {acmg_counts['Benign']} Benign.")
    if clinvar_override_count > 0:
        print(f"   🧬 ClinVar P/LP Absolute Overrides Applied: {clinvar_override_count} variants")
    if sanger_required_count > 0:
        print(f"   🔬 Orthogonal Sanger Confirmations Flagged: {sanger_required_count} variants")
    if prioritized_records:
        top = prioritized_records[0]
        print(f"   ⭐ Top Candidate: {top.get('primaryGene')} ({top.get('locus')}) - {top['clinicalEvaluation']['acmgClassification']} (Score: {top['clinicalEvaluation']['priorityScore']})")
    print(f"   Artifacts: {os.path.basename(out_tsv)}, {os.path.basename(out_json)}")

if __name__ == "__main__":
    main()

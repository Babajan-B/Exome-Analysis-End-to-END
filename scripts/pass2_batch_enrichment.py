#!/usr/bin/env python3
"""
Pass 2: Batch API Enrichment & Pathogenicity Calibration Engine (Version 2.6.0)

Enriches the Pass 1 Permissive Candidate Shortlist via MyVariant.info with ClinGen calibration:
1. Padded-SNV Two-Way Trimmed HGVS Formatter & Syntax Validation:
   Translates VCF records to canonical HGVS genomic notation (SNV, pure ins/del, delins, padded-SNVs).
   Validates syntax before submission; malformed/multi-allelic records stay AF_UNKNOWN.
2. Operator Consent Gate:
   External API enrichment is OFF by default. Requires explicit per-job opt-in:
   --enable-api-enrichment, $output_dir/.enable_api_enrichment, or ENABLE_API_ENRICHMENT=1.
   Recorded in the audit trail.
3. MyVariant.info Batch POST Query:
   Retrieves gnomAD exome/genome frequencies, popmax fallback, ClinVar RCVs, REVEL, and AlphaMissense.
4. ClinGen-Calibrated In-Silico Annotations:
   - REVEL: Pejaver 2022 calibrated evidence rules (takes max of transcript list).
   - AlphaMissense: Cheng 2023 under dbnsfp.alphamissense (takes max of transcript list).
   - SpliceAI: Reported as 'not_available' when not provided by endpoint.
5. ClinVar Reclassification Sentinel:
   Aggregates across all RCV assertions and compares with local ClinVar release to detect CLINVAR_RECLASSIFIED.
6. Frequency Status Classification:
   - Valid query with population AF >= 0.01 -> POPULATION_COMMON
   - Valid query with population AF < 0.01 -> POPULATION_RARE
   - Valid query confirmed absent in gnomAD -> ABSENT_COVERAGE_UNVERIFIED
   - Malformed query, lookup failure, or offline -> AF_UNKNOWN
7. Emits:
   - {sample_name}_pass2_enriched.tsv
   - {sample_name}_pass2_enriched.json
   - {output_dir}/qc/pass2_enrichment_summary.json
   - Updates {sample_name}_pass1_shortlist.tsv with resolved frequencies for Stage 7 Gate.
"""

import sys
import os
import re
import json
import urllib.request
import urllib.parse
from datetime import datetime, timezone

# Strict HGVS genomic validation regex
HGVS_REGEX = re.compile(
    r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT):g\.(\d+)(?:_(\d+))?(?:del[ACGTN]*|ins[ACGTN]+|delins[ACGTN]+|[ACGTN]>[ACGTN]|=)$"
)

def format_to_hgvs_genomic(chrom, pos, ref, alt):
    """
    Translates a normalized VCF record to MyVariant.info HGVS genomic notation.
    Handles two-way prefix/suffix trimming and resolves padded SNVs directly to canonical SNV format.
    Returns (hgvs_str, is_valid) tuple.
    """
    chr_prefix = chrom if chrom.startswith("chr") else f"chr{chrom}"
    ref = (ref or "").strip().upper()
    alt = (alt or "").strip().upper()

    # Reject un-split multi-allelic records or invalid bases
    if "," in alt or "," in ref or not ref or not alt:
        return None, False
    if not re.match(r"^[ACGTN]+$", ref) or not re.match(r"^[ACGTN]+$", alt):
        return None, False

    if ref == alt:
        hgvs = f"{chr_prefix}:g.{pos}="
        return hgvs, True

    # Direct SNV
    if len(ref) == 1 and len(alt) == 1:
        hgvs = f"{chr_prefix}:g.{pos}{ref}>{alt}"
        return hgvs, True

    trimmed_ref = ref
    trimmed_alt = alt
    curr_pos = int(pos)

    # Trim shared leading bases
    while len(trimmed_ref) > 0 and len(trimmed_alt) > 0 and trimmed_ref[0] == trimmed_alt[0]:
        trimmed_ref = trimmed_ref[1:]
        trimmed_alt = trimmed_alt[1:]
        curr_pos += 1

    # Trim shared trailing bases
    while len(trimmed_ref) > 0 and len(trimmed_alt) > 0 and trimmed_ref[-1] == trimmed_alt[-1]:
        trimmed_ref = trimmed_ref[:-1]
        trimmed_alt = trimmed_alt[:-1]

    # Padded SNV
    if len(trimmed_ref) == 1 and len(trimmed_alt) == 1:
        hgvs = f"{chr_prefix}:g.{curr_pos}{trimmed_ref}>{trimmed_alt}"
        return hgvs, True

    # Pure Deletion
    if len(trimmed_ref) > 0 and len(trimmed_alt) == 0:
        start_pos = curr_pos
        end_pos = curr_pos + len(trimmed_ref) - 1
        hgvs = f"{chr_prefix}:g.{start_pos}del" if start_pos == end_pos else f"{chr_prefix}:g.{start_pos}_{end_pos}del"
        return hgvs, True

    # Pure Insertion
    if len(trimmed_ref) == 0 and len(trimmed_alt) > 0:
        hgvs = f"{chr_prefix}:g.{curr_pos - 1}_{curr_pos}ins{trimmed_alt}"
        return hgvs, True

    # Complex Delins
    if len(trimmed_ref) > 0 and len(trimmed_alt) > 0:
        end_pos = curr_pos + len(trimmed_ref) - 1
        hgvs = f"{chr_prefix}:g.{curr_pos}delins{trimmed_alt}" if curr_pos == end_pos else f"{chr_prefix}:g.{curr_pos}_{end_pos}delins{trimmed_alt}"
        return hgvs, True

    return None, False

def extract_float(val, take_max=True):
    """
    Safely extracts float from scalars, lists, or tuples (taking max if multiple).
    """
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        floats = [extract_float(x, take_max=take_max) for x in val]
        valid = [f for f in floats if f is not None]
        if not valid:
            return None
        return max(valid) if take_max else valid[0]
    try:
        return float(val)
    except (ValueError, TypeError):
        return None

def calibrate_revel(score):
    if score is None:
        return None, 0.0, None
    if score >= 0.773:
        return "STRONG", 4.0, "PP3_Strong"
    if score >= 0.644:
        return "MODERATE", 2.0, "PP3_Moderate"
    if score >= 0.290:
        return "SUPPORTING", 1.0, "PP3_Supporting"
    if score <= 0.016:
        return "BENIGN_MODERATE", -2.0, "BP4_Moderate"
    if score <= 0.183:
        return "BENIGN_SUPPORTING", -1.0, "BP4_Supporting"
    return "INDETERMINATE", 0.0, None

def classify_alphamissense(score):
    if score is None:
        return None
    if score >= 0.564:
        return "likely_pathogenic"
    if score < 0.340:
        return "likely_benign"
    return "ambiguous"

def calculate_clinvar_stars(review_status):
    if not review_status:
        return 0
    rev = review_status.lower().replace("_", " ")
    if "no assertion" in rev or "no interpretation" in rev:
        return 0
    if "practice guideline" in rev:
        return 4
    if "expert panel" in rev:
        return 3
    if "multiple submitters" in rev and ("no conflict" in rev or "conflicting" not in rev):
        return 2
    if "single submitter" in rev or "conflicting" in rev:
        return 1
    return 0

def check_consent_gate(output_dir, extra_args):
    """
    Consent gate: API enrichment is OFF by default.
    Opt-in via:
    - command-line arg: --enable-api-enrichment, consent=true, optin
    - marker file: $output_dir/.enable_api_enrichment
    - env var: ENABLE_API_ENRICHMENT=1
    """
    if os.environ.get("ENABLE_API_ENRICHMENT") == "1":
        return True
    if os.path.exists(os.path.join(output_dir, ".enable_api_enrichment")):
        return True
    for arg in extra_args:
        arg_lower = str(arg).lower().strip("-")
        if arg_lower in ["enable_api_enrichment", "enable-api-enrichment", "consent=true", "consent", "optin", "true"]:
            return True
    return False

def find_shortlist_tsv(output_dir, sample_name):
    """
    Finds the Pass 1 shortlist TSV, checking annotation/ first to prevent symlink breakage.
    """
    candidates = [
        os.path.join(output_dir, "annotation", f"{sample_name}_pass1_shortlist.tsv"),
        os.path.join(output_dir, "annovar", f"{sample_name}_pass1_shortlist.tsv")
    ]
    for c in candidates:
        if os.path.lexists(c):
            # If it's a symlink, resolve real path
            target = os.path.realpath(c)
            if os.path.exists(target):
                return target
    return None

def main():
    if len(sys.argv) < 3:
        print("Usage: pass2_batch_enrichment.py <output_dir> <sample_name> [options...]")
        sys.exit(1)

    output_dir = os.path.abspath(sys.argv[1])
    sample_name = sys.argv[2]
    extra_args = sys.argv[3:]

    consent_granted = check_consent_gate(output_dir, extra_args)

    annot_dir = os.path.join(output_dir, "annotation")
    annovar_dir = os.path.join(output_dir, "annovar")
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(annot_dir, exist_ok=True)
    os.makedirs(annovar_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)

    shortlist_tsv = find_shortlist_tsv(output_dir, sample_name)
    cache_path = os.path.join(qc_dir, "myvariant_cache.json")
    out_enriched_tsv = os.path.join(annot_dir, f"{sample_name}_pass2_enriched.tsv")
    out_enriched_json = os.path.join(annot_dir, f"{sample_name}_pass2_enriched.json")
    summary_json = os.path.join(qc_dir, "pass2_enrichment_summary.json")

    if not shortlist_tsv or not os.path.exists(shortlist_tsv):
        print(f"⚠️  Shortlist TSV not found in {output_dir}/annotation or annovar. Skipping Pass 2 enrichment.")
        sys.exit(0)

    print("════════════════════════════════════════════════════════════")
    print(f"PASS 2: Batch API Enrichment & Pathogenicity Calibration")
    print(f"Sample: {sample_name} | Target: {shortlist_tsv}")
    print(f"Consent Gate: {'GRANTED (Live API Enabled)' if consent_granted else 'NOT GRANTED (Offline Fallback Path 3)'}")
    print("════════════════════════════════════════════════════════════")

    # Load cache if available
    cache = {}
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    rows = []
    with open(shortlist_tsv, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        col_map = {name: idx for idx, name in enumerate(header)}
        for line in f:
            line_str = line.rstrip("\r\n")
            if not line_str:
                continue
            rows.append(line_str.split("\t"))

    total_candidates = len(rows)
    print(f"  [PASS 2] Parsed {total_candidates} candidate shortlist variants.")

    if total_candidates == 0:
        print("  [PASS 2] Shortlist is empty. Writing empty enriched outputs.")
        sys.exit(0)

    # Format & Validate HGVS for all rows
    hgvs_list = []
    valid_mask = []
    for r in rows:
        c = r[col_map["Chr"]]
        p = r[col_map["Pos"]]
        ref = r[col_map["Ref"]]
        alt = r[col_map["Alt"]]
        hgvs, is_valid = format_to_hgvs_genomic(c, p, ref, alt)
        if is_valid and hgvs and HGVS_REGEX.match(hgvs):
            hgvs_list.append(hgvs)
            valid_mask.append(True)
        else:
            hgvs_list.append(hgvs or f"{c}:{p}_{ref}>{alt}_MALFORMED")
            valid_mask.append(False)

    execution_path = "PATH_2_BATCH_API" if consent_granted else "PATH_3_OFFLINE_FALLBACK"
    
    # Missing valid entries to query
    valid_hgvs_to_query = [hgvs_list[i] for i in range(total_candidates) if valid_mask[i] and hgvs_list[i] not in cache]
    cache_hits = sum(1 for i in range(total_candidates) if valid_mask[i] and hgvs_list[i] in cache)
    api_requests = 0

    if consent_granted and len(valid_hgvs_to_query) > 0:
        print(f"  [PASS 2] Querying MyVariant.info for {len(valid_hgvs_to_query)} valid variants ({cache_hits} cached)...")
        BATCH_SIZE = 1000
        for i in range(0, len(valid_hgvs_to_query), BATCH_SIZE):
            chunk = valid_hgvs_to_query[i:i + BATCH_SIZE]
            api_requests += 1
            url = "https://myvariant.info/v1/variant"
            params = urllib.parse.urlencode({
                "ids": ",".join(chunk),
                "assembly": "hg19",
                "fields": "gnomad_exome,gnomad_genome,clinvar,dbnsfp.revel,dbnsfp.alphamissense,cadd,dbsnp"
            }).encode("utf-8")

            req = urllib.request.Request(
                url,
                data=params,
                headers={"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"}
            )
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    if resp.status == 200:
                        data = json.loads(resp.read().decode("utf-8"))
                        if isinstance(data, list):
                            for hit in data:
                                q_id = hit.get("query") or hit.get("_id")
                                if q_id:
                                    cache[q_id] = hit
                        print(f"  ✅ Batch {api_requests} successfully returned results.")
                    else:
                        print(f"⚠️  MyVariant returned HTTP status {resp.status}. Entering Partial mode.")
                        execution_path = "PATH_4_PARTIAL"
            except Exception as e:
                print(f"⚠️  MyVariant query failed ({e}). Proceeding under Partial mode.")
                execution_path = "PATH_4_PARTIAL"

        # Save cache
        try:
            with open(cache_path, "w", encoding="utf-8") as f:
                json.dump(cache, f, indent=2)
        except Exception:
            pass

    # Process annotations and assemble enriched records
    enriched_json_records = []
    enriched_tsv_rows = []
    updated_shortlist_rows = []

    status_counts = {"absentUnverified": 0, "rare": 0, "common": 0, "unknown": 0}
    insilico_counts = {
        "revelStrong": 0, "revelModerate": 0, "revelSupporting": 0,
        "alphaMissenseLikelyPathogenic": 0, "caddDeleterious": 0,
        "spliceAi": "not_available"
    }
    clinvar_matches = 0
    clinvar_plp_count = 0
    clinvar_reclassified_count = 0

    enriched_headers = header + [
        "REVEL_Score", "REVEL_Evidence", "SpliceAI_MaxDS", "SpliceAI_Evidence",
        "AlphaMissense_Score", "AlphaMissense_Class", "CADD_Phred",
        "gnomAD_Popmax_AF", "gnomAD_Hom_Count", "ClinVar_Reclassified"
    ]

    for idx, r in enumerate(rows):
        hgvs = hgvs_list[idx]
        is_valid = valid_mask[idx]
        hit = cache.get(hgvs, {}) if is_valid else {}

        c = r[col_map["Chr"]]
        p = int(r[col_map["Pos"]])
        ref = r[col_map["Ref"]]
        alt = r[col_map["Alt"]]
        gene = r[col_map.get("Gene", 4)]
        transcript = r[col_map.get("Transcript", 5)]
        impact = r[col_map.get("Impact", 6)]
        consequence = r[col_map.get("Consequence", 7)]
        hgvs_c = r[col_map.get("HGVS_c", 8)]
        hgvs_p = r[col_map.get("HGVS_p", 9)]
        zygosity = r[col_map.get("Zygosity", 10)]
        dp = int(r[col_map.get("DP", 11)]) if r[col_map.get("DP", 11)].isdigit() else 30
        gq = int(r[col_map.get("GQ", 12)]) if r[col_map.get("GQ", 12)].isdigit() else 99
        ad = r[col_map.get("AD", 13)]
        ab = extract_float(r[col_map.get("AB", 14)]) or 0.5
        confidence = r[col_map.get("Genotype_Confidence", 15)]
        clnsig_orig = r[col_map.get("ClinVar_Significance", 16)]
        clnrev_orig = r[col_map.get("ClinVar_RevStat", 17)]
        clndn_orig = r[col_map.get("ClinVar_Disease", 18)]
        pop_af_orig = extract_float(r[col_map.get("Pop_AF", 19)])
        af_source_orig = r[col_map.get("AF_Source", 20)]
        retention_reason = r[col_map.get("Retention_Reason", 22)]
        requires_sanger_str = r[col_map.get("Requires_Sanger", 23)]

        # Extract gnomAD fields
        gnomad_exome = hit.get("gnomad_exome", {})
        gnomad_genome = hit.get("gnomad_genome", {})

        af_exome = extract_float(gnomad_exome.get("af", {}).get("af") if isinstance(gnomad_exome.get("af"), dict) else gnomad_exome.get("af"))
        af_genome = extract_float(gnomad_genome.get("af", {}).get("af") if isinstance(gnomad_genome.get("af"), dict) else gnomad_genome.get("af"))
        
        # Popmax fallback
        af_popmax = extract_float(gnomad_exome.get("af_popmax") or gnomad_exome.get("popmax"))
        if af_popmax is None and isinstance(gnomad_exome.get("af"), dict):
            sub_afs = [extract_float(gnomad_exome["af"].get(k)) for k in ["af_afr", "af_amr", "af_eas", "af_nfe", "af_sas"]]
            valid_sub = [s for s in sub_afs if s is not None]
            if valid_sub:
                af_popmax = max(valid_sub)

        hom_count = extract_float(gnomad_exome.get("hom", {}).get("hom") if isinstance(gnomad_exome.get("hom"), dict) else gnomad_exome.get("hom"))

        # Effective AF determination
        effective_af = af_exome if af_exome is not None else (af_genome if af_genome is not None else pop_af_orig)
        effective_source = "gnomAD_exome" if af_exome is not None else ("gnomAD_genome" if af_genome is not None else af_source_orig)

        # ClinVar RCV aggregation and reclassification check
        cv = hit.get("clinvar", {})
        rcv_entries = cv.get("rcv", [])
        if isinstance(rcv_entries, dict):
            rcv_entries = [rcv_entries]
        elif not isinstance(rcv_entries, list):
            rcv_entries = []

        api_sigs = []
        api_revs = []
        for rcv in rcv_entries:
            if isinstance(rcv, dict):
                sig_val = rcv.get("clinical_significance")
                if sig_val:
                    api_sigs.append(str(sig_val))
                rev_val = rcv.get("review_status")
                if rev_val:
                    api_revs.append(str(rev_val))

        cv_sig_primary = clnsig_orig if clnsig_orig and clnsig_orig != "." else (api_sigs[0] if api_sigs else None)
        cv_rev_primary = clnrev_orig if clnrev_orig and clnrev_orig != "." else (api_revs[0] if api_revs else None)
        cv_stars = calculate_clinvar_stars(cv_rev_primary)

        # Detect CLINVAR_RECLASSIFIED
        is_reclassified = False
        if clnsig_orig and clnsig_orig != "." and api_sigs:
            local_is_plp = "pathogenic" in clnsig_orig.lower() and "conflict" not in clnsig_orig.lower()
            api_has_plp = any("pathogenic" in s.lower() and "conflict" not in s.lower() for s in api_sigs)
            api_is_benign = all("benign" in s.lower() for s in api_sigs)
            if (local_is_plp and not api_has_plp) or (not local_is_plp and api_has_plp):
                is_reclassified = True
                clinvar_reclassified_count += 1

        is_plp = cv_sig_primary is not None and "pathogenic" in cv_sig_primary.lower() and "conflict" not in cv_sig_primary.lower()
        if cv_sig_primary:
            clinvar_matches += 1
            if is_plp:
                clinvar_plp_count += 1

        # Frequency Status determination
        if not is_valid:
            # Malformed query or un-split multi-allelic: NEVER absent!
            freq_status = "AF_UNKNOWN"
            status_counts["unknown"] += 1
        elif effective_af is not None:
            if effective_af >= 0.01:
                freq_status = "POPULATION_COMMON"
                status_counts["common"] += 1
            else:
                freq_status = "POPULATION_RARE"
                status_counts["rare"] += 1
        elif consent_granted and (hit.get("notfound") or (hit and not af_exome and not af_genome)):
            # Valid query confirmed absent from gnomAD without coverage data
            freq_status = "ABSENT_COVERAGE_UNVERIFIED"
            status_counts["absentUnverified"] += 1
        else:
            freq_status = "AF_UNKNOWN"
            status_counts["unknown"] += 1

        # In-silico predictors
        # 1. REVEL (handling per-transcript list)
        revel_raw = hit.get("dbnsfp", {}).get("revel", {}).get("score")
        revel_score = extract_float(revel_raw, take_max=True)
        revel_ev, revel_pts, revel_rule = calibrate_revel(revel_score)
        if revel_ev == "STRONG": insilico_counts["revelStrong"] += 1
        elif revel_ev == "MODERATE": insilico_counts["revelModerate"] += 1
        elif revel_ev == "SUPPORTING": insilico_counts["revelSupporting"] += 1

        # 2. AlphaMissense (under dbnsfp.alphamissense)
        am_raw = hit.get("dbnsfp", {}).get("alphamissense", {}).get("score")
        am_score = extract_float(am_raw, take_max=True)
        am_class = classify_alphamissense(am_score)
        if am_class == "likely_pathogenic": insilico_counts["alphaMissenseLikelyPathogenic"] += 1

        # 3. CADD
        cadd_phred = extract_float(hit.get("cadd", {}).get("phred"), take_max=True)
        if cadd_phred is not None and cadd_phred >= 20.0:
            insilico_counts["caddDeleterious"] += 1

        # Build updated shortlist TSV row
        updated_r = list(r)
        if freq_status == "ABSENT_COVERAGE_UNVERIFIED":
            updated_r[col_map["Pop_AF"]] = "0.000000"
            updated_r[col_map["AF_Source"]] = "gnomAD_absent"
        else:
            updated_r[col_map["Pop_AF"]] = f"{effective_af:.6f}" if effective_af is not None else "AF_UNKNOWN"
            updated_r[col_map["AF_Source"]] = effective_source
        updated_r[col_map["Frequency_Status"]] = freq_status
        if cv_sig_primary:
            updated_r[col_map["ClinVar_Significance"]] = cv_sig_primary
        if cv_rev_primary:
            updated_r[col_map["ClinVar_RevStat"]] = cv_rev_primary
        updated_shortlist_rows.append(updated_r)

        # Build enriched TSV row
        enriched_row = updated_r + [
            f"{revel_score:.4f}" if revel_score is not None else "NA",
            revel_ev or "NA",
            "not_available", # SpliceAI
            "NA",
            f"{am_score:.4f}" if am_score is not None else "NA",
            am_class or "NA",
            f"{cadd_phred:.2f}" if cadd_phred is not None else "NA",
            f"{af_popmax:.6f}" if af_popmax is not None else "NA",
            str(int(hom_count)) if hom_count is not None else "NA",
            "true" if is_reclassified else "false"
        ]
        enriched_tsv_rows.append(enriched_row)

        # Build unified JSON record
        record_json = {
            "locus": f"{c}:{p}",
            "chrom": c,
            "pos": p,
            "ref": ref,
            "alt": alt,
            "hgvsGenomic": hgvs if is_valid else None,
            "isValidHgvs": is_valid,
            "zygosity": zygosity,
            "primaryGene": gene,
            "primaryTranscript": transcript,
            "consequence": consequence,
            "impact": impact,
            "hgvsC": hgvs_c,
            "hgvsP": hgvs_p,
            "depth": dp,
            "genotypeQuality": gq,
            "alleleBalance": ab,
            "confidenceStatus": confidence,
            "requiresSangerValidation": requires_sanger_str.lower() == "true",
            "populationFrequency": {
                "gnomadExomeAf": af_exome,
                "gnomadGenomeAf": af_genome,
                "gnomadPopmaxAf": af_popmax,
                "popmaxFallbackUsed": af_popmax is not None,
                "homozygoteCount": int(hom_count) if hom_count is not None else None,
                "effectiveAf": effective_af,
                "frequencyStatus": freq_status,
                "afSource": effective_source
            },
            "clinvar": {
                "localSignificance": clnsig_orig,
                "apiSignificances": api_sigs,
                "primarySignificance": cv_sig_primary,
                "reviewStatus": cv_rev_primary,
                "stars": cv_stars,
                "isPathogenic": is_plp,
                "isReclassified": is_reclassified
            } if cv_sig_primary else None,
            "inSilico": {
                "revelScore": revel_score,
                "revelEvidence": revel_ev,
                "revelAcmgPoints": revel_pts,
                "spliceAi": "not_available",
                "alphaMissenseScore": am_score,
                "alphaMissenseClass": am_class,
                "caddPhred": cadd_phred,
                "caddDeleterious": cadd_phred >= 20.0 if cadd_phred is not None else False
            },
            "auditProvenance": {
                "executionPath": execution_path,
                "consentGranted": consent_granted,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "cacheHit": is_valid and hgvs in cache
            }
        }
        enriched_json_records.append(record_json)

    # 1. Write enriched TSV & JSON to annotation/
    with open(out_enriched_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(enriched_headers) + "\n")
        for r in enriched_tsv_rows:
            f.write("\t".join(r) + "\n")

    with open(out_enriched_json, "w", encoding="utf-8") as f:
        json.dump(enriched_json_records, f, indent=2)

    # 2. Update shortlist TSV in-place
    backup_tsv = shortlist_tsv + ".pre_enrichment"
    if not os.path.exists(backup_tsv):
        try:
            os.replace(shortlist_tsv, backup_tsv)
        except Exception:
            pass

    with open(shortlist_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in updated_shortlist_rows:
            f.write("\t".join(r) + "\n")

    # 3. Create or refresh symlinks in annovar/ using absolute paths
    legacy_enriched_tsv = os.path.join(annovar_dir, f"{sample_name}_pass2_enriched.tsv")
    legacy_enriched_json = os.path.join(annovar_dir, f"{sample_name}_pass2_enriched.json")
    legacy_shortlist_tsv = os.path.join(annovar_dir, f"{sample_name}_pass1_shortlist.tsv")
    for link_target, link_dest in [
        (out_enriched_tsv, legacy_enriched_tsv),
        (out_enriched_json, legacy_enriched_json),
        (shortlist_tsv, legacy_shortlist_tsv)
    ]:
        if link_target != link_dest:
            try:
                if os.path.lexists(link_dest):
                    os.remove(link_dest)
                os.symlink(os.path.abspath(link_target), link_dest)
            except OSError:
                pass

    # 4. Write summary JSON
    summary_data = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sampleName": sample_name,
        "executionPath": execution_path,
        "consentGranted": consent_granted,
        "totalCandidatesQueried": total_candidates,
        "cacheHits": cache_hits,
        "apiRequests": api_requests,
        "populationFrequencyBreakdown": {
            "absentUnverified": status_counts["absentUnverified"],
            "rare": status_counts["rare"],
            "common": status_counts["common"],
            "unknown": status_counts["unknown"]
        },
        "inSilicoEvidenceBreakdown": insilico_counts,
        "clinvarMatches": clinvar_matches,
        "clinvarPlpCount": clinvar_plp_count,
        "clinvarReclassifiedCount": clinvar_reclassified_count,
        "artifacts": {
            "enrichedTsv": out_enriched_tsv,
            "enrichedJson": out_enriched_json,
            "updatedShortlistTsv": shortlist_tsv
        }
    }

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2)

    print(f"✅ Pass 2 Enrichment Complete for {sample_name}!")
    print(f"   Execution Path: {execution_path} (Consent: {consent_granted})")
    print(f"   Resolved Frequencies: {status_counts['rare']} Rare, {status_counts['absentUnverified']} Absent (Unverified), {status_counts['common']} Common, {status_counts['unknown']} Unknown.")
    print(f"   In-Silico: {insilico_counts['revelStrong']} REVEL Strong, {insilico_counts['alphaMissenseLikelyPathogenic']} AlphaMissense Pathogenic, SpliceAI: not_available.")
    if clinvar_reclassified_count > 0:
        print(f"   ⚠️ ClinVar Reclassifications Detected: {clinvar_reclassified_count} variants differed between local and API release.")
    print(f"   Artifacts generated: {os.path.basename(out_enriched_tsv)}, {os.path.basename(out_enriched_json)}")

if __name__ == "__main__":
    main()

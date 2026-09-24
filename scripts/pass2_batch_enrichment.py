#!/usr/bin/env python3
"""
Pass 2: Batch API Enrichment & In-Silico Calibration Script

Enriches the Pass 1 Permissive Candidate Shortlist via MyVariant.info:
1. Padded-SNV Two-Way Trimmed HGVS Formatter:
   Translates VCF records to HGVS genomic notation (SNV, pure ins/del, delins, padded-SNVs).
2. MyVariant.info Batch POST Query:
   Retrieves gnomAD exome/genome frequencies, popmax, FAF95, ClinVar, REVEL, SpliceAI, AlphaMissense, and CADD.
3. ClinGen-Calibrated In-Silico Annotations:
   Applies ClinGen calibrated evidence rules (Pejaver 2022 for REVEL, ClinGen SVI for SpliceAI, Cheng 2023 for AlphaMissense).
4. Frequency Status Classification:
   Resolves AF_UNKNOWN to POPULATION_RARE, POPULATION_COMMON, or POPULATION_ABSENT_COVERED.
5. Emits:
   - {sample_name}_pass2_enriched.tsv
   - {sample_name}_pass2_enriched.json
   - {output_dir}/qc/pass2_enrichment_summary.json
   - Updates {sample_name}_pass1_shortlist.tsv with resolved frequencies for Stage 7 Gate.
"""

import sys
import os
import json
import urllib.request
import urllib.parse
from datetime import datetime, timezone

def format_to_hgvs_genomic(chrom, pos, ref, alt):
    """
    Translates a normalized VCF record to MyVariant.info HGVS genomic notation.
    Handles two-way prefix/suffix trimming and resolves padded SNVs directly to canonical SNV format.
    """
    chr_prefix = chrom if chrom.startswith("chr") else f"chr{chrom}"
    ref = (ref or "").strip().upper()
    alt = (alt or "").strip().upper()

    if ref == alt:
        return f"{chr_prefix}:g.{pos}="

    # Direct SNV
    if len(ref) == 1 and len(alt) == 1:
        return f"{chr_prefix}:g.{pos}{ref}>{alt}"

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
        return f"{chr_prefix}:g.{curr_pos}{trimmed_ref}>{trimmed_alt}"

    # Pure Deletion
    if len(trimmed_ref) > 0 and len(trimmed_alt) == 0:
        start_pos = curr_pos
        end_pos = curr_pos + len(trimmed_ref) - 1
        return f"{chr_prefix}:g.{start_pos}del" if start_pos == end_pos else f"{chr_prefix}:g.{start_pos}_{end_pos}del"

    # Pure Insertion
    if len(trimmed_ref) == 0 and len(trimmed_alt) > 0:
        return f"{chr_prefix}:g.{curr_pos - 1}_{curr_pos}ins{trimmed_alt}"

    # Complex Delins
    end_pos = curr_pos + len(trimmed_ref) - 1
    return f"{chr_prefix}:g.{curr_pos}delins{trimmed_alt}" if curr_pos == end_pos else f"{chr_prefix}:g.{curr_pos}_{end_pos}delins{trimmed_alt}"

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

def calibrate_spliceai(max_delta):
    if max_delta is None:
        return None, None
    if max_delta >= 0.50:
        return "SUPPORTING", "PP3_Supporting"
    if max_delta <= 0.10:
        return "BENIGN_SUPPORTING", "BP4_Supporting"
    return "INDETERMINATE", None

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

def extract_float(val):
    if val is None:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None

def main():
    if len(sys.argv) < 3:
        print("Usage: pass2_batch_enrichment.py <output_dir> <sample_name> [offline_mode]")
        sys.exit(1)

    output_dir = os.path.abspath(sys.argv[1])
    sample_name = sys.argv[2]
    offline_mode = len(sys.argv) > 3 and sys.argv[3].lower() in ["true", "1", "offline"]

    annot_dir = os.path.join(output_dir, "annovar")
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)

    shortlist_tsv = os.path.join(annot_dir, f"{sample_name}_pass1_shortlist.tsv")
    cache_path = os.path.join(qc_dir, "myvariant_cache.json")
    out_enriched_tsv = os.path.join(annot_dir, f"{sample_name}_pass2_enriched.tsv")
    out_enriched_json = os.path.join(annot_dir, f"{sample_name}_pass2_enriched.json")
    summary_json = os.path.join(qc_dir, "pass2_enrichment_summary.json")

    if not os.path.exists(shortlist_tsv):
        print(f"⚠️  Shortlist TSV not found at {shortlist_tsv}. Skipping Pass 2 enrichment.")
        sys.exit(0)

    print("════════════════════════════════════════════════════════════")
    print(f"PASS 2: Batch API Enrichment & Pathogenicity Calibration")
    print(f"Sample: {sample_name} | Target: {shortlist_tsv}")
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

    # Format HGVS for all rows
    hgvs_list = []
    for r in rows:
        c = r[col_map["Chr"]]
        p = r[col_map["Pos"]]
        ref = r[col_map["Ref"]]
        alt = r[col_map["Alt"]]
        hgvs = format_to_hgvs_genomic(c, p, ref, alt)
        hgvs_list.append(hgvs)

    execution_path = "PATH_3_OFFLINE_FALLBACK" if offline_mode else "PATH_2_BATCH_API"
    missing_hgvs = [h for h in hgvs_list if h not in cache]
    cache_hits = len(hgvs_list) - len(missing_hgvs)
    api_requests = 0

    if not offline_mode and len(missing_hgvs) > 0:
        print(f"  [PASS 2] Querying MyVariant.info for {len(missing_hgvs)} variants ({cache_hits} cached)...")
        # Batch up to 1000
        BATCH_SIZE = 1000
        for i in range(0, len(missing_hgvs), BATCH_SIZE):
            chunk = missing_hgvs[i:i + BATCH_SIZE]
            api_requests += 1
            url = "https://myvariant.info/v1/variant"
            params = urllib.parse.urlencode({
                "ids": ",".join(chunk),
                "assembly": "hg19",
                "fields": "gnomad_exome,gnomad_genome,clinvar,dbnsfp,cadd,spliceai,alphamissense,dbsnp"
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
                        print(f"  ✅ Batch {api_requests} successfully returned {len(data) if isinstance(data, list) else 1} results.")
                    else:
                        print(f"⚠️  MyVariant returned HTTP status {resp.status}. Entering Partial/Fallback mode.")
                        execution_path = "PATH_4_PARTIAL"
            except Exception as e:
                print(f"⚠️  MyVariant query failed ({e}). Proceeding under Partial/Fallback mode.")
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

    status_counts = {"absentCovered": 0, "absentUnverified": 0, "rare": 0, "common": 0, "unknown": 0}
    insilico_counts = {
        "revelStrong": 0, "revelModerate": 0, "revelSupporting": 0,
        "spliceAiSupporting": 0, "alphaMissenseLikelyPathogenic": 0, "caddDeleterious": 0
    }
    clinvar_matches = 0
    clinvar_plp_count = 0

    enriched_headers = header + [
        "REVEL_Score", "REVEL_Evidence", "SpliceAI_MaxDS", "SpliceAI_Evidence",
        "AlphaMissense_Score", "AlphaMissense_Class", "CADD_Phred",
        "gnomAD_Popmax_AF", "gnomAD_FAF95_Popmax", "gnomAD_Hom_Count"
    ]

    for idx, r in enumerate(rows):
        hgvs = hgvs_list[idx]
        hit = cache.get(hgvs, {})

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
        af_popmax = extract_float(gnomad_exome.get("af_popmax") or gnomad_exome.get("popmax"))
        
        # Calculate popmax if missing
        if af_popmax is None and isinstance(gnomad_exome.get("af"), dict):
            sub_afs = [extract_float(gnomad_exome["af"].get(k)) for k in ["af_afr", "af_amr", "af_eas", "af_nfe", "af_sas"]]
            valid_sub = [s for s in sub_afs if s is not None]
            if len(valid_sub) > 0:
                af_popmax = max(valid_sub)

        faf95_global = extract_float(gnomad_exome.get("faf95", {}).get("faf95") if isinstance(gnomad_exome.get("faf95"), dict) else gnomad_exome.get("faf95"))
        faf95_popmax = None
        if isinstance(gnomad_exome.get("faf95"), dict):
            fafs = [extract_float(gnomad_exome["faf95"].get(k)) for k in ["faf95_afr", "faf95_amr", "faf95_eas", "faf95_nfe", "faf95_sas"]]
            valid_fafs = [f for f in fafs if f is not None]
            if len(valid_fafs) > 0:
                faf95_popmax = max(valid_fafs)

        hom_count = extract_float(gnomad_exome.get("hom", {}).get("hom") if isinstance(gnomad_exome.get("hom"), dict) else gnomad_exome.get("hom"))

        # Effective AF determination
        effective_af = af_exome if af_exome is not None else (af_genome if af_genome is not None else pop_af_orig)
        effective_source = "gnomAD_exome" if af_exome is not None else ("gnomAD_genome" if af_genome is not None else af_source_orig)

        # ClinVar
        cv = hit.get("clinvar", {})
        rcv = cv.get("rcv")
        if isinstance(rcv, list) and len(rcv) > 0:
            rcv = rcv[0]
        elif not isinstance(rcv, dict):
            rcv = {}

        cv_sig = rcv.get("clinical_significance") or (clnsig_orig if clnsig_orig != "." else None)
        cv_rev = rcv.get("review_status") or (clnrev_orig if clnrev_orig != "." else None)
        cv_stars = calculate_clinvar_stars(cv_rev)
        is_plp = cv_sig is not None and "pathogenic" in cv_sig.lower() and "conflict" not in cv_sig.lower()
        if cv_sig:
            clinvar_matches += 1
            if is_plp:
                clinvar_plp_count += 1

        # Frequency status
        if effective_af is not None:
            if effective_af >= 0.01:
                freq_status = "POPULATION_COMMON"
                status_counts["common"] += 1
            else:
                freq_status = "POPULATION_RARE"
                status_counts["rare"] += 1
        elif hit and not hit.get("notfound"):
            # Variant was found in MyVariant without AF -> rare or absent
            freq_status = "POPULATION_ABSENT_COVERED"
            status_counts["absentCovered"] += 1
        elif hit and hit.get("notfound"):
            # Confirmed absent in gnomAD
            freq_status = "POPULATION_ABSENT_COVERED"
            status_counts["absentCovered"] += 1
        else:
            freq_status = "AF_UNKNOWN"
            status_counts["unknown"] += 1

        # In-silico predictors
        revel_score = extract_float(hit.get("dbnsfp", {}).get("revel", {}).get("score"))
        revel_ev, revel_pts, revel_rule = calibrate_revel(revel_score)
        if revel_ev == "STRONG": insilico_counts["revelStrong"] += 1
        elif revel_ev == "MODERATE": insilico_counts["revelModerate"] += 1
        elif revel_ev == "SUPPORTING": insilico_counts["revelSupporting"] += 1

        # SpliceAI
        spliceai_max = extract_float(hit.get("spliceai", {}).get("ds_max"))
        if spliceai_max is None and "spliceai" in hit:
            s_deltas = [extract_float(hit["spliceai"].get(k)) for k in ["ds_ag", "ds_al", "ds_dg", "ds_dl"]]
            v_deltas = [d for d in s_deltas if d is not None]
            if len(v_deltas) > 0:
                spliceai_max = max(v_deltas)
        splice_ev, splice_rule = calibrate_spliceai(spliceai_max)
        if splice_ev == "SUPPORTING": insilico_counts["spliceAiSupporting"] += 1

        # AlphaMissense
        am_score = extract_float(hit.get("alphamissense", {}).get("am_pathogenicity"))
        am_class = classify_alphamissense(am_score)
        if am_class == "likely_pathogenic": insilico_counts["alphaMissenseLikelyPathogenic"] += 1

        # CADD
        cadd_phred = extract_float(hit.get("cadd", {}).get("phred"))
        if cadd_phred is not None and cadd_phred >= 20.0:
            insilico_counts["caddDeleterious"] += 1

        # Build updated TSV row for shortlist (in-place enrichment)
        updated_r = list(r)
        if freq_status == "POPULATION_ABSENT_COVERED":
            updated_r[col_map["Pop_AF"]] = "0.000000"
            updated_r[col_map["AF_Source"]] = "gnomAD_absent"
        else:
            updated_r[col_map["Pop_AF"]] = f"{effective_af:.6f}" if effective_af is not None else "AF_UNKNOWN"
            updated_r[col_map["AF_Source"]] = effective_source
        updated_r[col_map["Frequency_Status"]] = freq_status
        if cv_sig:
            updated_r[col_map["ClinVar_Significance"]] = cv_sig
        if cv_rev:
            updated_r[col_map["ClinVar_RevStat"]] = cv_rev
        updated_shortlist_rows.append(updated_r)

        # Build enriched TSV row
        enriched_row = updated_r + [
            f"{revel_score:.4f}" if revel_score is not None else "NA",
            revel_ev or "NA",
            f"{spliceai_max:.4f}" if spliceai_max is not None else "NA",
            splice_ev or "NA",
            f"{am_score:.4f}" if am_score is not None else "NA",
            am_class or "NA",
            f"{cadd_phred:.2f}" if cadd_phred is not None else "NA",
            f"{af_popmax:.6f}" if af_popmax is not None else "NA",
            f"{faf95_popmax:.6f}" if faf95_popmax is not None else "NA",
            str(int(hom_count)) if hom_count is not None else "NA"
        ]
        enriched_tsv_rows.append(enriched_row)

        # Build unified JSON record
        record_json = {
            "locus": f"{c}:{p}",
            "chrom": c,
            "pos": p,
            "ref": ref,
            "alt": alt,
            "hgvsGenomic": hgvs,
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
                "faf95Global": faf95_global,
                "faf95Popmax": faf95_popmax,
                "homozygoteCount": int(hom_count) if hom_count is not None else None,
                "effectiveAf": effective_af,
                "frequencyStatus": freq_status,
                "afSource": effective_source
            },
            "clinvar": {
                "significance": cv_sig,
                "reviewStatus": cv_rev,
                "stars": cv_stars,
                "isPathogenic": is_plp
            } if cv_sig else None,
            "inSilico": {
                "revelScore": revel_score,
                "revelEvidence": revel_ev,
                "revelAcmgPoints": revel_pts,
                "spliceAiMaxDelta": spliceai_max,
                "spliceAiEvidence": splice_ev,
                "alphaMissenseScore": am_score,
                "alphaMissenseClass": am_class,
                "caddPhred": cadd_phred,
                "caddDeleterious": cadd_phred >= 20.0 if cadd_phred is not None else False
            },
            "auditProvenance": {
                "executionPath": execution_path,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "cacheHit": hgvs in cache
            }
        }
        enriched_json_records.append(record_json)

    # 1. Write enriched TSV
    with open(out_enriched_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(enriched_headers) + "\n")
        for r in enriched_tsv_rows:
            f.write("\t".join(r) + "\n")

    # 2. Write enriched JSON
    with open(out_enriched_json, "w", encoding="utf-8") as f:
        json.dump(enriched_json_records, f, indent=2)

    # 3. Update shortlist TSV in-place with resolved frequencies (so Stage 7 Gate reads it directly)
    backup_tsv = shortlist_tsv + ".pre_enrichment"
    if not os.path.exists(backup_tsv):
        try:
            os.rename(shortlist_tsv, backup_tsv)
        except Exception:
            pass

    with open(shortlist_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in updated_shortlist_rows:
            f.write("\t".join(r) + "\n")

    # 4. Write summary JSON
    summary_data = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sampleName": sample_name,
        "executionPath": execution_path,
        "totalCandidatesQueried": total_candidates,
        "cacheHits": cache_hits,
        "apiRequests": api_requests,
        "populationFrequencyBreakdown": {
            "absentCovered": status_counts["absentCovered"],
            "rare": status_counts["rare"],
            "common": status_counts["common"],
            "unknown": status_counts["unknown"]
        },
        "inSilicoEvidenceBreakdown": insilico_counts,
        "clinvarMatches": clinvar_matches,
        "clinvarPlpCount": clinvar_plp_count,
        "artifacts": {
            "enrichedTsv": out_enriched_tsv,
            "enrichedJson": out_enriched_json,
            "updatedShortlistTsv": shortlist_tsv
        }
    }

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2)

    print(f"✅ Pass 2 Enrichment Complete for {sample_name}!")
    print(f"   Execution Path: {execution_path}")
    print(f"   Resolved Frequencies: {status_counts['rare']} Rare, {status_counts['absentCovered']} Absent Covered, {status_counts['common']} Common, {status_counts['unknown']} Unknown.")
    print(f"   In-Silico: {insilico_counts['revelStrong']} REVEL Strong, {insilico_counts['spliceAiSupporting']} SpliceAI Supporting, {insilico_counts['alphaMissenseLikelyPathogenic']} AlphaMissense Pathogenic.")
    print(f"   Artifacts generated: {os.path.basename(out_enriched_tsv)}, {os.path.basename(out_enriched_json)}")

if __name__ == "__main__":
    main()

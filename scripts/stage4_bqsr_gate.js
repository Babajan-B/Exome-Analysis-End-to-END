#!/usr/bin/env node
/**
 * Inline Stage 4 BQSR Supervisor Gate & Reasoning Engine (Layer 3 Conductor)
 *
 * Invoked by ULTIMATE_MASTER_PIPELINE.sh directly after GATK BaseRecalibrator & ApplyBQSR.
 * Implements the 3-Tier BQSR & Covariate Drift Triage against qc.json:
 *   1. Tier 1: Clinical Grade Pass (Observations >= 500M, Drift <= 4.0 Phred) -> Exits 0
 *   2. Tier 2: Research Grade Qualified (Observations >= 100M, Drift <= 6.0 Phred) -> Exits 0 (FLAG_AND_CONTINUE)
 *   3. Tier 2: Research Bypass (Known-Sites missing, dedup.bam passed to caller) -> Exits 0 (FLAG_AND_CONTINUE)
 *   4. Tier 3: Floor Failure (Observations < 50M or Drift > 8.0 Phred) -> Exits 1 (Quality Halt with Operator Opinion Gate)
 *   5. Tool Crash / Corrupt Outputs -> Exits 2
 */
const fs = require("fs");
const path = require("path");

process.on("uncaughtException", (err) => {
  console.error("❌ [SUPERVISOR TOOL CRASH] Uncaught exception in Stage 4 BQSR gate:", err);
  process.exit(2);
});

const outputDir = process.argv[2];
const sampleName = process.argv[3] || "sample";
const referencePathArg = process.argv[4] || process.env.REFERENCE || null;

if (!outputDir) {
  console.error("Usage: node stage4_bqsr_gate.js <output_dir> [sample_name]");
  process.exit(2);
}

const recalTablePath = path.join(outputDir, "bqsr", "recal_data.table");
const recalBamPath = path.join(outputDir, "bqsr", "recal.bam");
const dedupBamPath = path.join(outputDir, "dedup", "dedup.bam");
const bypassedMarkerPath = path.join(outputDir, ".bqsr_bypassed");
const overrideMarkerPath = path.join(outputDir, ".override_bqsr_gate");
const auditTrailPath = path.join(outputDir, "audit_trail.json");
const haltReportPath = path.join(outputDir, "halt_report.json");
const reasoningPath = path.join(outputDir, "stage4_bqsr_reasoning.json");
const stage3ReasoningPath = path.join(outputDir, "stage3_dedup_reasoning.json");

// 1. Load qc.json (Single Source of Truth)
let qcPolicy = null;
const candidateQcPaths = [
  path.join(outputDir, "qc.json"),
  path.resolve(__dirname, "../../../../qc.json"),
  path.resolve(__dirname, "../../../qc.json"),
  path.resolve(process.cwd(), "qc.json"),
];
for (const p of candidateQcPaths) {
  try {
    if (fs.existsSync(p)) {
      qcPolicy = JSON.parse(fs.readFileSync(p, "utf8"));
      break;
    }
  } catch {}
}

// 2. Load BQSR Triage Engine
let bqsrTriageEngine = null;
const candidateEnginePaths = [
  path.resolve(__dirname, "../../../../lib/exome/bqsr-triage-engine.js"),
  path.resolve(__dirname, "../../../lib/exome/bqsr-triage-engine.js"),
  path.resolve(process.cwd(), "lib/exome/bqsr-triage-engine.js"),
];
for (const p of candidateEnginePaths) {
  try {
    if (fs.existsSync(p)) {
      bqsrTriageEngine = require(p);
      break;
    }
  } catch {}
}

if (!bqsrTriageEngine) {
  console.error("❌ [SUPERVISOR TOOL CRASH] Failed to load bqsr-triage-engine.js");
  process.exit(2);
}

// 3. Load Stage 3 Directives
let stage3Directives = null;
try {
  if (fs.existsSync(stage3ReasoningPath)) {
    const s3 = JSON.parse(fs.readFileSync(stage3ReasoningPath, "utf8"));
    stage3Directives = s3.downstreamDirectives || null;
  }
} catch {}

// 4. Check for Operator Override
const hasOverride = fs.existsSync(overrideMarkerPath);

// 5. Check if BQSR was bypassed
const isBypassed = fs.existsSync(bypassedMarkerPath);

// 6. Read pipeline known_sites.txt if available
const knownSitesFile = path.join(outputDir, "bqsr", "known_sites.txt");
let pipelineKnownSites = [];
if (fs.existsSync(knownSitesFile)) {
  try {
    pipelineKnownSites = fs
      .readFileSync(knownSitesFile, "utf8")
      .split(/\r?\n/)
      .map((s) => s.trim())
      .filter(Boolean)
      .map((p) => path.basename(p));
  } catch {}
}

// 7. Read reference build if available
const refBuildFile = path.join(outputDir, "bqsr", "ref_build.txt");
let refBuild = "hg19";
if (fs.existsSync(refBuildFile)) {
  try {
    refBuild = fs.readFileSync(refBuildFile, "utf8").trim() || "hg19";
  } catch {}
}

// 8. Verify Contig Format Parity (Bidirectional: BAM alignment_idxstats.txt vs Reference Genome)
let contigParity = {
  namingConvention: "unknown",
  hasChrPrefix: false,
  verified: false,
  bam: {
    convention: "unknown",
    hasChrPrefix: false,
    sampleContigs: [],
  },
  reference: {
    convention: "unknown",
    hasChrPrefix: false,
    source: null,
    sampleContigs: [],
  },
  parityStatus: "UNVERIFIED",
  isParityMatch: true,
};

// Side 1: Inspect BAM contigs from alignment_idxstats.txt
const idxstatsPath = path.join(outputDir, "aligned", "alignment_idxstats.txt");
if (fs.existsSync(idxstatsPath)) {
  try {
    const idxContent = fs.readFileSync(idxstatsPath, "utf8");
    const rawContigs = idxContent
      .split(/\r?\n/)
      .map((line) => line.split(/\s+/)[0])
      .filter((c) => c && c !== "*");
    if (rawContigs.length > 0) {
      contigParity.bam.sampleContigs = rawContigs.slice(0, 5);
      const hasChr = rawContigs.some((c) => c.startsWith("chr"));
      contigParity.bam.hasChrPrefix = hasChr;
      contigParity.bam.convention = hasChr ? "chr_prefixed" : "bare_numeric";
      // Backwards-compatible top-level keys
      contigParity.namingConvention = contigParity.bam.convention;
      contigParity.hasChrPrefix = hasChr;
      contigParity.verified = true;
    }
  } catch {}
}

// Side 2: Inspect Reference Genome contigs (.dict / .fai / FASTA)
function extractContigsFromFile(filePath) {
  if (!filePath || !fs.existsSync(filePath)) return null;
  try {
    const stat = fs.statSync(filePath);
    if (stat.size === 0) return null;
    const content = fs.readFileSync(filePath, "utf8");
    if (filePath.endsWith(".dict")) {
      const matches = [];
      const regex = /@SQ\s+SN:([^\s\t]+)/g;
      let m;
      while ((m = regex.exec(content)) !== null) {
        matches.push(m[1]);
        if (matches.length >= 10) break;
      }
      return matches.length > 0 ? matches : null;
    } else if (filePath.endsWith(".fai")) {
      const lines = content.split(/\r?\n/).filter(Boolean);
      const names = lines.map((l) => l.split(/\t/)[0]).filter((n) => n && n !== "*");
      return names.length > 0 ? names.slice(0, 10) : null;
    } else if (filePath.endsWith(".fa") || filePath.endsWith(".fasta")) {
      const matches = [];
      const regex = /^>([^\s\r\n]+)/gm;
      let m;
      while ((m = regex.exec(content)) !== null) {
        matches.push(m[1]);
        if (matches.length >= 10) break;
      }
      return matches.length > 0 ? matches : null;
    }
  } catch {}
  return null;
}

let refContigs = [];
let refSource = null;

if (referencePathArg) {
  const dictCand = referencePathArg.replace(/\.(fa|fasta)$/, ".dict");
  const faiCand = referencePathArg + ".fai";
  const faiCand2 = referencePathArg.replace(/\.(fa|fasta)$/, ".fai");

  for (const cand of [dictCand, faiCand, faiCand2, referencePathArg]) {
    const extracted = extractContigsFromFile(cand);
    if (extracted && extracted.length > 0) {
      refContigs = extracted;
      refSource = cand;
      break;
    }
  }
}

if (refContigs.length === 0) {
  const candidateRefDirs = [
    path.join(outputDir, "reference"),
    path.resolve(outputDir, "..", "reference"),
    path.resolve(outputDir, "../..", "reference"),
    path.resolve(__dirname, "../reference"),
    path.resolve(__dirname, "../../reference"),
  ];
  for (const dir of candidateRefDirs) {
    if (!fs.existsSync(dir)) continue;
    try {
      const files = fs.readdirSync(dir);
      const dictFile = files.find((f) => f.endsWith(".dict"));
      const faiFile = files.find((f) => f.endsWith(".fai"));
      const faFile = files.find((f) => f.endsWith(".fa") || f.endsWith(".fasta"));
      const best = dictFile || faiFile || faFile;
      if (best) {
        const fullPath = path.join(dir, best);
        const extracted = extractContigsFromFile(fullPath);
        if (extracted && extracted.length > 0) {
          refContigs = extracted;
          refSource = fullPath;
          break;
        }
      }
    } catch {}
  }
}

if (refContigs.length > 0) {
  contigParity.reference.sampleContigs = refContigs.slice(0, 5);
  contigParity.reference.source = refSource;
  const refHasChr = refContigs.some((c) => c.startsWith("chr"));
  contigParity.reference.hasChrPrefix = refHasChr;
  contigParity.reference.convention = refHasChr ? "chr_prefixed" : "bare_numeric";

  if (contigParity.bam.convention !== "unknown") {
    const match = contigParity.bam.hasChrPrefix === refHasChr;
    contigParity.isParityMatch = match;
    contigParity.parityStatus = match ? "PARITY_CONFIRMED" : "MISMATCH_DETECTED";
  } else {
    contigParity.parityStatus = "REFERENCE_ONLY_VERIFIED";
  }
} else {
  contigParity.parityStatus = contigParity.bam.convention !== "unknown" ? "BAM_ONLY_VERIFIED" : "UNVERIFIED";
  contigParity.isParityMatch = true;
}

let parsedMetrics = {
  totalObservations: 0,
  totalErrors: 0,
  empiricalQuality: 0,
  estimatedQuality: 0,
  meanQualityDrift: 0,
  readGroups: [],
  qualityMapping: [],
  knownSites: [],
  arguments: {},
  isValid: false,
};

if (!isBypassed) {
  // If not bypassed, recal_data.table must exist
  if (!fs.existsSync(recalTablePath)) {
    console.error(`❌ [SUPERVISOR TOOL CRASH] Missing recalibration table at: ${recalTablePath}`);
    process.exit(2);
  }

  let tableContent = "";
  try {
    tableContent = fs.readFileSync(recalTablePath, "utf8");
  } catch (err) {
    console.error("❌ [SUPERVISOR TOOL CRASH] Failed to read recal_data.table:", err.message);
    process.exit(2);
  }

  parsedMetrics = bqsrTriageEngine.parseRecalTable(tableContent);
  // If pipeline recorded known-sites explicitly, use that authoritative list
  if (pipelineKnownSites.length > 0) {
    parsedMetrics.knownSites = pipelineKnownSites;
  }

  if (!parsedMetrics.isValid && !hasOverride) {
    console.error("❌ [SUPERVISOR TOOL CRASH] recal_data.table was empty or contained 0 base observations.");
    process.exit(2);
  }

  // Verify recal.bam presence and basic size sanity (> 100 bytes)
  if (!fs.existsSync(recalBamPath) && !hasOverride) {
    console.error(`❌ [SUPERVISOR TOOL CRASH] Missing recalibrated BAM at: ${recalBamPath}`);
    process.exit(2);
  }
  if (fs.existsSync(recalBamPath)) {
    const bamStat = fs.statSync(recalBamPath);
    if (bamStat.size < 100 && !hasOverride) {
      console.error(`❌ [SUPERVISOR TOOL CRASH] recal.bam is truncated/empty (${bamStat.size} bytes).`);
      process.exit(2);
    }
  }
} else {
  // Bypassed mode: assign known-sites if recorded
  if (pipelineKnownSites.length > 0) {
    parsedMetrics.knownSites = pipelineKnownSites;
  }
}

// 9. Optional two-pass post-recalibration residual evaluation
let postRecalData = null;
const postRecalPath = path.join(outputDir, "bqsr", "post_recal_data.table");
if (fs.existsSync(postRecalPath)) {
  try {
    const postContent = fs.readFileSync(postRecalPath, "utf8");
    const parsedPost = bqsrTriageEngine.parseRecalTable(postContent);
    if (parsedPost.isValid) {
      postRecalData = {
        meanResidualDrift: parsedPost.meanQualityDrift,
        empiricalQuality: parsedPost.empiricalQuality,
        mode: "two_pass_evaluated",
      };
    }
  } catch {}
}

// 10. Evaluate BQSR Triage
const triage = bqsrTriageEngine.evaluateBqsrTriage(parsedMetrics, {
  sampleName,
  hasOverride,
  isBypassed,
  qcPolicy,
  postRecalResidualDrift: postRecalData ? postRecalData.meanResidualDrift : null,
});

// Enforce Contig Format Parity if mismatch detected between BAM and Reference Genome
const enforceContigParity = qcPolicy?.tiers?.tier_bqsr_recalibration_qc?.contig_format_parity?.enforce_chr_prefix_match !== false;
if (enforceContigParity && !contigParity.isParityMatch && !hasOverride) {
  triage.exitCode = 1;
  triage.tier = "TIER_3_BQSR_HALT";
  if (triage.classification) {
    triage.classification.contigParity = "MISMATCH_DETECTED";
  }
  triage.operatorExplanation = {
    rootCause: "CONTIG_FORMAT_PARITY_MISMATCH",
    details: `Contig naming convention mismatch between BAM (${contigParity.bam.convention}: ${contigParity.bam.sampleContigs.join(", ")}) and Reference Genome (${contigParity.reference.convention}: ${contigParity.reference.sampleContigs.join(", ")}).`,
    sequencerPhysics: "GATK requires exact string parity between BAM sequence headers and reference genome FASTA dictionary. Mismatched contig prefixes (e.g. 'chr1' vs '1') cause silent masking failure of known polymorphic sites or fatal coordinate misalignment.",
    suggestedAction: "Re-header BAM to match reference dictionary or re-align to matching reference build.",
  };
  triage.supervisorThought = `🛑 [Supervisor Cognitive Halt] Sample ${sampleName} failed Contig Format Parity: BAM uses ${contigParity.bam.convention} while Reference uses ${contigParity.reference.convention}.`;
}

// 11. Write stage4_bqsr_reasoning.json
const reasoningPayload = {
  stage: "stage_4_base_quality_score_recalibration",
  sampleName,
  timestamp: new Date().toISOString(),
  tier: triage.tier,
  exitCode: triage.exitCode,
  classification: triage.classification,
  metrics: {
    totalObservations: parsedMetrics.totalObservations,
    totalErrors: parsedMetrics.totalErrors,
    empiricalQuality: parsedMetrics.empiricalQuality,
    estimatedQuality: parsedMetrics.estimatedQuality,
    meanQualityDrift: parsedMetrics.meanQualityDrift,
    readGroupsCount: parsedMetrics.readGroups.length,
    knownSites: parsedMetrics.knownSites,
    isBypassed,
    floorSwitch: triage.floorSwitch,
  },
  floorSwitch: triage.floorSwitch,
  readGroups: parsedMetrics.readGroups,
  qualityMapping: parsedMetrics.qualityMapping.slice(0, 30), // top calibration bins
  contigParity,
  postRecalData: postRecalData || { mode: "single_pass", status: "NOT_EVALUATED_SINGLE_PASS" },
  refBuild,
  operatorExplanation: triage.operatorExplanation,
  downstreamDirectives: triage.downstreamDirectives,
  stage3DirectivesReceived: stage3Directives,
  recommendation: triage.recommendation,
  supervisorThought: triage.supervisorThought,
};

try {
  fs.writeFileSync(reasoningPath, JSON.stringify(reasoningPayload, null, 2));
} catch (err) {
  console.error("⚠️ [SUPERVISOR] Warning: Failed to write stage4_bqsr_reasoning.json:", err.message);
}

// 8. Append to audit_trail.json
try {
  let auditTrail = [];
  if (fs.existsSync(auditTrailPath)) {
    try { auditTrail = JSON.parse(fs.readFileSync(auditTrailPath, "utf8")); } catch {}
  }
  auditTrail.push({
    timestamp: new Date().toISOString(),
    layer: 3,
    agent: "Cognitive Supervisor (Stage 4 BQSR)",
    tier: triage.tier,
    exitCode: triage.exitCode,
    totalObservations: parsedMetrics.totalObservations,
    meanQualityDrift: parsedMetrics.meanQualityDrift,
    message: triage.supervisorThought,
  });
  fs.writeFileSync(auditTrailPath, JSON.stringify(auditTrail, null, 2));
} catch (err) {
  console.error("⚠️ [SUPERVISOR] Warning: Failed to append to audit_trail.json:", err.message);
}

// 9. Handle Exit Conditions
if (triage.exitCode === 0) {
  if (hasOverride) {
    try {
      const appliedPath = overrideMarkerPath + ".applied";
      fs.renameSync(overrideMarkerPath, appliedPath);
    } catch {}
    console.log(`✅ [SUPERVISOR OVERRIDE] Stage 4 BQSR floor failure explicitly overridden by operator.`);
  } else if (isBypassed) {
    console.log(`ℹ️  [SUPERVISOR BYPASS] ${triage.supervisorThought}`);
  } else {
    console.log(`✅ [SUPERVISOR PASS] ${triage.supervisorThought}`);
  }
  process.exit(0);
} else if (triage.exitCode === 1) {
  // Quality Floor Failure: Write halt_report.json
  const haltPayload = {
    halted: true,
    timestamp: new Date().toISOString(),
    phase: "BQSR",
    stage: "stage_4_base_quality_score_recalibration",
    sampleName,
    exitCode: 1,
    faultClass: triage.operatorExplanation?.rootCause || "BQSR_QUALITY_FLOOR_FAILURE",
    diagnostics: {
      metrics: reasoningPayload.metrics,
      classification: triage.classification,
    },
    operatorExplanation: triage.operatorExplanation,
    remediationOptions: [
      {
        action: "ABORT",
        description: "Abort pipeline run. Re-sequence or inspect alignment to reference build.",
        command: "rm -rf " + outputDir,
      },
      {
        action: "OVERRIDE",
        description: "Engage Operator Override (High-Risk Research Protocol). Proceeds to HaplotypeCaller.",
        command: `touch "${overrideMarkerPath}" && restart pipeline`,
      },
    ],
  };

  try {
    fs.writeFileSync(haltReportPath, JSON.stringify(haltPayload, null, 2));
  } catch (err) {
    console.error("⚠️ [SUPERVISOR] Warning: Failed to write halt_report.json:", err.message);
  }

  console.error(`🛑 [SUPERVISOR HALT] ${triage.supervisorThought}`);
  process.exit(1);
} else {
  console.error(`❌ [SUPERVISOR CRASH] Tool failure during BQSR processing.`);
  process.exit(2);
}

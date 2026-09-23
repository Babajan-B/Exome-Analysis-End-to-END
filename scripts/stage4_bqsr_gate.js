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

// 8. Verify Contig Format Parity from alignment_idxstats.txt
let contigParity = {
  namingConvention: "unknown",
  hasChrPrefix: false,
  verified: false,
};
const idxstatsPath = path.join(outputDir, "aligned", "alignment_idxstats.txt");
if (fs.existsSync(idxstatsPath)) {
  try {
    const idxContent = fs.readFileSync(idxstatsPath, "utf8");
    const firstContig = idxContent.split(/\r?\n/)[0]?.split(/\s+/)[0] || "";
    if (firstContig.startsWith("chr")) {
      contigParity = { namingConvention: "chr_prefixed", hasChrPrefix: true, verified: true };
    } else if (/^[0-9XYM]/.test(firstContig)) {
      contigParity = { namingConvention: "bare_numeric", hasChrPrefix: false, verified: true };
    }
  } catch {}
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
});

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
  },
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

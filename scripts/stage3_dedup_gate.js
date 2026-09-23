#!/usr/bin/env node
/**
 * Inline Stage 3 Deduplication Supervisor Gate & Reasoning Engine (Layer 3 Conductor)
 *
 * Invoked by ULTIMATE_MASTER_PIPELINE.sh directly after GATK/Picard MarkDuplicates.
 * Implements the 3-Tier Duplication & Library Complexity Triage against qc.json:
 *   1. Tier 1: Clinical Grade Pass (Duplication <= 12%, Optical <= 2.0%, Library >= 50M) -> Exits 0
 *   2. Tier 2: Research Grade Qualified (Duplication <= 22%, Optical <= 4.5%) -> Exits 0 (FLAG_AND_CONTINUE)
 *   3. Tier 3: Floor Failure (Duplication >= 25% or Library < 10M) -> Exits 1 (Quality Halt with Operator Opinion Gate)
 *   4. Tool Crash / Corrupt Outputs -> Exits 2
 */
const fs = require("fs");
const path = require("path");

process.on("uncaughtException", (err) => {
  console.error("❌ [SUPERVISOR TOOL CRASH] Uncaught exception in Stage 3 Deduplication gate:", err);
  process.exit(2);
});

const outputDir = process.argv[2];
const sampleName = process.argv[3] || "sample";

if (!outputDir) {
  console.error("Usage: node stage3_dedup_gate.js <output_dir> [sample_name]");
  process.exit(2);
}

let metricsPath = path.join(outputDir, "dedup", "metrics.txt");
if (!fs.existsSync(metricsPath)) {
  const altCandidates = [
    path.join(outputDir, "dedup", "marked_dup_metrics.txt"),
    path.join(outputDir, "dedup", "dedup_metrics.txt"),
    path.join(outputDir, "metrics.txt"),
  ];
  for (const alt of altCandidates) {
    if (fs.existsSync(alt)) {
      metricsPath = alt;
      break;
    }
  }
}
const auditTrailPath = path.join(outputDir, "audit_trail.json");
const haltReportPath = path.join(outputDir, "halt_report.json");
const reasoningPath = path.join(outputDir, "stage3_dedup_reasoning.json");
const stage2ReasoningPath = path.join(outputDir, "stage2_alignment_reasoning.json");
const overrideMarkerPath = path.join(outputDir, ".override_dedup_gate");

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

// 2. Load Deduplication Triage Engine
let dedupTriageEngine = null;
const candidateEnginePaths = [
  path.resolve(__dirname, "../../../../lib/exome/dedup-triage-engine.js"),
  path.resolve(__dirname, "../../../lib/exome/dedup-triage-engine.js"),
  path.resolve(process.cwd(), "lib/exome/dedup-triage-engine.js"),
];
for (const p of candidateEnginePaths) {
  try {
    if (fs.existsSync(p)) {
      dedupTriageEngine = require(p);
      break;
    }
  } catch {}
}

if (!dedupTriageEngine) {
  console.error("❌ [SUPERVISOR TOOL CRASH] Failed to load dedup-triage-engine.js");
  process.exit(2);
}

// 3. Load Stage 2 Directives
let stage2Directives = null;
try {
  if (fs.existsSync(stage2ReasoningPath)) {
    const s2 = JSON.parse(fs.readFileSync(stage2ReasoningPath, "utf8"));
    stage2Directives = s2.downstreamDirectives || null;
  }
} catch {}

// 4. Check for Operator Override
const hasOverride = fs.existsSync(overrideMarkerPath);

// 5. Read Picard metrics from disk
let metricsContent = "";
try {
  if (fs.existsSync(metricsPath)) {
    metricsContent = fs.readFileSync(metricsPath, "utf8");
  }
} catch {}

// Parse metrics
let metrics = null;
if (dedupTriageEngine) {
  metrics = dedupTriageEngine.parsePicardMetrics(metricsContent);
}

// Evaluate triage
let result = null;
if (dedupTriageEngine && metrics) {
  result = dedupTriageEngine.evaluateDedupTriage(metrics, {
    sampleName,
    hasOverride,
    stage2Directives,
    qcPolicy,
  });
} else {
  result = {
    tier: "TOOL_CRASH",
    action: "HALT_TOOL_CRASH",
    exitCode: 2,
    sampleName,
    metrics: metrics || {
      library: "unknown",
      unpairedReadsExamined: 0,
      readPairsExamined: 0,
      secondaryOrSupplementaryRds: 0,
      unmappedReads: 0,
      unpairedReadDuplicates: 0,
      readPairDuplicates: 0,
      readPairOpticalDuplicates: 0,
      percentDuplication: 0,
      estimatedLibrarySize: 0,
      opticalDuplicateRate: 0,
      pcrDuplicateRate: 0,
    },
    supervisorThought: "Deduplication triage engine missing or output unparseable.",
    confidenceScore: 0.5,
    provider: "deterministic",
    generatedAt: new Date().toISOString(),
  };
}

// 6. Write Reasoning Artifact
try {
  fs.writeFileSync(reasoningPath, JSON.stringify(result, null, 2));
} catch (e) {
  console.error("Warning: Failed to write stage3_dedup_reasoning.json:", e.message);
}

// 7. Update Audit Trail
try {
  let auditTrail = [];
  if (fs.existsSync(auditTrailPath)) {
    auditTrail = JSON.parse(fs.readFileSync(auditTrailPath, "utf8"));
  }
  auditTrail.push({
    timestamp: new Date().toISOString(),
    stage: "Stage 3: Duplicate Marking",
    layer: 3,
    agent: "Supervisor Conductor (Deduplication Gate)",
    tier: result.tier,
    action: result.action,
    exitCode: result.exitCode,
    percentDuplication: result.metrics.percentDuplication,
    opticalDuplicateRate: result.metrics.opticalDuplicateRate,
    estimatedLibrarySize: result.metrics.estimatedLibrarySize,
    thought: result.supervisorThought,
  });
  fs.writeFileSync(auditTrailPath, JSON.stringify(auditTrail, null, 2));
} catch {}

// 8. If Floor Rejection (Halt), emit Halt Report for Web Dashboard
if (result.exitCode === 1) {
  try {
    const haltReport = {
      jobId: sampleName,
      sampleName,
      halted: true,
      phase: "Deduplication",
      stage: "dedup",
      timestamp: new Date().toISOString(),
      faultClass: "HIGH_DUPLICATION_LIBRARY_EXHAUSTION",
      summary: `Stage 3 Deduplication failed clinical quality floor (${result.metrics.percentDuplication.toFixed(1)}% duplication).`,
      rootCause: result.operatorExplanation?.rootCause || "Excessive duplication exceeding 25% rejection floor.",
      sequencerPhysics: result.operatorExplanation?.sequencerPhysics || "",
      downstreamRisks: result.operatorExplanation?.downstreamRisks || [],
      remediationOptions: result.operatorExplanation?.remediationOptions || [],
      supervisorThought: result.supervisorThought,
    };
    fs.writeFileSync(haltReportPath, JSON.stringify(haltReport, null, 2));
  } catch {}
} else if (result.exitCode === 0) {
  // Clear any previous halt report if this run passed or was overridden
  try {
    if (fs.existsSync(haltReportPath)) {
      const hr = JSON.parse(fs.readFileSync(haltReportPath, "utf8"));
      if (hr?.phase === "Deduplication") {
        fs.unlinkSync(haltReportPath);
      }
    }
  } catch {}
  // If override was applied, archive marker
  if (hasOverride) {
    try {
      fs.renameSync(overrideMarkerPath, overrideMarkerPath + ".applied");
    } catch {}
  }
}

// 9. Console output
console.log("");
console.log("┌────────────────────────────────────────────────────────────────────────┐");
console.log(`│   STAGE 3 SUPERVISOR GATE: ${result.tier.padEnd(43)} │`);
console.log("├────────────────────────────────────────────────────────────────────────┤");
console.log(`│ Sample: ${sampleName.padEnd(62)} │`);
console.log(`│ Duplication Rate:     ${(result.metrics.percentDuplication.toFixed(2) + "%").padEnd(48)} │`);
console.log(`│ Optical Duplicate:    ${(result.metrics.opticalDuplicateRate.toFixed(2) + "%").padEnd(48)} │`);
console.log(`│ Est. Library Size:    ${((result.metrics.estimatedLibrarySize / 1e6).toFixed(1) + "M fragments").padEnd(48)} │`);
console.log(`│ Read Pairs Examined:  ${result.metrics.readPairsExamined.toLocaleString().padEnd(48)} │`);
console.log(`│ Action:               ${result.action.padEnd(48)} │`);
console.log(`│ Exit Code:            ${String(result.exitCode).padEnd(48)} │`);
console.log("└────────────────────────────────────────────────────────────────────────┘");
console.log("");
console.log("Cognitive Rationale:");
console.log(result.supervisorThought);
console.log("");

process.exit(result.exitCode);

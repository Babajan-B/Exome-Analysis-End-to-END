#!/usr/bin/env node
/**
 * Inline Stage 2 Alignment Supervisor Gate & Reasoning Engine (Layer 3 Conductor)
 *
 * Invoked by ULTIMATE_MASTER_PIPELINE.sh directly after BWA-MEM alignment & metrics generation.
 * Implements the 3-Tier Alignment Quality Triage against qc.json (tier_2_alignment_complexity_qc):
 *   1. Tier 1: Clinical Grade Pass (Mapping >= 98%, Properly Paired >= 95%, Chimeric < 1.5%) -> Exits 0
 *   2. Tier 2: Research Grade Qualified (Mapping >= 95%, Properly Paired >= 90%) -> Exits 0 (FLAG_AND_CONTINUE)
 *   3. Tier 3: Floor Failure (Mapping < 90% or Properly Paired < 85%) -> Exits 1 (Quality Halt with Operator Opinion Gate)
 *   4. Tool Crash / Corrupt Outputs -> Exits 2
 */
const fs = require("fs");
const path = require("path");

process.on("uncaughtException", (err) => {
  console.error("❌ [SUPERVISOR TOOL CRASH] Uncaught exception in Stage 2 Alignment gate:", err);
  process.exit(2);
});

const outputDir = process.argv[2];
const sampleName = process.argv[3] || "sample";

if (!outputDir) {
  console.error("Usage: node stage2_align_gate.js <output_dir> [sample_name]");
  process.exit(2);
}

const flagstatPath = path.join(outputDir, "aligned", "alignment_flagstat.txt");
const statsPath = path.join(outputDir, "aligned", "alignment_stats.txt");
const auditTrailPath = path.join(outputDir, "audit_trail.json");
const haltReportPath = path.join(outputDir, "halt_report.json");
const reasoningPath = path.join(outputDir, "stage2_alignment_reasoning.json");
const stage1ReasoningPath = path.join(outputDir, "supervisor_reasoning.json");
const overrideMarkerPath = path.join(outputDir, ".override_align_gate");

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

// 2. Load Alignment Triage Engine
let alignTriageEngine = null;
const candidateEnginePaths = [
  path.resolve(__dirname, "../../../../lib/exome/alignment-triage-engine.js"),
  path.resolve(__dirname, "../../../lib/exome/alignment-triage-engine.js"),
  path.resolve(process.cwd(), "lib/exome/alignment-triage-engine.js"),
];
for (const p of candidateEnginePaths) {
  try {
    if (fs.existsSync(p)) {
      alignTriageEngine = require(p);
      break;
    }
  } catch {}
}

if (!alignTriageEngine) {
  console.error("❌ [SUPERVISOR TOOL CRASH] Failed to load alignment-triage-engine.js");
  process.exit(2);
}

// 3. Load Stage 1 Directives
let stage1Directives = null;
try {
  if (fs.existsSync(stage1ReasoningPath)) {
    const s1 = JSON.parse(fs.readFileSync(stage1ReasoningPath, "utf8"));
    stage1Directives = s1.downstreamDirectives || null;
  }
} catch {}

// 4. Check for Operator Override
const hasOverride = fs.existsSync(overrideMarkerPath);

// 5. Read alignment metrics from disk
let flagstatContent = "";
let statsContent = "";

try {
  if (fs.existsSync(flagstatPath)) {
    flagstatContent = fs.readFileSync(flagstatPath, "utf8");
  }
} catch {}

try {
  if (fs.existsSync(statsPath)) {
    statsContent = fs.readFileSync(statsPath, "utf8");
  }
} catch {}

// Parse metrics
let metrics = null;
if (alignTriageEngine) {
  const fMetrics = alignTriageEngine.parseFlagstat(flagstatContent);
  const sMetrics = alignTriageEngine.parseStats(statsContent);
  metrics = { ...fMetrics, ...sMetrics };
}

// Evaluate triage
let result = null;
if (alignTriageEngine && metrics) {
  result = alignTriageEngine.evaluateAlignmentTriage(metrics, {
    sampleName,
    hasOverride,
    stage1Directives,
    qcPolicy,
  });
} else {
  // Fallback if engine missing
  result = {
    tier: "TOOL_CRASH",
    action: "HALT_TOOL_CRASH",
    exitCode: 2,
    sampleName,
    metrics: metrics || { totalReads: 0, mappedReads: 0, mappingPct: 0, properlyPairedReads: 0, properlyPairedPct: 0, singletons: 0, singletonsPct: 0, secondary: 0, supplementary: 0, duplicates: 0, chimericSplitPct: 0 },
    supervisorThought: "Alignment triage engine missing or output unparseable.",
    confidenceScore: 0.5,
    provider: "deterministic",
    generatedAt: new Date().toISOString(),
  };
}

// 6. Write Reasoning Artifact
try {
  fs.writeFileSync(reasoningPath, JSON.stringify(result, null, 2));
} catch (e) {
  console.error("Warning: Failed to write stage2_alignment_reasoning.json:", e.message);
}

// 7. Update Audit Trail
try {
  let auditTrail = [];
  if (fs.existsSync(auditTrailPath)) {
    auditTrail = JSON.parse(fs.readFileSync(auditTrailPath, "utf8"));
  }
  auditTrail.push({
    timestamp: new Date().toISOString(),
    stage: "Stage 2: Read Alignment",
    layer: 3,
    agent: "Supervisor Conductor (Alignment Gate)",
    tier: result.tier,
    action: result.action,
    exitCode: result.exitCode,
    mappingPct: result.metrics.mappingPct,
    properlyPairedPct: result.metrics.properlyPairedPct,
    thought: result.supervisorThought,
  });
  fs.writeFileSync(auditTrailPath, JSON.stringify(auditTrail, null, 2));
} catch {}

// 8. Console Output
console.log("");
console.log("═══════════════════════════════════════════════════════════════");
console.log(" 🧠 [SUPERVISOR AGENT] STAGE 2 ALIGNMENT COGNITIVE GATE");
console.log("═══════════════════════════════════════════════════════════════");
console.log(` Sample:       ${sampleName}`);
console.log(` Mapping Rate: ${result.metrics.mappingPct.toFixed(2)}% (${result.metrics.mappedReads.toLocaleString()} / ${result.metrics.totalReads.toLocaleString()} reads)`);
console.log(` Paired Rate:  ${result.metrics.properlyPairedPct.toFixed(2)}% (${result.metrics.properlyPairedReads.toLocaleString()} reads)`);
if (result.metrics.insertSizeAverage) {
  console.log(` Insert Size:  ${result.metrics.insertSizeAverage.toFixed(0)} bp ± ${result.metrics.insertSizeStandardDeviation ? result.metrics.insertSizeStandardDeviation.toFixed(0) : 0} bp`);
}
console.log(` Chimeric:     ${result.metrics.chimericSplitPct.toFixed(2)}%`);
console.log("───────────────────────────────────────────────────────────────");

if (result.tier === "TIER_1_CLINICAL_ALIGNMENT") {
  console.log(" 🟢 VERDICT: TIER 1 — CLINICAL GRADE ALIGNMENT PASS");
  console.log(`    ${result.supervisorThought}`);
  console.log("    Direct handoff to MarkDuplicates authorized.");
  console.log("");
  process.exit(0);
}

if (result.tier === "TIER_2_RESEARCH_ALIGNMENT") {
  console.log(" 🟡 VERDICT: TIER 2 — RESEARCH GRADE QUALIFIED (FLAG_AND_CONTINUE)");
  console.log(`    ${result.supervisorThought}`);
  console.log("    Proceeding to MarkDuplicates with research advisory tag.");
  console.log("");
  process.exit(0);
}

if (result.tier === "TIER_3_OPERATOR_OVERRIDDEN") {
  console.log(" ⚠️  VERDICT: TIER 3 — OPERATOR OVERRIDE ACTIVE");
  console.log(`    ${result.supervisorThought}`);
  try {
    fs.renameSync(overrideMarkerPath, overrideMarkerPath + ".applied");
  } catch {}
  process.exit(0);
}

if (result.tier === "TOOL_CRASH") {
  console.log(" 💥 VERDICT: TOOL CRASH / CORRUPTED ALIGNMENT OUTPUT");
  console.log(`    ${result.supervisorThought}`);
  try {
    fs.writeFileSync(haltReportPath, JSON.stringify({
      halted: true,
      phase: "Alignment",
      faultClass: "TOOL_ERROR",
      severity: "critical",
      timestamp: new Date().toISOString(),
      sampleName,
      metrics: result.metrics,
      rootCause: "BWA-MEM or samtools crashed, run was terminated, or outputs were corrupted.",
    }, null, 2));
  } catch {}
  process.exit(2);
}

// Tier 3 Quality Floor Failure -> HALT
console.log(" 🔴 VERDICT: TIER 3 — ALIGNMENT QUALITY REJECTION FLOOR FAILURE");
console.log(`    ${result.supervisorThought}`);
console.log("");
console.log(" ═══════════════════════════════════════════════════════════════");
console.log(" 🛑 HUMAN-IN-THE-LOOP OPERATOR OPINION GATE INVOKED");
console.log(" ═══════════════════════════════════════════════════════════════");
if (result.operatorExplanation) {
  console.log(` [Root Cause Diagnosis]  ${result.operatorExplanation.rootCause}`);
  console.log(` [Sequencer Physics]     ${result.operatorExplanation.sequencerPhysics}`);
  console.log(" [Projected Downstream Risks]:");
  result.operatorExplanation.downstreamRisks.forEach((r, idx) => console.log(`   ${idx + 1}. ${r}`));
}

try {
  fs.writeFileSync(haltReportPath, JSON.stringify({
    halted: true,
    phase: "Alignment",
    faultClass: "ALIGNMENT_QUALITY_FAILURE",
    severity: "critical",
    timestamp: new Date().toISOString(),
    sampleName,
    metrics: result.metrics,
    operatorExplanation: result.operatorExplanation,
    supervisorThought: result.supervisorThought,
  }, null, 2));
} catch {}

process.exit(1);

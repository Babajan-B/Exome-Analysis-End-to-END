#!/usr/bin/env node
/**
 * Inline Stage 1 QC Supervisor Gate & Reasoning Engine (Layer 3 Conductor)
 *
 * Invoked by ULTIMATE_MASTER_PIPELINE.sh directly after fastp trimming / profiling.
 * Implements the 3-Tier Quality Triage & Operator Opinion Gate:
 *   1. Tier 1: Pristine Quality (Phred >= 30, Q30 >= 85%, Adapters < 0.5%)
 *      -> Direct to Alignment. Bypasses quality trimming, preserves full 150bp read length. Exits 0.
 *   2. Tier 2: Moderate Quality (Phred 20 - 30)
 *      -> Branch A (Pre-trim): Autonomous trimming loop & re-check. Prescribes fastp flags. Exits 42.
 *      -> Branch B (Post-trim): Considered & Qualified as Research Grade (FLAG_AND_CONTINUE). Exits 0.
 *   3. Tier 3: Low Quality / Floor Failure (Phred < 20, Q30 < 75%)
 *      -> Halts execution before Alignment. Exits 1.
 *      -> Emits deep diagnostic explanation + Human-in-the-Loop Operator Opinion Gate.
 */
const fs = require("fs");
const path = require("path");

process.on("uncaughtException", (err) => {
  console.error("❌ [SUPERVISOR TOOL CRASH] Uncaught exception in Stage 1 QC gate:", err);
  process.exit(2);
});

const outputDir = process.argv[2];
const sampleName = process.argv[3] || "sample";

if (!outputDir) {
  console.error("Usage: node stage1_qc_gate.js <output_dir> [sample_name]");
  process.exit(2);
}

const fastpJsonPath = path.join(outputDir, "trimmed", "fastp_report.json");
const auditTrailPath = path.join(outputDir, "audit_trail.json");
const haltReportPath = path.join(outputDir, "halt_report.json");
const reasoningPath = path.join(outputDir, "supervisor_reasoning.json");
const retryFile = path.join(outputDir, "trimmed", ".retry_count");
const remediationFlagsFile = path.join(outputDir, "trimmed", ".remediation_flags");

// Load policy standards from qc.json (Single Source of Truth)
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

// Load unified QC Triage Engine (Single Source of Truth for tier logic & thought strings)
let qcTriageEngine = null;
const candidateEnginePaths = [
  path.resolve(__dirname, "../../../../lib/exome/qc-triage-engine.js"),
  path.resolve(__dirname, "../../../lib/exome/qc-triage-engine.js"),
  path.resolve(process.cwd(), "lib/exome/qc-triage-engine.js"),
];
for (const p of candidateEnginePaths) {
  try {
    if (fs.existsSync(p)) {
      qcTriageEngine = require(p);
      break;
    }
  } catch {}
}

if (!qcTriageEngine) {
  console.error("❌ [SUPERVISOR TOOL CRASH] Failed to load qc-triage-engine.js");
  process.exit(2);
}

const rawReadQc = qcPolicy?.tiers?.tier_1_raw_read_qc || {};
const Q30_CLINICAL = rawReadQc.Q30_percentage?.clinical_grade ?? 85.0;
const Q30_RESEARCH = rawReadQc.Q30_percentage?.research_grade ?? 78.0;
const Q30_FLOOR = rawReadQc.Q30_percentage?.rejection_floor ?? 75.0;
const ADAPTER_CLINICAL = rawReadQc.adapter_percentage?.clinical_grade ?? 0.5;
const ADAPTER_MAX = rawReadQc.adapter_percentage?.research_grade ?? 2.0;
const ADAPTER_FLOOR = rawReadQc.adapter_percentage?.rejection_floor ?? 3.0;
const PHRED_CLINICAL = 30.0;
const PHRED_FLOOR = 20.0;

// 1. Load audit trail
let auditTrail = [];
try {
  if (fs.existsSync(auditTrailPath)) {
    auditTrail = JSON.parse(fs.readFileSync(auditTrailPath, "utf8"));
  }
} catch {
  auditTrail = [];
}

function appendAudit(event) {
  event.timestamp = new Date().toISOString();
  auditTrail.push(event);
  try {
    fs.writeFileSync(auditTrailPath, JSON.stringify(auditTrail, null, 2));
  } catch (err) {
    console.error("Failed to write audit trail:", err.message);
  }
}

// 2. Verify fastp output exists
if (!fs.existsSync(fastpJsonPath)) {
  console.error(`\n❌ [SUPERVISOR] Tool Crash / File Missing: fastp_report.json missing in ${outputDir}/trimmed.`);
  const toolCrashReasoning = {
    tier: "TOOL_EXECUTION_FAILURE",
    action: "HALT_TOOL_CRASH",
    exitCode: 2,
    sampleName,
    faultClass: "TOOL_CRASH",
    title: "fastp Execution Aborted / Report Missing",
    rootCause: "fastp did not generate fastp_report.json. The process may have run out of memory (OOM), experienced an unhandled exception, or aborted unexpectedly.",
    evidence: `Expected report file missing at ${fastpJsonPath}`,
    supervisorThought: `[Supervisor Alert: Tool Failure] fastp failed to generate an execution report for ${sampleName}. This is a tool or infrastructure failure, not a biological quality failure. Compute halted.`,
    downstreamDirectives: { targetStage: "none", bwaFlags: [], softClipRisk: "HIGH", inDelArtifactRisk: "HIGH", coverageRisk: "HIGH", watchpoints: ["Check fastp stdout/stderr and system logs for OOM or binary crash."] },
    confidenceScore: 1.0,
    provider: "deterministic",
    generatedAt: new Date().toISOString(),
  };
  try { fs.writeFileSync(reasoningPath, JSON.stringify(toolCrashReasoning, null, 2)); } catch {}
  appendAudit({
    layer: 3,
    agent: "Supervisor Conductor",
    stage: "QC",
    status: "failed",
    message: "Critical: fastp_report.json not produced. Tool crash or system abort.",
  });
  process.exit(2);
}

let report;
try {
  report = JSON.parse(fs.readFileSync(fastpJsonPath, "utf8"));
} catch (err) {
  console.error(`\n❌ [SUPERVISOR] Corrupted fastp JSON:`, err.message);
  const jsonCorruptReasoning = {
    tier: "TOOL_EXECUTION_FAILURE",
    action: "HALT_TOOL_CRASH",
    exitCode: 2,
    sampleName,
    faultClass: "TOOL_CRASH",
    title: "fastp Report JSON Corrupt",
    rootCause: `Corrupted JSON payload in fastp_report.json: ${err.message}`,
    evidence: `File exists but could not be parsed as valid JSON.`,
    supervisorThought: `[Supervisor Alert: Tool Failure] fastp_report.json for ${sampleName} is truncated or malformed. Tool execution did not complete cleanly.`,
    downstreamDirectives: { targetStage: "none", bwaFlags: [], softClipRisk: "HIGH", inDelArtifactRisk: "HIGH", coverageRisk: "HIGH", watchpoints: ["Re-run fastp after inspecting disk IOPS and process logs."] },
    confidenceScore: 1.0,
    provider: "deterministic",
    generatedAt: new Date().toISOString(),
  };
  try { fs.writeFileSync(reasoningPath, JSON.stringify(jsonCorruptReasoning, null, 2)); } catch {}
  process.exit(2);
}

// 3. Extract metrics
const after = report.summary?.after_filtering || {};
const before = report.summary?.before_filtering || {};
const adapter = report.adapter_cutting || {};

const q30 = (after.q30_rate != null ? after.q30_rate : before.q30_rate || 0) * 100;
const q30Before = (before.q30_rate || 0) * 100;
const totalBefore = before.total_reads || 0;
const totalAfter = after.total_reads || totalBefore;
const retainedPct = totalBefore > 0 ? (totalAfter / totalBefore) * 100 : 100;
const adapterTrimmed = adapter.adapter_trimmed_reads || 0;
const adapterPct = totalBefore > 0 ? (adapterTrimmed / totalBefore) * 100 : 0;

// Extract quality curves
const r1Curve = report.read1_before_filtering?.quality_curves?.mean || report.read1_after_filtering?.quality_curves?.mean;
const r2Curve = report.read2_before_filtering?.quality_curves?.mean || report.read2_after_filtering?.quality_curves?.mean;

// Calculate or infer mean Phred using shared engine
const meanPhred = qcTriageEngine.calculateMeanPhred(r1Curve, r2Curve, q30);

// Read retry count
let retryCount = 0;
try {
  if (fs.existsSync(retryFile)) {
    retryCount = parseInt(fs.readFileSync(retryFile, "utf8").trim(), 10) || 0;
  }
} catch {
  retryCount = 0;
}
const isPostTrim = retryCount > 0;

console.log("\n═══════════════════════════════════════════════════════════════");
console.log(" 🧠 [SUPERVISOR AGENT] COGNITIVE REASONING & QUALITY TRIAGE");
console.log("═══════════════════════════════════════════════════════════════");
console.log(` Sample:       ${sampleName}`);
console.log(` Mean Phred:   ${meanPhred.toFixed(1)}`);
console.log(` Q30 Quality:  ${q30.toFixed(2)}% (Pre-filter: ${q30Before.toFixed(2)}%)`);
console.log(` Adapters:     ${adapterPct.toFixed(2)}% (${adapterTrimmed.toLocaleString()} reads)`);
console.log(` Retained:     ${retainedPct.toFixed(1)}% (${totalAfter.toLocaleString()} / ${totalBefore.toLocaleString()})`);
console.log(` Cycle / Pass: ${isPostTrim ? `Post-Trim (Attempt ${retryCount})` : "Initial Evaluation"}`);
console.log("───────────────────────────────────────────────────────────────");

const now = new Date().toISOString();

// Check for Operator Override
const overrideFile = path.join(outputDir, ".override_qc_gate");
if (fs.existsSync(overrideFile)) {
  console.log("\n ⚠️  [SUPERVISOR] Human-in-the-Loop Operator Override Active.");
  console.log("    User choice: [Override & Force Run Anyway]");
  console.log("    Authorizing handoff to Alignment Agent under HIGH-RISK RESEARCH protocol.\n");

  const overrideReasoning = {
    tier: "TIER_3_OPERATOR_OVERRIDDEN",
    action: "FLAG_AND_CONTINUE",
    exitCode: 0,
    sampleName,
    meanPhred,
    q30Pct: q30,
    q30BeforePct: q30Before,
    adapterPct,
    retainedPct,
    retryCount,
    supervisorThought:
      `[Supervisor Operator Override] Sample ${sampleName} failed quality baseline (Mean Phred ${meanPhred.toFixed(1)}, Q30 ${q30.toFixed(1)}%), ` +
      `but the Human Operator explicitly exercised override authority. Advancing to Alignment with strict warning flags.`,
    downstreamDirectives: {
      targetStage: "align",
      bwaFlags: ["-M", "-Y"],
      softClipRisk: "HIGH",
      inDelArtifactRisk: "HIGH",
      coverageRisk: "HIGH",
      watchpoints: [
        "CRITICAL: Run proceeded via Operator Override. Downstream variant calling will contain high false-positive rate.",
      ],
    },
    confidenceScore: 0.5,
    provider: "deterministic",
    generatedAt: now,
  };
  fs.writeFileSync(reasoningPath, JSON.stringify(overrideReasoning, null, 2));

  appendAudit({
    layer: 3,
    agent: "Supervisor Conductor",
    stage: "QC",
    status: "warning",
    tier: "TIER_3_OPERATOR_OVERRIDDEN",
    message: `Supervisor Handoff Authorized via Operator Override. Sample quality is sub-threshold (Phred ${meanPhred.toFixed(1)}).`,
  });

  // Consume the override marker file so it does not permanently poison future fresh runs
  try {
    const appliedFile = path.join(outputDir, ".override_qc_gate.applied");
    fs.renameSync(overrideFile, appliedFile);
  } catch {}
  // Clean up retry flags
  try { if (fs.existsSync(retryFile)) fs.unlinkSync(retryFile); } catch {}
  try { if (fs.existsSync(remediationFlagsFile)) fs.unlinkSync(remediationFlagsFile); } catch {}

  process.exit(0);
}

// ─────────────────────────────────────────────────────────────
// Execute Unified Cognitive QC Triage Engine
// ─────────────────────────────────────────────────────────────
const triageResult = qcTriageEngine.evaluateQcTriage({
  sampleName,
  meanPhred,
  q30,
  q30Before,
  adapterPct,
  retainedPct,
  retryCount,
  isPostTrim,
  qcPolicy,
  now,
});

fs.writeFileSync(reasoningPath, JSON.stringify(triageResult, null, 2));

// ─────────────────────────────────────────────────────────────
// TIER 1: Pristine High Quality (Phred >= 30, Q30 >= 85%, Adapters < 0.5%)
// ─────────────────────────────────────────────────────────────
if (triageResult.tier === "TIER_1_DIRECT_ALIGNMENT") {
  console.log(" 🟢 VERDICT: TIER 1 — DIRECT TO ALIGNMENT (Phred >= 30)");
  console.log("    " + triageResult.supervisorThought);
  console.log("    Direct handoff to Alignment Agent authorized.\n");

  appendAudit({
    layer: 3,
    agent: "Supervisor Conductor",
    stage: "QC",
    status: "passed",
    tier: "TIER_1_DIRECT_ALIGNMENT",
    metrics: { meanPhred, q30, adapterPct, retainedPct },
    message: `Supervisor Approved (Tier 1 Direct Alignment): Mean Phred ${meanPhred.toFixed(1)} >= ${PHRED_CLINICAL}, Q30 ${q30.toFixed(1)}% >= ${Q30_CLINICAL}%. Trimming bypassed to preserve read length.`,
  });

  // Clean up retry flags if any
  try { if (fs.existsSync(retryFile)) fs.unlinkSync(retryFile); } catch {}
  try { if (fs.existsSync(remediationFlagsFile)) fs.unlinkSync(remediationFlagsFile); } catch {}

  process.exit(0);
}

// ─────────────────────────────────────────────────────────────
// TIER 2: Moderate Quality (Phred 20 - 30) — Branch A: Autonomous Trimming Loop
// ─────────────────────────────────────────────────────────────
if (triageResult.tier === "TIER_2_TRIM_AND_RECHECK") {
  const trimFlags = (triageResult.prescribedTrimFlags || []).join(" ");
  fs.writeFileSync(retryFile, String(triageResult.retryCount));
  fs.writeFileSync(remediationFlagsFile, trimFlags);

  console.log(" 🟡 VERDICT: TIER 2 — AUTONOMOUS TRIMMING & RE-CHECK (Phred 20 - 30)");
  console.log("    " + triageResult.supervisorThought);
  console.log(`    Prescribed Trimming Flags: ${trimFlags}`);
  console.log("    Signaling pipeline worker to execute trimming and loop back (Code 42)...\n");

  appendAudit({
    layer: 3,
    agent: "Supervisor Conductor",
    stage: "QC",
    status: "info",
    action: "AUTONOMOUS_TRIMMING_TRIGGERED",
    tier: "TIER_2_TRIM_AND_RECHECK",
    prescribedParams: trimFlags,
    message: `Supervisor Triggered Trimming (Tier 2): Mean Phred ${meanPhred.toFixed(1)} in [${PHRED_FLOOR}, ${PHRED_CLINICAL}]. Prescribing sliding window & adapter removal.`,
  });

  process.exit(42);
}

// ─────────────────────────────────────────────────────────────
// TIER 2: Moderate Quality (Phred 20 - 30) — Branch B: Research Grade Qualified
// ─────────────────────────────────────────────────────────────
if (triageResult.tier === "TIER_2_POST_TRIM_RESEARCH_QUALIFIED") {
  console.log(" 🟢 VERDICT: TIER 2 — POST-TRIM CONSIDERED & APPROVED (Research Grade)");
  console.log("    " + triageResult.supervisorThought);
  console.log("    Authorizing handoff to Alignment Agent with downstream advisory.\n");

  appendAudit({
    layer: 3,
    agent: "Supervisor Conductor",
    stage: "QC",
    status: "warning",
    tier: "TIER_2_POST_TRIM_RESEARCH_QUALIFIED",
    grade: "RESEARCH_GRADE_QUALIFIED",
    behavior: "FLAG_AND_CONTINUE",
    metrics: { meanPhred, q30, adapterPct, retainedPct },
    message: `Supervisor Approved (Tier 2 Research Grade): Mean Phred ${meanPhred.toFixed(1)} in [${PHRED_FLOOR}, ${PHRED_CLINICAL}]. Considered & advancing with downstream watchpoint.`,
  });

  // Clean up retry flags
  try { if (fs.existsSync(retryFile)) fs.unlinkSync(retryFile); } catch {}
  try { if (fs.existsSync(remediationFlagsFile)) fs.unlinkSync(remediationFlagsFile); } catch {}

  process.exit(0);
}

// ─────────────────────────────────────────────────────────────
// TIER 3: Low Quality / Floor Failure (Phred < 20, Q30 < 75%)
// ─────────────────────────────────────────────────────────────
const haltReport = {
  halted: true,
  stage: "Quality Control",
  faultClass: "POOR_LIBRARY_QUALITY",
  title: "Sub-Threshold Quality Below Rejection Floor (Phred < 20)",
  rootCause: triageResult.operatorExplanation?.rootCause,
  evidence: triageResult.operatorExplanation?.evidence,
  severity: "critical",
  operatorExplanation: triageResult.operatorExplanation,
  remediation: {
    kind: "user_opinion",
    label: "Operator Opinion Required",
    detail: triageResult.operatorExplanation?.summary,
  },
  generatedAt: now,
};

fs.writeFileSync(haltReportPath, JSON.stringify(haltReport, null, 2));

console.log("\n 🔴 VERDICT: TIER 3 — CRITICAL QUALITY HALT (Phred < 20)");
console.log("    " + triageResult.supervisorThought);
console.log("\n 📋 OPERATOR EXPLANATION & RISK AUDIT:");
console.log(`    • Root Cause:        ${triageResult.operatorExplanation?.rootCause}`);
console.log(`    • Sequencer Physics: ${triageResult.operatorExplanation?.sequencerPhysics}`);
console.log("    • Downstream Risks:");
(triageResult.operatorExplanation?.downstreamRisks || []).forEach((r) => console.log(`       - ${r}`));
console.log("\n ⏸️  Awaiting Operator Opinion: [Abort Analysis] vs. [Override & Force Run]\n");

appendAudit({
  layer: 4,
  agent: "Reviewer Agent",
  stage: "QC",
  status: "failed",
  tier: "TIER_3_HALT_USER_OPINION",
  faultClass: "POOR_LIBRARY_QUALITY",
  metrics: { meanPhred, q30, adapterPct, retainedPct },
  message: `Supervisor Halted (Tier 3 Critical Quality): Mean Phred ${meanPhred.toFixed(1)} < 20. Operator opinion required.`,
});

// Clean up retry flags on halt so future fresh runs start at cycle 0
try { if (fs.existsSync(retryFile)) fs.unlinkSync(retryFile); } catch {}
try { if (fs.existsSync(remediationFlagsFile)) fs.unlinkSync(remediationFlagsFile); } catch {}

process.exit(1);

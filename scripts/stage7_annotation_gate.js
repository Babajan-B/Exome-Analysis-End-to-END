#!/usr/bin/env node
/**
 * Stage 7 Variant Functional Annotation & Integrity Supervisor Gate
 *
 * Implements:
 * 1. Fast streaming VCF metric parsing (ANN consequence completeness, rsID matching, population AF)
 * 2. 3-Tier Quality Triage against Single Source of Truth (`qc.json`)
 * 3. Multi-Engine Registry detection (snpEff, ANNOVAR, VEP)
 * 4. Human-in-the-Loop Operator Opinion Gate on rejection floor failure (.override_annotation_gate)
 * 5. Downstream directives and candidate locus handoff for Stage 8 ACMG Interpretation
 */

const fs = require("fs");
const path = require("path");

const outputDir = process.argv[2];
const sampleName = process.argv[3] || (outputDir ? path.basename(outputDir) : "Sample");

if (!outputDir) {
  console.error("Usage: node stage7_annotation_gate.js <output_dir> [sample_name]");
  process.exit(1);
}

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

// 2. Load unified Annotation Triage Engine
let annotationTriageEngine = null;
const candidateEnginePaths = [
  path.resolve(__dirname, "../../../../lib/exome/annotation-triage-engine.js"),
  path.resolve(__dirname, "../../../lib/exome/annotation-triage-engine.js"),
  path.resolve(process.cwd(), "lib/exome/annotation-triage-engine.js"),
];
for (const p of candidateEnginePaths) {
  try {
    if (fs.existsSync(p)) {
      annotationTriageEngine = require(p);
      break;
    }
  } catch {}
}

if (!annotationTriageEngine) {
  console.error("❌ [SUPERVISOR TOOL CRASH] Failed to load annotation-triage-engine.js");
  process.exit(2);
}

// 3. Define artifact paths
const annovarDir = path.join(outputDir, "annovar");
const snpeffDir = path.join(annovarDir, "snpeff");
const qcDir = path.join(outputDir, "qc");
if (!fs.existsSync(qcDir)) {
  try { fs.mkdirSync(qcDir, { recursive: true }); } catch {}
}

const overrideMarkerPath = path.join(outputDir, ".override_annotation_gate");
const appliedMarkerPath = path.join(outputDir, ".override_annotation_gate.applied");
const reasoningPath = path.join(qcDir, "stage7_annotation_reasoning.json");
const rootReasoningPath = path.join(outputDir, "stage7_annotation_reasoning.json");
const auditTrailPath = path.join(outputDir, "audit_trail.json");
const haltReportPath = path.join(outputDir, "halt_report.json");

// 4. Locate annotated VCF
const candidateVcfs = [
  path.join(snpeffDir, `${sampleName}_snpEff_annotated.vcf`),
  path.join(snpeffDir, `${sampleName}_snpEff_annotated.vcf.gz`),
  path.join(outputDir, "annotated", `${sampleName}_annotated.vcf`),
  path.join(annovarDir, `annotated_${sampleName}.vcf`),
];

let activeVcf = "";
for (const c of candidateVcfs) {
  if (fs.existsSync(c)) {
    activeVcf = c;
    break;
  }
}

// Fallback search in snpeff directory
if (!activeVcf && fs.existsSync(snpeffDir)) {
  const files = fs.readdirSync(snpeffDir);
  const found = files.find(f => f.endsWith(".vcf") || f.endsWith(".vcf.gz"));
  if (found) activeVcf = path.join(snpeffDir, found);
}

const hasOverride = fs.existsSync(overrideMarkerPath);

if (!activeVcf && !hasOverride) {
  console.error(`❌ [SUPERVISOR TOOL CRASH] Missing annotated VCF outputs in ${snpeffDir} or ${annovarDir}`);
  const haltPayload = {
    halted: true,
    timestamp: new Date().toISOString(),
    phase: "Annotation",
    stage: "stage_8_annotation_integrity_and_clinical_critic",
    sampleName,
    exitCode: 1,
    faultClass: "MISSING_ANNOTATION_OUTPUTS",
    diagnostics: {
      message: `No annotated VCF found for sample ${sampleName}. Checked: ${candidateVcfs.join(", ")}`
    }
  };
  try { fs.writeFileSync(haltReportPath, JSON.stringify(haltPayload, null, 2)); } catch {}
  process.exit(1);
}

// 5. Parse annotated VCF metrics
let metrics = {
  totalRecords: 0,
  annotatedWithGeneCount: 0,
  annotatedWithConsequenceCount: 0,
  consequenceCompletenessPct: 0,
  rsIdCount: 0,
  rsIdMatchRatePct: 0,
  populationAfCount: 0,
  populationAfRatePct: 0,
  impactBreakdown: { HIGH: 0, MODERATE: 0, LOW: 0, MODIFIER: 0, UNKNOWN: 0 },
  consequenceBreakdown: {},
  topGenes: [],
  highImpactCandidateCount: 0,
  highImpactCandidates: []
};

if (activeVcf && fs.existsSync(activeVcf)) {
  try {
    metrics = annotationTriageEngine.parseAnnotatedVcf(activeVcf);
  } catch (err) {
    console.error(`⚠️ [SUPERVISOR] Failed to parse annotated VCF: ${err.message}`);
  }
}

// 6. Detect available annotator backends
const annotatorEngines = annotationTriageEngine.detectAvailableAnnotators();

// 7. Evaluate 3-tier triage
const triage = annotationTriageEngine.evaluateAnnotationTriage({
  sampleName,
  outputDir,
  normLogPath: path.join(qcDir, "bcftools_norm.log"),
  normStatusPath: path.join(qcDir, "norm_status.json"),
  pass1SummaryPath: path.join(qcDir, "pass1_filter_summary.json"),
  haltReportPath,
  metrics,
  qcPolicy
});

// Override handling for Stage 7's own floor failure
if (hasOverride && triage.exitCode === 1 && !triage.upstreamHalt) {
  triage.tier = "TIER_2_RESEARCH_QUALIFIED";
  triage.exitCode = 0;
  triage.thought = `[SUPERVISOR OVERRIDE] Stage 7 Annotation Rejection Floor overridden by operator. Downstream analysis flagged as Research Grade.`;
}

// 8. Write reasoning artifacts
const reasoningPayload = {
  timestamp: new Date().toISOString(),
  sampleName,
  stage: "stage_8_annotation_integrity_and_clinical_critic",
  stageId: 7,
  stageName: "Variant Functional Annotation",
  tier: triage.tier,
  exitCode: triage.exitCode,
  activeVcf: path.relative(outputDir, activeVcf || ""),
  metricsSummary: triage.metricsSummary,
  classification: triage.classification,
  annotatorEngines,
  directivesForStage8: triage.directivesForStage8,
  operatorExplanation: triage.operatorExplanation,
  thought: triage.thought
};

try {
  fs.writeFileSync(reasoningPath, JSON.stringify(reasoningPayload, null, 2));
  fs.writeFileSync(rootReasoningPath, JSON.stringify(reasoningPayload, null, 2));
} catch (err) {
  console.error("⚠️ [SUPERVISOR] Warning: Failed to write stage7_annotation_reasoning.json:", err.message);
}

// 9. Append to audit_trail.json
try {
  let auditTrail = [];
  if (fs.existsSync(auditTrailPath)) {
    try { auditTrail = JSON.parse(fs.readFileSync(auditTrailPath, "utf8")); } catch {}
  }
  auditTrail.push({
    timestamp: new Date().toISOString(),
    layer: 3,
    agent: "Cognitive Supervisor (Stage 7 Variant Functional Annotation)",
    tier: triage.tier,
    exitCode: triage.exitCode,
    totalAnnotatedVariants: metrics.totalRecords,
    consequenceCompleteness: `${metrics.consequenceCompletenessPct}%`,
    rsIdMatchRate: `${metrics.rsIdMatchRatePct}%`,
    highImpactVariants: metrics.impactBreakdown?.HIGH || 0,
    message: triage.thought
  });
  fs.writeFileSync(auditTrailPath, JSON.stringify(auditTrail, null, 2));
} catch (err) {
  console.error("⚠️ [SUPERVISOR] Warning: Failed to append to audit_trail.json:", err.message);
}

// 10. Exit
if (triage.exitCode === 0) {
  if (hasOverride) {
    try { fs.renameSync(overrideMarkerPath, appliedMarkerPath); } catch {}
    console.log(`✅ [SUPERVISOR OVERRIDE] Stage 7 Annotation floor failure overridden by operator.`);
  } else if (triage.tier === "TIER_2_RESEARCH_QUALIFIED") {
    console.log(`ℹ️  [SUPERVISOR RESEARCH] ${triage.thought}`);
  } else {
    console.log(`✅ [SUPERVISOR PASS] ${triage.thought}`);
  }
  process.exit(0);
} else {
  console.log(`🛑 [SUPERVISOR HALT] ${triage.thought}`);
  process.exit(1);
}

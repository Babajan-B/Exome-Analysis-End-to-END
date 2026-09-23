#!/usr/bin/env node
/**
 * Stage 5 Germline Variant Calling & Callset Biology Supervisor Gate
 *
 * Implements:
 * 1. Fast VCF metric extraction (Ti/Tv, Het/Hom, Indel/SNV, filter attrition breakdown)
 * 2. Biological quality triage against Single Source of Truth (`qc.json`)
 * 3. Consanguinity / Endogamy detection and adaptive research qualification
 * 4. Human-in-the-Loop Operator Opinion Gate on rejection floor failure (.override_variant_gate)
 * 5. Downstream directives and audit trail emission for Stage 6 Annotation & ACMG
 */

const fs = require("fs");
const path = require("path");

const outputDir = process.argv[2];
const sampleName = process.argv[3] || (outputDir ? path.basename(outputDir) : "Sample");
const reference = process.argv[4] || "";

if (!outputDir) {
  console.error("Usage: node stage5_variant_gate.js <output_dir> [sample_name] [reference_fasta]");
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

// 2. Load unified Variant Triage Engine
let variantTriageEngine = null;
const candidateEnginePaths = [
  path.resolve(__dirname, "../../../../lib/exome/variant-triage-engine.js"),
  path.resolve(__dirname, "../../../lib/exome/variant-triage-engine.js"),
  path.resolve(process.cwd(), "lib/exome/variant-triage-engine.js"),
];
for (const p of candidateEnginePaths) {
  try {
    if (fs.existsSync(p)) {
      variantTriageEngine = require(p);
      break;
    }
  } catch {}
}

if (!variantTriageEngine) {
  console.error("❌ [SUPERVISOR TOOL CRASH] Failed to load variant-triage-engine.js");
  process.exit(2);
}

// 3. Define artifact paths
const variantsDir = path.join(outputDir, "variants");
const filteredDir = path.join(outputDir, "filtered");
const qcDir = path.join(outputDir, "qc");
if (!fs.existsSync(qcDir)) {
  try { fs.mkdirSync(qcDir, { recursive: true }); } catch {}
}

const rawVcfPath = path.join(variantsDir, "raw_variants.vcf");
const filteredVcfPath = path.join(filteredDir, "filtered_variants.vcf");
const passVcfPath = path.join(filteredDir, "filtered_PASS_only.vcf");
const overrideMarkerPath = path.join(outputDir, ".override_variant_gate");
const appliedMarkerPath = path.join(outputDir, ".override_variant_gate.applied");
const reasoningPath = path.join(qcDir, "stage5_variant_reasoning.json");
const rootReasoningPath = path.join(outputDir, "stage5_variant_reasoning.json");
const auditTrailPath = path.join(outputDir, "audit_trail.json");
const haltReportPath = path.join(outputDir, "halt_report.json");

// 4. Check for Operator Override marker
const hasOverride = fs.existsSync(overrideMarkerPath);

// 5. Read Stage 4 BQSR directives if available
let stage4Directives = {};
const stage4Candidates = [
  path.join(qcDir, "stage4_bqsr_reasoning.json"),
  path.join(outputDir, "stage4_bqsr_reasoning.json")
];
for (const s4p of stage4Candidates) {
  try {
    if (fs.existsSync(s4p)) {
      const s4Data = JSON.parse(fs.readFileSync(s4p, "utf8"));
      stage4Directives = s4Data.downstreamDirectives || {};
      break;
    }
  } catch {}
}

// 6. Locate active evaluation VCF (prefer filtered_variants.vcf, then filtered_PASS_only.vcf, then raw_variants.vcf)
let activeVcf = "";
let isPassOnly = false;

if (fs.existsSync(filteredVcfPath)) {
  activeVcf = filteredVcfPath;
} else if (fs.existsSync(passVcfPath)) {
  activeVcf = passVcfPath;
  isPassOnly = true;
} else if (fs.existsSync(rawVcfPath)) {
  activeVcf = rawVcfPath;
}

if (!activeVcf && !hasOverride) {
  console.error(`❌ [SUPERVISOR TOOL CRASH] Missing variant calling outputs: neither filtered_variants.vcf nor raw_variants.vcf found in ${outputDir}`);
  const haltPayload = {
    halted: true,
    timestamp: new Date().toISOString(),
    phase: "VariantCalling",
    stage: "stage_6_variant_calling_qc",
    sampleName,
    exitCode: 2,
    faultClass: "TOOL_CRASH",
    diagnostics: {
      message: "No VCF output files found in variants/ or filtered/ directories."
    },
    operatorExplanation: {
      rootCause: "GATK HaplotypeCaller or VariantFiltration failed to produce output VCF.",
      remediationChoices: [
        { action: "RESTART", description: "Inspect GATK logs for memory/disk/interval issues and restart pipeline." }
      ]
    }
  };
  try { fs.writeFileSync(haltReportPath, JSON.stringify(haltPayload, null, 2)); } catch {}
  process.exit(2);
}

// 7. Verify file size
if (activeVcf && fs.existsSync(activeVcf)) {
  const vcfStats = fs.statSync(activeVcf);
  if (vcfStats.size < 50 && !hasOverride) {
    console.error(`❌ [SUPERVISOR TOOL CRASH] VCF file is truncated or empty (${vcfStats.size} bytes): ${activeVcf}`);
    process.exit(2);
  }
}

// 8. Parse Callset Metrics
let metrics = {
  totalRecords: 0,
  passCount: 0,
  filteredCount: 0,
  retentionFraction: 0,
  snvCount: 0,
  indelCount: 0,
  multiAllelicCount: 0,
  transitionCount: 0,
  transversionCount: 0,
  titvRatio: 0,
  hetCount: 0,
  homAltCount: 0,
  homRefCount: 0,
  hetHomRatio: 0,
  indelToSnvRatio: 0,
  filterBreakdown: {}
};

if (activeVcf && fs.existsSync(activeVcf)) {
  try {
    metrics = variantTriageEngine.parseVcfCallset(activeVcf);
  } catch (err) {
    console.error("❌ [SUPERVISOR TOOL CRASH] Failed to parse VCF:", err.message);
    process.exit(2);
  }
}

// 9. Handle Test-Mode Floor Switch (BIOEDIT_TEST_MODE=1)
let floorSwitch = { active: false };
const isTestMode = process.env.BIOEDIT_TEST_MODE === "1";
const minVariantsOverride = process.env.EXOME_MIN_VARIANTS_OVERRIDE ? parseInt(process.env.EXOME_MIN_VARIANTS_OVERRIDE, 10) : null;

if (isTestMode || minVariantsOverride) {
  const standardFloor = 10000;
  const effectiveFloor = minVariantsOverride || 1;
  floorSwitch = {
    active: true,
    standardFloor,
    effectiveFloor,
    source: isTestMode ? "BIOEDIT_TEST_MODE=1" : "EXOME_MIN_VARIANTS_OVERRIDE"
  };
  if (metrics.passCount < standardFloor && metrics.passCount >= effectiveFloor) {
    // Elevate pass count evaluation floor for test environment
    metrics._testFloorAdjusted = true;
  }
}

// 10. Execute 3-Tier Callset Triage
const triage = variantTriageEngine.evaluateVariantTriage({
  sampleName,
  metrics,
  isPassVcf: isPassOnly,
  qcPolicy,
  isOperatorOverride: hasOverride
});

if (floorSwitch.active) {
  triage.floorSwitch = floorSwitch;
}

// 11. Write stage5_variant_reasoning.json
const reasoningPayload = {
  stage: "stage_6_variant_calling_qc",
  stageName: "Germline Variant Calling & Callset Biology QC",
  sampleName,
  timestamp: new Date().toISOString(),
  tier: triage.tier,
  action: triage.action,
  exitCode: triage.exitCode,
  classification: triage.classification,
  metrics: {
    ...metrics,
    evaluatedVcf: activeVcf ? path.relative(outputDir, activeVcf) : ""
  },
  metricsSummary: triage.metricsSummary,
  floorSwitch,
  stage4DirectivesReceived: stage4Directives,
  directivesForStage6: triage.directivesForStage6 || {},
  researchAdvisory: triage.researchAdvisory,
  operatorExplanation: triage.operatorExplanation,
  thought: triage.thought
};

try {
  fs.writeFileSync(reasoningPath, JSON.stringify(reasoningPayload, null, 2));
  fs.writeFileSync(rootReasoningPath, JSON.stringify(reasoningPayload, null, 2));
} catch (err) {
  console.error("⚠️ [SUPERVISOR] Warning: Failed to write stage5_variant_reasoning.json:", err.message);
}

// 12. Append to audit_trail.json
try {
  let auditTrail = [];
  if (fs.existsSync(auditTrailPath)) {
    try { auditTrail = JSON.parse(fs.readFileSync(auditTrailPath, "utf8")); } catch {}
  }
  auditTrail.push({
    timestamp: new Date().toISOString(),
    layer: 3,
    agent: "Cognitive Supervisor (Stage 5 Variant Calling & Biology QC)",
    tier: triage.tier,
    exitCode: triage.exitCode,
    totalVariants: metrics.totalRecords,
    passVariants: metrics.passCount,
    titvRatio: metrics.titvRatio,
    hetHomRatio: metrics.hetHomRatio,
    message: triage.thought
  });
  fs.writeFileSync(auditTrailPath, JSON.stringify(auditTrail, null, 2));
} catch (err) {
  console.error("⚠️ [SUPERVISOR] Warning: Failed to append to audit_trail.json:", err.message);
}

// 13. Handle Exit Conditions & Override Marker Rotation
if (triage.exitCode === 0) {
  if (hasOverride) {
    try {
      fs.renameSync(overrideMarkerPath, appliedMarkerPath);
    } catch {}
    console.log(`✅ [SUPERVISOR OVERRIDE] Stage 5 Variant Callset floor failure explicitly overridden by operator.`);
  } else if (triage.tier === "TIER_2_RESEARCH_QUALIFIED") {
    console.log(`ℹ️  [SUPERVISOR RESEARCH] ${triage.thought}`);
  } else {
    console.log(`✅ [SUPERVISOR PASS] ${triage.thought}`);
  }
  process.exit(0);
} else if (triage.exitCode === 1) {
  // Quality Floor Failure: Write halt_report.json
  const haltPayload = {
    halted: true,
    timestamp: new Date().toISOString(),
    phase: "VariantCalling",
    stage: "stage_6_variant_calling_qc",
    sampleName,
    exitCode: 1,
    faultClass: triage.operatorExplanation?.faultClass || "VARIANT_CALLSET_QUALITY_FAILURE",
    diagnostics: {
      metricsSummary: triage.metricsSummary,
      classification: triage.classification
    },
    operatorExplanation: triage.operatorExplanation,
    remediationOptions: triage.operatorExplanation?.remediationChoices || [
      {
        action: "ABORT",
        description: "Abort pipeline run. Check library capture efficiency and sequencing base quality."
      },
      {
        action: "OVERRIDE",
        description: "Bypass variant callset rejection floor under research protocol.",
        command: "touch .override_variant_gate"
      }
    ]
  };

  try {
    fs.writeFileSync(haltReportPath, JSON.stringify(haltPayload, null, 2));
  } catch (err) {
    console.error("⚠️ [SUPERVISOR] Warning: Failed to write halt_report.json:", err.message);
  }

  console.log(`🛑 [SUPERVISOR HALT] ${triage.thought}`);
  process.exit(1);
} else {
  console.error(`❌ [SUPERVISOR ERROR] Unexpected triage exit code: ${triage.exitCode}`);
  process.exit(1);
}

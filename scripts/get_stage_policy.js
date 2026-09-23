#!/usr/bin/env node
/**
 * CLI Helper: Extract Stage Policy or Metric Threshold from SSoT
 *
 * Usage:
 *   node get_stage_policy.js 5
 *   node get_stage_policy.js bqsr empirical_quality_floor clinical_limit
 *   node get_stage_policy.js 1 q30_percentage rejection_floor
 */

const path = require("path");
let policyLoader = null;

const candidatePaths = [
  path.resolve(__dirname, "../../../lib/exome/policy-loader.js"),
  path.resolve(__dirname, "../../../../lib/exome/policy-loader.js"),
  path.resolve(process.cwd(), "lib/exome/policy-loader.js")
];

for (const p of candidatePaths) {
  try {
    policyLoader = require(p);
    break;
  } catch {}
}

if (!policyLoader) {
  console.error("❌ Failed to load policy-loader.js");
  process.exit(1);
}

const stageArg = process.argv[2];
const metricArg = process.argv[3];
const fieldArg = process.argv[4];

if (!stageArg) {
  console.error("Usage: node get_stage_policy.js <stage_num|name> [metric_name] [clinical_limit|research_limit|rejection_floor]");
  process.exit(1);
}

try {
  const slice = policyLoader.getStagePolicy(stageArg);
  if (!metricArg) {
    console.log(JSON.stringify(slice, null, 2));
    process.exit(0);
  }

  const m = slice.metrics[metricArg] || slice.raw_stage[metricArg];
  if (!m) {
    console.error(`Metric '${metricArg}' not found in stage ${stageArg}`);
    process.exit(1);
  }

  if (!fieldArg) {
    console.log(JSON.stringify(m, null, 2));
  } else {
    console.log(m[fieldArg] !== undefined ? m[fieldArg] : "");
  }
} catch (err) {
  console.error(`❌ ${err.message}`);
  process.exit(1);
}

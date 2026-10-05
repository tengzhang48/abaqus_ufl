import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";
import { validateReport } from "./load-livebench-parser.mjs";

const manifest = JSON.parse(await readFile(new URL("../../tools/livebench_cases.json", import.meta.url), "utf8"));

function validReport() {
  return {
    schema_version: 1, started_at: "2026-10-05T12:00:00+00:00", finished_at: "2026-10-05T12:01:00+00:00",
    duration_seconds: 60, status: "passed", source: { commit: "a".repeat(40), dirty: false, repository: "https://github.com/tengzhang48/abaqus_ufl", ref: null },
    run_url: null, environment: { python: "3.12", platform: "Linux", gfortran: "GNU Fortran", packages: { numpy: "2.4" } },
    cases: Object.entries(manifest).map(([id, item]) => ({
      id, title: item.title, target: item.target, status: "passed",
      checks: item.scripts.map((script) => ({ script, label: "Check", status: "passed", exit_code: 0, duration_seconds: 1, metrics: {},
        log: `logs/${id}-${script.replace(/\.py$/, "")}.log`, source_path: `${item.directory}/${script}` })),
    })),
  };
}

test("a complete run has six bundles plus the FE mesh checks", () => {
  const result = validateReport(validReport(), manifest);
  assert.equal(result.cases.length, 7);
  assert.equal(result.cases.at(-1).target, "FE");
});

test("missing compiled checks, duplicate cases, and nonfinite metrics cannot show a pass", () => {
  for (const mutate of [
    (report) => { report.cases[0].checks.pop(); },
    (report) => { report.cases.push(report.cases[0]); },
    (report) => { report.cases[0].checks[0].metrics.error = NaN; },
    (report) => { report.cases = []; },
    (report) => { report.cases[0].checks[0].status = "failed"; },
    (report) => { report.cases[0].checks[0].log = "https://example.com/log"; },
    (report) => { report.cases[0].checks[0].exit_code = 1; },
  ]) {
    const report = validReport();
    mutate(report);
    assert.throws(() => validateReport(report, manifest));
  }
});

test("a real failed check stays failed, while selected subsets remain explicit", () => {
  const report = validReport();
  report.cases[0].checks[0].status = "timed_out";
  report.cases[0].checks[0].exit_code = -9;
  report.cases[0].status = "failed";
  report.status = "failed";
  assert.equal(validateReport(report, manifest).status, "failed");
  const subset = validReport();
  subset.cases = subset.cases.slice(0, 1);
  assert.equal(validateReport(subset, manifest).cases.length, 1);
});

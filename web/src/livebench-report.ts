export type BenchmarkCheck = {
  script: string;
  label: string;
  source_path: string;
  status: "passed" | "failed" | "timed_out";
  exit_code: number | null;
  duration_seconds: number;
  metrics: Record<string, number>;
  log: string;
};

export type BenchmarkCase = {
  id: string;
  title: string;
  target: "UMAT" | "UEL" | "FE";
  status: "passed" | "failed";
  checks: BenchmarkCheck[];
};

export type BenchmarkReport = {
  schema_version: 1;
  started_at: string;
  finished_at: string;
  duration_seconds: number;
  status: "passed" | "failed";
  source: { commit: string | null; dirty: boolean | null; repository: string; ref: string | null };
  run_url: string | null;
  environment: { python: string; platform: string; gfortran: string | null; packages: Record<string, string | null> };
  cases: BenchmarkCase[];
};

type Manifest = Record<string, { title: string; target: string; directory: string; scripts: string[] }>;

// Treat unavailable, incomplete, and contradictory records explicitly. Never
// infer a successful run from an HTTP response or an empty list of checks.
export function validateReport(value: unknown, manifest: Manifest): BenchmarkReport {
  const report = value as BenchmarkReport;
  const fail = () => { throw new Error("The benchmark report is incomplete or invalid."); };
  const finite = (number: unknown) => typeof number === "number" && Number.isFinite(number) && number >= 0;
  if (!report || report.schema_version !== 1 || !["passed", "failed"].includes(report.status)) fail();
  if (!Number.isFinite(Date.parse(report.started_at)) || !Number.isFinite(Date.parse(report.finished_at))
      || Date.parse(report.finished_at) < Date.parse(report.started_at) || !finite(report.duration_seconds)) fail();
  if (!report.source || report.source.repository !== "https://github.com/tengzhang48/abaqus_ufl"
      || !(report.source.commit === null || /^[a-f0-9]{40,64}$/.test(report.source.commit))
      || !(report.source.dirty === null || typeof report.source.dirty === "boolean")) fail();
  if (!(report.run_url === null || /^https:\/\/github\.com\/[\w.-]+\/[\w.-]+\/actions\/runs\/\d+$/.test(report.run_url))) fail();
  if (!report.environment || typeof report.environment.python !== "string"
      || typeof report.environment.platform !== "string" || !report.environment.packages
      || !(report.environment.gfortran === null || typeof report.environment.gfortran === "string")) fail();
  for (const version of Object.values(report.environment.packages)) {
    if (!(version === null || typeof version === "string")) fail();
  }
  if (!Array.isArray(report.cases) || report.cases.length === 0) fail();
  const seen = new Set<string>();
  for (const item of report.cases) {
    const expected = manifest[item.id];
    if (!expected || seen.has(item.id) || item.title !== expected.title || item.target !== expected.target) fail();
    seen.add(item.id);
    if (!Array.isArray(item.checks) || item.checks.length !== expected.scripts.length) fail();
    for (const [index, check] of item.checks.entries()) {
      if (check.script !== expected.scripts[index] || typeof check.label !== "string"
          || check.source_path !== `${expected.directory}/${check.script}`
          || check.log !== `logs/${item.id}-${check.script.replace(/\.py$/, "")}.log`
          || !["passed", "failed", "timed_out"].includes(check.status)
          || !finite(check.duration_seconds)
          || !(check.exit_code === null || Number.isInteger(check.exit_code))) fail();
      if (!check.metrics || typeof check.metrics !== "object" || Array.isArray(check.metrics)) fail();
      for (const [name, metric] of Object.entries(check.metrics)) {
        if (!/^[a-z][a-z0-9_]*$/.test(name) || typeof metric !== "number" || !Number.isFinite(metric)) fail();
      }
      if (check.status === "passed" && check.exit_code !== 0) fail();
    }
    const passed = item.checks.every((check) => check.status === "passed");
    if (item.status !== (passed ? "passed" : "failed")) fail();
  }
  if (report.status !== (report.cases.every((item) => item.status === "passed") ? "passed" : "failed")) fail();
  return report;
}

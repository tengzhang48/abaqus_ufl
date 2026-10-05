export type CurvePlot = {
  kind: "curve"; id: string; title: string; description: string;
  x_label: string; y_label: string; x: number[]; x_scale: "linear" | "log"; y_scale: "linear" | "log";
  series: { label: string; role: "computed" | "reference" | "guide"; values: number[] }[];
  comparison?: { max_abs_error: number; rtol: number; atol: number };
  figure?: string; difference_figure?: string;
};

export type MeshPlot = {
  kind: "mesh"; id: string; title: string; description: string;
  nodes: number[][]; elements: number[][]; values: number[]; displacements?: number[][];
  figure?: string;
};

export type BenchmarkPlot = CurvePlot | MeshPlot;

export type BenchmarkCheck = {
  script: string;
  label: string;
  source_path: string;
  status: "passed" | "failed" | "timed_out";
  exit_code: number | null;
  duration_seconds: number;
  metrics: Record<string, number>;
  log: string;
  plots?: BenchmarkPlot[];
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
      if (check.plots !== undefined) {
        if (!Array.isArray(check.plots)) fail();
        const ids = new Set<string>();
        for (const plot of check.plots) {
          validatePlot(plot);
          if (plot.figure !== undefined && plot.figure !== `figures/${item.id}/${plot.id}.svg`) fail();
          if (plot.kind === "curve" && plot.difference_figure !== undefined
              && plot.difference_figure !== `figures/${item.id}/${plot.id}-difference.svg`) fail();
          if (ids.has(plot.id)) fail();
          ids.add(plot.id);
        }
      }
      if (check.status === "passed" && check.exit_code !== 0) fail();
    }
    const passed = item.checks.every((check) => check.status === "passed");
    if (item.status !== (passed ? "passed" : "failed")) fail();
  }
  if (report.status !== (report.cases.every((item) => item.status === "passed") ? "passed" : "failed")) fail();
  return report;
}

export function validatePlot(value: unknown): BenchmarkPlot {
  const plot = value as BenchmarkPlot;
  const fail = () => { throw new Error("The benchmark plot data is incomplete or invalid."); };
  const array = (v: unknown, length?: number): v is number[] => Array.isArray(v)
    && v.length > 0 && (length === undefined || v.length === length)
    && v.every((n) => typeof n === "number" && Number.isFinite(n));
  if (!plot || !/^[a-z][a-z0-9-]*$/.test(plot.id) || typeof plot.title !== "string"
      || typeof plot.description !== "string") fail();
  if (plot.kind === "curve") {
    if (!array(plot.x) || plot.x.length < 2 || typeof plot.x_label !== "string" || typeof plot.y_label !== "string"
        || !["linear", "log"].includes(plot.x_scale) || !["linear", "log"].includes(plot.y_scale)
        || plot.x.some((x, i) => (i > 0 && x <= plot.x[i - 1]) || (plot.x_scale === "log" && x <= 0))
        || !Array.isArray(plot.series) || plot.series.length === 0) fail();
    for (const series of plot.series) {
      if (typeof series.label !== "string" || !["computed", "reference", "guide"].includes(series.role)
          || !array(series.values, plot.x.length) || (plot.y_scale === "log" && series.values.some((y) => y <= 0))) fail();
    }
    if (!plot.series.some((s) => s.role === "computed")) fail();
    if (plot.comparison !== undefined) {
      const computed = plot.series.filter((s) => s.role === "computed");
      const reference = plot.series.filter((s) => s.role === "reference");
      if (computed.length !== 1 || reference.length !== 1) fail();
      const c = plot.comparison;
      if (![c.max_abs_error, c.rtol, c.atol].every((n) => typeof n === "number" && Number.isFinite(n) && n >= 0)) fail();
      const errors = computed[0].values.map((n, i) => Math.abs(n - reference[0].values[i]));
      if (Math.abs(Math.max(...errors) - c.max_abs_error) > 8 * Number.EPSILON * Math.max(1e-300, c.max_abs_error)
          || errors.some((e, i) => e > c.atol + c.rtol * Math.abs(reference[0].values[i]))) fail();
    }
  } else if (plot.kind === "mesh") {
    if (!Array.isArray(plot.nodes) || !array(plot.values, plot.nodes.length)
        || plot.nodes.some((node) => !array(node, 2)) || !Array.isArray(plot.elements) || !plot.elements.length
        || plot.elements.some((element) => !Array.isArray(element) || element.length !== 4
          || new Set(element).size !== 4 || element.some((i) => !Number.isInteger(i) || i < 0 || i >= plot.nodes.length))) fail();
    if (plot.displacements !== undefined && (!Array.isArray(plot.displacements)
        || plot.displacements.length !== plot.nodes.length || plot.displacements.some((u) => !array(u, 2)))) fail();
  } else fail();
  return plot;
}

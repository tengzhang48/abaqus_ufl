import assert from "node:assert/strict";
import { access, readFile } from "node:fs/promises";
import { validateReport } from "./load-livebench-parser.mjs";

const base = new URL("../public/livebench/", import.meta.url);
const manifest = JSON.parse(await readFile(new URL("../../tools/livebench_cases.json", import.meta.url), "utf8"));
const report = validateReport(JSON.parse(await readFile(new URL("report.json", base), "utf8")), manifest);
assert.deepEqual(report.cases.map((item) => item.id).sort(), Object.keys(manifest).sort(), "publish all seven benchmark cases");
assert.equal(report.source.dirty, false, "a published Actions report must come from a clean checkout");
assert.match(report.source.commit, /^[a-f0-9]{40}$/);
for (const item of report.cases) {
  for (const check of item.checks) await access(new URL(check.log, base));
  if (item.status === "passed") assert.deepEqual(item.checks.flatMap((check) => check.plots ?? []).map((plot) => plot.id), manifest[item.id].plots, `record every numerical plot for ${item.id}`);
  for (const plot of item.checks.flatMap((check) => check.plots ?? [])) {
    if (item.status === "passed") assert.ok(plot.figure, `export a standalone figure for ${item.id}/${plot.id}`);
    if (plot.figure) await access(new URL(plot.figure, base));
    if (plot.kind === "curve" && plot.comparison && item.status === "passed") assert.ok(plot.difference_figure);
    if (plot.kind === "curve" && plot.difference_figure) await access(new URL(plot.difference_figure, base));
  }
}
console.log(`Checked a complete ${report.status} report and its ${report.cases.flatMap((item) => item.checks).length} logs.`);

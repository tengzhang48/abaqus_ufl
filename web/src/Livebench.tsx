import { useEffect, useState } from "react";
import manifest from "../../tools/livebench_cases.json";
import { validateReport } from "./livebench-report";
import type { BenchmarkReport } from "./livebench-report";

const repository = "https://github.com/tengzhang48/abaqus_ufl";
const resultBase = `${import.meta.env.BASE_URL}livebench/`;

export default function Livebench() {
  const [report, setReport] = useState<BenchmarkReport | null>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(true);
  const [refresh, setRefresh] = useState(0);
  const [target, setTarget] = useState("all");
  const [status, setStatus] = useState("all");

  useEffect(() => {
    const abort = new AbortController();
    setLoading(true);
    setError("");
    setReport(null);
    fetch(`${resultBase}report.json`, { cache: "no-store", signal: abort.signal })
      .then((response) => {
        if (!response.ok) throw new Error("A published benchmark report is not available.");
        return response.json();
      })
      .then((value: unknown) => setReport(validateReport(value, manifest)))
      .catch((reason: Error) => { if (!abort.signal.aborted) setError(reason.message); })
      .finally(() => { if (!abort.signal.aborted) setLoading(false); });
    return () => abort.abort();
  }, [refresh]);

  const visible = report?.cases.filter((item) => (target === "all" || item.target === target)
    && (status === "all" || item.status === status)) ?? [];
  const passed = report?.cases.filter((item) => item.status === "passed").length ?? 0;
  const checks = report?.cases.flatMap((item) => item.checks) ?? [];

  return (
    <section className="section livebench-section" id="livebench">
      <div className="shell">
        <div className="section-heading split-heading">
          <div><p className="eyebrow">Runnable verification · no Abaqus license</p><h2>Livebench</h2></div>
          <p>Six material and element bundles, plus mesh solves with the optional serial FE runtime. The generated Fortran runs through f2py; this page shows the recorded results.</p>
        </div>
        <div className="bench-actions">
          <a className="button bench-button" href={`${repository}/blob/main/docs/LIVEBENCH.md`}>Run these benchmarks ↗</a>
          <a href={`${repository}/actions/workflows/pages.yml`}>GitHub Actions ↗</a>
          <button type="button" onClick={() => setRefresh((value) => value + 1)} disabled={loading}>Refresh results</button>
        </div>
        <div aria-live="polite">
          {loading && <p className="bench-unavailable">Loading the latest published run…</p>}
          {error && <p className="bench-unavailable">Results unavailable. {error} Open GitHub Actions for the latest run.</p>}
        </div>
        {report && (
          <>
            <div className="bench-summary">
              <div><strong>{passed}/{report.cases.length}</strong><span>cases passed</span></div>
              <div><strong>{checks.filter((item) => item.status === "passed").length}/{checks.length}</strong><span>checks passed</span></div>
              <div><strong>{report.duration_seconds.toFixed(1)}s</strong><span>run including builds</span></div>
              <div><strong className={`bench-status ${report.status}`}>{report.status === "passed" ? "Passed" : "Failed"}</strong><span>recorded run status</span></div>
            </div>
            <div className="bench-run">
              <p>Finished <time dateTime={report.finished_at}>{new Date(report.finished_at).toLocaleString(undefined, { timeZone: "UTC" })} UTC</time> · source {report.source.commit ? <a href={`${repository}/commit/${report.source.commit}`}><code>{report.source.commit.slice(0, 8)}</code></a> : "unversioned"}{report.source.dirty ? " + local changes" : ""}</p>
              <div>{report.run_url && <a href={report.run_url}>Open this run ↗</a>}<a href={`${resultBase}report.json`} download>Download JSON ↓</a></div>
            </div>
            {report.cases.length !== Object.keys(manifest).length && <p className="bench-unavailable">This report covers a selected subset of the seven benchmarks.</p>}
            <div className="bench-filters">
              <div className="bench-filter"><label htmlFor="bench-type">Type</label><select id="bench-type" value={target} onChange={(event) => setTarget(event.target.value)}><option value="all">All types</option><option>UMAT</option><option>UEL</option><option value="FE">Serial FE</option></select></div>
              <div className="bench-filter"><label htmlFor="bench-result">Result</label><select id="bench-result" value={status} onChange={(event) => setStatus(event.target.value)}><option value="all">All results</option><option value="passed">Passed</option><option value="failed">Failed</option></select></div>
              <span>{visible.length} of {report.cases.length} cases</span>
            </div>
            <div className="bench-cases">
              {visible.map((item) => (
                <details className="bench-case" key={item.id}>
                  <summary><span className="bench-case-title">{item.title}<small>{item.target}</small></span><span className={`bench-status ${item.status}`}>{item.status === "passed" ? "Passed" : "Failed"}</span><span className="bench-time">{item.checks.reduce((sum, check) => sum + check.duration_seconds, 0).toFixed(1)}s</span><span className="bench-expand" aria-hidden="true">＋</span></summary>
                  <div className="bench-detail">
                    {item.checks.map((check) => (
                      <div className="bench-check" key={check.script}>
                        <div className="bench-check-heading"><strong>{check.label}</strong><span className={`bench-status ${check.status}`}>{check.status.replace("_", " ")}</span><a href={`${resultBase}${check.log}`}>Log ↗</a><a href={`${repository}/blob/${report.source.commit || "main"}/${check.source_path}`}>Check source ↗</a></div>
                        {Object.keys(check.metrics).length > 0 && <dl className="bench-metrics">{Object.entries(check.metrics).map(([name, value]) => <div key={name}><dt>{name.replaceAll("_", " ")}</dt><dd>{value.toExponential(4)}</dd></div>)}</dl>}
                      </div>
                    ))}
                  </div>
                </details>
              ))}
              {visible.length === 0 && <p>No cases match these filters.</p>}
            </div>
            <details className="bench-environment"><summary>Run environment</summary><p>Python {report.environment.python} · {report.environment.platform}</p><p>{report.environment.gfortran || "Fortran compiler unavailable"}</p><p>{Object.entries(report.environment.packages).map(([name, version]) => `${name} ${version || "unavailable"}`).join(" · ")}</p></details>
          </>
        )}
        <p className="bench-scope">These material, element, and small-mesh checks complement the paper's evidence records. Paper-scale Abaqus analyses and figure reproductions are documented separately. Timings include compilation and depend on the displayed environment.</p>
        <pre className="bench-command"><code>pip install -e ".[dev]"{"\n"}python tools/run_livebench.py</code></pre>
      </div>
    </section>
  );
}

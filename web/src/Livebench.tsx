import { useEffect, useState } from "react";
import manifest from "../../tools/livebench_cases.json";
import { validateReport } from "./livebench-report";
import type { BenchmarkReport } from "./livebench-report";
import Plot from "./BenchmarkPlot";

const repository = "https://github.com/tengzhang48/abaqus_ufl";
const resultBase = `${import.meta.env.BASE_URL}livebench/`;
const explanations: Record<
  string,
  { title: string; setup: string; reference: string }
> = {
  thermal_bending_bvp: {
    title: "Thermal bending of a clamped strip",
    setup: "Watch a strip bend as heat diffuses between hot and cold baths. One end is clamped; the remaining surfaces are traction-free. Explore 20 solved time steps, temperature, and the free-tip history.",
    reference: "An independent Fourier series checks transient temperature. A slender-beam approximation checks deflection within its stated tolerance; mesh refinement, heat dissipation, and clamp reaction balances provide separate checks.",
  },
  heated_plate_bvp: {
    title: "Transient heating of a 2D plate",
    setup: "A varying temperature bath heats the top of an initially cold plate while the other edges stay cold. Follow the solved 2D field over 20 time steps, inspect stored heat, and compare mesh and time refinement.",
    reference: "Independent continuous-time and backward-Euler Fourier solutions distinguish time error from spatial error. Refinement checks second-order space and first-order time accuracy, with heat balance at every accepted step.",
  },
  serial_fe: {
    title: "Heat diffusion on a mesh",
    setup:
      "Compare a sinusoidal temperature field before and after a diffusion step, inspect the solved mesh, and check spatial refinement.",
    reference:
      "A closed-form Q1 consistent-mass eigenmode checks the discrete solve. Refinement approaches the continuum backward-Euler solution at second order.",
  },
  neo_hookean_umat: {
    title: "Neo-Hookean elasticity",
    setup:
      "Stretch and shear a material point. The shear path also exposes the nonlinear normal-stress response.",
    reference:
      "Hand-derived Cauchy stress formulas, plus an independent finite-difference check of the compiled material tangent.",
  },
  ogden_umat: {
    title: "Ogden elasticity",
    setup:
      "An isochoric stretch sweep exercises the spectral material response with repeated transverse eigenvalues.",
    reference:
      "Principal-stretch powers give the stress directly, without an eigensolver. Separate checks cover rotated states and repeated-spectrum tangents.",
  },
  small_strain_j2_umat: {
    title: "Yield and plastic hardening",
    setup:
      "Follow a shear loading path through elastic yield and into plastic flow. Inspect both stress and accumulated plastic strain.",
    reference:
      "The exact proportional-loading J2 consistency solution with linear isotropic hardening checks every increment and the returned material state.",
  },
  small_strain_viscoelastic_umat: {
    title: "Viscoelastic relaxation",
    setup:
      "Apply a step in shear, hold the strain, and follow the decay toward the equilibrium stress.",
    reference:
      "The exact backward-Euler recurrence checks the compiled history. A separate continuous curve shows the time-step approximation.",
  },
  scalar_diffusion_uel: {
    title: "Thermal expansion and storage",
    setup:
      "A constrained Quad4 couples temperature to mechanical forces. A uniform heating step checks the integrated heat-storage balance.",
    reference:
      "Closed-form shape-gradient integrals check thermal loading, while conservation gives the exact signed storage rate.",
  },
  thermo_mechanics_quad8: {
    title: "Heat flux under deformation",
    setup:
      "Stretch a mixed-order Quad8 while holding a linear referential temperature gradient. The thermal load changes with deformation.",
    reference:
      "The C⁻¹ pull-back gives a thermal-load ratio of 1/λ², checked directly against the compiled element.",
  },
};

export default function Livebench() {
  const [report, setReport] = useState<BenchmarkReport | null>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(true);
  const [refresh, setRefresh] = useState(0);
  const [caseId, setCaseId] = useState("thermal_bending_bvp");
  const [plotId, setPlotId] = useState("");

  useEffect(() => {
    const abort = new AbortController();
    setLoading(true);
    setError("");
    setReport(null);
    fetch(`${resultBase}report.json`, {
      cache: "no-store",
      signal: abort.signal,
    })
      .then((response) => {
        if (!response.ok)
          throw new Error("A published benchmark report is not available.");
        return response.json();
      })
      .then((value: unknown) => setReport(validateReport(value, manifest)))
      .catch((reason: Error) => {
        if (!abort.signal.aborted) setError(reason.message);
      })
      .finally(() => {
        if (!abort.signal.aborted) setLoading(false);
      });
    return () => abort.abort();
  }, [refresh]);

  const selected =
    report?.cases.find((item) => item.id === caseId) ?? report?.cases[0];
  const plots = selected?.checks.flatMap((check) => check.plots ?? []) ?? [];
  const plot = plots.find((item) => item.id === plotId) ?? plots[0];
  const info = selected ? explanations[selected.id] : null;
  const passed =
    report?.cases.filter((item) => item.status === "passed").length ?? 0;
  const checks = report?.cases.flatMap((item) => item.checks) ?? [];
  const ordered = report
    ? [...report.cases].sort((a, b) => {
        const order = ["thermal_bending_bvp", "heated_plate_bvp", "serial_fe"];
        const rank = (id: string) => order.includes(id) ? order.indexOf(id) : order.length;
        return rank(a.id) - rank(b.id);
      })
    : [];

  return (
    <section className="section livebench-section" id="livebench">
      <div className="shell">
        <div className="bench-intro">
          <div>
            <p className="eyebrow">Mechanics · numerical comparisons</p>
            <h2>Livebench</h2>
          </div>
          <p>
            Explore boundary-value simulations with prescribed loads and
            supports, transient fields, and deformation. Compare solved results
            with independent references and inspect the material and element checks.
          </p>
        </div>
        {loading && (
          <p className="bench-unavailable" role="status">
            Loading the latest recorded benchmark run…
          </p>
        )}
        {error && (
          <div className="bench-unavailable" role="alert">
            <strong>Results unavailable.</strong>
            <p>{error}</p>
            <a href={`${repository}/actions/workflows/pages.yml`}>
              Open the latest runs ↗
            </a>
            <button type="button" onClick={() => setRefresh((n) => n + 1)}>
              Try again
            </button>
          </div>
        )}
        {report && selected && info && (
          <>
            <div className={`bench-record ${report.status}`}>
              <span className="run-state">
                <i />
                {report.status === "passed"
                  ? "All recorded cases passed"
                  : `${report.cases.length - passed} recorded case(s) failed`}
              </span>
              <span>
                {report.cases.length} cases · {checks.length} checks
              </span>
              <time dateTime={report.finished_at}>
                {new Date(report.finished_at).toLocaleDateString(undefined, {
                  timeZone: "UTC",
                  year: "numeric",
                  month: "short",
                  day: "numeric",
                })}
              </time>
              <div>
                {report.run_url && <a href={report.run_url}>Recorded run ↗</a>}
                <button type="button" onClick={() => setRefresh((n) => n + 1)}>
                  Refresh
                </button>
              </div>
            </div>
            {report.cases.length !== Object.keys(manifest).length && (
              <p className="bench-unavailable">
                This recorded run covers {report.cases.length} of the {Object.keys(manifest).length} current cases.
              </p>
            )}
            <div className="bench-explorer">
              <aside className="bench-sidebar" aria-label="Benchmark cases">
                <p className="bench-sidebar-label">Choose a problem</p>
                <div className="bench-case-nav">
                  {ordered.map((item) => (
                    <button
                      type="button"
                      key={item.id}
                      aria-pressed={selected.id === item.id}
                      onClick={() => {
                        setCaseId(item.id);
                        setPlotId("");
                      }}
                    >
                      <span className="case-kind">
                        {item.id.endsWith("_bvp") ? "Boundary-value simulation" : item.target === "FE" ? "Mesh check" : item.target}
                      </span>
                      <strong>{explanations[item.id].title}</strong>
                      <span
                        className={`case-dot ${item.status}`}
                        aria-label={item.status}
                      />
                    </button>
                  ))}
                </div>
                <div className="bench-mobile-case">
                  <label htmlFor="bench-case">Problem</label>
                  <select
                    id="bench-case"
                    value={selected.id}
                    onChange={(event) => {
                      setCaseId(event.target.value);
                      setPlotId("");
                    }}
                  >
                    {ordered.map((item) => (
                      <option key={item.id} value={item.id}>
                        {explanations[item.id].title}
                        {item.status === "failed" ? " · failed" : ""}
                      </option>
                    ))}
                  </select>
                </div>
                <p className="bench-sidebar-note">
                  Computed with f2py.
                  <br />
                  No Abaqus installation needed.
                </p>
              </aside>
              <article className="bench-workspace" aria-label={info.title}>
                <div className="bench-problem-heading">
                  <div>
                    <p className="eyebrow">
                      {selected.target === "FE"
                        ? selected.id.endsWith("_bvp") ? "Transient boundary-value simulation" : "Serial finite-element check"
                        : `${selected.target} · generated Fortran`}
                    </p>
                    <h3>{info.title}</h3>
                  </div>
                  <span className={`bench-status ${selected.status}`}>
                    {selected.status === "passed" ? "Passed" : "Failed"}
                  </span>
                </div>
                <p className="bench-problem-summary">{info.setup}</p>
                {plot?.setup && (
                  <div className="simulation-setup">
                    <p>{plot.setup.domain}</p>
                    <dl className="simulation-facts">
                      <div><dt>Mesh</dt><dd>{plot.setup.elements} elements · {plot.setup.nodes} nodes</dd></div>
                      <div><dt>Time</dt><dd>{plot.setup.steps} increments · Δt = {plot.setup.time_step}</dd></div>
                      <div><dt>Solved DOFs</dt><dd>{plot.setup.free_temperature_dofs} temperature · {plot.setup.free_displacement_dofs} displacement</dd></div>
                    </dl>
                    <strong>Boundary conditions</strong>
                    <ul>{plot.setup.boundary_conditions.map((condition) => <li key={condition}>{condition}</li>)}</ul>
                    <details>
                      <summary>Equations, initial state & material</summary>
                      {plot.setup.equations.map((equation) => <p key={equation}>{equation}</p>)}
                      <p>{plot.setup.initial_condition}</p>
                      <p>{plot.setup.properties.map((p) => `${p.name} = ${p.value}`).join(" · ")}</p>
                    </details>
                  </div>
                )}
                {plot ? (
                  <>
                    <div className="bench-view-select">
                      <label htmlFor="bench-view">View</label>
                      <select
                        id="bench-view"
                        value={plot.id}
                        onChange={(event) => setPlotId(event.target.value)}
                      >
                        {plots.map((item) => (
                          <option key={item.id} value={item.id}>
                            {item.title}
                          </option>
                        ))}
                      </select>
                    </div>
                    <Plot
                      key={`${selected.id}-${plot.id}`}
                      plot={plot}
                      caseId={selected.id}
                    />
                    <p className="bench-plot-description">{plot.description}</p>
                  </>
                ) : (
                  <div className="bench-unavailable">
                    {selected.status === "failed"
                      ? "This check failed; a verified plot is unavailable. Read the diagnostics below."
                      : "Response data was not recorded in this older run. The numerical checks are available below."}
                  </div>
                )}
                <div className="bench-reference">
                  <span>Reference & acceptance</span>
                  <p>{info.reference}</p>
                </div>
                <details className="bench-diagnostics">
                  <summary>Numerical checks & source</summary>
                  {selected.checks.map((check) => (
                    <div className="bench-check" key={check.script}>
                      <div className="bench-check-heading">
                        <strong>{check.label}</strong>
                        <span className={`bench-status ${check.status}`}>
                          {check.status.replace("_", " ")}
                        </span>
                        <a href={`${resultBase}${check.log}`}>Log ↗</a>
                        <a
                          href={`${repository}/blob/${report.source.commit || "main"}/${check.source_path}`}
                        >
                          Source ↗
                        </a>
                      </div>
                      {Object.keys(check.metrics).length > 0 && (
                        <dl className="bench-metrics">
                          {Object.entries(check.metrics).map(
                            ([name, value]) => (
                              <div key={name}>
                                <dt>{name.replaceAll("_", " ")}</dt>
                                <dd>{value.toExponential(4)}</dd>
                              </div>
                            ),
                          )}
                        </dl>
                      )}
                    </div>
                  ))}
                </details>
              </article>
            </div>
            <div className="bench-footer">
              <details>
                <summary>Reproduce this case</summary>
                <p>
                  Clone the repository, install a Fortran compiler and the
                  development dependencies, then run this case:
                </p>
                <pre>
                  <code>
                    git clone {repository}.git{"\n"}cd abaqus_ufl{"\n"}pip
                    install -e ".[dev]"{"\n"}python tools/run_livebench.py
                    --case {selected.id}
                  </code>
                </pre>
                <a href={`${repository}/blob/main/docs/LIVEBENCH.md`}>
                  Setup & all benchmarks ↗
                </a>
              </details>
              <details>
                <summary>Run provenance & environment</summary>
                <p>
                  Finished{" "}
                  {new Date(report.finished_at).toLocaleString(undefined, {
                    timeZone: "UTC",
                  })}{" "}
                  UTC · {report.duration_seconds.toFixed(1)}s including builds
                </p>
                <p>
                  Source{" "}
                  {report.source.commit ? (
                    <a href={`${repository}/commit/${report.source.commit}`}>
                      <code>{report.source.commit.slice(0, 8)}</code>
                    </a>
                  ) : (
                    "unversioned"
                  )}
                  {report.source.dirty ? " + local changes" : ""} · Python{" "}
                  {report.environment.python}
                </p>
                <p>
                  {report.environment.gfortran ||
                    "Fortran compiler unavailable"}
                </p>
                <p>
                  {Object.entries(report.environment.packages)
                    .map(
                      ([name, version]) =>
                        `${name} ${version || "unavailable"}`,
                    )
                    .join(" · ")}
                </p>
                <a href={`${resultBase}report.json`} download>
                  Full report & plot data ↓
                </a>
              </details>
            </div>
          </>
        )}
        <p className="bench-scope">
          The plots show recorded f2py computations; browser controls inspect
          those data. These material, element, and small-mesh benchmarks
          complement the paper's separately documented Abaqus analyses.
        </p>
      </div>
    </section>
  );
}

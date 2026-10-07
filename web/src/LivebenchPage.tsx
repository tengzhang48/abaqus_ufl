import { Fragment, useEffect, useState } from "react";
import type { AnchorHTMLAttributes } from "react";
import manifest from "../../tools/livebench_cases.json";
import Plot from "./BenchmarkPlot";
import SiteHeader, { livebenchUrl } from "./SiteHeader";
import BenchmarkPreview from "./BenchmarkPreview";
import { modelCards } from "./livebench-cards";
import { validateReport } from "./livebench-report";
import type { BenchmarkReport, SimulationSetup } from "./livebench-report";
import {
  repository,
  modelCode,
  steps,
  supportCases,
  walkthroughs,
} from "./livebench-content";
import type { Equation, Walkthrough } from "./livebench-content";

const resultBase = `${import.meta.env.BASE_URL}livebench/`;
const caseUrl = (id: string, step?: string) =>
  `${livebenchUrl}?case=${encodeURIComponent(id)}${step ? `&step=${encodeURIComponent(step)}` : ""}`;

function PageLink({
  href,
  onClick,
  ...props
}: AnchorHTMLAttributes<HTMLAnchorElement>) {
  return (
    <a
      {...props}
      href={href}
      onClick={(event) => {
        onClick?.(event);
        if (
          !href ||
          event.defaultPrevented ||
          event.button !== 0 ||
          event.metaKey ||
          event.ctrlKey ||
          event.shiftKey ||
          event.altKey
        )
          return;
        if (
          new URL(href, window.location.href).pathname !==
          window.location.pathname
        )
          return;
        event.preventDefault();
        window.history.pushState(null, "", href);
        window.dispatchEvent(new PopStateEvent("popstate"));
        window.scrollTo({ top: 0 });
      }}
    />
  );
}

function ProblemSketch({ id }: { id: string }) {
  const plate = id === "heated_plate_bvp";
  const shear = id === "ogden_bvp" || id === "plasticity_bvp";
  return (
    <svg
      className="problem-sketch"
      viewBox="0 0 360 200"
      role="img"
      aria-label={
        shear
          ? "Square shear block with fixed bottom, driven top and traction-free sides"
          : plate
            ? "Unit square with a heated top edge and three cold edges"
            : "Undeformed strip with a clamped left end, hot top and cold bottom"
      }
    >
      {shear ? (
        <>
          <rect
            x="110"
            y="38"
            width="140"
            height="130"
            fill="#f2f6f4"
            stroke="#a4bcb2"
          />
          <line
            x1="110"
            y1="168"
            x2="250"
            y2="168"
            stroke="#4b655b"
            strokeWidth="4"
          />
          {[110, 130, 150, 170, 190, 210, 230, 250].map((x) => (
            <line
              key={x}
              x1={x - 8}
              y1="178"
              x2={x}
              y2="168"
              stroke="#4b655b"
            />
          ))}
          <line
            x1="110"
            y1="38"
            x2="250"
            y2="38"
            stroke="#bd593a"
            strokeWidth="4"
          />
          <path
            d="M 150 22 H 222 l -8 -5 m 8 5 l -8 5"
            fill="none"
            stroke="#bd593a"
            strokeWidth="2"
          />
          <text x="180" y="110" textAnchor="middle">
            Traction-free sides
          </text>
          <text x="180" y="198" textAnchor="middle">
            Bottom fixed · top driven
          </text>
        </>
      ) : plate ? (
        <>
          <rect
            x="110"
            y="35"
            width="140"
            height="140"
            fill="#f2f6f4"
            stroke="#407d95"
            strokeWidth="3"
          />
          <line
            x1="110"
            y1="35"
            x2="250"
            y2="35"
            stroke="#bd593a"
            strokeWidth="4"
          />
          <text x="180" y="20" textAnchor="middle">
            Varying temperature bath
          </text>
          <text x="180" y="107" textAnchor="middle">
            Initially T = 0
          </text>
          <text x="180" y="195" textAnchor="middle">
            Three edges held at T = 0
          </text>
        </>
      ) : (
        <>
          <rect
            x="32"
            y="75"
            width="295"
            height="45"
            fill="#f2f6f4"
            stroke="#a4bcb2"
          />
          <line
            x1="32"
            y1="75"
            x2="327"
            y2="75"
            stroke="#bd593a"
            strokeWidth="4"
          />
          <line
            x1="32"
            y1="120"
            x2="327"
            y2="120"
            stroke="#407d95"
            strokeWidth="4"
          />
          <line
            x1="32"
            y1="65"
            x2="32"
            y2="130"
            stroke="#4b655b"
            strokeWidth="4"
          />
          {[65, 78, 91, 104, 117, 130].map((y) => (
            <line key={y} x1="21" y1={y + 8} x2="32" y2={y} stroke="#4b655b" />
          ))}
          <text x="180" y="55" textAnchor="middle">
            Hot bath
          </text>
          <text x="180" y="145" textAnchor="middle">
            Cold bath
          </text>
          <text x="32" y="180">
            Clamped
          </text>
          <text x="250" y="180">
            Free to bend
          </text>
        </>
      )}
    </svg>
  );
}

function Equations({ items }: { items: Equation[] }) {
  return (
    <div className="walkthrough-equations">
      {items.map((item) => (
        <article className="equation-card" key={item.title}>
          <h3>{item.title}</h3>
          <p className="equation-expression">{item.expression}</p>
          <p>{item.explanation}</p>
        </article>
      ))}
    </div>
  );
}

function SetupFacts({ id, setup }: { id: string; setup?: SimulationSetup }) {
  const dofs = [
    [setup?.free_temperature_dofs, "temperature"],
    [setup?.free_displacement_dofs, "displacement"],
  ].filter(([count]) => count);
  return setup ? (
    <dl className="simulation-facts">
      <div>
        <dt>Mesh</dt>
        <dd>
          {setup.elements} elements · {setup.nodes} nodes
        </dd>
      </div>
      <div>
        <dt>{setup.free_temperature_dofs ? "Time" : "Load steps"}</dt>
        <dd>
          {setup.steps} increments · {id === "ogden_bvp" ? "Δs" : "Δt"} ={" "}
          {Number(setup.time_step.toPrecision(4))}
        </dd>
      </div>
      <div>
        <dt>Solved DOFs</dt>
        <dd>
          {dofs.map(([count, field]) => `${count} ${field}`).join(" · ")}
        </dd>
      </div>
    </dl>
  ) : null;
}

function CaseResults({
  id,
  report,
  loading,
  error,
  refresh,
}: {
  id: string;
  report: BenchmarkReport | null;
  loading: boolean;
  error: string;
  refresh: () => void;
}) {
  const [plotId, setPlotId] = useState("");
  const selected = report?.cases.find((item) => item.id === id);
  const plots = selected?.checks.flatMap((check) => check.plots ?? []) ?? [];
  const plot = plots.find((item) => item.id === plotId) ?? plots[0];
  const reference = walkthroughs[id]?.reference ?? supportCases[id]?.reference;
  const sourceLink = (path: string) =>
    `${repository}/blob/${report?.source.commit || "main"}/${path}`;
  return (
    <>
      <h2>
        {walkthroughs[id]
          ? "Run the generated model"
          : "Recorded check results"}
      </h2>
      {walkthroughs[id] && (
        <p>
          The generated{" "}
          {id === "plasticity_bvp"
            ? "UMAT is compiled and called at every integration point"
            : "UEL is compiled and executed for every element evaluation"}
          . A local verification driver applies the boundary conditions and
          solves each accepted increment. The controls below explore the
          recorded solution.
        </p>
      )}
      {loading && (
        <p className="bench-unavailable" role="status">
          Loading the recorded run…
        </p>
      )}
      {error && (
        <div className="bench-unavailable" role="alert">
          <strong>Results unavailable.</strong>
          <p>{error}</p>
          <button type="button" onClick={refresh}>
            Try again
          </button>{" "}
          <a href={`${repository}/actions/workflows/pages.yml`}>
            Open the workflow ↗
          </a>
        </div>
      )}
      {report && !selected && (
        <p className="bench-unavailable">
          This recorded run does not include this case.
        </p>
      )}
      {selected && (
        <>
          <div className="walkthrough-run">
            <span className={`bench-status ${selected.status}`}>
              {selected.status === "passed"
                ? "Recorded checks passed"
                : "Recorded checks failed"}
            </span>
            <time dateTime={report!.finished_at}>
              {new Date(report!.finished_at).toLocaleDateString(undefined, {
                timeZone: "UTC",
                year: "numeric",
                month: "short",
                day: "numeric",
              })}
            </time>
            <button type="button" onClick={refresh}>
              Refresh results
            </button>
          </div>
          <SetupFacts id={id} setup={plot?.setup} />
          {plot ? (
            <>
              <div className="bench-view-select">
                <label htmlFor="bench-view">Result</label>
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
              <Plot key={`${id}-${plot.id}`} plot={plot} caseId={id} />
              <p className="bench-plot-description">{plot.description}</p>
            </>
          ) : (
            <p className="bench-unavailable">
              {selected.status === "failed"
                ? "This run failed. Inspect the diagnostics below."
                : "This older run contains checks without recorded plot data."}
            </p>
          )}
          <div className="bench-reference">
            <span>What checks this result?</span>
            <p>{reference}</p>
          </div>
          <details className="bench-diagnostics">
            <summary>Numerical checks, tolerances & source</summary>
            {selected.checks.map((check) => (
              <div className="bench-check" key={check.script}>
                <div className="bench-check-heading">
                  <strong>{check.label}</strong>
                  <span className={`bench-status ${check.status}`}>
                    {check.status.replace("_", " ")}
                  </span>
                  <a href={`${resultBase}${check.log}`}>Log ↗</a>
                  <a href={sourceLink(check.source_path)}>Source ↗</a>
                </div>
                {Object.keys(check.metrics).length > 0 && (
                  <dl className="bench-metrics">
                    {Object.entries(check.metrics).map(([name, value]) => (
                      <div key={name}>
                        <dt>{name.replaceAll("_", " ")}</dt>
                        <dd>{value.toExponential(4)}</dd>
                      </div>
                    ))}
                  </dl>
                )}
              </div>
            ))}
          </details>
          <div className="bench-footer">
            <details>
              <summary>Reproduce this case</summary>
              <p>
                Install a Fortran compiler and the development dependencies,
                then run:
              </p>
              <pre>
                <code>{`git clone ${repository}.git\ncd abaqus_ufl\npip install -e ".[dev]"\npython tools/run_livebench.py --case ${id}`}</code>
              </pre>
              <a href={sourceLink("docs/LIVEBENCH.md")}>
                Setup instructions ↗
              </a>
            </details>
            <details>
              <summary>Recorded source & environment</summary>
              <p>
                Source{" "}
                {report!.source.commit ? (
                  <a href={`${repository}/commit/${report!.source.commit}`}>
                    <code>{report!.source.commit.slice(0, 8)}</code>
                  </a>
                ) : (
                  "unversioned"
                )}
                {report!.source.dirty ? " + local changes" : ""} · Python{" "}
                {report!.environment.python}
              </p>
              <p>
                {report!.environment.gfortran || "Fortran compiler unavailable"}
              </p>
              <p>
                {Object.entries(report!.environment.packages)
                  .map(
                    ([name, version]) => `${name} ${version || "unavailable"}`,
                  )
                  .join(" · ")}
              </p>
              {report!.run_url && (
                <p>
                  <a href={report!.run_url}>Recorded workflow run ↗</a>
                </p>
              )}
              <a href={`${resultBase}report.json`} download>
                Full report & plot data ↓
              </a>
            </details>
          </div>
        </>
      )}
    </>
  );
}

function WalkthroughStep({
  id,
  model,
  step,
  sourceLink,
}: {
  id: string;
  model: Walkthrough;
  step: string;
  sourceLink: (path: string) => string;
}) {
  const code = modelCode[id];
  if (step === "equations")
    return (
      <>
        <h2>Define the physical problem</h2>
        <p>{model.fields}</p>
        <Equations items={model.equations} />
        <div className="model-setup">
          <div>
            <h3>Domain & boundary conditions</h3>
            <p>{model.domain}</p>
            <ul>
              {model.boundaries.map((item) => (
                <li key={item}>{item}</li>
              ))}
            </ul>
          </div>
          <ProblemSketch id={id} />
        </div>
        <p className="walkthrough-note">
          <strong>Material parameters.</strong> {model.properties}
        </p>
      </>
    );
  if (step === "weak-form")
    return (
      <>
        <h2>Turn the balances into residuals</h2>
        <p>{model.testFunctions}</p>
        <Equations items={model.weakForms} />
        <article className="equation-card">
          <h3>Choose the element fields</h3>
          <p className="equation-expression">
            {model.constitutive ? "uₕ = Σₐ Nₐuₐ" : "Tₕ = Σₐ NₐTₐ, uₕ = Σₐ Nₐuₐ"}
          </p>
          <p>
            Bilinear Quad4 shape functions interpolate the nodal fields. The
            same basis supplies their test functions. Assemble the residuals and
            solve the free rows R = 0 with tangent Kₐᵦ = ∂Rₐ/∂Uᵦ.
          </p>
        </article>
        <div className="walkthrough-note">
          <h3>
            {model.constitutive
              ? "Connect the material law to element equilibrium"
              : "How the declaration represents the thermal weak form"}
          </h3>
          {model.constitutive === "ogden" ? (
            <p>
              The displacement-only declaration's <code>momentum_equation</code>{" "}
              returns first Piola stress. The generator integrates its
              contraction with the reference test-function gradient and
              differentiates the element residual.
            </p>
          ) : model.constitutive === "plasticity" ? (
            <p>
              The UMAT defines a local constitutive update. The verification
              host separately integrates <code>Bᵀσ</code> and{" "}
              <code>BᵀDDSDDE B</code>. Accepted integration-point history stays
              fixed throughout every Newton and line-search trial.
            </p>
          ) : (
            <p>
              <code>transport_equation</code> returns{" "}
              <code>(storage, flux)</code>. The generator assembles{" "}
              <code>storage × θ − flux · Grad θ</code>. Returning{" "}
              <code>flux = −k × grad_T</code> produces the positive diffusion
              term shown above.
            </p>
          )}
        </div>
        <a
          className="text-link"
          href={sourceLink(
            model.constitutive === "ogden"
              ? "examples/ogden_umat/README.md"
              : model.constitutive === "plasticity"
                ? "examples/small_strain_j2_umat/README.md"
                : "examples/scalar_diffusion_uel/README.md",
          )}
        >
          Read the model's conventions ↗
        </a>
      </>
    );
  if (step === "python")
    return (
      <>
        <h2>Declare the model in Python</h2>
        <p>{code.pythonIntro}</p>
        {code.python.map((block) => (
          <Fragment key={block.title}>
            <h3>{block.title}</h3>
            <p className="source-caption">
              <a href={sourceLink(block.path)}>{block.path} ↗</a> · excerpt
              from the shipped source
            </p>
            <pre className="walkthrough-code">
              <code>{block.code}</code>
            </pre>
          </Fragment>
        ))}
        {code.pythonNote && (
          <p className="walkthrough-note">{code.pythonNote}</p>
        )}
        {code.pythonMore && (
          <details className="walkthrough-source">
            <summary>{code.pythonMore.title}</summary>
            <p className="source-caption">
              <a href={sourceLink(code.pythonMore.path)}>
                {code.pythonMore.path} ↗
              </a>
            </p>
            <pre className="walkthrough-code">
              <code>{code.pythonMore.code}</code>
            </pre>
          </details>
        )}
        <details className="walkthrough-source">
          <summary>{code.boundaries.title}</summary>
          <p>
            The simulation script prescribes boundary DOFs at the end of each
            increment. The remaining nodal fields are solved.
          </p>
          <p className="source-caption">
            <a href={sourceLink(code.boundaries.path)}>
              {code.boundaries.path} ↗
            </a>
          </p>
          <pre className="walkthrough-code">
            <code>{code.boundaries.code}</code>
          </pre>
        </details>
      </>
    );
  return (
    <>
      <h2>{code.generationTitle}</h2>
      <p>{code.generationSummary}</p>
      <p className="source-caption">
        <a href={sourceLink(code.fortranPath)}>{code.fortranPath} ↗</a> ·{" "}
        {code.fortranCaption}, excerpt from the generated file
      </p>
      <pre className="walkthrough-code">
        <code>{code.fortranExcerpt}</code>
      </pre>
      <ol className="generation-path">
        <li>
          <strong>Generate</strong>
          <span>
            The Python declaration becomes the material law, residual and
            complex-step tangent in Fortran.
          </span>
        </li>
        <li>
          <strong>Compile</strong>
          <span>
            For these runs, f2py compiles the file with a small driver (
            <a href={sourceLink(code.compilePath)}>{code.compilePath} ↗</a>).
            The same file is a standard Abaqus user subroutine.
          </span>
        </li>
        <li>
          <strong>Execute</strong>
          <span>
            The runtime supplies properties, coordinates, current values and
            accepted old values.
          </span>
        </li>
      </ol>
      <details className="walkthrough-source">
        <summary>{code.interfaceTitle}</summary>
        <pre className="walkthrough-code">
          <code>{code.interface}</code>
        </pre>
        <a href={sourceLink(code.fortranPath)}>
          Open the complete generated Fortran ↗
        </a>
      </details>
      <p className="walkthrough-note">{code.runtimeNote}</p>
      <a
        className="text-link"
        href={sourceLink("docs/BOUNDARY_VALUE_SIMULATIONS.md")}
      >
        Numerical setup & acceptance gates ↗
      </a>
    </>
  );
}

function Catalog({ report }: { report: BenchmarkReport | null }) {
  const status = (id: string) =>
    report?.cases.find((item) => item.id === id)?.status;
  return (
    <>
      <p className="eyebrow">Livebench · worked models</p>
      <h1>Follow a model from equations to results</h1>
      <p className="walkthrough-lede">
        Choose a problem, then trace its equations, weak form, Python
        declaration and generated subroutine. The last step shows the computed
        solution and the checks behind it.
      </p>
      <div className="walkthrough-catalog">
        {Object.entries(walkthroughs).map(([id, model]) => (
          <PageLink className="walkthrough-card" href={caseUrl(id)} key={id}>
            <BenchmarkPreview
              caseId={id}
              plot={modelCards.find((card) => card.id === id)!.plot}
              alt={modelCards.find((card) => card.id === id)!.alt}
            />
            <div>
              <span className="case-kind">{model.kind}</span>
              <h2>{model.title}</h2>
              <p>{model.summary}</p>
              <span className="text-link">Start with the equations →</span>
              {status(id) && (
                <span className={`catalog-status bench-status ${status(id)}`}>
                  {status(id) === "passed"
                    ? "Recorded checks passed"
                    : "Recorded checks failed"}
                </span>
              )}
            </div>
          </PageLink>
        ))}
      </div>
      <details className="support-checks">
        <summary>
          Material & element checks{" "}
          <span>{Object.keys(supportCases).length} supporting cases</span>
        </summary>
        <p>
          These focused material-point, element and mesh checks support the
          generation pipeline. Open a case to inspect its numerical evidence.
        </p>
        <div className="support-grid">
          {Object.entries(supportCases).map(([id, info]) => (
            <PageLink href={caseUrl(id)} key={id}>
              <h3>{info.title}</h3>
              <p>{info.summary}</p>
              <span className="text-link">Open the check →</span>
            </PageLink>
          ))}
        </div>
      </details>
      <p className="walkthrough-note">
        The simulation results are recorded mesh runs using compiled generated
        subroutines, with independent quantitative checks. The separately
        documented paper examples contain the available Abaqus analysis
        evidence.
      </p>
    </>
  );
}

export default function LivebenchPage() {
  const [search, setSearch] = useState(window.location.search);
  const [report, setReport] = useState<BenchmarkReport | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [refresh, setRefresh] = useState(0);
  const params = new URLSearchParams(search),
    id = params.get("case") || "";
  const model = walkthroughs[id],
    support = supportCases[id];
  const step = steps.find((item) => item.id === params.get("step")) ?? steps[0];
  const stepIndex = steps.indexOf(step);
  const sourceLink = (path: string) =>
    `${repository}/blob/${report?.source.commit || "main"}/${path}`;
  useEffect(() => {
    const update = () => setSearch(window.location.search);
    window.addEventListener("popstate", update);
    return () => window.removeEventListener("popstate", update);
  }, []);
  useEffect(() => {
    document.title = `${model?.title || support?.title || "Livebench"}${model ? ` · ${step.title}` : ""} — abaqus_ufl`;
  }, [model, support, step]);
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
          throw new Error("A recorded benchmark report is not available.");
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
  const resultProps = {
    id,
    report,
    loading,
    error,
    refresh: () => setRefresh((n) => n + 1),
  };
  return (
    <>
      <SiteHeader livebench />
      <main id="main" className="livebench-section walkthrough-page">
        <div className="shell">
          <nav className="walkthrough-breadcrumb" aria-label="Breadcrumb">
            <a href={import.meta.env.BASE_URL}>Home</a>
            <span aria-hidden="true">/</span>
            <PageLink href={livebenchUrl}>Livebench</PageLink>
            {id && (
              <>
                <span aria-hidden="true">/</span>
                <span>{model?.title || support?.title || "Unknown case"}</span>
              </>
            )}
          </nav>
          {!id ? (
            <Catalog report={report} />
          ) : !model && !support ? (
            <>
              <h1>Case not found</h1>
              <p>
                Choose a model from the{" "}
                <PageLink href={livebenchUrl}>livebench index</PageLink>.
              </p>
            </>
          ) : (
            <>
              <p className="eyebrow">
                {model?.kind || "Material & element evidence"}
              </p>
              <h1>{model?.title || support.title}</h1>
              <p className="walkthrough-lede">
                {model?.summary || support.summary}
              </p>
              {model ? (
                <div className="walkthrough-layout">
                  <nav
                    className="walkthrough-steps"
                    aria-label="Model walkthrough"
                  >
                    {steps.map((item, index) => (
                      <PageLink
                        key={item.id}
                        href={caseUrl(id, item.id)}
                        aria-current={step.id === item.id ? "step" : undefined}
                      >
                        <span className="walkthrough-step-number">
                          0{index + 1}
                        </span>
                        <span>
                          <strong>{item.title}</strong>
                          <small>{item.detail}</small>
                        </span>
                      </PageLink>
                    ))}
                  </nav>
                  <article
                    className="walkthrough-body"
                    aria-label={step.title}
                    key={`${id}-${step.id}`}
                  >
                    {step.id === "simulation" ? (
                      <>
                        <CaseResults {...resultProps} />
                        <p className="walkthrough-note">{model.scope}</p>
                      </>
                    ) : (
                      <WalkthroughStep
                        id={id}
                        model={model}
                        step={step.id}
                        sourceLink={sourceLink}
                      />
                    )}
                    <nav
                      className="walkthrough-next"
                      aria-label="Walkthrough pagination"
                    >
                      {stepIndex > 0 ? (
                        <PageLink href={caseUrl(id, steps[stepIndex - 1].id)}>
                          ← {steps[stepIndex - 1].title}
                        </PageLink>
                      ) : (
                        <PageLink href={livebenchUrl}>← All models</PageLink>
                      )}
                      {stepIndex < steps.length - 1 ? (
                        <PageLink href={caseUrl(id, steps[stepIndex + 1].id)}>
                          Next: {steps[stepIndex + 1].title} →
                        </PageLink>
                      ) : (
                        <PageLink href={livebenchUrl}>
                          Explore another model →
                        </PageLink>
                      )}
                    </nav>
                  </article>
                </div>
              ) : (
                <article className="walkthrough-body support-results" key={id}>
                  <CaseResults {...resultProps} />
                </article>
              )}
            </>
          )}
        </div>
      </main>
      <footer className="walkthrough-footer">
        <a href={import.meta.env.BASE_URL}>abaqus_ufl project</a>
        <a href={repository}>Source on GitHub ↗</a>
      </footer>
    </>
  );
}

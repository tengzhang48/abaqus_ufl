import siteData from "../site-data.json";
import type { MouseEvent } from "react";

import tetDisplacement from "../../paper_examples/stabilized_tet4/figure/Tet4_u_mag.png";
import tetTheta from "../../paper_examples/stabilized_tet4/figure/Tet4_NT11.png";
import pasta0 from "../../paper_examples/morphing_hex8/figure/pasta_t0.png";
import pasta45 from "../../paper_examples/morphing_hex8/figure/pasta_t45.png";
import pasta90 from "../../paper_examples/morphing_hex8/figure/pasta_t90.png";
import pasta150 from "../../paper_examples/morphing_hex8/figure/pasta_t150.png";
import pasta200 from "../../paper_examples/morphing_hex8/figure/pasta_t200.png";
import pasta360 from "../../paper_examples/morphing_hex8/figure/pasta_t360.png";
import gel0 from "../../paper_examples/gel_bilayer/figure/gel_bilayer_00min.png";
import gel30 from "../../paper_examples/gel_bilayer/figure/gel_bilayer_30min.png";
import gel60 from "../../paper_examples/gel_bilayer/figure/gel_bilayer_1h.png";
import gel360 from "../../paper_examples/gel_bilayer/figure/gel_bilayer_6h.png";

const repository = siteData.repository.url;

type FigureAsset = {
  src: string;
  path: string;
  alt: string;
  caption: string;
};

const figureAssets: Record<string, FigureAsset[]> = {
  tet4: [
    {
      src: tetDisplacement,
      path: siteData.paperEvidence[0].figurePaths[0],
      alt: "Abaqus rendering of displacement magnitude in the compressed stabilized Tet4 block",
      caption: "Displacement magnitude at full follower pressure",
    },
    {
      src: tetTheta,
      path: siteData.paperEvidence[0].figurePaths[1],
      alt: "Abaqus rendering of the stabilized Tet4 theta field in the compressed block",
      caption: "θ̃ field from the same Abaqus run",
    },
  ],
  pasta: [
    { src: pasta0, path: siteData.paperEvidence[1].figurePaths[0], alt: "Frame 1 of six from the grooved gel sheet morphing sequence", caption: "Sequence frame 1" },
    { src: pasta45, path: siteData.paperEvidence[1].figurePaths[1], alt: "Frame 2 of six from the grooved gel sheet morphing sequence", caption: "Sequence frame 2" },
    { src: pasta90, path: siteData.paperEvidence[1].figurePaths[2], alt: "Frame 3 of six from the grooved gel sheet morphing sequence", caption: "Sequence frame 3" },
    { src: pasta150, path: siteData.paperEvidence[1].figurePaths[3], alt: "Frame 4 of six from the grooved gel sheet morphing sequence", caption: "Sequence frame 4" },
    { src: pasta200, path: siteData.paperEvidence[1].figurePaths[4], alt: "Frame 5 of six from the grooved gel sheet morphing sequence", caption: "Sequence frame 5" },
    { src: pasta360, path: siteData.paperEvidence[1].figurePaths[5], alt: "Frame 6 of six from the grooved gel sheet morphing sequence", caption: "Sequence frame 6" },
  ],
  gel: [
    { src: gel0, path: siteData.paperEvidence[2].figurePaths[0], alt: "Gel bilayer in its initial flat state", caption: "0 min" },
    { src: gel30, path: siteData.paperEvidence[2].figurePaths[1], alt: "Gel bilayer bending after 30 minutes", caption: "30 min" },
    { src: gel60, path: siteData.paperEvidence[2].figurePaths[2], alt: "Gel bilayer bending after one hour", caption: "1 h" },
    { src: gel360, path: siteData.paperEvidence[2].figurePaths[3], alt: "Gel bilayer bent configuration after six hours", caption: "6 h" },
  ],
  corrosion: [],
};

function repoFile(path: string) {
  return `${repository}/blob/main/${path}`;
}

function Arrow() {
  return <span aria-hidden="true">↗</span>;
}

function closeMobileMenu(event: MouseEvent<HTMLAnchorElement>) {
  event.currentTarget.closest("details")?.removeAttribute("open");
}

function App() {
  return (
    <>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="abaqus ufl home">
          <span className="brand-mark" aria-hidden="true">a/u</span>
          <span>abaqus_<strong>ufl</strong></span>
        </a>
        <nav className="desktop-nav" aria-label="Project website">
          <a href="#how-it-works">How it works</a>
          <a href="#scope">Scope</a>
          <a href="#examples">Examples</a>
          <a href="#paper-evidence">Manuscript examples</a>
          <a href="#status">Status</a>
        </nav>
        <a className="header-link" href={repository}>GitHub <Arrow /></a>
        <details className="mobile-menu">
          <summary aria-label="Open navigation"><span aria-hidden="true">Menu</span></summary>
          <nav aria-label="Mobile project website">
            <a href="#how-it-works" onClick={closeMobileMenu}>How it works</a>
            <a href="#scope" onClick={closeMobileMenu}>Scope</a>
            <a href="#examples" onClick={closeMobileMenu}>Examples</a>
            <a href="#paper-evidence" onClick={closeMobileMenu}>Manuscript examples</a>
            <a href="#status" onClick={closeMobileMenu}>Status</a>
            <a href={repository} onClick={closeMobileMenu}>GitHub</a>
          </nav>
        </details>
      </header>

      <main id="main">
        <section className="hero" id="top">
          <div className="hero-grid shell">
            <div className="hero-copy">
              <p className="eyebrow light">Open research software · v{siteData.repository.version}</p>
              <h1>Generate Abaqus UMAT and UEL source from Python models.</h1>
              <p className="hero-lede">
                <code>abaqus_ufl</code> generates self-contained fixed-form Fortran UMATs and UELs from supported Python models. The Python model and generated source show the fields, constitutive response, tangent construction, and Abaqus interface.
              </p>
              <div className="hero-actions">
                <a className="button primary" href={`${repository}#try-a-complete-local-workflow`}>Run an example <Arrow /></a>
                <a className="button ghost" href="#paper-evidence">View manuscript examples</a>
              </div>
              <p className="hero-boundary">
                The package supports a defined set of models. It does not implement FEniCS UFL, compile arbitrary weak forms, or create a complete Abaqus analysis.
              </p>
            </div>

            <aside className="declaration-card" aria-label="Example Python declaration">
              <div className="codebar">
                <span>neo_hookean.py</span>
                <span className="code-status">Python → UMAT</span>
              </div>
              <pre><code><span className="kw">class</span> <span className="type">NeoHookean</span>(au.Material):{"\n"}  props = dict(G=<span className="num">0.5</span>, K=<span className="num">50.0</span>){"\n\n"}  <span className="kw">def</span> <span className="fn">stress_PK1</span>(self, F):{"\n"}    J = det(F){"\n"}    FinvT = inv(F).T{"\n"}    <span className="kw">return</span> (self.G * (F - FinvT) +{"\n"}      self.K * log(J) * FinvT){"\n\n"}model = NeoHookean(){"\n"}model.verify(){"\n"}au.generate_umat(model, <span className="str">"model.for"</span>)</code></pre>
              <div className="code-footer">
                <span><i className="signal green" /> material and tangent checks</span>
                <span><i className="signal amber" /> application validation depends on the example</span>
              </div>
            </aside>
          </div>
        </section>

        <section className="workflow-section section" id="how-it-works">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">Generation and verification</p>
                <h2>From a Python model to an Abaqus analysis</h2>
              </div>
              <p>The Python model, compiled subroutine, and Abaqus analysis are checked separately because they can fail in different ways.</p>
            </div>
            <ol className="workflow" aria-label="abaqus ufl workflow">
              {siteData.workflow.map((step) => (
                <li key={step.number}>
                  <span className="step-number">{step.number}</span>
                  <h3>{step.title}</h3>
                  <p>{step.detail}</p>
                </li>
              ))}
            </ol>
          </div>
        </section>

        <section className="section scope-section" id="scope">
          <div className="shell">
            <div className="section-heading centered">
              <p className="eyebrow">Current generation scope</p>
              <h2>Supported UMAT and UEL generation</h2>
              <p>The current API covers selected material-response and coupled-element patterns. Users still define the mesh, steps, loads, boundary conditions, units, and solver controls in the Abaqus model.</p>
            </div>
            <div className="scope-grid">
              {siteData.scopes.map((scope) => (
                <article className="scope-card" key={scope.target}>
                  <div className="scope-title">
                    <span>{scope.target}</span>
                    <h3>{scope.title}</h3>
                  </div>
                  <p>{scope.description}</p>
                  <ul>
                    {scope.included.map((item) => <li key={item}>{item}</li>)}
                  </ul>
                  <div className="boundary compact">
                    <strong>Limits</strong>
                    <p>{scope.boundary}</p>
                  </div>
                </article>
              ))}
            </div>
            <p className="scope-note">
              The declaration style is inspired by FEniCS UFL, but <code>abaqus_ufl</code> neither depends on nor implements UFL. <a href={repoFile("CREDITS.md")}>See credits and references <Arrow /></a>
            </p>
          </div>
        </section>

        <section className="section examples-section" id="examples">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">Included examples</p>
                <h2>Four UMAT and two UEL examples</h2>
              </div>
              <p>Each directory contains the Python model, a reference check, generated source, compiled checks, and any Abaqus result available for that example.</p>
            </div>
            <div className="examples-grid">
              {siteData.examples.map((example, index) => (
                <article className="example-card" key={example.title}>
                  <div className="example-topline">
                    <span className={`target ${example.target.toLowerCase()}`}>{example.target}</span>
                    <span className="example-index">0{index + 1}</span>
                  </div>
                  <h3>{example.title}</h3>
                  <p>{example.summary}</p>
                  <dl>
                    <dt>Checks included</dt>
                    <dd>{example.evidence}</dd>
                    <dt>Limits</dt>
                    <dd>{example.boundary}</dd>
                  </dl>
                  <a className="text-link" href={repoFile(example.path)}>Open example <Arrow /></a>
                </article>
              ))}
            </div>
          </div>
        </section>

        <section className="section paper-section" id="paper-evidence">
          <div className="shell">
            <div className="section-heading paper-heading">
              <p className="eyebrow light">Submitted manuscript</p>
              <h2>Examples from <em>Making coupled-field Abaqus user elements simple</em></h2>
              <p>These four example packages contain the available declarations, generated sources, selected Abaqus decks, reduced results, and figure inputs. Each package states which checks were repeated from a fresh clone.</p>
            </div>

            {siteData.paperEvidence.map((paper, paperIndex) => {
              const figures = figureAssets[paper.id];
              return (
                <article className={`paper-case ${paper.id}`} key={paper.id}>
                  <div className="paper-copy">
                    <div className="paper-number">0{paperIndex + 1}</div>
                    <p className="eyebrow">{paper.eyebrow}</p>
                    <h3>{paper.title}</h3>
                    <p className="paper-summary">{paper.summary}</p>
                    <div className="result-line"><strong>Available result</strong><p>{paper.result}</p></div>
                    <div className="result-line"><strong>Fresh-clone check</strong><p>{paper.freshCheck}</p></div>
                    <div className="boundary"><strong>Limits</strong><p>{paper.boundary}</p></div>
                    <a className="text-link" href={repoFile(paper.path)}>Read the package record <Arrow /></a>
                  </div>
                  {figures.length > 0 ? (
                    <div className={`figure-grid figure-grid-${figures.length}`}>
                      {figures.map((figure) => (
                        <figure key={figure.path}>
                          <a href={repoFile(figure.path)} aria-label={`Open source figure: ${figure.caption}`}>
                            <img src={figure.src} alt={figure.alt} loading="lazy" />
                          </a>
                          <figcaption>{figure.caption} <a href={repoFile(figure.path)}>source</a></figcaption>
                        </figure>
                      ))}
                    </div>
                  ) : (
                    <aside className="provenance-panel" aria-label="Corrosion artifact provenance boundary">
                      <span className="provenance-mark" aria-hidden="true">P</span>
                      <div>
                        <p className="eyebrow light">Text record only</p>
                        <h4>Figure not republished on this site</h4>
                        <p>The comparison mesh came from a third-party distribution whose complete BSD license notice has not been located. The corrosion figure is therefore omitted from this site.</p>
                        <a className="text-link" href={repoFile("CREDITS.md")}>Read the provenance record <Arrow /></a>
                      </div>
                    </aside>
                  )}
                </article>
              );
            })}

            <div className="paper-index-link">
              <p>The repository keeps the submitted source and the current generated version separately when they differ.</p>
              <a className="button light-button" href={repoFile("paper_examples/README.md")}>Open the paper-package index <Arrow /></a>
            </div>
          </div>
        </section>

        <section className="section status-section" id="status">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">Fresh-clone record · {siteData.validation.date}</p>
                <h2>Fresh-clone checks</h2>
              </div>
              <p>The check began at revision <code>{siteData.validation.sourceRevision}</code>; corrections and the report are in <code>{siteData.validation.reportRevision}</code>. Abaqus checks used Abaqus/Standard 2022, Intel Fortran 19.1.1.217, and the documented Python environment.</p>
            </div>
            <div className="metric-grid">
              {siteData.validation.metrics.map((metric) => (
                <div className="metric" key={metric.label}>
                  <strong>{metric.value}</strong>
                  <span>{metric.label}</span>
                </div>
              ))}
            </div>
            <div className="status-grid">
              <article>
                <span className="status-kicker complete">Completed solver checks</span>
                <h3>Small verification analyses</h3>
                <ul>{siteData.validation.completeSolves.map((item) => <li key={item}>{item}</li>)}</ul>
              </article>
              <article>
                <span className="status-kicker datacheck">Datacheck only</span>
                <h3>Manuscript input and compilation checks</h3>
                <ul>{siteData.validation.datachecks.map((item) => <li key={item}>{item}</li>)}</ul>
              </article>
            </div>
            <div className="boundary status-boundary">
              <strong>What an Abaqus datacheck establishes</strong>
              <p>{siteData.validation.boundary}</p>
              <a className="text-link" href={repoFile(siteData.validation.sourcePath)}>Read the complete validation report <Arrow /></a>
            </div>
          </div>
        </section>

        <section className="section start-section" id="start">
          <div className="shell start-grid">
            <div>
              <p className="eyebrow light">Neo-Hookean UMAT example</p>
              <h2>Generate and verify the example locally</h2>
              <p>These commands generate the UMAT, run its closed-form reference checks, compile it with gfortran and f2py, and compare its output with the Python calculation. Python 3.8+, NumPy, and SymPy are required; compiled checks also require gfortran, f2py, Meson, and Ninja.</p>
              <div className="hero-actions">
                <a className="button primary" href={repoFile("docs/API_USAGE.md")}>Read the API guide <Arrow /></a>
                <a className="button ghost" href={repoFile("HOWTO_ADD_AN_EXAMPLE.md")}>Example contract</a>
              </div>
            </div>
            <pre className="install"><code><span>$</span> git clone {repository}.git{"\n"}<span>$</span> cd abaqus_ufl{"\n"}<span>$</span> pip install -e <b>".[dev]"</b>{"\n"}<span>$</span> cd examples/neo_hookean_umat{"\n"}<span>$</span> python build.py{"\n"}<span>$</span> python check_reference.py{"\n"}<span>$</span> python check_compiled.py</code></pre>
          </div>
        </section>
      </main>

      <footer>
        <div className="shell footer-grid">
          <div className="footer-brand">
            <span className="brand-mark" aria-hidden="true">a/u</span>
            <div><strong>abaqus_ufl</strong><p>Project-authored source is MIT licensed. Third-party materials and items with unresolved redistribution terms are documented separately.</p></div>
          </div>
          <nav aria-label="Documentation">
            <strong>Documentation</strong>
            <a href={repoFile("docs/README.md")}>Documentation index</a>
            <a href={repoFile("docs/API_USAGE.md")}>API usage</a>
            <a href={repoFile("docs/theory.md")}>Theory & conventions</a>
            <a href={repoFile("docs/lessons/README.md")}>Lessons learned</a>
          </nav>
          <nav aria-label="Project records">
            <strong>Project records</strong>
            <a href={repoFile("CITATION.cff")}>Citation</a>
            <a href={repoFile("CREDITS.md")}>Credits & provenance</a>
            <a href={repoFile("LICENSE")}>MIT license</a>
            <a href={`${import.meta.env.BASE_URL}THIRD_PARTY_NOTICES.txt`}>Website third-party notices</a>
            <a href={siteData.repository.contact}>Contact Teng Zhang</a>
          </nav>
        </div>
        <div className="shell footer-bottom">
          <span>v{siteData.repository.version} · authored and maintained by {siteData.repository.author}</span>
          <a href={repository}>View source on GitHub <Arrow /></a>
        </div>
      </footer>
    </>
  );
}

export default App;

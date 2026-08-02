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
      caption: "Retained θ̃ field from the same completed run",
    },
  ],
  pasta: [
    { src: pasta0, path: siteData.paperEvidence[1].figurePaths[0], alt: "First retained export of the grooved gel sheet morphing sequence", caption: "Retained export 1" },
    { src: pasta45, path: siteData.paperEvidence[1].figurePaths[1], alt: "Second retained export of the grooved gel sheet morphing sequence", caption: "Retained export 2" },
    { src: pasta90, path: siteData.paperEvidence[1].figurePaths[2], alt: "Third retained export of the grooved gel sheet morphing sequence", caption: "Retained export 3" },
    { src: pasta150, path: siteData.paperEvidence[1].figurePaths[3], alt: "Fourth retained export of the grooved gel sheet morphing sequence", caption: "Retained export 4" },
    { src: pasta200, path: siteData.paperEvidence[1].figurePaths[4], alt: "Fifth retained export of the grooved gel sheet morphing sequence", caption: "Retained export 5" },
    { src: pasta360, path: siteData.paperEvidence[1].figurePaths[5], alt: "Sixth retained export of the grooved gel sheet morphing sequence", caption: "Retained export 6" },
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
          <a href="#paper-evidence">Paper evidence</a>
          <a href="#status">Status</a>
        </nav>
        <a className="header-link" href={repository}>GitHub <Arrow /></a>
        <details className="mobile-menu">
          <summary aria-label="Open navigation"><span aria-hidden="true">Menu</span></summary>
          <nav aria-label="Mobile project website">
            <a href="#how-it-works" onClick={closeMobileMenu}>How it works</a>
            <a href="#scope" onClick={closeMobileMenu}>Scope</a>
            <a href="#examples" onClick={closeMobileMenu}>Examples</a>
            <a href="#paper-evidence" onClick={closeMobileMenu}>Paper evidence</a>
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
              <h1>From Python declarations to inspectable Abaqus user subroutines.</h1>
              <p className="hero-lede">
                <code>abaqus_ufl</code> generates self-contained Fortran for supported UMAT and UEL models—while keeping field definitions, constitutive choices, tangent construction, and interface conventions open to review.
              </p>
              <div className="hero-actions">
                <a className="button primary" href={`${repository}#try-a-complete-local-workflow`}>Try the local workflow <Arrow /></a>
                <a className="button ghost" href="#paper-evidence">Inspect the evidence</a>
              </div>
              <p className="hero-boundary">
                Not UFL-compatible, not a general variational compiler, and not a generator for complete Abaqus analyses.
              </p>
            </div>

            <aside className="declaration-card" aria-label="Example Python declaration">
              <div className="codebar">
                <span>neo_hookean.py</span>
                <span className="code-status">Python → UMAT</span>
              </div>
              <pre><code><span className="kw">class</span> <span className="type">NeoHookean</span>(au.Material):{"\n"}  props = dict(G=<span className="num">0.5</span>, K=<span className="num">50.0</span>){"\n\n"}  <span className="kw">def</span> <span className="fn">stress_PK1</span>(self, F):{"\n"}    J = det(F){"\n"}    FinvT = inv(F).T{"\n"}    <span className="kw">return</span> (self.G * (F - FinvT) +{"\n"}      self.K * log(J) * FinvT){"\n\n"}model = NeoHookean(){"\n"}model.verify(){"\n"}au.generate_umat(model, <span className="str">"model.for"</span>)</code></pre>
              <div className="code-footer">
                <span><i className="signal green" /> declaration check</span>
                <span><i className="signal amber" /> scientific validation is example-owned</span>
              </div>
            </aside>
          </div>
        </section>

        <section className="workflow-section section" id="how-it-works">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">A visible generation path</p>
                <h2>Keep the declaration short. Keep the generated boundary inspectable.</h2>
              </div>
              <p>Each gate catches a different class of mistake. Passing an earlier gate does not replace assembled checks, solver execution, or comparison with physical evidence.</p>
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
              <h2>Two targets, with different responsibilities</h2>
              <p>The public API supports selected material-response and coupled-element patterns. The surrounding Abaqus model stays explicit and user-owned.</p>
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
                    <strong>Boundary</strong>
                    <p>{scope.boundary}</p>
                  </div>
                </article>
              ))}
            </div>
            <p className="scope-note">
              The declaration style is inspired by FEniCS UFL, but <code>abaqus_ufl</code> neither depends on nor implements UFL. <a href={repoFile("CREDITS.md")}>Read the scientific lineage <Arrow /></a>
            </p>
          </div>
        </section>

        <section className="section examples-section" id="examples">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">Six release examples</p>
                <h2>Small enough to audit, complete enough to exercise the boundary</h2>
              </div>
              <p>These are the public allowlist—not a capability inventory for a larger research archive. Each bundle records its own oracle, generated source, compiled path, and execution boundary.</p>
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
                    <dt>Retained evidence</dt>
                    <dd>{example.evidence}</dd>
                    <dt>Boundary</dt>
                    <dd>{example.boundary}</dd>
                  </dl>
                  <a className="text-link" href={repoFile(example.path)}>Open example record <Arrow /></a>
                </article>
              ))}
            </div>
          </div>
        </section>

        <section className="section paper-section" id="paper-evidence">
          <div className="shell">
            <div className="section-heading paper-heading">
              <p className="eyebrow light">Manuscript-scale evidence</p>
              <h2>Results first. Source, provenance, and limits beside them.</h2>
              <p>The four packages accompany <em>Making coupled-field Abaqus user elements simple</em>. They retain declarations, sources, selected decks, reduced records, and figure inputs at different evidence levels.</p>
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
                    <div className="result-line"><strong>Retained result</strong><p>{paper.result}</p></div>
                    <div className="result-line"><strong>Fresh-clone check</strong><p>{paper.freshCheck}</p></div>
                    <div className="boundary"><strong>Claim boundary</strong><p>{paper.boundary}</p></div>
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
                        <p>The public repository records a derived-mesh provenance item that must be resolved before further distribution. This website therefore does not copy the corrosion image into its Pages artifact.</p>
                        <a className="text-link" href={repoFile("CREDITS.md")}>Read credits and redistribution boundaries <Arrow /></a>
                      </div>
                    </aside>
                  )}
                </article>
              );
            })}

            <div className="paper-index-link">
              <p>Exact submitted sources and current generator outputs are retained separately where they differ.</p>
              <a className="button light-button" href={repoFile("paper_examples/README.md")}>Open the paper-package index <Arrow /></a>
            </div>
          </div>
        </section>

        <section className="section status-section" id="status">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">Fresh-clone record · {siteData.validation.date}</p>
                <h2>What was actually exercised</h2>
              </div>
              <p>A fresh clone started from <code>{siteData.validation.sourceRevision}</code>; corrections and this report are retained in <code>{siteData.validation.reportRevision}</code>. Checks used Abaqus/Standard 2022, Intel Fortran 19.1.1.217, and the documented Python environment.</p>
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
                <h3>Small analyses run through Abaqus</h3>
                <ul>{siteData.validation.completeSolves.map((item) => <li key={item}>{item}</li>)}</ul>
              </article>
              <article>
                <span className="status-kicker datacheck">Datacheck only</span>
                <h3>Paper-scale setup and compilation checks</h3>
                <ul>{siteData.validation.datachecks.map((item) => <li key={item}>{item}</li>)}</ul>
              </article>
            </div>
            <div className="boundary status-boundary">
              <strong>How to read this record</strong>
              <p>{siteData.validation.boundary}</p>
              <a className="text-link" href={repoFile(siteData.validation.sourcePath)}>Read the complete validation report <Arrow /></a>
            </div>
          </div>
        </section>

        <section className="section start-section" id="start">
          <div className="shell start-grid">
            <div>
              <p className="eyebrow light">Start with one inspectable material</p>
              <h2>Generate locally. Review the Fortran. Add evidence for your claim.</h2>
              <p>Python-only generation requires Python 3.8+, NumPy, and SymPy. Compiled example checks additionally use gfortran, f2py, Meson, and Ninja. Abaqus is not required by the package CI.</p>
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
            <div><strong>abaqus_ufl</strong><p>Project-authored source is MIT licensed. Third-party and provenance-pending artifacts retain their separately documented terms.</p></div>
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

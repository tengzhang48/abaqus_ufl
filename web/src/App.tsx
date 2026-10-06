import siteData from "../site-data.json";
import type { MouseEvent } from "react";
import { useEffect } from "react";
import publication from "../../docs/publication.json";
import Livebench from "./Livebench";
import PaperExamples from "./PaperExamples";
import type { FigureAsset } from "./PaperExamples";

import tetDisplacement from "./assets/figures/Tet4_u_mag.webp";
import tetTheta from "./assets/figures/Tet4_NT11.webp";
import pasta0 from "./assets/figures/pasta_t0.webp";
import pasta45 from "./assets/figures/pasta_t45.webp";
import pasta90 from "./assets/figures/pasta_t90.webp";
import pasta150 from "./assets/figures/pasta_t150.webp";
import pasta200 from "./assets/figures/pasta_t200.webp";
import pasta360 from "./assets/figures/pasta_t360.webp";
import gel0 from "./assets/figures/gel_bilayer_00min.webp";
import gel30 from "./assets/figures/gel_bilayer_30min.webp";
import gel60 from "./assets/figures/gel_bilayer_1h.webp";
import gel360 from "./assets/figures/gel_bilayer_6h.webp";

const repository = siteData.repository.url;

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
  useEffect(() => {
    const hash = window.location.hash.slice(1);
    const anchor = hash === "examples" ? "livebench" : hash;
    if (!anchor) return;
    const frame = requestAnimationFrame(() => document.getElementById(anchor)?.scrollIntoView({ behavior: "instant" }));
    return () => cancelAnimationFrame(frame);
  }, []);
  return (
    <>
      <header className="site-header">
        <a className="brand" href="#top" aria-label="abaqus ufl home">
          <span className="brand-mark" aria-hidden="true">a/u</span>
          <span>abaqus_<strong>ufl</strong></span>
        </a>
        <nav className="desktop-nav" aria-label="Project website">
          <a href="#livebench">Livebench</a>
          <a href="#how-it-works">How it works</a>
          <a href="#paper-evidence">Paper examples</a>
          <a href="#status">Status</a>
        </nav>
        <a className="header-link" href={repository}>GitHub <Arrow /></a>
        <details className="mobile-menu">
          <summary aria-label="Open navigation"><span aria-hidden="true">Menu</span></summary>
          <nav aria-label="Mobile project website">
            <a href="#livebench" onClick={closeMobileMenu}>Livebench</a>
            <a href="#how-it-works" onClick={closeMobileMenu}>How it works</a>
            <a href="#paper-evidence" onClick={closeMobileMenu}>Paper examples</a>
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
                Define a material or coupled-field element in Python. Generate Fortran, check its numerical response, and use it in Abaqus.
              </p>
              <div className="hero-actions">
                <a className="button primary" href="#livebench">Open the livebench <span aria-hidden="true">↓</span></a>
                <a className="button ghost" href={publication.url}>Read the published paper <Arrow /></a>
              </div>
              <p className="publication-link">Published in <em>{publication.journal}</em> {publication.volume} ({publication.year}), {publication.articleNumber} · <a href={publication.url}>DOI <Arrow /></a></p>
              <p className="hero-note">
                Run the benchmarks with f2py, without an Abaqus installation.
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

        <Livebench />

        <section className="workflow-section section" id="how-it-works">
          <div className="shell">
            <div className="workflow-heading">
              <div>
                <p className="eyebrow">How it works</p>
                <h2>Python → Fortran → Abaqus</h2>
              </div>
              <a className="text-link" href={repoFile("docs/API_USAGE.md")}>API guide <Arrow /></a>
            </div>
            <ol className="workflow" aria-label="abaqus ufl workflow">
              {siteData.workflow.map((step) => (
                <li key={step.number}>
                  <span className="step-number">{step.number}</span>
                  <div>
                    <h3>{step.title}</h3>
                    <p>{step.detail}</p>
                  </div>
                </li>
              ))}
            </ol>
            <details className="generation-details" id="scope">
              <summary>Supported formulations and limits</summary>
              <div className="generation-detail-body">
                <div className="formulation-grid">
                  {siteData.scopes.map((scope) => (
                    <article key={scope.target}>
                      <h3><code>{scope.target}</code> · {scope.title}</h3>
                      <p>{scope.description}</p>
                      <ul>
                        {scope.included.map((item) => <li key={item}>{item}</li>)}
                      </ul>
                      <p>{scope.boundary}</p>
                    </article>
                  ))}
                </div>
                <p>Generation produces self-contained fixed-form Fortran. Python checks, compiled subroutine comparisons, and Abaqus analyses provide separate evidence; each livebench case records which checks are available.</p>
                <p>The declaration style is inspired by FEniCS UFL, but <code>abaqus_ufl</code> neither depends on nor implements UFL. <a href={repoFile("CREDITS.md")}>Credits and references <Arrow /></a></p>
              </div>
            </details>
          </div>
        </section>

        <PaperExamples figures={figureAssets} />

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
            <a href={publication.url}>Published paper</a>
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

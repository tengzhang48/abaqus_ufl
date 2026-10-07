import { useEffect } from "react";
import siteData from "../site-data.json";
import publication from "../../docs/publication.json";
import Livebench from "./Livebench";
import BenchmarkPreview from "./BenchmarkPreview";
import SiteHeader, { livebenchUrl } from "./SiteHeader";
import tetDisplacement from "../../paper_examples/stabilized_tet4/figure/Tet4_u_mag.png";
import pastaFinal from "../../paper_examples/morphing_hex8/figure/pasta_t360.png";
import gelFinal from "../../paper_examples/gel_bilayer/figure/gel_bilayer_6h.png";

const repository = siteData.repository.url;
const repoFile = (path: string) => `${repository}/blob/main/${path}`;
const paperImages = [
  {
    src: tetDisplacement,
    path: siteData.paperEvidence[0].figurePaths[0],
    alt: "Historical Abaqus displacement field in the compressed stabilized Tet4 block",
    caption: "Displacement at full follower pressure",
  },
  {
    src: pastaFinal,
    path: siteData.paperEvidence[1].figurePaths[5],
    alt: "Sixth retained Abaqus frame of a diffusion-driven grooved gel sheet",
    caption: "Retained morphing sequence · frame 6",
  },
  {
    src: gelFinal,
    path: siteData.paperEvidence[2].figurePaths[3],
    alt: "Historical Abaqus image of the gel bilayer after six hours",
    caption: "Retained gel-bilayer export · 6 h",
  },
];

function Arrow() {
  return <span aria-hidden="true">↗</span>;
}

export default function App() {
  useEffect(() => {
    const anchor = window.location.hash.slice(1);
    if (!anchor) return;
    const frame = requestAnimationFrame(() =>
      document.getElementById(anchor)?.scrollIntoView({ behavior: "instant" }),
    );
    return () => cancelAnimationFrame(frame);
  }, []);
  return (
    <>
      <SiteHeader />
      <main id="main" className="software-home">
        <section className="hero home-hero" id="top">
          <div className="hero-grid shell">
            <div className="hero-copy">
              <p className="eyebrow light">
                Open research software · v{siteData.repository.version}
              </p>
              <h1>
                From Python
                <br />
                to Abaqus Fortran.
              </h1>
              <p className="hero-lede">
                Write material laws and coupled-field declarations in Python.
                Generate inspectable UMAT and UEL Fortran, with consistent
                tangents and reproducible checks.
              </p>
              <div className="hero-actions">
                <a className="button primary" href="#livebench">
                  Explore the livebench <span aria-hidden="true">→</span>
                </a>
                <a
                  className="button ghost"
                  href={repoFile("docs/API_USAGE.md")}
                >
                  Read the documentation <Arrow />
                </a>
              </div>
              <div className="hero-facts" aria-label="Package features">
                <span>UMAT + UEL</span>
                <span>Self-contained Fortran</span>
                <span>MIT source</span>
              </div>
              <p className="hero-paper">
                Published in{" "}
                <a href={publication.url}>
                  <em>{publication.journal}</em> <Arrow />
                </a>
              </p>
            </div>
            <a
              className="hero-model"
              href={`${livebenchUrl}?case=ogden_bvp`}
              aria-label="Open the Ogden equation-to-simulation walkthrough"
            >
              <div className="hero-model-heading">
                <span>Ogden elasticity</span>
                <span>Python → UEL → simulation</span>
              </div>
              <BenchmarkPreview
                caseId="ogden_bvp"
                plot="deformation"
                alt="Recorded finite-strain deformation of a sheared Ogden block, colored by displacement magnitude"
                eager
              />
              <div className="hero-model-footer">
                <span>Finite shear · actual deformation ×1</span>
                <strong>Follow the model →</strong>
              </div>
            </a>
          </div>
        </section>

        <Livebench />

        <section className="section home-workflow" id="how-it-works">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">A readable path to generated code</p>
                <h2>Keep the physics in view.</h2>
              </div>
              <p>
                The livebench connects each mathematical step to the source and
                numerical result. Start with a worked model, then use the same
                public API for your own.
              </p>
            </div>
            <ol
              className="home-flow"
              aria-label="Equation-to-simulation workflow"
            >
              {[
                [
                  "Equations",
                  "Define the physics, fields, and boundary conditions.",
                ],
                ["Weak form", "Write the test functions and residuals."],
                ["Python", "Declare a supported material or coupled problem."],
                ["Fortran", "Generate and compile the Abaqus subroutine."],
                ["Simulation", "Inspect mesh results and independent checks."],
              ].map(([title, detail], index) => (
                <li key={title}>
                  <span className="step-number">0{index + 1}</span>
                  <h3>{title}</h3>
                  <p>{detail}</p>
                </li>
              ))}
            </ol>
          </div>
        </section>

        <section className="section home-scope" id="scope">
          <div className="shell home-scope-grid">
            <div>
              <p className="eyebrow">Two generation targets</p>
              <h2>
                Material laws.
                <br />
                Coupled elements.
              </h2>
              <p>
                The API covers selected finite- and small-strain material laws
                and coupled-field element patterns. Your Abaqus model supplies
                the mesh, analysis steps, loads, and solver settings.
              </p>
              <a className="text-link" href={repoFile("docs/API_USAGE.md")}>
                Check the supported API <Arrow />
              </a>
            </div>
            <div className="home-targets">
              {siteData.scopes.map((scope) => (
                <article key={scope.target}>
                  <span className={`target ${scope.target.toLowerCase()}`}>
                    {scope.target}
                  </span>
                  <h3>{scope.title}</h3>
                  <p>{scope.description}</p>
                  <details>
                    <summary>Supported patterns and scope</summary>
                    <ul>
                      {scope.included.map((item) => (
                        <li key={item}>{item}</li>
                      ))}
                    </ul>
                    <p>{scope.boundary}</p>
                  </details>
                </article>
              ))}
            </div>
          </div>
        </section>

        <section className="section home-paper" id="paper-evidence">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow light">Published research</p>
                <h2>From declarations to deforming structures.</h2>
              </div>
              <p>
                Selected historical Abaqus results from{" "}
                <a href={publication.url}>
                  <em>{publication.title}</em>
                </a>
                . Each package records its source, retained evidence, and
                reproduction status.
              </p>
            </div>
            <div className="research-gallery">
              {siteData.paperEvidence.slice(0, 3).map((paper, index) => {
                const figure = paperImages[index];
                return (
                  <article className="research-card" key={paper.id}>
                    <figure>
                      <a href={repoFile(figure.path)}>
                        <img
                          src={figure.src}
                          alt={figure.alt}
                          loading="lazy"
                          width="640"
                          height="360"
                        />
                      </a>
                      <figcaption>{figure.caption}</figcaption>
                    </figure>
                    <div className="research-copy">
                      <p className="eyebrow light">{paper.eyebrow}</p>
                      <h3>{paper.title}</h3>
                      <p>{paper.summary}</p>
                      <a className="text-link" href={repoFile(paper.path)}>
                        Read the source and result record <Arrow />
                      </a>
                      <details>
                        <summary>Evidence and limits</summary>
                        <p>{paper.result}</p>
                        <p>
                          <strong>Fresh-clone check:</strong> {paper.freshCheck}
                        </p>
                        <p>{paper.boundary}</p>
                      </details>
                    </div>
                  </article>
                );
              })}
            </div>
            <p className="research-index">
              The{" "}
              <a href={repoFile("paper_examples/README.md")}>
                paper-package index <Arrow />
              </a>{" "}
              also includes the{" "}
              <a href={repoFile(siteData.paperEvidence[3].path)}>
                text record for the corrosion comparison
              </a>{" "}
              and its figure-provenance limits.
            </p>
          </div>
        </section>

        <section className="section home-evidence" id="status">
          <div className="shell">
            <div className="section-heading split-heading">
              <div>
                <p className="eyebrow">Inspect the evidence</p>
                <h2>Checks you can follow and repeat.</h2>
              </div>
              <p>
                Model checks, compiled subroutine calls, and Abaqus runs have
                separate records. The livebench links each recorded result to
                its source revision and execution log.
              </p>
            </div>
            <div className="evidence-links">
              <a href={livebenchUrl}>
                <span>01 · Livebench</span>
                <h3>Numerical references</h3>
                <p>
                  Boundary-value results, mesh refinement, equilibrium, and
                  material history.
                </p>
                <strong>Inspect the checks →</strong>
              </a>
              <a href={repoFile("HOWTO_ADD_AN_EXAMPLE.md")}>
                <span>02 · Development</span>
                <h3>A repeatable model pipeline</h3>
                <p>
                  Python verification, deterministic generation, compiled
                  execution, and independent references.
                </p>
                <strong>
                  Read the example contract <Arrow />
                </strong>
              </a>
              <a href={repoFile(siteData.validation.sourcePath)}>
                <span>03 · Abaqus record</span>
                <h3>Solver evidence and scope</h3>
                <p>
                  The {siteData.validation.date} fresh-clone record
                  distinguishes completed analyses from datachecks.
                </p>
                <strong>
                  Read the validation record <Arrow />
                </strong>
              </a>
            </div>
          </div>
        </section>

        <section className="section start-section home-start" id="start">
          <div className="shell start-grid">
            <div>
              <p className="eyebrow light">Start with a working model</p>
              <h2>Generate your first UMAT.</h2>
              <p>
                Clone the package and run the Neo-Hookean example. Python,
                NumPy, and SymPy are required; compiled checks also use
                gfortran, f2py, Meson, and Ninja.
              </p>
              <div className="hero-actions">
                <a
                  className="button primary"
                  href={repoFile("docs/API_USAGE.md")}
                >
                  Open the API guide <Arrow />
                </a>
                <a className="button ghost" href={repository}>
                  Get the source <Arrow />
                </a>
              </div>
            </div>
            <pre className="install">
              <code>
                <span>$</span> git clone {repository}.git{"\n"}
                <span>$</span> cd abaqus_ufl{"\n"}
                <span>$</span> pip install -e <b>".[dev]"</b>
                {"\n"}
                <span>$</span> cd examples/neo_hookean_umat{"\n"}
                <span>$</span> python build.py{"\n"}
                <span>$</span> python check_reference.py{"\n"}
                <span>$</span> python check_compiled.py
              </code>
            </pre>
          </div>
        </section>
      </main>
      <footer>
        <div className="shell footer-grid">
          <div className="footer-brand">
            <span className="brand-mark" aria-hidden="true">
              a/u
            </span>
            <div>
              <strong>abaqus_ufl</strong>
              <p>
                Open research software by {siteData.repository.author}.
                Project-authored source is MIT licensed; third-party materials
                have separate provenance records.
              </p>
            </div>
          </div>
          <nav aria-label="Documentation">
            <strong>Documentation</strong>
            <a href={repoFile("docs/README.md")}>Documentation index</a>
            <a href={repoFile("docs/API_USAGE.md")}>API usage</a>
            <a href={repoFile("docs/theory.md")}>Theory and conventions</a>
            <a href={livebenchUrl}>Livebench</a>
          </nav>
          <nav aria-label="Project records">
            <strong>Project records</strong>
            <a href={publication.url}>Published paper</a>
            <a href={repoFile("CITATION.cff")}>Citation</a>
            <a href={repoFile("CREDITS.md")}>Credits and provenance</a>
            <a href={repoFile("LICENSE")}>MIT license</a>
            <a href={`${import.meta.env.BASE_URL}THIRD_PARTY_NOTICES.txt`}>
              Website third-party notices
            </a>
            <a href={siteData.repository.contact}>Contact Teng Zhang</a>
          </nav>
        </div>
        <div className="shell footer-bottom">
          <span>
            v{siteData.repository.version} · {siteData.repository.author}
          </span>
          <a href={repository}>
            View source on GitHub <Arrow />
          </a>
        </div>
      </footer>
    </>
  );
}

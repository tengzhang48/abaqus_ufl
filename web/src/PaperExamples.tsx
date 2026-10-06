import { useState } from "react";
import siteData from "../site-data.json";
import publication from "../../docs/publication.json";

export type FigureAsset = {
  src: string;
  path: string;
  alt: string;
  caption: string;
};

const repository = siteData.repository.url;
const repoFile = (path: string) => `${repository}/blob/main/${path}`;
const checkDate = new Date(
  `${siteData.validation.date}T00:00:00Z`,
).toLocaleDateString("en-GB", {
  day: "numeric",
  month: "long",
  year: "numeric",
  timeZone: "UTC",
});

function PaperExamples({
  figures,
}: {
  figures: Record<string, FigureAsset[]>;
}) {
  const [selectedId, setSelectedId] = useState(siteData.paperEvidence[0].id);
  const paper = siteData.paperEvidence.find((item) => item.id === selectedId)!;
  const available = figures[paper.id];
  const lead = available[paper.id === "tet4" ? 0 : available.length - 1];
  const remaining = available.filter((figure) => figure !== lead);

  function renderFigure(figure: FigureAsset) {
    return (
      <figure key={figure.path}>
        <a
          href={repoFile(figure.path)}
          aria-label={`Open source figure: ${figure.caption}`}
        >
          <img
            src={figure.src}
            alt={figure.alt}
            width="960"
            height="695"
            loading="lazy"
          />
        </a>
        <figcaption>
          {figure.caption} <a href={repoFile(figure.path)}>source ↗</a>
        </figcaption>
      </figure>
    );
  }

  return (
    <section className="section paper-section" id="paper-evidence">
      <div className="shell">
        <div className="section-heading paper-heading">
          <p className="eyebrow light">
            Published paper · {publication.journal}
          </p>
          <h2>Paper examples</h2>
          <p>
            <a className="text-link" href={publication.url}>
              {publication.title} ↗
            </a>
          </p>
          <p className="paper-validation" id="status">
            Abaqus checks, {checkDate}:{" "}
            {siteData.validation.completeSolves.length} small analyses
            completed; {siteData.validation.datachecks.length} paper inputs
            passed datacheck.{" "}
            <a href={repoFile(siteData.validation.sourcePath)}>
              Validation record ↗
            </a>
          </p>
        </div>

        <div className="paper-explorer">
          <aside className="paper-sidebar" aria-label="Paper example picker">
            <p className="paper-sidebar-label">Choose an example</p>
            <div className="paper-case-nav">
              {siteData.paperEvidence.map((item) => (
                <button
                  type="button"
                  key={item.id}
                  aria-pressed={selectedId === item.id}
                  onClick={() => setSelectedId(item.id)}
                >
                  <span>{item.eyebrow}</span>
                  <strong>{item.title}</strong>
                </button>
              ))}
            </div>
            <div className="paper-mobile-case">
              <label htmlFor="paper-case">Example</label>
              <select
                id="paper-case"
                value={selectedId}
                onChange={(event) => setSelectedId(event.target.value)}
              >
                {siteData.paperEvidence.map((item) => (
                  <option key={item.id} value={item.id}>
                    {item.title}
                  </option>
                ))}
              </select>
            </div>
          </aside>

          <article
            className="paper-case"
            key={paper.id}
            aria-label={paper.title}
          >
            <div className="paper-copy">
              <p className="eyebrow light">{paper.eyebrow}</p>
              <h3>{paper.title}</h3>
              <p className="paper-summary">{paper.summary}</p>
              <div className="result-line">
                <strong>Available result</strong>
                <p>{paper.result}</p>
              </div>
              <details className="paper-checks">
                <summary>Checks and limits</summary>
                <div className="result-line">
                  <strong>Check recorded on {checkDate}</strong>
                  <p>{paper.freshCheck}</p>
                </div>
                <div className="boundary">
                  <strong>Limits</strong>
                  <p>{paper.boundary}</p>
                </div>
              </details>
              <a className="text-link" href={repoFile(paper.path)}>
                Read the package record ↗
              </a>
            </div>

            {lead ? (
              <div className="paper-figures">
                <div className="figure-grid figure-grid-1">
                  {renderFigure(lead)}
                </div>
                {remaining.length > 0 && (
                  <details className="paper-frames">
                    <summary>
                      {paper.id === "tet4" ? "More figures" : "More frames"} (
                      {remaining.length})
                    </summary>
                    <div className="figure-grid figure-grid-2">
                      {remaining.map(renderFigure)}
                    </div>
                  </details>
                )}
              </div>
            ) : (
              <aside
                className="provenance-panel"
                aria-label="Corrosion artifact provenance boundary"
              >
                <div>
                  <p className="eyebrow light">Text record only</p>
                  <h4>Figure not republished on this site</h4>
                  <p>
                    The comparison mesh came from a third-party distribution
                    whose complete BSD license notice has not been located. The
                    corrosion figure is therefore omitted from this site.
                  </p>
                  <a className="text-link" href={repoFile("CREDITS.md")}>
                    Read the provenance record ↗
                  </a>
                </div>
              </aside>
            )}
          </article>
        </div>

        <div className="paper-index-link">
          <p>
            Submitted and current generated versions are retained when they
            differ.
          </p>
          <a className="text-link" href={repoFile("paper_examples/README.md")}>
            All paper packages ↗
          </a>
        </div>
      </div>
    </section>
  );
}

export default PaperExamples;

import BenchmarkPreview from "./BenchmarkPreview";
import { modelCards } from "./livebench-cards";
import { livebenchUrl } from "./SiteHeader";

export default function Livebench() {
  return (
    <section className="section livebench-section home-models" id="livebench">
      <div className="shell">
        <div className="section-heading split-heading">
          <div>
            <p className="eyebrow">The livebench</p>
            <h2>See the equations become a simulation.</h2>
          </div>
          <p>
            Choose a model. Follow its equations, weak form, Python declaration,
            generated Fortran, and checked mesh results on a dedicated page.
          </p>
        </div>
        <div className="model-gallery">
          {modelCards.map((model) => (
            <a
              className="model-card"
              key={model.id}
              href={`${livebenchUrl}?case=${model.id}`}
            >
              <BenchmarkPreview
                caseId={model.id}
                plot={model.plot}
                alt={model.alt}
              />
              <div className="model-card-copy">
                <p className="eyebrow">{model.kind}</p>
                <h3>{model.title}</h3>
                <p>{model.summary}</p>
                <span className="model-caption">{model.caption}</span>
                <span className="text-link">
                  Follow this model <span aria-hidden="true">→</span>
                </span>
              </div>
            </a>
          ))}
        </div>
        <p className="gallery-note">
          Images show recorded numerical results. Each walkthrough includes
          boundary conditions, references, and reproducible checks.{" "}
          <a href={livebenchUrl}>View the livebench index →</a>
        </p>
      </div>
    </section>
  );
}

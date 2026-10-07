import { livebenchUrl } from "./SiteHeader";

export default function Livebench() {
  return (
    <section className="section livebench-section" id="livebench">
      <div className="shell">
        <div className="section-heading split-heading">
          <div>
            <p className="eyebrow">Explore a worked model</p>
            <h2>From equations to simulation</h2>
          </div>
          <p>
            See how a physical model becomes a weak form, a Python declaration,
            and a generated Abaqus user element. Each example opens its own
            walkthrough.
          </p>
        </div>
        <ol
          className="livebench-preview-flow"
          aria-label="Livebench learning path"
        >
          <li>
            <span>01</span> Equations & boundaries
          </li>
          <li>
            <span>02</span> Weak form & Python
          </li>
          <li>
            <span>03</span> Generated code & results
          </li>
        </ol>
        <div className="livebench-preview-cases">
          <a href={`${livebenchUrl}?case=heated_plate_bvp`}>
            <span className="case-kind">Heat transfer · Quad4 UEL</span>
            <h3>A plate heated from one edge</h3>
            <p>
              Follow the heat equation through its weak form to a transient 2D
              temperature field.
            </p>
            <span className="text-link">Start the walkthrough →</span>
          </a>
          <a href={`${livebenchUrl}?case=thermal_bending_bvp`}>
            <span className="case-kind">Thermomechanics · Quad4 UEL</span>
            <h3>A strip that bends as it heats</h3>
            <p>
              Connect heat diffusion, thermal stress, and a clamped boundary to
              the solved deformation.
            </p>
            <span className="text-link">Start the walkthrough →</span>
          </a>
        </div>
        <a className="button livebench-open" href={livebenchUrl}>
          Open the livebench →
        </a>
      </div>
    </section>
  );
}

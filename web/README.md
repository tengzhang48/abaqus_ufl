# Public website

This directory is the maintainable source for the `abaqus_ufl` GitHub Pages
website. The deployed site is a static React/TypeScript build; it does not run
Abaqus or execute generated subroutines in a visitor's browser.

The livebench reads `public/livebench/report.json` and presents recorded
two genuine transient boundary-value simulations, stress/strain histories,
relaxation, thermal-element comparisons, and FE profiles, mesh fields, and
convergence. The initial view is a thermally bending strip with a clamped end
and otherwise traction-free surfaces; a heated plate provides an independent
2D transport BVP. Visitors inspect the boundary conditions and free field DOF
counts, play recorded accepted time steps, select a physical problem
and view, overlay reference curves, inspect samples, and download CSV/SVG
data. Scientific plots use the actual f2py results; browser controls do not
recompute them. Metrics and logs are in each case's details.
Standalone SVGs are rendered by Matplotlib from the recorded arrays and
shipped with the report; mesh SVG downloads explicitly identify the final
frame. Frame CSV and full-history CSV downloads use the same full-precision
time/nodal data. Playback uses fixed geometry and color scales, and boundary
overlays preserve the prescribed edge node identities.

The Python runner and browser share `../tools/livebench_cases.json`.
The parser rejects missing checks, duplicate cases, nonfinite metrics or
plot arrays, inconsistent agreement claims, and invalid mesh connectivity.
Older reports without plots show an explicit data-unavailable message.
The publication gate requires every declared plot for successful cases.
An unavailable report stays unavailable.
Published article metadata comes from `../docs/publication.json` and is checked
against `../CITATION.cff` during the build.

## Evidence model

- `site-data.json` holds the public copy and links used by the interface.
- `scripts/check-evidence.mjs` checks that linked example records and figure
  inputs exist and that headline validation facts remain supported by
  `docs/ABAQUS_VALIDATION_2026-07-30.md`.
- Selected tracked figure assets with documented provenance are imported from the
  `paper_examples/` tree at build time. Each rendered figure links back to its
  exact source path in GitHub. The corrosion figure is intentionally not
  bundled while its derived-mesh notice and redistribution status remain open.
- Manuscript-package limitations remain beside the related visual, including
  the external gel mesh seed and the open corrosion-mesh provenance item.

The website deliberately distinguishes generated-source checks, direct
compiled calls, completed Abaqus solves, datachecks, retained historical
results, and scientific validation. Do not collapse those categories into a
single “validated” label.

## Local development

Use Node 22 or newer:

```bash
cd web
npm ci
npm run dev
```

Run the same source, evidence, type, and build checks used by Pages:

```bash
npm run check
```

To preview real local results:

```bash
cd ..
python tools/run_livebench.py
mkdir -p web/public/livebench
cp -R benchmark-results/. web/public/livebench/
cd web
npm run dev
```

Local runs can include uncommitted changes, which are labeled in the interface.
The Actions publication gate additionally requires a clean checkout and all
all nine cases. The publication gate additionally requires the BVP setup
and accepted time-history arrays. The parser rejects nonfinite, mismatched or
nonmonotonic frames and inconsistent final-frame/mesh/setup data. Generated
reports and logs are ignored by Git; the workflow
delivers them as build artifacts rather than committing changing output.

The production base path is `/abaqus_ufl/`; set `VITE_BASE_PATH=/` for a
root-mounted local production build if needed.

## Deployment

`.github/workflows/pages.yml` first runs the f2py livebench, retains its report
and logs, builds this directory with those results, and publishes `web/dist`
from the upstream `main` branch. Pull requests and forks retain a downloadable
website preview and benchmark artifact without deploying the upstream site.
GitHub Pages must use **GitHub Actions** as its source in the repository settings.

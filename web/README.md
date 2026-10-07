# Public website

This directory is the maintainable source for the `abaqus_ufl` GitHub Pages
website. The deployed site is a static React/TypeScript build; it does not run
Abaqus or execute generated subroutines in a visitor's browser.

The homepage contains a compact introduction linking to the dedicated
`/abaqus_ufl/livebench/` page. Choose a model to follow **Equations → Weak
form → Python → Generated Fortran → Simulation**, one stage at a time.
The four worked boundary-value models are a heated plate, a clamped thermal
bending strip, a finite-shear Ogden block, and a J2 loading/unloading block.
The homepage places four clickable model images directly below a concise
package introduction, with direct documentation and getting-started actions.
Detailed research evidence stays in linked records and expandable summaries.
Source excerpts are imported from the public Python and
Fortran files rather than maintained as separate example implementations.
Seven supporting material, element and mesh checks have separate result pages.

The simulation stage reads `public/livebench/report.json`. Visitors play
recorded accepted steps, overlay references, inspect convergence, and download
CSV/SVG data. FEM views show the field and optional element edges, with no
default node dots. Opening **Inspect nodal values** reveals the nodal data
and a small cross at the selected node. Scientific plots use the actual f2py
results; browser controls do not recompute them. Metrics and logs are in
each case's details.
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
- The homepage loads 960 px WebP copies of three tracked `paper_examples/`
  figures from `src/assets/figures/`. Each rendered figure links back to its
  original PNG in GitHub. The corrosion figure is intentionally not bundled
  while its derived-mesh notice and redistribution status remain open.
- Manuscript-package limitations remain beside the related visual, including
  the external gel mesh seed and the open corrosion-mesh provenance item.

The WebP copies were made once with Pillow, which is not a build dependency.
To regenerate them from the repository root:

```bash
python - <<'PY'
from pathlib import Path
from PIL import Image

output = Path("web/src/assets/figures")
for source in ["paper_examples/stabilized_tet4/figure/Tet4_u_mag.png",
               "paper_examples/morphing_hex8/figure/pasta_t360.png",
               "paper_examples/gel_bilayer/figure/gel_bilayer_6h.png"]:
    with Image.open(source) as image:
        image.thumbnail((960, 960), Image.Resampling.LANCZOS)
        image.save(output / (Path(source).stem + ".webp"),
                   "WEBP", quality=82, method=6)
PY
```

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
The Actions publication gate additionally requires a clean checkout and
all eleven cases. The publication gate additionally requires the BVP setup
and accepted time-history arrays. The parser rejects nonfinite, mismatched or
nonmonotonic frames and inconsistent final-frame/mesh/setup data. Generated
reports and logs are ignored by Git; the workflow
delivers them as build artifacts rather than committing changing output.

The production base path is `/abaqus_ufl/`; set `VITE_BASE_PATH=/` for a
root-mounted local production build if needed.

Vite builds two HTML entries: `index.html` and `livebench/index.html`.
Case and stage links use query parameters, for example
`/abaqus_ufl/livebench/?case=heated_plate_bvp&step=weak-form`. Direct visits,
refresh, new tabs, and browser back/forward work on static GitHub Pages without
a server rewrite. Report files remain alongside the livebench HTML entry.

## Deployment

`.github/workflows/pages.yml` first runs the f2py livebench, retains its report
and logs, builds this directory with those results, and publishes `web/dist`
from the upstream `main` branch. Pull requests and forks retain a downloadable
website preview and benchmark artifact without deploying the upstream site.
GitHub Pages must use **GitHub Actions** as its source in the repository settings.

# Public website

This directory is the maintainable source for the `abaqus_ufl` GitHub Pages
website. The deployed site is a static React/TypeScript build; it does not run
Abaqus or execute generated subroutines in a visitor's browser.

The livebench reads `public/livebench/report.json` and presents recorded
stress/strain histories, relaxation, thermal-element comparisons, and FE
profiles, mesh fields, and convergence. Visitors select a physical problem
and view, overlay reference curves, inspect samples, and download CSV/SVG
data. Scientific plots use the actual f2py results; browser controls do not
recompute them. Each UMAT/UEL case includes its example guide and model/source
links. The case's details retain the example's verification evidence and Abaqus
coverage alongside the numerical metrics and logs.
Standalone SVGs are rendered by Matplotlib from the recorded arrays and
shipped with the report; CSV downloads use the same full-precision data.

The livebench follows the introduction, ahead of the compact three-step
Python → Fortran → Abaqus workflow. Supported formulations, generation limits,
and the UFL relationship are in an expandable panel within that workflow.
The `#scope` link still reaches that panel; `#examples` reaches the livebench.

Paper examples use their own picker, with one lead figure at a time. Remaining
frames and the exact dated checks and limits are expandable. The `#status`
anchor reaches the July 2026 Abaqus record above the picker, linked to the
complete report; historical automated-test counts are not presented as current.

Response plots show reference bands beneath computed markers, with matching
legend swatches and round axis ticks. The continuum backward-Euler guide starts
hidden in the FE profiles and can be enabled in the legend. The Difference view,
full-precision downloads, recorded numerical arrays, and acceptance checks keep
their original data.

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

The browser loads 960 px WebP copies from `src/assets/figures/`, with intrinsic
dimensions declared to reserve layout space. Source links and evidence paths
still point to the original PNGs in `paper_examples/`. The copies were made
once using the existing Pillow installation; Pillow is not a build dependency.
To regenerate them from the repository root:

```bash
python - <<'PY'
import json
from pathlib import Path
from PIL import Image

data = json.loads(Path("web/site-data.json").read_text())
output = Path("web/src/assets/figures")
output.mkdir(parents=True, exist_ok=True)
for case in data["paperEvidence"]:
    for source in case["figurePaths"]:
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

`public/social-preview.png` is the 1200 × 630 link-preview card. It combines
project text with the Tet4 displacement web copy; the source figure remains
`paper_examples/stabilized_tet4/figure/Tet4_u_mag.png`. Canonical and Open Graph
metadata point to the public GitHub Pages URL.

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
seven cases. Generated reports and logs are ignored by Git; the workflow
delivers them as build artifacts rather than committing changing output.

The production base path is `/abaqus_ufl/`; set `VITE_BASE_PATH=/` for a
root-mounted local production build if needed.

## Deployment

`.github/workflows/pages.yml` first runs the f2py livebench, retains its report
and logs, builds this directory with those results, and publishes `web/dist`
from the upstream `main` branch. Pull requests and forks retain a downloadable
website preview and benchmark artifact without deploying the upstream site.
GitHub Pages must use **GitHub Actions** as its source in the repository settings.

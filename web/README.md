# Public website

This directory is the maintainable source for the `abaqus_ufl` GitHub Pages
website. The deployed site is a static React/TypeScript build; it does not run
Abaqus or execute generated subroutines in a visitor's browser.

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

The production base path is `/abaqus_ufl/`; set `VITE_BASE_PATH=/` for a
root-mounted local production build if needed.

## Deployment

`.github/workflows/pages.yml` builds this directory from `main` and publishes
only `web/dist`. GitHub Pages must use **GitHub Actions** as its source in the
repository settings.

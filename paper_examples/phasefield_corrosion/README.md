# Phase-field stress-corrosion cracking (Sec. 3.1, Fig. 3)

Code-to-code benchmark against the reference implementation of Cui, Ma,
and Martinez-Paneda (JMPS 2021).

- `phasefield_corrosion_cui.py` — the declaration and generation
  pipeline (three fields u/phi/c on a reduced-integration Quad8;
  J2 plasticity with damage; kappa_r and the tangent-block suppression
  are declaration constants). Running it verifies the declaration and
  generates the main UEL; `generate_diag()` produces the diagonal
  comparison variant used for Fig. 3.
- `generated/` — the current generated sources for both variants.
- `archived_submitted/` — the exact sources associated with the
  accepted comparison runs, preserved unmodified.
- `abaqus/SCC_abaqus_ufl_Tdummy_diag.inp` — our comparison deck
  (mesh derived from the attributed reference distribution; see the open
  license/provenance boundary below).
- `reference_data/extracted_results_original/` — reduced results
  extracted from the reference UEL's run (pit-depth history and final
  fields), the basis of the Fig. 3 comparison.
- `figure/` — the Fig. 3 script and its image inputs.

## Generation and Abaqus check

Run `python phasefield_corrosion_cui.py` from this directory. It verifies the
material tangents and writes both current sources under `generated/`.

On 2026-07-30, both sources regenerated and passed fixed-form gfortran
compilation. The tracked full deck and freshly generated diagonal UEL also
completed an Abaqus/Standard 2022 datacheck with zero errors. The full
corrosion analysis was not repeated in that fresh-clone validation.

Reference distribution (not redistributed): obtain from the authors and
verify SHA-256:
- `PhaseFieldSCC.f`:
  `555f4434272534ef1aaf96e0eea3cd669ace4d68dd67c39f0be31866c47d0e57`
- `SCC.inp`:
  `d51a1008f5297558550008a490f003131dcd6189a0f39f36d51bb8b0bfc7a56e`

## License and provenance boundary

The reference UEL header says that the code is distributed under a BSD
license. A 2026-08-02 audit of upstream Git commit
`cd1fb320a90ada8ebb7a9437254549a0d181a0e0`, its complete public history, and
the separate Oxford Mechanics of Materials Lab download found no exact BSD
variant, license text, or copyright notice. The hashes above match those
audited distributions, but the missing terms must not be inferred.

The comparison deck retains project changes and attribution, but its mesh has
reference-distribution lineage. The public website therefore gives this case
a text-only record and does not copy its figure into the Pages artifact. Before
republishing the deck or mesh-derived output elsewhere, obtain the exact notice
from the authors or replace the mesh with a project-authored one. See
[`CREDITS.md`](../../CREDITS.md) for the project-wide record.

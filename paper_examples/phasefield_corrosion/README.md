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
  accepted comparison runs (job 52411708 era), preserved unmodified.
- `abaqus/SCC_abaqus_ufl_Tdummy_diag.inp` — our comparison deck
  (mesh derived from the BSD-licensed reference distribution, with
  attribution).
- `reference_data/extracted_results_original/` — reduced results
  extracted from the reference UEL's run (pit-depth history and final
  fields), the basis of the Fig. 3 comparison.
- `figure/` — the Fig. 3 script and its image inputs.

Reference distribution (not redistributed): obtain from the authors and
verify SHA-256:
- `PhaseFieldSCC.f`:
  `555f4434272534ef1aaf96e0eea3cd669ace4d68dd67c39f0be31866c47d0e57`
- `SCC.inp`:
  `d51a1008f5297558550008a490f003131dcd6189a0f39f36d51bb8b0bfc7a56e`

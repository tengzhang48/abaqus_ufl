# Pasta Morphing With the Generated Local-Pressure Hex8 UEL

This folder is a portable copy of the exact Abaqus input files and user
subroutine used for the completed three-dimensional pasta-morphing simulation
in the EML manuscript. It contains no scheduler-specific submission script and
no generated Abaqus result files.

The files were copied without modification from the completed Abaqus/Standard
2022 run on 2026-07-20.

## Files

| File | Purpose |
|---|---|
| `Pasta_W15_G5_H20_B10_L127_T50_M10_V_abaqus_ufl_hex8.inp` | Main Abaqus input deck |
| `Pasta_W15_G5_H20_B10_L127_T50_M10_V_edummy.inp` | C3D8 visualization-element connectivity included by the main deck |
| `ElasticGelProps_mm.inp` | Nine UEL material properties included by the main deck |
| `chester_anand_local_pressure_hex8_pasta.for` | Exact generated UEL and UVARM source submitted with the completed run |

Keep all four files in the same directory. The main deck uses relative
`*Include` references to the visualization connectivity and material-property
files.

The large ODB, solver logs, and cluster submission script are intentionally not
included. Reduced result data and the manuscript figure inputs are retained
under `EML_paper/figures/`.

## Formulation

The real analysis mesh contains 54,000 eight-node U3 elements. Each node has
the active Abaqus degrees of freedom

```text
1, 2, 3, 11
```

which the UEL interprets as

```text
u1, u2, u3, mu
```

The global unknowns are displacement `u` and chemical potential `mu`. One
constant pressure `p` is solved locally in each Hex8 element, stored in
`SVARS(1)`, and statically condensed from the global element equations.

The polymer volume fraction is derived from

```text
Je  = exp(p/K)
phi = phi0 * Je / det(F)
```

The UEL has 32 global element DOFs and uses eight integration points. The deck
declares `Variables=17`: one pressure state plus eight polymer-fraction and
eight pressure diagnostic slots.

## Visualization Fields

The duplicate `elDummy` C3D8 mesh shares the physical nodes with the U3 mesh
and makes the deformed geometry and integration-point diagnostics visible in
the ODB.

| ODB field | Meaning |
|---|---|
| `U` | Nodal displacement in mm |
| `NT11` | Nodal chemical potential `mu`; it is not physical temperature |
| `UVARM1` | Polymer volume fraction `phi` |
| `UVARM2` | Condensed element-local pressure `p` |
| `LE` | Logarithmic strain reported by the duplicate C3D8 mesh |

`UVARM1` and `UVARM2` are transferred from each U3 element to its corresponding
C3D8 visualization element. The submitted source uses
`LOCALP_UVARM_OFFSET=1000000`, matching dummy element labels
`1000001` through `1054000`.

For visualization, display `elDummy` and hide the U3 analysis elements and the
unconnected `extraElement` C3D8T control element. The solved mesh is a
half-thickness model with `ZSYMM` on `z=0`; mirror about the XY plane when a
full-thickness view is required.

## Abaqus Run

First run a serial datacheck from this directory:

```bash
abaqus \
  job=pasta_morphing_hex8_datacheck \
  input=Pasta_W15_G5_H20_B10_L127_T50_M10_V_abaqus_ufl_hex8.inp \
  user=chester_anand_local_pressure_hex8_pasta.for \
  datacheck \
  cpus=1 \
  interactive \
  ask_delete=OFF
```

The packaged files were rechecked with Abaqus/Standard 2022 on 2026-07-28.
The datacheck completed in 54 s with zero errors, 13 `.dat` warnings, and one
`.msg` warning. The warnings include the two intentionally unconnected regions,
the dummy coupled-thermal/contact packaging, and the unsupported request for
element output directly on user elements.

After the datacheck passes, run the analysis:

```bash
abaqus \
  job=pasta_morphing_hex8 \
  input=Pasta_W15_G5_H20_B10_L127_T50_M10_V_abaqus_ufl_hex8.inp \
  user=chester_anand_local_pressure_hex8_pasta.for \
  cpus=32 \
  mp_mode=threads \
  interactive \
  ask_delete=OFF
```

Adjust `cpus` to the available hardware and Abaqus license. Keep
`mp_mode=threads` when `UVARM1` and `UVARM2` are required. The UEL-to-UVARM
transfer uses process-local Fortran module storage; a multi-rank MPI run does
not guarantee that a visualization element can read data written by its
corresponding U3 element on another rank.

This is a large nonlinear run with 245,888 equations. It should normally be
executed through the site's scheduler, but the scheduler wrapper is deliberately
left outside this portable example.

## Completed-Run Record

The retained files produced a successful Abaqus/Standard 2022 analysis with:

```text
Step:                  Swell_Deswell
Step time:             360 s
Parallel configuration: 1 MPI rank x 32 threads
Increments:            110
Automatic cutbacks:   0
Equilibrium iterations: 318
Error messages:       0
Wallclock time:       1762 s
Final max |U|:        12.5329 mm
```

The solver reported 53 analysis warnings, including 52 negative-eigenvalue
warnings and no numerical-problem warnings. The final frame contained finite,
fully populated diagnostics:

```text
NT11:   -0.818389 to approximately 0
UVARM1:  0.167997 to 0.633727
UVARM2: -2424.60 to 175.156
```

## Exact-Source Provenance

SHA-256 checksums of the completed-run inputs are:

```text
f85968a58da1b2a99631376a4276ab6381f5c3fc88e8b15eb8a9f00963fa5b7b  Pasta_W15_G5_H20_B10_L127_T50_M10_V_abaqus_ufl_hex8.inp
e323f3ccddcb5b43a7eef012414c3012cea7d4b34d91afc87330984a9b10ef7b  Pasta_W15_G5_H20_B10_L127_T50_M10_V_edummy.inp
32db8a591a6f4d8af128d3ecc08613d389e51a7ec15b49c1525af3846e1c46f6  ElasticGelProps_mm.inp
e3d788b287717faf2482e64365996ec31d0f6f5efcc17bd80c524f6bb2a8c780  chester_anand_local_pressure_hex8_pasta.for
```

The submitted UEL is intentionally retained rather than replaced by
`examples/gel_chester_anand/chester_anand_local_pressure_hex8.for`. The current
generator source has subsequently changed, including the default UVARM element
offset and shared tensor-helper template. A current-generator rerun must update
the UEL and dummy-element numbering together and must be treated as a separate
validation run.

The inherited heading in the main input deck says "Plane strain swell induced
bending." The retained model itself is three-dimensional; the old heading was
left unchanged to preserve the exact submitted deck.

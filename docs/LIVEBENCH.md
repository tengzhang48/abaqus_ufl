# Livebench

The livebench runs two transient boundary-value simulations, six public
UMAT/UEL verification bundles, and a serial FE verification case through
Python and the actual generated Fortran subroutines
using f2py. It requires no Abaqus installation, license, or commercial solver
service.

## Run locally

From a repository checkout with Python 3.8+, gfortran, Meson, and Ninja:

```bash
pip install -e ".[dev]"
python tools/run_livebench.py
python tools/run_livebench.py --case neo_hookean_umat --output benchmark-results/neo
python tools/run_livebench.py --case heated_plate_bvp --case thermal_bending_bvp
```

Repeat `--case` to select several cases. `--timeout` sets a positive time
limit in seconds per check (default 300); a timeout also stops the check's
compiler processes on Linux. There are no toolchain-related skips: a
missing compiler or build backend makes the compiled gate fail.

## Evidence and results

Each of the six example cases runs its existing `check_reference.py` and
`check_compiled.py`.
The Quad4 and Quad8 UELs also run `check_assembled.py`. These checks retain
their numerical tolerances and deliberately broken controls. The compiled
checks regenerate into temporary directories, compare with the committed
Fortran, compile, call through f2py, and check output and tangent/state
contracts. They do not overwrite the committed generated sources.

The serial verification case runs `tools/check_fe_runtime.py` through the optional
`abaqus_ufl.fe` runtime. It compiles the Quad4 UEL once and reuses the module
for an affine mechanics patch and three refined transient diffusion meshes.
Both tests have independent quantitative closed-form oracles. The case
manifest is `tools/livebench_cases.json`, shared by the runner and website.

The two BVP cases use the same generated Quad4 with genuine boundary
conditions and solve every free field DOF at each of 20 time increments:

- `check_heated_plate.py`: a 24 × 24 plate, a spatially varying hot top
  boundary, and three cold edges; 529 free thermal DOFs. Independent Fourier
  solutions separate space and time error. Mesh studies use 6², 12² and 24²
  elements; time studies use four step sizes on a fixed 12² mesh.
- `check_thermal_bending.py`: a 48 × 8 plane-strain strip clamped at its left
  end, with zero traction elsewhere and opposite thermal baths; 343 free
  thermal and 864 free displacement DOFs. It checks the thermal Fourier
  solution, thermal dissipation, clamp force/moment balance, refinement, and
  a small-strain slender-beam approximation with explicit tolerances.

See [the equations, oracles, and acceptance gates](BOUNDARY_VALUE_SIMULATIONS.md).
CI includes deliberately broken compiled controls for lost nodal history,
delayed thermal boundary values, prescribed strip displacements, and disabled
thermal expansion; each must fail its relevant gate.

CI runs the complete test suite on Python 3.10 and 3.12. The livebench job
also runs the bundles on Python 3.12 to produce the website's numerical
report. Compiled coverage remains independent of the website workflow.

The livebench records 28 plots and fields during these checks:

| Case | Numerical comparisons |
| --- | --- |
| Neo-Hookean | Constrained stretch, shear stress, and normal stress in shear versus closed forms |
| Ogden | Isochoric stretch versus direct principal-stretch powers, including repeated eigenvalues |
| J2 | Shear stress and plastic strain versus the exact consistency solution at every increment |
| Viscoelastic | Stress relaxation versus the exact discrete recurrence, with the continuous solution shown separately |
| Quad4 | Thermal nodal loading and integrated heat storage versus closed forms |
| Quad8 | Deformed thermal-load ratio versus the exact C⁻¹ pull-back factor 1/λ² |
| Serial FE | Three temperature profiles, the temperature field, spatial convergence, and an affine displacement patch |
| Heated plate BVP | Full temperature history, center response, final profile, heat content, and space/time convergence |
| Thermal bending BVP | Full displacement/temperature histories, tip deflection, through-thickness temperature, centerline, and accepted equilibrium residuals |

Optional `record` callbacks capture the existing state histories and mesh
solutions. Additional material/element sweeps reuse the module already
compiled by the check. Exporting plots does not build a second copy of the
Fortran module or repeat the mesh solves. Every response-curve sample is
gated against its independent reference; all numerical data remain finite.
Model-unit labels and setup parameters accompany the plots.

The runner saves:

- `report.json`: every required check's result, returned numerical metrics,
  response curves and nodal mesh fields,
  wall time, source commit and dirty-tree flag, Python/compiler/package
  versions, and GitHub Actions run URL when available;
- `logs/`: the complete output and traceback of each check.
- `figures/`: standalone Matplotlib SVG figures for every view and each
  response-curve discrepancy, including the setup and source revision.

Metric values come from the check functions' return values, not rounded
console output. Checks that return no metrics still retain their assertions
and logs. Timings include generation, compilation, and figure export where applicable and
are diagnostics for the recorded environment, not comparative solver-speed
claims. Exit status is zero only when all requested checks pass.

## GitHub Actions and live page

The **Livebench and website** workflow runs the full suite on
pushes to `main`, pull requests, and manual dispatch. Anyone can run the
benchmarks locally or in their own fork; manual dispatch in this repository
requires repository write access. The downloadable `livebench-results`
artifact includes the report and logs, including a failed run's report.

For the upstream repository, the workflow builds the existing website with
that run's report and publishes it to GitHub Pages on `main`. A numerical
failure remains a failing workflow and can still be published as a failed
result. A setup or website-build failure leaves the previous published
report in place; the displayed timestamp and commit identify that older
run. Pull requests and forks do not deploy the upstream website.

The boundary-value simulations lead the website's problem picker. Each shows
its domain, prescribed/natural boundary conditions, initial state, properties,
mesh, step size, and actual free DOF counts from the solve. Time controls play
or inspect recorded accepted frames. Geometry and color scales remain fixed
during playback; displayed deformation is magnified by ten. Boundary overlays
identify the prescribed edges. Frame CSV and complete-history CSV downloads
retain time and authoritative nodal fields; standalone field SVGs explicitly
show the final frame.

Choose a physical problem, then select its response, field, or convergence
view. Curves overlay the computed samples and the independent reference;
the difference view exposes discrepancies that are hidden by overlapping
lines. Pointer and keyboard controls inspect individual samples or nodes.
CSV downloads retain the original numerical precision. Downloadable SVG
figures are generated with Matplotlib from the same recorded arrays.
The source, metrics, logs, and run environment remain available
in the details.

Plot data are optional additions to the version-1 report format. Older
reports remain readable and explicitly state that response data were not
recorded. The publication gate requires all current plot IDs for successful
cases. Invalid arrays, nonfinite data, inconsistent agreement claims, or
invalid mesh connectivity are rejected. Without an available report the
page shows an unavailable state. Browsers inspect recorded results; they
do not run Fortran. To preview local data, run:

```bash
python tools/run_livebench.py --output web/public/livebench
cd web
npm run dev
```

## Scientific scope

The associated paper is Teng Zhang, *Making coupled-field Abaqus user
elements simple*, **Extreme Mechanics Letters 89** (2026), 102530,
[doi:10.1016/j.eml.2026.102530](https://doi.org/10.1016/j.eml.2026.102530).
Publication metadata is recorded in `publication.json`.

Reference-model agreement and compiled element execution have separate
roles. A passing livebench is not a new Abaqus solve, complete reproduction
of the four paper figures, or physical validation. The BVPs are license-free
global FE simulations; the bending model uses one-way heat-to-mechanics
coupling and a quasi-static plane-strain mechanical solve. See the individual
example records and `paper_examples/README.md` for those evidence levels.

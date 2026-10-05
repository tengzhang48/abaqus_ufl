# Livebench

The livebench runs the six public UMAT/UEL verification bundles and a serial
FE mesh benchmark through Python and the actual generated Fortran subroutines
using f2py. It requires no Abaqus installation, license, or commercial solver
service.

## Run locally

From a repository checkout with Python 3.8+, gfortran, Meson, and Ninja:

```bash
pip install -e ".[dev]"
python tools/run_livebench.py
python tools/run_livebench.py --case neo_hookean_umat --output benchmark-results/neo
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

The seventh case runs `tools/check_fe_runtime.py` through the optional
`abaqus_ufl.fe` runtime. It compiles the Quad4 UEL once and reuses the module
for an affine mechanics patch and three refined transient diffusion meshes.
Both tests have independent quantitative closed-form oracles. The case
manifest is `tools/livebench_cases.json`, shared by the runner and website.

The Python 3.12 CI suite delegates these same compiled bundle and mesh gates
to the livebench job, avoiding a second execution on that Python version.
Python 3.10 retains them for compatibility coverage. Independent generator,
compiler-helper, manifest, and failure-path tests remain in both CI jobs.

The runner saves:

- `report.json`: every required check's result, returned numerical metrics,
  wall time, source commit and dirty-tree flag, Python/compiler/package
  versions, and GitHub Actions run URL when available;
- `logs/`: the complete output and traceback of each check.

Metric values come from the check functions' return values, not rounded
console output. Checks that return no metrics still retain their assertions
and logs. Timings include generation and compilation where applicable and
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

The page permits filtering results by case type and status and downloading the
JSON and logs. Without an available report it shows an unavailable state,
never an inferred pass. Browsers inspect recorded results; they do not run
Fortran. To inspect a local report on the site, copy its directory into
`web/public/livebench/`, then follow `web/README.md`.

## Scientific scope

The associated paper is Teng Zhang, *Making coupled-field Abaqus user
elements simple*, **Extreme Mechanics Letters 89** (2026), 102530,
[doi:10.1016/j.eml.2026.102530](https://doi.org/10.1016/j.eml.2026.102530).
Publication metadata is recorded in `publication.json`.

Reference-model agreement and compiled element execution have separate
roles. A passing livebench is not a new Abaqus solve, complete reproduction
of the four paper figures, or physical validation. See the individual
example records and `paper_examples/README.md` for those evidence levels.

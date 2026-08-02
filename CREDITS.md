# Scientific lineage, reference implementations, and credits

`abaqus_ufl` is authored and maintained by Teng Zhang. Its development drew
on published mechanics, shared research codes, teaching software, benchmark
data, and discussions with other researchers. This record distinguishes those
roles; it does not use “contributor” as a blanket term for people whose theory,
code, data, or advice entered the project in different ways.

## Interface inspiration

The package's declaration style is conceptually inspired by the
[Unified Form Language (UFL)](https://doi.org/10.1145/2566630) and the wider
FEniCS ecosystem. `abaqus_ufl` neither depends on nor implements UFL, and its
Python declarations are not source-compatible with FEniCS UFL.

## Shared and published user-subroutine references

### Chester–Di Leo–Anand gel implementation

The coupled gel theory and original supplemental Abaqus UEL/decks developed by
Shawn A. Chester, Claudio V. Di Leo, and Lallit Anand were important
scientific and finite-element references:

> Shawn A. Chester, Claudio V. Di Leo, and Lallit Anand, “A finite element
> implementation of a coupled diffusion-deformation theory for elastomeric
> gels,” *International Journal of Solids and Structures* 52 (2015), 1–18.
> [doi:10.1016/j.ijsolstr.2014.08.015](https://doi.org/10.1016/j.ijsolstr.2014.08.015)

The retained original package did not include an explicit software license.
The original UEL source is therefore not redistributed here.

The mixed-order bilayer example includes our deck for the swell-induced
bending problem of Chester, Di Leo, and Anand. The deck was written by this
project; its mesh discretization follows their supplemental example, with
attribution. Their original supplemental files are not redistributed.
Regenerating the accepted deck with `build.py` requires a separately obtained
copy of their supplemental input as the mesh seed. This provenance is recorded
in
[`paper_examples/gel_bilayer/README.md`](paper_examples/gel_bilayer/README.md)
and in the generated deck's internal header.

The project also studied the STATEV, initialization, and time-stepping
patterns in Shawn Chester's historically shared `umat_simpleFeFp.for`. The
retained project record does not contain a source URL, revision, or license,
so this credit does not claim that the UMAT source is part of this repository.

### Datta–Nguyen hydrogel UEL

The modular implementation by Bibekananda Datta and Thao D. Nguyen was studied
as a hydrogel UEL reference, including total-Lagrangian organization,
centroid-based F-bar treatment, and local-solver practices:

> Bibekananda Datta and Thao D. Nguyen, “A finite element model and Abaqus user
> element (UEL) implementation of hydrogel chemo-mechanics,” Zenodo (2025).
> [doi:10.5281/zenodo.15725220](https://doi.org/10.5281/zenodo.15725220)

Their source package states a BSD 3-Clause license and its documentation a
CC BY-NC-SA 4.0 license. Their implementation is not bundled here, and this
project does not claim a wholesale port of their source.

### Cui–Ma–Martínez-Pañeda corrosion UEL

The corrosion example is a code-to-code comparison against the reference
formulation and UEL by Chuanjie Cui, Rujin Ma, and Emilio
Martínez-Pañeda:

> Chuanjie Cui, Rujin Ma, and Emilio Martínez-Pañeda, “A phase field
> formulation for dissolution-driven stress corrosion cracking,” *Journal of
> the Mechanics and Physics of Solids* 147 (2021), 104254.
> [doi:10.1016/j.jmps.2020.104254](https://doi.org/10.1016/j.jmps.2020.104254)

The authors' reference UEL header states that the code is distributed under a
BSD license. The original UEL and deck are not redistributed here. Their exact
identities are recorded by SHA-256 in
[`paper_examples/phasefield_corrosion/README.md`](paper_examples/phasefield_corrosion/README.md).
The public comparison deck derives its mesh from that distribution and retains
attribution. The source repository is
[`ChuanjieCui/Phase-field-modelling-of-corrosion`](https://github.com/ChuanjieCui/Phase-field-modelling-of-corrosion).

An audit on 2026-08-02 checked upstream Git commit
`cd1fb320a90ada8ebb7a9437254549a0d181a0e0`, its complete public history, and
the separate Oxford Mechanics of Materials Lab download. The recorded source
hashes match those distributions, but neither distribution includes the exact
BSD variant, license text, or copyright notice. The project therefore does not
infer missing terms: the public website does not bundle the corrosion image,
and the derived-mesh redistribution question remains open. Before republishing
that deck or its mesh-derived output elsewhere, obtain an exact notice from the
authors or replace the mesh with a project-authored one.

The original package's `VisualMesh.m`, which is not included here, credits
E. Martínez-Pañeda and M. Muñiz-Calvente and requests citation of the
Abaqus2Matlab work by G. Papazafeiropoulos, M. Muñiz-Calvente, and
E. Martínez-Pañeda:
[doi:10.1016/j.advengsoft.2017.01.006](https://doi.org/10.1016/j.advengsoft.2017.01.006).

## Solver lineage

The internal `feacheap` research solver is based on Professor Allan Bower's
EN234_FEA teaching code at Brown University. The project added solver back
ends, VTK/output support, build tooling, and an `abaqus_ufl` integration
layer. `feacheap` was used as a non-Abaqus execution host for
Abaqus-compatible subroutines during development.

`feacheap` is not included in this public repository. The retained internal
copy does not record a pinned upstream revision or a complete redistribution
license for EN234_FEA, and its optional UMFPACK wrapper refers to a license
file absent from that tree. Those items must be resolved before any future
redistribution. The MIT license of `abaqus_ufl` must not be interpreted as a
license for Bower's original code or SuiteSparse/UMFPACK.

The retained UMFPACK Fortran wrapper carries the copyright notice of Timothy
A. Davis and interfaces with
[SuiteSparse/UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html).
The absent referenced license file is why neither the wrapper nor `feacheap`
is treated here as redistribution-cleared.

## Published formulations and benchmarks

Example-specific READMEs cite publications that supply formulations, published
problems, or numerical benchmarks. In particular:

- The Tet4 formulation and block-compression benchmark follow Guglielmo
  Scovazzi, Rubén Zorrilla, and Riccardo Rossi,
  [doi:10.1016/j.cma.2023.116076](https://doi.org/10.1016/j.cma.2023.116076).
  This is publication/formulation attribution, not a claim of source-code
  reuse.
- The pressure-based gel formulation, its reparameterization, the grooved-sheet
  geometry, and the solvent-exposure setting follow the morphing model of Ye
  Tao and coauthors,
  [doi:10.1126/sciadv.abf4098](https://doi.org/10.1126/sciadv.abf4098).
  The generated pressure-based UEL is a later implementation and is not the
  original Chester–Di Leo–Anand UEL.

The private research repository contains additional research ports. Each must
retain its own publication, source, version, modification, comparison, and
license record; omission from this public subset does not erase that lineage
or define the full project's capability.

## Discussions and research support

Teng Zhang acknowledges:

- Professor Lallit Anand, Massachusetts Institute of Technology, for many
  inspiring discussions on multiphysics modeling and simulation;
- Professor Allan Bower, Brown University, for his course on computational
  solid and structural mechanics and for sharing the `feacheap`/EN234_FEA
  code; and
- Professor Kenichi Soga and Yaobin Yang, University of California, Berkeley,
  for helpful discussions.

The project also builds on Teng Zhang's earlier gel and morphing research with
Professor Lining Yao's Morphing Matter Lab.

## AI-assisted development

Generative-AI systems—including Anthropic's Claude; Fable; OpenAI's ChatGPT
and Codex; Google's Gemini; Moonshot's Kimi; DeepSeek; Qwen; and Zhipu's
GLM—assisted at recorded stages with discussion of theory, software
development, validation, review, and documentation. These labels record tools
used; they do not transfer scientific authorship or publication responsibility
from Teng Zhang.

## How to add a new reference source

Before a shared UEL, UMAT, deck, mesh, dataset, or derived artifact enters a
public example, record:

1. authors and citation;
2. stable URL, DOI, version, commit, or source hash;
3. the exact role: inspiration, code reuse, adaptation, execution oracle,
   benchmark data, or comparison only;
4. modifications made locally;
5. the license and required notices; and
6. whether the original or a derivative is redistributed.

When any item is unknown, do not add the artifact or broaden its distribution.
Record the provenance gap explicitly rather than inferring permission. If an
artifact is already public, stop copying it to new surfaces and either obtain
the missing terms or replace or quarantine it according to the documented
risk.

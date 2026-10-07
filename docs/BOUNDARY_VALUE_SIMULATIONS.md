# Livebench boundary-value simulations

These four examples solve global finite-element boundary-value problems using
generated public UMAT/UEL models. Python assembles element residuals/tangents,
and SciPy solves the free nodal DOFs. Thermal and Ogden cases execute a
generated Quad4 UEL at every element evaluation; the J2 case calls its generated
UMAT at every integration point in a case-local Quad4 host. Prescribed
boundary values are applied at the end of an increment; accepted old nodal
values remain fixed throughout Newton and line search. No analytic response
is substituted for a solved field.

```bash
pip install -e ".[dev]"
python tools/run_livebench.py --case heated_plate_bvp --case thermal_bending_bvp
python tools/run_livebench.py --case ogden_bvp --case plasticity_bvp
```

All examples use consistent model units. They are license-free verification
simulations, separate from the paper's Abaqus evidence records.

## Heated plate

On `0 ≤ X,Y ≤ 1`, solve `ρcp ∂T/∂t = k ΔT`, initially `T=0`.
The top edge has `T=sin(πX) g(t)`, with `g(t)=1−exp(−2t)`;
the other three edges have `T=0`. Displacements are fixed throughout this
heat-transfer problem. Properties are the public model defaults, including
`k=0.5, ρcp=1`.

The displayed 24 × 24 mesh has 625 nodes and **529 solved interior thermal
DOFs**. Twenty backward-Euler steps of `Δt=0.05` reach `t=1`. Every accepted
nodal field is retained.

An independent continuum solution uses the harmonic lifting
`s(Y)=sinh(πY)/sinh(π)`, with sine coefficients
`cₙ=2(−1)ⁿ⁺¹ n/[π(1+n²)]`. For diffusivity `D=k/ρcp` and
`λₙ=Dπ²(1+n²)`, write

`T=sin(πX)[g(t)s(Y)+Σ bₙ(t) sin(nπY)]`.

The exact mode obeys `bₙ′+λₙ bₙ=−cₙ g′`, starting at zero.
The continuous-time formula is
`bₙ=−cₙ ω[exp(−ωt)−exp(−λₙt)]/(λₙ−ω)`, with its analytic
resonant limit when `λₙ=ω`. A separate continuum backward-Euler recurrence is

`bₙʲ=[bₙʲ⁻¹−cₙ(gʲ−gʲ⁻¹)]/[1+Δt λₙ]`.

Both references use 256 modes. An independent 4096-mode review changed the
recorded reference values by less than `7e-8`.

The publication gates include:

- final maximum spatial error against the continuum BE series ≤ `6e-4`;
- final maximum error against the continuous-time series ≤ `1.5e-3`;
- second-order spatial refinement, measured order between 1.8 and 2.2, on
  6², 12² and 24² meshes at fixed `Δt=0.05`;
- first-order temporal refinement, measured order between 0.8 and 1.2,
  using successive FE solution differences at `Δt=0.2,0.1,0.05,0.025` on
  the same 12² mesh;
- the center temperature, final profile, and integrated stored heat against
  their independent continuum BE series; and
- nonzero interior heating, exact prescribed values, free residual, supplied
  heat balance, and the temperature-norm dissipation identity at every step.

Successive time-step differences cancel the spatial bias. Comparing only
against the continuous-time series can reach the spatial error floor and
misidentify the time integration order.

## Thermal bending strip

On `0 ≤ X ≤ 8, 0 ≤ Y ≤ 1`, solve quasi-static plane-strain mechanics
`Div P=0` and the same heat equation. The left end is fully clamped:
`u₁=u₂=0`. All other mechanical edges have zero traction. The top and bottom
thermal baths are `T=±(1−exp(−4t))`; the left and right end faces are
insulated. Initially `u=T=0`.

The generated constitutive law is

`P=G(F−F⁻ᵀ)+K log(J)F⁻ᵀ−KαTF⁻ᵀ`.

Properties are `G=1,K=2,α=1e-3,k=0.5,ρcp=1`. Here `K` is the coefficient of
the logarithmic volumetric term, and `α` is the thermal-pressure coupling
coefficient; it is not independently labelled the usual thermal-expansion
coefficient. The selected `K/G` avoids a nearly incompressible regime in this
small benchmark.
Twenty `Δt=0.05` steps solve all free thermal and displacement fields together.
The displayed 48 × 8 mesh has 441 nodes, **343 free thermal DOFs**, and
**864 free displacement DOFs**. A 24 × 4 run provides a refinement check.

Heat transfer drives thermal expansion. The current public material has no
mechanical heating or deformation-dependent transport, so this is one-way
thermal coupling, with finite-strain mechanical equilibrium at every step.

For the independent temperature solution, replace the plate lifting by
`s(Y)=2Y/H−1`, coefficients by `cₙ=−2[1+(−1)ⁿ]/(nπ)`, and decay rates by
`λₙ=D(nπ/H)²`, and apply the same exact/BE mode equations. The maximum error
across all retained fine-mesh temperatures must be ≤ `2e-3` and improve by
more than a factor of three under mesh refinement.

The small-strain material limit is
`σ=2G ε+K tr(ε)I−KαTI`. Eliminating transverse stress in plane strain gives
`β=Kα/[2(G+K)]`. A slender, free-bending cross-section has

`κ(t)=−12β/H³ ∫(Y−H/2)T(Y,t) dY`, and `v_tip≈κL²/2`.

The thermal moment is integrated directly from the independent Fourier
coefficients. The beam reference is an **approximation**, because full end
clamping introduces an end effect and the solve uses finite strain. The
final tip must agree within 3%; every plotted transient tip sample must agree
within 8% plus `1e-7` absolute tolerance. Fine-mesh tip agreement must improve
over the coarse mesh. The centerline beam overlay is a guide.

Other gates check actual free DOF counts, boundary values, free residuals,
nonzero bending, positive thermal diffusion dissipation, hot/cold supplied
heat signs, and clamp force/moment balance. Finite-strain moment balance uses
deformed nodal coordinates.

## Balances and deliberately broken controls

The independent Q1 reference-domain heat integral is
`H=ρcp Σ area_e mean(T_e)`. At each accepted step,
`(H_new−H_old)/Δt=Σ constrained R_T`.

Independent Q1 quadrature also checks the BE temperature-norm identity:

`ρcp/(2Δt)(||T_new||²−||T_old||²+||T_new−T_old||²)`
`+k||Grad T_new||²=Σ constrained T_new R_T`.

This is a thermal diffusion identity, not full thermomechanical energy
conservation. The antisymmetric strip has zero net heat content, so a summed
heat balance alone is insufficient; individual bath signs, positive
dissipation, and the quantitative temperature oracle remain required.

[`test_boundary_value_problems.py`](../tests/test_boundary_value_problems.py)
reuses one compiled UEL to prove that losing old nodal temperature, delaying
the heating boundary, prescribing all strip displacements, and disabling
compiled thermal expansion fail their relevant gates. Nonfinite histories
are also rejected.

## Ogden finite-shear block

The unit-square plane-strain block has fixed bottom displacement, top
`u=(0.6s,0)` for quasi-static loading parameter `0≤s≤1`, and traction-free
sides. The public [`OgdenOneTerm`](../examples/ogden_umat/build.py) law is
reused in a displacement-only [`OgdenShearProblem`](../tools/check_ogden_bvp.py)
verification UEL. Its energy and equilibrium are

`W=(2μ/α²)(Σᵢ λ̄ᵢᵅ−3)+(K/2)(ln J)²`,
`λ̄ᵢ=J⁻¹ᐟ³λᵢ`, `P=∂W/∂F`, `Div_X P=0`.

The weak form is `∫Ω₀ Grad_X v:P dΩ₀=0`, with tests vanishing at both grips.
The three-dimensional stretch spectrum includes the plane-strain direction.
Properties `μ=1,α=3.5,K=10` explicitly override the shipped material's
`K=100`; this displacement-Q1 example demonstrates moderate compressibility,
not near-incompressible performance. The original material example ships a
UMAT; the BVP separately generates an eight-DOF Quad4 UEL from the same law.

Twelve increments on a 16 × 16 mesh retain 13 frames, 289 nodes and **510 free
displacement DOFs**. Side freedom produces a non-affine deformation. An
independently written real-SVD principal-stretch stress/energy calculation and
rectangular-Q1 integration audit every free residual and grip reaction.
Centered differences of independently integrated energy with respect to
virtual grip displacement check reaction work. Gates also check positive
`J`, force and current-coordinate moment balance, a closed-form material shear
patch, and the eigenvalue-free `α=2` limit at actual nonisochoric mesh states.

The 4/8/16 refinement sequence must reduce successive grip-reaction changes
by at least a factor of 1.5, with the final change below 3%, and decrease
independent stored energy. The observed final reaction change is about 2.44%.
This is a refinement-change criterion, not an exact continuum error or a claim
that corner stress peaks have converged. The affine material shear patch is
a constitutive oracle, not the exact solution of the free-sided block.

[`test_ogden_bvp.py`](../tests/test_ogden_bvp.py) includes wrong-exponent,
all-prescribed displacement, delayed-grip and nonfinite-history controls.
Rigid-rotation tests require zero force and rotated-tangent objectivity.
Those tests exposed a shared real-value eigensolver defect at roundoff-size
near-isotropic spectra. Python and generated Fortran now preserve the existing
small-scale fallback and also select it using a relative near-isotropic guard.
Direct generated UMAT checks cover zero-stress rigid rotation and its tangent.

## J2 loading and unloading block

The unit-square plane-strain block has fixed bottom displacement, free sides,
and top `u₂=0`, with `u₁=0.02t` on `0≤t≤1` and `u₁=0.02(2−t)` on `1<t≤2`.
This is rate-independent quasi-static pseudo-time. The public
[`SmallStrainJ2`](../examples/small_strain_j2_umat/build.py) UMAT uses
`G=10,λ=20,σy=0.1,H=5`, a radial return and equivalent plastic strain
`STATEV(1)`. The weak residual is `∫Ω ε(v):σ_new dΩ=0`; the case-local
[`J2Quad4Host`](../tools/plasticity_bvp_common.py) integrates `Bᵀσ` and
`BᵀDDSDDE B` at four Gauss points per element.

The host calls the full 3D UMAT (`NDI=3,NSHR=3,NTENS=6`) using engineering
strain ordering `[εxx,εyy,0,γxy,0,0]` and retains `σzz`. The generated UMAT
already converts compression-positive Python tensors to tension-positive
Abaqus arrays. The host selects active tangent indices `[0,1,3]`; it does not
apply another sign flip or a plane-stress condensation.

Twenty increments on the 12 × 12 mesh retain 21 frames, 169 nodes, **286 free
displacement DOFs**, and 576 Gauss-point histories. Every Newton/line-search
trial begins from copies of accepted stress, strain and plastic strain. Only
converged equilibrium commits all states atomically. Rejected trials, failed
solves and UMAT cutbacks cannot change accepted arrays. This example does
not add generic material-state ownership to `abaqus_ufl.fe`.

An independent scalar/Voigt radial-return implementation evolves its own
accepted history and separately assembles nodal forces. Elastic-normal and
proportional-shear closed forms check stress, state and engineering-shear
tangent conventions. The cycle explicitly reaches elastic loading, plastic
flow, elastic unloading and reverse plasticity before returning to zero grip
displacement. It checks yield consistency, irreversible plastic strain,
prescribed-boundary force/reference-coordinate moment balance, and free
equilibrium. The 4/8/12 peak-reaction refinement changes must decrease by more
than a factor of two, with the last change below 2.5% (observed about 1.95%).

The independent exact endpoint-work identity checks each accepted increment:

`Σ_boundary R_new·Δu = ∫Ω [Δ(W_el+H ep²/2)+σy Δep`
`+Δσ:C⁻¹:Δσ/2+H(Δep)²/2] dΩ`,

where `W_el=p²/(2K_bulk)+s:s/(4G)` and `K_bulk=λ+2G/3`.
Physical plastic dissipation is `σy Δep`; the stress-increment and hardening
increment terms are nonnegative endpoint numerical dissipation. The two are
recorded separately. The tolerance is `1e-10` in consistent model units.

[`test_plasticity_bvp.py`](../tests/test_plasticity_bvp.py) checks fixed-history
assembled tangents, rejected-trial/failed-solve/cutback ownership, lost-state,
all-prescribed displacement, shear-convention and cross-tangent controls.
Each deliberately broken implementation must fail its numerical gate.

## Output and website

In thermal cases, recorded temperature is nodal DOF 3; displacement is nodal
DOFs 1 and 2. Ogden and J2 cases have only the two displacement DOFs.
Node identities and Quad4 connectivity are preserved. Displacement magnitude
is computed directly from the solved displacement components. The site
colors the Q1 interpolation of those nodal arrays; it does not rerun the
solver or create additional physical results. Thermal deformation display
uses ×10; Ogden and J2 use ×1. The report declares the display scale,
with fixed geometry and color scales across the history.

Each BVP has a dedicated walkthrough connecting equations, weak form, the
actual Python declaration and generated Fortran to its recorded simulation.
FEM views show the field and element edges without default node markers;
opening nodal inspection adds a small cross at the selected node.
The final stage exposes every accepted time frame, prescribed-edge overlays,
setup, material properties, refinement, reactions/balance diagnostics, and
full-precision frame/history CSV downloads. Standalone Matplotlib field SVGs
show the final frame. The report parser and publication gate reject invalid
time histories, missing setup, bad boundary node identities, nonfinite data,
and inconsistent final-frame fallbacks.

Plastic-strain colors are a visualization-only reference-area-weighted average
of incident element Gauss-weighted means. Raw `STATEV(1)` arrays, Gauss-point
coordinates, weights, element/point identity and every accepted frame accompany
that plot. Quantitative histories use raw points, not projected colors. The
report parser checks the raw coverage and reproduces the nodal projection.
Hysteresis curves preserve chronological sample order explicitly, including
the returning displacement branch. Entrance PNGs and final-frame SVGs are
rendered with Matplotlib from the same recorded arrays, without FEM node dots.

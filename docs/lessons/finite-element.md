# Lessons: Finite-Element Implementation

A collection of implementation traps encountered when writing finite-element
kernels — particularly finite-strain user elements and user materials.
Each entry states the trap, why it happens, the fix, and a generic snippet.

---

### Initialize F to identity, not zeros (plane strain)

**Trap.** In 2D plane strain the deformation gradient `F` is still a 3×3
tensor with `F₃₃ = 1`. If you allocate `F` as zeros and fill only the 2×2
in-plane block, then `F₃₃ = 0`, so `det(F) = 0` and `inv(F)` returns NaN.

**Fix.** Always initialize `F` (and `F_old`) to the 3×3 identity before
accumulating the in-plane displacement gradient:

```fortran
      CALL eye33d(F)
      DO a = 1, nNode
        DO i = 1, ndim    ! ndim = 2
          DO j = 1, ndim
            F(i,j) = F(i,j) + u_node(i,a)*dsh(a,j)
          END DO
        END DO
      END DO
```

The same applies to any tensor that must carry an out-of-plane "1" in
reduced dimensions.

---

### One Jacobian for all fields (mixed-degree elements)

**Trap.** In a mixed element — e.g. quadratic (8-node) geometry with a
linear (4-node) pressure field — the lower-degree shape-function
derivatives must be mapped with the *same* isoparametric Jacobian as the
element geometry. Computing a separate Jacobian from only the corner nodes
(`coords_corner`) gives wrong physical gradients on curved elements.

**Fix.** Compute the Jacobian once from the full element geometry, then
reuse its inverse for every field:

```fortran
C     Geometry Jacobian from the full element (once)
      CALL map_grad_2d(dshxi8, coords, 8, dsh8, detJ, Jinv, stat)
C     Apply the same Jacobian to lower-degree fields
      CALL apply_jinv(dshxi4, Jinv, 4, dsh4)
```

---

### Element-local DOF indexing: 0-based front-end vs 1-based arrays

**Trap.** The element residual and stiffness are dimensioned by the number
of element DOFs (in an Abaqus UEL, `RHS` and `AMATRX` of size `NDOFEL`,
1-based). A Python or symbolic front-end typically hands you 0-based DOF
maps, and global DOF labels are a different numbering entirely. Indexing
element arrays with 0-based or global indices crashes or corrupts assembly.

**Fix.** Build element-local DOF (`edof`) arrays during DOF parsing; the
running counter already gives the correct 1-based, element-local position:

```fortran
      idx = 0
      DO a = 1, NNODE
        DO i = 1, ndim
          idx = idx + 1
          u_node(i, a) = U(idx)
          edof_u(i, a) = idx    ! 1-based element-local index
        END DO
      END DO
```

---

### F-bar: contract the pressure-like stress over in-plane indices only

**Trap.** In the F-bar tangent correction, the pressure-like stress
`Q_iJ = Σ A_iJmN F̄_mN` must sum over in-plane indices only
(`m, N = 1..ndim`), **not** `1..3`. In plane strain `F̄₃₃ = 1` is a
constant, so `dF̄₃₃/dα = 0` and the `(3,3)` component must be excluded from
the chain-rule correction.

**Fix.** All F-bar correction loops run `1, ndim`, never `1, 3`.

---

### F-bar: the volumetric exponent is dimension-dependent

**Trap.** The F-bar scaling is `α = (J̄/J)^(1/ndim)`. Hardcoding `0.5`
happens to be right in 2D but is wrong in 3D.

**Fix.** Compute the exponent from the dimension:

```fortran
      alpha = (Jbar/J) ** (1.0d0/DBLE(ndim))
```

---

### The push-forward of dP/dF is not the Truesdell tangent (Jaumann-tangent trap)

**Trap.** A routine that builds a spatial elasticity tensor by pushing
`∂P/∂F` forward through `F` on two indices,

```
AFF_ijkl = (∂P_iJ / ∂F_kL) · F_jJ · F_lL
```

*looks* like it computes the Truesdell tangent `c` divided by `J`. It does
not. Because `P = F · S`, the chain rule gives

```
∂P_iJ / ∂F_kL = δ_ik · S_LJ + F_iI · ∂S_IJ / ∂F_kL
```

and after contracting with `F_jJ F_lL` the first piece becomes
`δ_ik · τ_jl`. Dividing by `J`:

```
AFF_ijkl / J  =  c^Truesdell_ijkl  +  δ_ik · σ_jl       (★)
```

The extra `δ_ik σ_jl` is **geometric, not material** — it is present for
every constitutive law built through this path. If you assume
`AFF/J = c` and apply the standard Bonet Jaumann correction
`+ (1/2)(σ ⊗̄ I + I ⊗̄ σ)` directly to `AFF/J`, the resulting tangent is
wrong by exactly `δ_ik σ_jl`.

**Why it hides.** The error is linear in stress: it is machine-zero at
`σ = 0`, so it is invisible to `F = I` and stress-free rotation tests; it
is only ~1e-3 in a uniaxial elastic bar test (masked by the Newton
tolerance); but it reaches ~40% in any finite-stress multi-axial
deformation. A tangent with this bug can converge in one iteration per step
on a uniaxial bar test indefinitely — a textbook false positive.

**Fix.** Make identity (★) explicit in the docstring of any push-forward
routine that starts from `∂P/∂F`. The correct consistent tangent for the
Jaumann rate of Kirchhoff stress is:

```
DDSDDE_ijkl = AFF_ijkl / J
            - δ_ik · σ_jl                                    (bookkeeping)
            + (1/2)(σ_ik δ_jl + σ_jk δ_il
                  + σ_il δ_jk + σ_jl δ_ik)                   (Jaumann)
```

The bookkeeping line is **non-optional** — it is the only difference
between "works for everything we tested" and "actually correct."

**How to test it.** Three tests, in order, each cheap. Neither of the first
two is sufficient on its own; the third is mandatory before claiming a
finite-strain tangent works:

1. **F = I, stress-free.** Confirms the tensor assembly itself is correct.
   Catches typos. Does not exercise the Jaumann correction.
2. **Stress-free rotation (F = R, σ = 0).** Confirms objectivity of the
   stress update. Does not exercise the Jaumann correction.
3. **Finite stress, multi-axial, non-aligned D.** A simple-shear
   neo-Hookean at `γ ≳ 0.3`, or a biaxial state with `σ_11 ≠ σ_22` and a
   stretching `D` whose eigenvectors are not aligned with `σ`. With the
   correct tangent, Newton converges quadratically (1–2 iter/step to
   1e-10). With the buggy tangent, expect linear convergence or failure on
   the same input.

A uniaxial bar test is **not discriminating**: it has `tr(D) ≈ 0`
(quasi-incompressibility) and small `|σ|`, which puts the bug in the null
space of the test.

**Numerical reference (NumPy reproducer).** Independent of any FE run, the
formula can be unit-tested in ~50 lines by building the rate-defined
reference and comparing it to what the tangent predicts:

```python
# Compute reference: tau^J / J directly from the rate definition
tau_dot     = AFF : L + tau · L^T          # material rate of Kirchhoff
tau_jaumann = tau_dot - W·tau + tau·W      # Jaumann rate of Kirchhoff
ref         = tau_jaumann / J              # what DDSDDE must produce

# Compare to the generated tangent contracted with D
err = DDSDDE_generated : D - ref
assert np.linalg.norm(err) / np.linalg.norm(ref) < 1e-9
```

Wire this into CI for any code generator that emits finite-strain tangents.
It catches the whole class of "tangent looks right at F=I but is silently
wrong at finite stress" bugs in milliseconds.

---

### Module-scope scratch arrays are a thread-safety hazard under OpenMP

**Trap.** Element-utility modules often declare per-element scratch arrays —
shape functions, spatial derivatives, volume-averaged derivatives,
integration points and weights, the element Jacobian — as **module-scope
`SAVE` variables**. This is fine while assembly is serial: one element is
in flight at a time.

The moment the element loop is parallelized with OpenMP, every thread
writes to the *same* shared memory. Threads stomp on each other's
intermediate results, producing garbage. A common failure mode is a run
that **exits cleanly, with exit code 0, without ever reaching the
analysis** — no error message, just a silent stop after the mesh summary.

**How it presents.** The bug is triggered only above whatever element-count
threshold turns the parallel loop on, so small meshes pass and larger ones
fail silently. Crucially, a **debug build often omits `-fopenmp`** so that
bounds-checking can run without races — which means the loop goes serial
and the bug is *invisible in debug*. Only the optimized, OpenMP-enabled
build hits it.

**Fix.** Treat every module-scope scratch variable used by an element
routine as a thread-safety hazard the instant any caller goes parallel.
Either pass scratch through the call signature (slower, always safe) or mark
it `THREADPRIVATE` so each thread gets a private copy:

```fortran
      real (prec), save :: shape_functions_3D(20)
      !$OMP THREADPRIVATE(shape_functions_3D)
```

One directive per `SAVE` variable (or per group). Keep each directive under
the 132-character free-form line limit — split rather than concatenate.

**Two supporting habits.**

- **Do not trust a green debug build to validate parallel correctness.** If
  the debug configuration drops `-fopenmp`, OpenMP-only bugs cannot appear
  there. Sanity-check the *release* build with `OMP_NUM_THREADS=1` to
  separate parallel races from algorithmic errors.
- **Do not stop at the first plausible compiler warning.** A concurrent
  `-fcheck=all` run may surface a real but unrelated defect (e.g. a
  string-length mismatch on an assumed-length dummy argument). Fixing it
  clears the warning without touching the silent exit. The discriminator
  here: the silent exit reproduces under a *release* build that never prints
  the warning — a sign that the warning is a bystander, not the root cause.

Resolving the race also removes the artificial ceiling on problem size:
once the parallel loop is correct, mesh refinements that were previously
unreachable run, and performance comparisons become meaningful at scale.

# Symbolic Tangent Design

## Current Scope

This document describes two related but distinct ideas:

1. **Current implementation:** whole-method symbolic tangent generation with
   complex-step verification and fallback. This is implemented in
   `abaqus_ufl/core/symbolic_tangent.py` as `SymbolicTangent`.
2. **Future research direction:** true subexpression-level hybrid
   symbolic/complex-step differentiation, where only a hard node such as an
   eigensolve or return map is differentiated by local complex-step and then
   chained analytically. This is not implemented yet.

Complex-step remains the safe default. Symbolic tangents are a performance
opt-in for closed-form materials where SymPy succeeds and the generated
tangent verifies against the complex-step oracle.

## Motivation

Complex-step (CS) tangent computation costs ~10× more than an analytical
tangent: for a 3×3 deformation gradient, 9 perturbation directions ×
~3× cost per complex evaluation = ~27× overhead vs. a single real
evaluation. For many models (NeoHookean, gel, Mooney-Rivlin), every
term in dP/dF has a known closed-form derivative — the 10× CS overhead
is entirely unnecessary. For models with return mapping or
eigendecomposition, some terms require CS, but not all.

The implemented symbolic approach differentiates a whole material method with
SymPy and emits pure real Fortran for the tangent. The future hybrid approach
would decompose the tangent into individual terms, compute each term using the
cheapest method, and use CS only at hard nodes.

## Relationship to Existing Frameworks

This is the same computation-graph approach used in FEniCSx (UFL/FFCx)
and AceGen (Korelc & Wriggers, "Automation of Finite Element Methods,"
Springer, 2016):

| Framework | Graph approach | At hard nodes |
|-----------|---------------|---------------|
| FEniCSx | Full symbolic AD | Fails (forbidden) |
| AceGen | Forward-mode AD | Manual special-casing |
| **abaqus_ufl current** | **Whole-method symbolic where possible** | **Full CS fallback** |
| **abaqus_ufl future** | **Subexpression symbolic + local CS** | **Local CS fallback** |

FEniCSx walks the expression DAG and applies symbolic AD at every node.
When it hits eigendecomposition, d(λᵢ)/dC contains 1/(λᵢ−λⱼ), which
is singular at repeated eigenvalues. FEniCSx has no fallback — the
model simply cannot be expressed. AceGen handles hard nodes through
manual implementation (Daleckii-Krein formula for eigenvalue
derivatives), requiring case-by-case treatment.

Current contribution: the generator can classify closed-form methods as
symbolic, emit pure real Fortran, and verify against CS. Future contribution:
local-CS hard-node isolation.

## Future Core Idea: Per-Node Classification

The following per-node examples describe the future hybrid research path, not
the current `SymbolicTangent` implementation.

Example — gel stress:

```
P = μ₀F + (λ₀ ln(Jₑ) - μ₀) F⁻ᵀ + mixing(φ) · J · F⁻ᵀ

dP/dF = μ₀ I₄                           ← trivial
      + λ₀ (F⁻ᵀ ⊗ F⁻ᵀ)                ← outer product, symbolic
      + (λ₀ lnJₑ - μ₀) d(F⁻ᵀ)/dF       ← Layer 2 identity
      + d(mixing)/dF · J · F⁻ᵀ + ...    ← chain rule, all symbolic
```

Every node is symbolic → 0 CS perturbations needed.

Example — J2 plasticity:

```
F → Fe_trial → σ_trial → [return mapping: Δγ] → σ → P
      ↑            ↑             ↑                ↑     ↑
    known       known      HARD NODE           known  known

dP/dF = dP/dσ · dσ/dΔγ · dΔγ/dσ_trial · dσ_trial/dFe · dFe/dF
         sym      sym      CS (6 perturbs)    sym          sym
```

Only the return mapping needs CS. And σ_trial is symmetric, so only
6 perturbations (not 9 for full F). Each perturbation evaluates only
the return mapping, not the full material function.

Example — Ogden (eigendecomposition):

```
F → C → [eig(C): λᵢ, Nᵢ] → ψ(λᵢ) → S → P
                  ↑
              HARD NODE

dP/dF = dP/dS · dS/dψ · dψ/dλ · dλ/dC · dC/dF
         sym     sym     sym    CS (6)    sym: 2F·δ
```

Only eigendecomposition needs CS. 6 perturbations of symmetric C.

### Cost comparison

| Model | Current CS | Symbolic or future hybrid | Savings |
|-------|-----------|--------|---------|
| NeoHookean | 9 complex mat calls | 0 (all symbolic) | 10× |
| Gel 3-field | 13 complex mat calls | 0 (all symbolic) | 10× |
| Mooney-Rivlin | 9 complex mat calls | 0 (all symbolic) | 10× |
| J2 plasticity | 9 complex mat calls | 6 complex return-map calls | ~4× |
| Ogden | 9 complex mat calls | 6 complex eig calls | ~4× |
| Hencky | 9 complex mat calls | 6 complex eig calls | ~4× |

## Current Detection: Run with SymPy, See Whether the Method Works

No AST parsing needed. Feed SymPy symbolic matrices into the user's
material function. If the function completes, the whole method is symbolic.
If SymPy raises TypeError, the current implementation falls back to CS for that
method. Hard-node identification is future work.

```python
F_sym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'F_{i+1}{j+1}'))
try:
    P_sym = material.stress_PK1(F_sym)
    mode = 'symbolic'       # whole function is symbolic
except TypeError as e:
    mode = 'cs'             # current behavior
    # Future: identify a hard node and apply local CS there.
```

This works because:
- SymPy matrices support +, -, *, /, det, inv, transpose, log, exp
- SymPy CANNOT do: eig() on symbolic matrix, for loops with convergence,
  if-branches that compare symbolic values
- Future work may use the failure point as a hard-node hint.

### User override

```python
class NeoHookean(au.Material):
    props = dict(mu=1.0, lam=10.0)
    tangent = 'symbolic'   # 'symbolic', 'cs', or 'auto' (default)

    def stress_PK1(self, F):
        ...
```

`'auto'` runs the SymPy probe. `'symbolic'` forces symbolic (fails
with clear error if SymPy can't handle it). `'cs'` forces full
complex-step (current behavior, always correct).

## Architecture: Three Layers

### Layer 1: Tensor Primitives (existing)

Already in `tensor_ops.for` — both real and complex versions:

```
det33, inv33, transpose33, matmul33, trace33,
eye33, sym33, dev33, outer33, contract44_33
```

No changes needed. These are called by both symbolic and CS paths.

### Layer 2: Derivative Identity Library (new)

Pre-coded, verified Fortran subroutines for standard tensor calculus
results. Each is 5–15 lines. Written once, tested against CS, reused
by every model that needs them.

```fortran
C     ================================================================
C     Layer 2: Derivative identities — DOUBLE PRECISION only
C     (CS path uses its own complex perturbation; these are for
C      the symbolic tangent path which is pure real)
C     ================================================================

C     dJ/dF: derivative of det(F) w.r.t. F
C     Result: dJdF(i,J) = J * Finv(J,i)
      SUBROUTINE ddetdF33(F, J, dJdF)
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3), J
      DOUBLE PRECISION, INTENT(OUT) :: dJdF(3,3)
      DOUBLE PRECISION :: Finv(3,3)
      CALL inv33(F, Finv)
      DO i = 1, 3
        DO J1 = 1, 3
          dJdF(i,J1) = J * Finv(J1,i)
        END DO
      END DO
      END SUBROUTINE

C     d(F^{-T})/dF: 4th-order tensor
C     Result: T4(i,J,k,L) = -Finv(J,k) * Finv(L,i)
      SUBROUTINE dFinvTdF(Finv, T4)
      DOUBLE PRECISION, INTENT(IN)  :: Finv(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: T4(3,3,3,3)
      DO i = 1, 3
        DO J1 = 1, 3
          DO k = 1, 3
            DO L = 1, 3
              T4(i,J1,k,L) = -Finv(J1,k) * Finv(L,i)
            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE

C     d(lnJ)/dF: derivative of log(det(F)) w.r.t. F
C     Result: dlnJdF(i,J) = Finv(J,i)  (= F^{-T})
      SUBROUTINE dlndetdF33(F, dlnJdF)
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: dlnJdF(3,3)
      DOUBLE PRECISION :: Finv(3,3)
      CALL inv33(F, Finv)
      DO i = 1, 3
        DO J1 = 1, 3
          dlnJdF(i,J1) = Finv(J1,i)
        END DO
      END DO
      END SUBROUTINE

C     d(C^{-1})/dC for symmetric C:
C     Result: T4(I,J,K,L) = -0.5*(Cinv(I,K)*Cinv(J,L)
C                                 +Cinv(I,L)*Cinv(J,K))
      SUBROUTINE dCinvdC(Cinv, T4)
      DOUBLE PRECISION, INTENT(IN)  :: Cinv(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: T4(3,3,3,3)
      DO I1 = 1, 3
        DO J1 = 1, 3
          DO K = 1, 3
            DO L = 1, 3
              T4(I1,J1,K,L) = -0.5d0*(Cinv(I1,K)*Cinv(J1,L)
     &                               +Cinv(I1,L)*Cinv(J1,K))
            END DO
          END DO
        END DO
      END DO
      END SUBROUTINE

C     d(I1)/dF where I1 = tr(F^T F) = F:F
C     Result: dI1dF(i,J) = 2 * F(i,J)
      SUBROUTINE dI1dF(F, result)
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: result(3,3)
      result = 2.0d0 * F
      END SUBROUTINE

C     d(I2)/dF where I2 = 0.5*(I1^2 - tr(C^2))
C     Result: dI2dF(i,J) = 2*(I1*F(i,J) - (F C)(i,J))
C     (For Mooney-Rivlin and similar)
      SUBROUTINE dI2dF(F, result)
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: result(3,3)
      DOUBLE PRECISION :: C(3,3), FC(3,3), I1
      CALL matmul33T1(F, F, C)       ! C = F^T F
      CALL matmul33(F, C, FC)        ! FC = F C
      I1 = C(1,1) + C(2,2) + C(3,3)
      result = 2.0d0 * (I1 * F - FC)
      END SUBROUTINE
```

**Full Layer 2 catalog** (to be implemented incrementally):

| Subroutine | Identity | Used by |
|------------|----------|---------|
| `ddetdF33` | dJ/dF = J F^{-T} | All models |
| `dFinvTdF` | d(F^{-T})/dF_{kL} | NeoHookean, any with F^{-T} |
| `dlndetdF33` | d(ln J)/dF = F^{-T} | NeoHookean, compressible |
| `dCinvdC` | d(C^{-1})/dC | Flux tangents |
| `dI1dF` | d(tr C)/dF = 2F | Mooney-Rivlin, Yeoh |
| `dI2dF` | d(I2)/dF | Mooney-Rivlin |
| `dPdF_NeoHookean` | Full 4th-order tangent | NeoHookean (composed) |
| `dPdF_MooneyRivlin` | Full 4th-order tangent | Mooney-Rivlin |
| `pushforward_4` | dP/dF → c (spatial) | UMAT push-forward |
| `jaumann_corr` | σ⊗I + I⊗σ correction | UMAT Jaumann |

The idea: each model-specific tangent (like `dPdF_NeoHookean`) is a
short subroutine that calls the generic identities. The generator can
either emit a call to an existing library tangent, or compose a custom
one from the building blocks.

### Layer 3: Model-Specific Generated Code (thin)

For a NeoHookean UMAT, the symbolic tangent generator produces:

```fortran
      SUBROUTINE neohookean_tangent(F, props, P, dPdF)
C     Generated by abaqus_ufl (symbolic mode)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3), props(2)
      DOUBLE PRECISION, INTENT(OUT) :: P(3,3), dPdF(3,3,3,3)
      DOUBLE PRECISION :: mu, lam, J, Finv(3,3), FinvT(3,3)
      DOUBLE PRECISION :: dlnJdF(3,3), dFinvT(3,3,3,3)

      mu  = props(1)
      lam = props(2)

C     Kinematics
      CALL det33(F, J)
      CALL inv33(F, Finv)
      CALL transpose33(Finv, FinvT)

C     PK1 stress: P = mu*F + (lam*ln(J) - mu)*F^{-T}
      P = mu * F + (lam * LOG(J) - mu) * FinvT

C     Tangent: dP/dF
C     = mu * I_4 + lam * F^{-T} ⊗ F^{-T}
C       + (lam*ln(J) - mu) * d(F^{-T})/dF
      CALL dlndetdF33(F, dlnJdF)       ! = F^{-T}
      CALL dFinvTdF(Finv, dFinvT)       ! 4th-order

      DO i = 1, 3
        DO J1 = 1, 3
          DO k = 1, 3
            DO L = 1, 3
C             mu * delta_ik * delta_JL  (identity 4th-order)
              dPdF(i,J1,k,L) = 0.0d0
              IF (i.EQ.k .AND. J1.EQ.L) THEN
                dPdF(i,J1,k,L) = mu
              END IF
C             + lam * FinvT_iJ * FinvT_kL  (outer product)
              dPdF(i,J1,k,L) = dPdF(i,J1,k,L)
     &          + lam * FinvT(i,J1) * FinvT(k,L)
C             + (lam*lnJ - mu) * dFinvT_iJ_kL
              dPdF(i,J1,k,L) = dPdF(i,J1,k,L)
     &          + (lam*LOG(J) - mu) * dFinvT(i,J1,k,L)
            END DO
          END DO
        END DO
      END DO

      END SUBROUTINE
```

This is ~40 lines of readable Fortran with no complex arithmetic. Compare
to the CS version which would need 9 complex material evaluations.

## Generator Pipeline

```
User Python Material class
         │
         ▼
    Run stress_PK1(F_sym) with SymPy symbolic F
         │
         ├── succeeds ──► mode = 'symbolic'
         │                  │
         │                  ▼
         │              SymPy diff(P, F) for all tangent blocks
         │                  │
         │                  ▼
         │              CSE + fcode → pure real Fortran
         │              (Layer 2 calls where patterns match)
         │
         ├── TypeError at eig/logm ──► future mode = 'hybrid'
         │                  │
         │                  ▼
         │              Symbolic chain rule for all nodes
         │              EXCEPT the hard node
         │                  │
         │                  ▼
         │              Local CS loop around hard node only
         │              (perturb its inputs, not full F)
         │                  │
         │                  ▼
         │              Emit: symbolic Fortran + local CS block
         │              + chain rule assembly
         │
         ├── TypeError at for/if ──► mode = 'cs'
         │                  │
         │                  ▼
         │              Full CS (current behavior, unchanged)
         │
         ▼
    .for output file
    (always verified against full CS oracle)
```

### Per-node vs per-perturbation-direction for UEL

The UEL has multiple tangent blocks. When CS perturbs one input
(e.g., F₁₁), ALL outputs come out simultaneously — dP/dF₁₁,
drp/dF₁₁, dflux/dF₁₁, dstorage/dF₁₁. So if ANY output needs CS
w.r.t. F, you do the F perturbations and get all F-derivative blocks
for free.

The savings come from perturbation directions where ALL dependent
outputs can be computed symbolically. For the gel model, every block
is closed-form → 0 perturbation directions needed. For a hypothetical
model with closed-form stress but iterative concentration solver,
the F perturbations might still be needed for dstorage/dF, but
dflux/dgrad_mu (2 directions) could be fully symbolic.

## SymPy Integration

### What SymPy does

1. User's `stress_PK1` is extracted as a symbolic expression
2. SymPy differentiates: `dP/dF`, and for UMAT, push-forward derivatives
3. SymPy's `cse()` simplifies the expression (common subexpression elimination)
4. Pattern matching recognizes known identities → library calls
5. Remaining terms → direct Fortran arithmetic via `fcode()`

### What SymPy does NOT do

- No AST translation of the full material function (that's the existing
  CS path's job)
- No runtime evaluation — SymPy runs at code generation time only
- No complex arithmetic — the symbolic path is pure real Fortran

### SymPy tensor representation

```python
import sympy as sp
from sympy import Matrix, symbols, log, sqrt, diff

# F as a 3x3 symbolic matrix
F = Matrix(3, 3, lambda i, j: sp.Symbol(f'F_{i+1}{j+1}'))

# Material parameters
mu, lam = sp.symbols('mu lam', positive=True)

# Kinematics
J = F.det()
Finv = F.inv()
FinvT = Finv.T

# PK1 stress
P = mu * F + (lam * sp.log(J) - mu) * FinvT

# Tangent: differentiate each P component w.r.t. each F component
dPdF = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3)
for i in range(3):
    for J1 in range(3):
        for k in range(3):
            for L in range(3):
                dPdF[i, J1, k, L] = sp.diff(P[i, J1],
                                             F[k, L])

# CSE to simplify
replacements, simplified = sp.cse(dPdF)

# Fortran code generation
from sympy.printing.fcode import fcode
for var, expr in replacements:
    print(f"      {var} = {fcode(expr)}")
```

### Pattern matching (optional optimization)

Before running CSE, scan the symbolic tangent for patterns:

```python
# Known patterns to match:
patterns = {
    'J * Finv.T':         'CALL ddetdF33(F, J, dJdF)',
    '-Finv.T @ Finv.T':   'CALL dFinvTdF(Finv, dFinvT4)',
    '2 * F':              'CALL dI1dF(F, dI1)',
}
```

If a subexpression matches, replace it with a library call and remove
those terms from the CSE input. This produces shorter, more readable
generated code.

The generator can skip pattern matching entirely and just use CSE + fcode.
The output is correct but verbose. Pattern matching is an optimization.

## Verification Protocol

Every generated tangent — symbolic or CS — is verified the same way:

### 1. Python-level verification (in `verify()`)

```python
def verify(self):
    F_test = random_F()  # nonzero deformation

    if self.tangent == 'symbolic':
        P_sym, dPdF_sym = self._symbolic_tangent(F_test)
        P_cs,  dPdF_cs  = self._cs_tangent(F_test)
        assert np.allclose(P_sym, P_cs, atol=1e-12)
        assert np.allclose(dPdF_sym, dPdF_cs, atol=1e-10)
        print("Symbolic tangent verified against CS")

    # Also check CS against finite differences
    dPdF_fd = self._fd_tangent(F_test, eps=1e-6)
    assert np.allclose(dPdF_cs, dPdF_fd, rtol=1e-4)
```

### 2. Fortran-level verification (single-element stiffness check)

Run a single-element stiffness check with the generated .for file: the
assembled AMATRX is compared against a finite-difference perturbation of
the RHS.

### 3. Cross-verification between modes

For any model that supports both symbolic and CS, generate both .for
files and compare outputs on the same test problem. They should match
to machine precision.

## Implementation Sequence

### Phase A: Layer 2 library (foundation)

Write and test the derivative identity subroutines. Add to
`templates/tangent_identities.for`. Verify each against CS.

Deliverable: `tangent_identities.for` with 8–10 identities, each
tested against CS to 1e-12.

### Phase B: SymPy probe + NeoHookean symbolic (proof of concept)

Implement `symbolic_tangent.py` with:
1. `analyze()` — run material with SymPy symbolic F, classify mode
2. `symbolic_tangent()` — SymPy diff + CSE
3. `generate_fortran_tangent()` — emit Fortran via fcode
4. `verify()` — cross-check against full CS

Test on NeoHookean: analyze → symbolic → generate → verify.

Deliverable: Python module that auto-classifies NeoHookean as
'symbolic' and generates compilable pure-real Fortran tangent.

### Phase C: Gel model full symbolic tangent

Extend to multi-field: all 12 tangent blocks for the gel model
(dP/dF, dP/dp, dflux/dgrad_mu, etc.). All closed-form, all symbolic.

Deliverable: All 12 gel tangent blocks generated as pure-real Fortran,
verified against CS.

### Future Phase D: Hybrid mode for Ogden/Hencky

Implement the future hybrid path: symbolic chain rule plus local CS around
eigendecomposition. This is a separate research project, not a near-term
generator refactor.

Deliverable: Ogden model with 6 CS perturbations of symmetric C
(not 9 of full F), chain-ruled with symbolic dC/dF.

### Phase E: Integration with generators

Wire `SymbolicTangent` into `umat_gen.py` and `uel_gen.py`.
Add `tangent='auto'/'symbolic'/'cs'` user option.

## What This Does NOT Change

- The CS engine is untouched — it remains the default and the oracle
- The UEL architecture is untouched — still total Lagrangian, PK1
- The Material class interface is untouched — users write the same Python
- The WeakForm is untouched
- Complex arithmetic in tensor_ops.for is untouched (still needed for CS)
- The AST translator for Fortran code generation is untouched
- Internal variable handling is untouched (CS only)

## File Organization

```
abaqus_ufl/
  core/
    tensor.py              (unchanged)
    cs_engine.py           (unchanged)
    sympy_tangent.py       (NEW — SymPy differentiation + codegen)
    ast_classifier.py      (NEW — symbolic/CS classification)
  generators/
    umat_gen.py            (modified — adds symbolic tangent path)
    uel_gen.py             (modified — future, per-block classification)
  templates/
    tensor_ops.for         (unchanged)
    tangent_identities.for (NEW — Layer 2 derivative library)
    cs_tangent_engine.for  (unchanged)
    umat_wrapper.for       (unchanged)
  examples/
    NeoHookean/
      neohookean.py
      neohookean_symbolic.for  (generated, pure real)
      neohookean_cs.for        (generated, complex-step)
```

## Phasing and Dependencies

Symbolic tangents are a performance enhancement layered on the validated
complex-step engine, which remains the default and the oracle.

Phases A–C only add new files, so they can proceed independently of the rest
of the generator. Phase D (the hybrid local-CS path) is a separate research
contribution and should not be conflated with the current whole-method
symbolic tangent implementation. Phase E wires the generator integration
together and depends on the symbolic pipeline being validated against the
complex-step oracle first.

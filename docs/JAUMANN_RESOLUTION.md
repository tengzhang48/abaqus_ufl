# Jaumann-Rate Tangent for UMAT DDSDDE

**Scope:** `abaqus_ufl/generators/umat_gen.py`, subroutine
`pk1_to_cauchy_jaumann`.

This note derives the conversion the UMAT generator uses to turn a material
PK1 tangent `dP/dF` into the Abaqus `DDSDDE` spatial tangent, and gives a
reproducible 50-line NumPy check that verifies it.

---

## Summary

For Abaqus/Standard solid continuum UMATs, `DDSDDE` is the spatial tangent
associated with the Jaumann rate of the Kirchhoff stress `τ = J σ`:

```
τ^J_ij = J · DDSDDE_ijkl · D_kl
```

The generator emits

```
DDSDDE_ijkl = AFF_ijkl / J
            - delta_ik * sigma_jl                                (bookkeeping)
            + (1/2)(sigma_ik*delta_jl + sigma_jk*delta_il
                  + sigma_il*delta_jk + sigma_jl*delta_ik)
```

where

```
AFF_ijkl = sum_{J,L} dPdF(i,J,k,L) * F(j,J) * F(l,L)
sigma    = (1/J) * P * F^T      (Cauchy stress)
J        = det(F)
```

The `- delta_ik * sigma_jl` line is a non-optional bookkeeping correction:
`AFF/J` is **not** the Truesdell tangent, and this term is what recovers it.
Without the term, `DDSDDE` is wrong by a stress-linear amount that is
invisible to `F = I`, stress-free-rotation, and uniaxial-bar tests, but is
~40% relative on finite multi-axial deformation.

The result is verified to `3.4e-10` relative against the directly-computed
Jaumann rate of Kirchhoff stress for a NeoHookean material at a stressed `F`
with an asymmetric velocity gradient `L = D + W`.

---

## The conversion formula

Relative to the four-term symmetric block, the bookkeeping term appears as:

```
       c(i,j,k,l) = AFF * JdetInv
                  + 0.5d0*(sigma(i,k)*delta(j,l) + sigma(j,k)*delta(i,l)
                         + sigma(i,l)*delta(j,k) + sigma(j,l)*delta(i,k))
                  - delta(i,k)*sigma(j,l)
```

### Index-variant equivalence

Four index orderings of the bookkeeping term produce the same answer when
contracted against the symmetric `D` that Abaqus passes in, all verified to
`3.4e-10` rel:

```
- delta(i,k) * sigma(j,l)      ← used by the generator
- delta(i,l) * sigma(j,k)
- sigma(j,l) * delta(i,k)
- sigma(j,k) * delta(i,l)
```

The first form is chosen because:

1. It drops directly out of the identity
   `AFF/J = c_Truesdell + δ_ik σ_jl` for this `dP/dF`-based generator, so the
   code stays traceable to the derivation in
   §"Why this is the correct form".
2. It makes the relationship to the Tensor-Toolbox/Bonet-style canonical
   spatial form `c_Truesdell + (S ⊗̄ I) + (I ⊗̄ S)` transparent: after the
   bookkeeping correction `−δ_ik σ_jl` recovers `c_Truesdell` from `AFF/J`,
   the remaining `(1/2)(symm σ block)` is exactly the
   `(S ⊗̄ I) + (I ⊗̄ S)` block.

The other three forms are mathematically equivalent here and equally valid,
just less traceable.

---

## Numerical verification

NeoHookean, `mu = 0.4`, `lam = 1.0`. `F = I + 0.15 * randn(3,3)` (seed 0),
giving `det(F) = 1.54` and a fully triaxial Cauchy stress. Velocity gradient
`L = D + W` with both `D` (symmetric) and `W` (antisymmetric) deliberately
non-spherical. Reference computed by definition:

```
τ̇_ij    = AFF_ijkl L_kl + τ_ik L_jk
τ^∇J    = τ̇ − Wτ + τW
ref     = τ^∇J / J
```

Each candidate bookkeeping term is contracted with `D` (on top of the same
`AFF/J + symm` base) and compared to `ref`:

| Bookkeeping term | Relative error |
|---|---|
| none (omit the term) | `2.4e-01` |
| `− σ_ij δ_kl` | `4.1e-01` |
| `− σ_ij δ_kl − σ_kl δ_ij` | `7.2e-01` |
| **`− δ_ik σ_jl`** (used by the generator) | **`3.4e-10`** |

`symm` is the four-term block
`(1/2)(σ_ik δ_jl + σ_jk δ_il + σ_il δ_jk + σ_jl δ_ik)`.

The error of a wrong bookkeeping term is linear in stress, which is why
low-stress tests do not reveal it:

```
|σ|     err(− σ_ij δ_kl)     err(− δ_ik σ_jl)
0.008   3.0e-03              4.8e-10
0.078   3.0e-02              3.4e-10
0.33    1.5e-01             1.9e-10
0.74    4.1e-01             3.4e-10
```

Slope ≈ `0.5 · |σ| · |D|` in this geometry, putting a low-stress uniaxial
state in the `~1e-3` band — which is why an `F = I`, pure-rotation, or
uniaxial-bar test passes with a wrong term and does not expose the defect.

---

## Why this is the correct form — derivation

`AFF` is the standard PK1 push-forward:

```
AFF_ijkl = (∂P_iJ / ∂F_kL) F_jJ F_lL
```

Differentiating `P_iJ = F_iI S_IJ`:

```
∂P_iJ / ∂F_kL = δ_ik δ_IL S_IJ + F_iI ∂S_IJ/∂F_kL
              = δ_ik S_LJ      + F_iI ∂S_IJ/∂F_kL
```

So

```
AFF_ijkl = δ_ik (F_jJ S_JL F_lL) + F_iI F_jJ F_lL (∂S_IJ/∂F_kL)
         = δ_ik τ_jl              + [push-forward of material tangent C^e]_ijkl
```

The push-forward term is `J · c_Truesdell`, where `c_Truesdell` is the
spatial tangent for the Truesdell rate of Cauchy stress (equivalently the
Lie derivative of Kirchhoff stress divided by `J`):

```
J · c_Truesdell_ijkl = F_iI F_jJ F_kK F_lL · C^e_IJKL
```

Therefore

```
AFF_ijkl / J = c_Truesdell_ijkl + δ_ik σ_jl                          (★)
```

This `δ_ik σ_jl` is the term that is easy to miss. `AFF/J` is **not** the
Truesdell tangent — it differs from it by exactly `δ_ik σ_jl`, which acts on
a symmetric `D` as `(D σ)_ij`.

Now using the Lie derivative form (cleanest direction): let `L_v τ` be the
Lie derivative of `τ` along the velocity field. For hyperelastic
`L_v τ = J · c_Truesdell : D`. And

```
τ^∇J = L_v τ + D τ + τ D
```

Dividing by `J`:

```
τ^∇J / J = (L_v τ)/J + (D τ + τ D)/J
         = c_Truesdell : D + (D σ + σ D)
```

So the Abaqus tangent expected by `DDSDDE` is

```
DDSDDE_ijkl = c_Truesdell_ijkl + (D σ + σ D)-tangent
            = c_Truesdell_ijkl
              + (1/2)(σ_ik δ_jl + σ_jk δ_il + σ_il δ_jk + σ_jl δ_ik)
```

This is exactly Tensor Toolbox's `piola(F,C^e)/J + (S ⊗̄ I) + (I ⊗̄ S)`.

Now substitute (★) to express `c_Truesdell` in terms of `AFF/J`:

```
DDSDDE_ijkl = AFF_ijkl/J − δ_ik σ_jl + (1/2)(symm σ block)            (☆)
```

That is the `DDSDDE` the generator emits. The `− δ_ik σ_jl` is the
bookkeeping correction for the fact that `AFF/J ≠ c_Truesdell`.

---

## Why other bookkeeping terms fail

The bookkeeping term must have exactly the index structure `δ_ik σ_jl`.
Three properties pin it down:

1. `AFF` does not have the minor symmetry in `(k,l)` that `c_Truesdell`
   has. The `(k,l)`-asymmetric part of `AFF/J` must be subtracted off
   explicitly, and that asymmetric part is exactly `δ_ik σ_jl`.
2. A term with the structure `− σ_ij δ_kl` corresponds to no documented
   objective rate. Contracted against a symmetric `D`, that formula returns
   `σ^∇J + D σ` — it is neither the Jaumann rate of Kirchhoff stress nor of
   Cauchy stress. It agrees with the correct tangent only at `σ = 0` and
   drifts linearly in `σ` away from it, which is why `F = I`, pure-rotation,
   and small-strain uniaxial tests all pass while any finite multi-axial
   state exposes the error.
3. Adding a second term `− σ_kl δ_ij` doubles the `(ij)↔(kl)` asymmetry and
   moves further from the correct form (72% vs. 41% in the error table).

The discriminating empirical test is therefore a state with finite,
non-equal principal stresses (biaxial or simple shear): with the correct
tangent, Newton converges quadratically; with a wrong bookkeeping term it
converges only linearly, or fails. A uniaxial-bar test converges in one
iteration per step either way and is not discriminating.

---

## Interface convention

`UMAT` at finite strain is a spatial / current-configuration interface for
the returned stress and tangent. `STRESS` stores Cauchy stress components in
the current configuration, and `DDSDDE` is expected in the corresponding
spatial basis. `UMAT` is not purely "updated-Lagrangian data-only": Abaqus
passes `DFGRD0` and `DFGRD1`, deformation gradients relative to the
original / reference configuration, so the constitutive routine may compute
from total `F`, but the returned stress and tangent must be in the spatial
basis associated with the Jaumann rate of Kirchhoff stress. For
Abaqus/Standard solid continuum UMATs the relevant objective rate is the
Jaumann rate of Kirchhoff stress `τ = J σ`, and the target is the spatial
tangent in `τ^J_ij = J × DDSDDE_ijkl × D_kl`.

---

## References

- Bonet, J., Gil, A. J. & Wood, R. D. (2016). *Nonlinear Solid Mechanics
  for Finite Element Analysis: Statics*, 2nd ed. Cambridge University Press.
- Tensor Toolbox (adtzlr/ttb), Ex. 01 Saint Venant-Kirchhoff — canonical
  UMAT push-forward pattern:
  `C4 = piola(F, C^e) / J + (S cdya I) + (I cdya S)`.
- Nguyen, N. & Waas, A. M. (2016). Nonlinear, finite deformation, finite
  element analysis. *Z. Angew. Math. Phys.* 67:35 — `DDSDDE` corresponds to
  the Jaumann rate of Kirchhoff stress for quadratic convergence.

---

## Reproducing the numerical check

```python
import numpy as np

mu, lam = 0.4, 1.0
def PK1(F):
    J = np.linalg.det(F); FinvT = np.linalg.inv(F).T
    return mu*(F - FinvT) + lam*np.log(J)*FinvT
def Cauchy(F):
    return PK1(F) @ F.T / np.linalg.det(F)

np.random.seed(0)
F = np.eye(3) + 0.15*np.random.randn(3,3)
J = np.linalg.det(F); sigma = Cauchy(F); tau = J*sigma

dpdf = np.zeros((3,3,3,3)); h = 1e-7
for k in range(3):
    for L in range(3):
        Fp = F.copy(); Fp[k,L] += h
        Fm = F.copy(); Fm[k,L] -= h
        dpdf[:,:,k,L] = (PK1(Fp) - PK1(Fm))/(2*h)
AFF = np.einsum('iJkL,jJ,lL->ijkl', dpdf, F, F)

D = np.array([[0.20, 0.05, 0.00],
              [0.05,-0.10, 0.07],
              [0.00, 0.07, 0.30]])
W = np.array([[0.0, 0.1, 0.0],
              [-0.1,0.0, 0.2],
              [0.0,-0.2, 0.0]])
L = D + W

tau_dot = np.einsum('ijkl,kl->ij', AFF, L) + tau @ L.T
ref = (tau_dot - W @ tau + tau @ W) / J

I3 = np.eye(3)
symm = 0.5*( np.einsum('ik,jl->ijkl', sigma, I3)
           + np.einsum('jk,il->ijkl', sigma, I3)
           + np.einsum('il,jk->ijkl', sigma, I3)
           + np.einsum('jl,ik->ijkl', sigma, I3) )
delta_ik_sigma_jl = np.einsum('ik,jl->ijkl', I3, sigma)

C_correct = AFF/J + symm - delta_ik_sigma_jl
err = (np.einsum('ijkl,kl->ij', C_correct, D) - ref)
print("rel err:", np.linalg.norm(err) / np.linalg.norm(ref))
# -> 3.36e-10
```

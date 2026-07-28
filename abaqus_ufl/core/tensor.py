"""
Tensor algebra for constitutive models.

Every function in this module:
  1. Works with real NumPy arrays and the infinitesimal complex-step
     perturbations used by verification/code generation
  2. Has a known Fortran counterpart for code generation

Most algebraic primitives also accept arbitrary complex values. ``eig`` has a
narrower constitutive contract: real symmetric state plus an infinitesimal
complex-step perturbation, not a general finite-complex eigensolver.

Fortran mapping (used by the code generator):
    det(A)    -> det33z(Az)
    inv(A)    -> inv33z(Az, Ainvz)
    trace(A)  -> A(1,1)+A(2,2)+A(3,3)
    sym(A)    -> 0.5*(A + A^T)
    dev(A)    -> A - (1/3)*trace(A)*I
    eye(n)    -> eye33z(Az)
    exp(x)    -> EXP(xz)              [scalar only]
    log(x)    -> LOG(xz)              [scalar only]
    sqrt(x)   -> SQRT(xz)             [scalar only]
    normalize(v) -> v / sqrt(v^T v)   [3-vector, no conjugation]
    dyad(a,b) -> a(i)*b(j) loop
    cross(a,b)-> explicit formula
    sqrtm(A)  -> sqrtm33z(A, S)       [3x3 symmetric matrix]
    logm(A)   -> logm33z(A, L)        [3x3 symmetric matrix]
    expm(A)   -> expm33z(A, E)        [3x3 symmetric matrix]
    polar(F)  -> polar33z(F, R, U)    [polar decomposition]

Matrix operations @ and .T are native Python/Fortran and don't need
explicit mapping (@ -> matmul33z, .T -> transpose33z).

Backend status:
  - Python sqrtm/logm/expm use robust iterative algorithms with fixed
    iteration counts. This is the material-point oracle path.
  - Generated Fortran defaults to fixed-count iterative helpers for
    sqrtm/logm/expm/polar. The eig-based matrix-function backend remains an
    explicit compatibility/debugging option.
  - Explicit eig calls use the trigonometric real-value path and a
    rotation-invariant, first-order complex-step fallback in both Python and
    generated Fortran.

The backend split is intentional for now and is tracked in the
matrix-function design notes.  f2py material tests compare generated Fortran
against the Python oracle at representative states.

See docs/matrix_functions_design.md for details.

IMPORTANT: Users must ONLY use functions from this module in their
Material methods. Using raw numpy (np.linalg.det, np.abs, etc.)
will work for verify() but will FAIL at code generation time.
"""

import numpy as np


# =====================================================================
# Scalar functions (complex-safe intrinsics)
# =====================================================================

def exp(x):
    """Exponential. Fortran: CDEXP(xz)."""
    return np.exp(x)


def log(x):
    """Natural logarithm. Fortran: CDLOG(xz)."""
    return np.log(x)


def sqrt(x):
    """Square root. Fortran: CDSQRT(xz)."""
    return np.sqrt(x)


# =====================================================================
# 3x3 matrix operations
# =====================================================================

def det(A):
    """
    Determinant of 3x3 matrix. Works for real and complex.

    Fortran: det33z(Az) or det33d(Ad).
    """
    return (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
          - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
          + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))


def inv(A):
    """
    Inverse of 3x3 matrix via cofactor formula. Works for complex.

    Fortran: inv33z(Az, Ainvz).
    """
    d = det(A)
    Ai = np.empty_like(A)
    Ai[0, 0] = (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) / d
    Ai[0, 1] = (A[0, 2] * A[2, 1] - A[0, 1] * A[2, 2]) / d
    Ai[0, 2] = (A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]) / d
    Ai[1, 0] = (A[1, 2] * A[2, 0] - A[1, 0] * A[2, 2]) / d
    Ai[1, 1] = (A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]) / d
    Ai[1, 2] = (A[0, 2] * A[1, 0] - A[0, 0] * A[1, 2]) / d
    Ai[2, 0] = (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]) / d
    Ai[2, 1] = (A[0, 1] * A[2, 0] - A[0, 0] * A[2, 1]) / d
    Ai[2, 2] = (A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) / d
    return Ai


def trace(A):
    """Trace of 3x3 matrix. Fortran: A(1,1)+A(2,2)+A(3,3)."""
    return A[0, 0] + A[1, 1] + A[2, 2]


def sym(A):
    """Symmetric part: (A + A^T) / 2."""
    return 0.5 * (A + A.T)


def dev(A):
    """Deviatoric part: A - (1/3)*tr(A)*I. Returns 3x3."""
    tr = trace(A)
    return A - (tr / 3.0) * eye(3, dtype=A.dtype)


def skew(A):
    """Skew-symmetric part: (A - A^T) / 2."""
    return 0.5 * (A - A.T)


def eye(n=3, dtype=float):
    """Identity matrix. Fortran: eye33z or eye33d."""
    return np.eye(n, dtype=dtype)


def dyad(a, b):
    """Dyadic product a ⊗ b. Returns (n, n) from two n-vectors."""
    return np.outer(a, b)


def sym3(a11, a22, a33, a12, a13, a23):
    """Symmetric 3x3 tensor from its 6 unique entries. Fortran: sym3z.

    Returns ``[[a11, a12, a13], [a12, a22, a23], [a13, a23, a33]]``.
    Useful for constant tensors (base fabrics, anisotropy directors,
    projection operators) that cannot be expressed via vector literals.
    """
    dtype = np.result_type(a11, a22, a33, a12, a13, a23, float)
    T = np.empty((3, 3), dtype=dtype)
    T[0, 0] = a11
    T[1, 1] = a22
    T[2, 2] = a33
    T[0, 1] = a12
    T[1, 0] = a12
    T[0, 2] = a13
    T[2, 0] = a13
    T[1, 2] = a23
    T[2, 1] = a23
    return T


def cross(a, b):
    """Cross product of 3-vectors."""
    c = np.empty_like(a)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c


def normalize(v):
    """Return a unit 3-vector using the CS-safe bilinear norm.

    The normalization is

    ``v / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])``.

    There is deliberately no conjugation, ``abs()``, or ``.real`` projection:
    around a nonzero real vector this expression is locally holomorphic and
    therefore preserves complex-step derivatives.  The input must be nonzero.

    This primitive is intended for constitutive formulas that require a unit
    direction, such as slip-system construction from a column returned by
    :func:`eig`.  ``eig`` does not guarantee the scale or normalization of its
    columns.  Matrix-function reconstruction must therefore use
    ``V @ diag(f(lam)) @ inv(V)``; call ``normalize`` only when the constitutive
    theory itself requires a unit direction.

    Fortran: an inline 3-vector normalization emitted by the generator.
    """
    scale = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return v / scale


# =====================================================================
# Invariants (convenience)
# =====================================================================

def I1(A):
    """First invariant = trace(A)."""
    return trace(A)


def I2(A):
    """Second invariant = 0.5*(tr(A)^2 - tr(A^2))."""
    trA = trace(A)
    trA2 = trace(A @ A)
    return 0.5 * (trA * trA - trA2)


# =====================================================================
# Eigendecomposition of 3x3 symmetric matrix
# =====================================================================

def eig(A):
    """
    Eigenvalues and eigenvectors of a 3x3 symmetric matrix.

    Distinct real spectra use the trigonometric solution for the depressed
    cubic. Repeated/quasi-repeated spectra and complex-step calls use a real
    symmetric principal basis followed by a perturbation fallback.

    Eigenvector column scale is an implementation detail: columns are not
    guaranteed to be unit or orthogonal. Use ``inv(V)`` (not ``V.T``) for
    reconstruction: ``f(A) = V @ diag(f(lam)) @ inv(V)``. If a constitutive
    equation needs a unit direction, normalize that column explicitly.

    Within the quasi-repeated grouping band, only eigenspace-invariant
    reconstructions are supported. Derivatives of an arbitrarily selected
    individual eigenvector or rank-one projector are gauge-dependent,
    ill-conditioned, and outside this API contract.

    Args:
        A: 3x3 symmetric real state, optionally carrying an infinitesimal
           complex-step perturbation. Arbitrary finite-complex matrix
           eigendecomposition is outside this constitutive DSL contract.

    Returns:
        lam: length-3 array of eigenvalues (sorted by real part)
        V:   3x3 array whose columns correspond to ``lam``. Column scale and
             normalization are unspecified.

    Fortran: eig33z(A, lam, V)
    """
    A = np.asarray(A, dtype=np.complex128)

    # A repeated spectrum is not a coordinate-dependent event.  During a
    # complex-step call, first diagonalize the real symmetric state, rotate the
    # full complex perturbation into that basis, and reuse the established
    # quasi-degenerate treatment there.  The former guard only recognized
    # diagonal repeated states and failed after a real orthogonal rotation.
    if np.max(np.abs(A.imag)) > 0.0:
        return _eig_rotated_fallback(A)

    # Invariants
    p1 = A[0, 1] ** 2 + A[0, 2] ** 2 + A[1, 2] ** 2
    q = (A[0, 0] + A[1, 1] + A[2, 2]) / 3.0  # trace / 3

    # B = A - q*I
    B = A.copy()
    B[0, 0] -= q
    B[1, 1] -= q
    B[2, 2] -= q

    p2 = B[0, 0] ** 2 + B[1, 1] ** 2 + B[2, 2] ** 2 + 2.0 * p1
    p = np.sqrt(p2 / 6.0)

    # Guard: A is a multiple of identity (possibly plus a tiny — e.g.
    # complex-step — deviation). The trigonometric path is singular
    # here; use the perturbation-theory fallback, which keeps the
    # deviation's contribution to eigenvalues AND eigenvectors to first
    # order. The old fallback (lam = q, V = I) silently zeroed CS
    # derivatives of spectral functions at C = I (audit finding H2).
    if abs(p.real) < 1e-30:
        return _eig_rotated_fallback(A)

    # For nearly-diagonal matrices the trigonometric formula is unstable
    # when two eigenvalues are nearly equal (|r| ≈ 1).  In that regime
    # tiny rounding errors in det(C) are amplified by the singular
    # derivative dλ/dr and produce spurious splitting.  Fall back to
    # perturbation theory about the diagonal (keeps the off-diagonal CS
    # parts that the old diag/V=I fallback discarded — finding H2).
    p2_abs = abs(p2.real)
    if p2_abs > 0.0 and abs(p1.real) < 1e-14 * p2_abs:
        return _eig_rotated_fallback(A)

    # C = B / p
    C = B / p

    # r = det(C) / 2
    det_C = (C[0, 0] * (C[1, 1] * C[2, 2] - C[1, 2] * C[2, 1])
             - C[0, 1] * (C[1, 0] * C[2, 2] - C[1, 2] * C[2, 0])
             + C[0, 2] * (C[1, 0] * C[2, 1] - C[1, 1] * C[2, 0]))
    r = det_C / 2.0

    # Repeated and quasi-repeated spectra are also rotation invariant on the
    # pure real value path. Cross-product columns become singular when the
    # repeated eigenspace is not coordinate-aligned, so reuse the real-basis
    # fallback before evaluating the cubic angle.
    if abs(1.0 - abs(r.real)) <= 1.0e-10:
        return _eig_rotated_fallback(A)

    # Clamp r to [-1, 1] to avoid arccos domain issues.
    # For complex-step inputs we clamp the real part only; the imag
    # derivative is independent and unaffected by the clamp.
    r_val = max(-1.0, min(1.0, r.real))
    r = complex(r_val, r.imag) if hasattr(r, 'imag') else r_val

    # Trigonometric solution for three eigenvalues
    phi = np.arccos(r) / 3.0
    lam = np.empty(3, dtype=np.complex128)
    lam[0] = q + 2.0 * p * np.cos(phi)
    lam[1] = q + 2.0 * p * np.cos(phi + 2.0 * np.pi / 3.0)
    lam[2] = q + 2.0 * p * np.cos(phi + 4.0 * np.pi / 3.0)

    # Sort by real part for consistent ordering
    idx = np.argsort(lam.real)
    lam = lam[idx]

    # Compute eigenvectors
    V = np.zeros((3, 3), dtype=np.complex128)
    I = np.eye(3, dtype=np.complex128)
    for k in range(3):
        V[:, k] = _eigvec(A - lam[k] * I, k)

    return lam, V


def _eig_rotated_fallback(A):
    """Rotation-invariant repeated-spectrum/complex-step eig fallback.

    ``A.real`` supplies a real orthonormal principal basis.  In that basis the
    full complex matrix is near diagonal, so ``_eig_near_diagonal`` can resolve
    repeated and quasi-repeated blocks without assuming that the original
    coordinate axes are principal axes. For complex-step input this is a
    first-order construction. Pure real, well-separated value calls stay on
    ``eig``'s trigonometric path; repeated/near-repeated values use this path.
    """
    _, Q = np.linalg.eigh(A.real)
    A_principal = Q.T @ A @ Q
    lam, V_principal = _eig_near_diagonal(A_principal)
    return lam, Q.astype(np.complex128) @ V_principal


def _eig_near_diagonal(A):
    """Eigendecomposition of a (near-)diagonal complex-symmetric 3x3 by
    quasi-degenerate perturbation theory.

    Used by eig()'s degeneracy guards, where the trigonometric path is
    singular. A = D + E with D = diag(A) and E small (real rounding
    and/or a complex-step perturbation). Diagonal entries are grouped by
    closeness; within each group the derivative-carrying imaginary-symmetric
    deviation is diagonalized during a complex-step call (the real-symmetric
    deviation is used for a pure value call; both use numpy.linalg.eigh with no
    recursion into eig), and across groups the eigenvectors receive the
    standard first-order correction

        V[i, a] += (e_i . A v0_a) / (lam_a - D_i).

    This preserves the imaginary (derivative-carrying) part of E in both
    eigenvalues and eigenvectors, so spectral reconstructions
    V f(L) inv(V) recover the correct divided-difference / f' CS
    derivatives at degenerate states (C = I at the first increment,
    diagonal C at uniaxial/biaxial states).
    """
    D = np.diag(A).copy()
    E = A - np.diag(D)
    scale = float(np.max(np.abs(A))) + 1e-300

    # Group (near-)degenerate diagonal entries, walking in sorted order.
    order = np.argsort(D.real)
    groups = [[order[0]]]
    for a in order[1:]:
        if abs(D[a].real - D[groups[-1][-1]].real) <= 1e-6 * scale:
            groups[-1].append(a)
        else:
            groups.append([a])

    lam = np.zeros(3, dtype=np.complex128)
    V = np.zeros((3, 3), dtype=np.complex128)
    for g in groups:
        idx = np.array(g)
        sub = A[np.ix_(idx, idx)]
        n = len(g)
        if n == 1:
            lam[idx[0]] = D[idx[0]]
            V[idx[0], idx[0]] = 1.0
            continue
        # Deviation from the group mean. A CS call must resolve its imaginary
        # direction inside this intentionally quasi-degenerate block; a pure
        # value call uses the real-symmetric deviation.
        dev = sub - (np.trace(sub) / n) * np.eye(n)
        # During a complex-step call, the imaginary block is the derivative
        # direction.  Resolve it even when a small real eigenvalue split is
        # numerically larger; treating every member of this intentionally
        # quasi-degenerate group in the real basis would erase the off-diagonal
        # derivative.  Purely real value calls still use the real deviation.
        M = dev.imag if np.max(np.abs(dev.imag)) > 0.0 else dev.real
        if np.max(np.abs(M)) == 0.0:
            Q = np.eye(n)
        else:
            _, Q = np.linalg.eigh(M)
        sub_rot = Q.T @ sub @ Q
        for a in range(n):
            lam[idx[a]] = sub_rot[a, a]
            V[idx, idx[a]] = Q[:, a]

    # First-order cross-group eigenvector corrections.
    for g in groups:
        others = [i for i in range(3) if i not in g]
        if not others:
            continue
        for a in g:
            v0 = V[:, a]
            for i in others:
                gap = lam[a] - D[i]
                V[i, a] += (A[i, :] @ v0) / gap

    # Sort by real part for consistent ordering (matches the main path).
    idx = np.argsort(lam.real)
    return lam[idx], V[:, idx]


def _eigvec(M, idx):
    """
    Compute eigenvector for singular matrix M = A - lambda*I.

    This helper is used only on ``eig``'s distinct-real trigonometric path and
    returns the raw (unnormalized) cross-product vector. The public ``eig``
    contract does not promise this scaling because other paths choose a real
    symmetric basis. Reconstruction uses ``inv(V)``, so column scaling cancels.
    A constitutive consumer that requires a unit direction should call
    ``normalize()``, which uses the locally holomorphic bilinear norm (not a
    conjugated/Hermitian norm).

    Falls back to coordinate basis for degenerate cases.
    """
    r0, r1, r2 = M[0, :], M[1, :], M[2, :]
    candidates = [_cross3(r0, r1), _cross3(r0, r2), _cross3(r1, r2)]

    # Pick candidate with largest magnitude.
    best = candidates[0]
    best_nsq = best[0]*best[0] + best[1]*best[1] + best[2]*best[2]
    for c in candidates[1:]:
        nsq = c[0]*c[0] + c[1]*c[1] + c[2]*c[2]
        if abs(nsq) > abs(best_nsq):
            best = c
            best_nsq = nsq

    # This path only receives pure-real distinct spectra. A fixed absolute
    # threshold is not scale invariant: for a small but otherwise ordinary
    # rotated tensor the cross products scale as ||A||^2 and their squared
    # norms as ||A||^4. Fall back only when the cross products are exactly
    # zero; sufficiently small spectra are already routed through the rotated
    # fallback before reaching this helper.
    if abs(best_nsq) == 0.0:
        v = np.zeros(3, dtype=np.complex128)
        v[idx] = 1.0
        return v

    return best


def _cross3(a, b):
    """Cross product of 3-vectors (complex-safe)."""
    c = np.empty(3, dtype=np.complex128)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c


# =====================================================================
# Matrix functions (iterative, fixed counts, CS-safe)
# =====================================================================

def sqrtm(A):
    """
    Matrix square root of 3x3 symmetric positive definite matrix.

    Denman-Beavers iteration with a fixed iteration count.

    Fortran: sqrtm33z(A, S)
    """
    A = np.asarray(A, dtype=np.complex128)
    Y = A.copy()
    Z = np.eye(3, dtype=np.complex128)
    for _ in range(20):
        Zinv = inv(Z)
        Yinv = inv(Y)
        Y = 0.5 * (Y + Zinv)
        Z = 0.5 * (Z + Yinv)
    return Y


def logm(A):
    """
    Matrix logarithm of 3x3 symmetric positive definite matrix.

    Inverse scaling-and-squaring followed by a fixed-order Taylor
    expansion.  Gives the Hencky strain: E_hencky = 0.5 * logm(C).

    Fortran: logm33z(A, L)
    """
    A = np.asarray(A, dtype=np.complex128)
    I_mat = np.eye(3, dtype=np.complex128)

    # Repeated square roots (fixed count)
    B = A.copy()
    for _ in range(6):
        Y = B.copy()
        Z = I_mat.copy()
        for _ in range(20):
            Zinv = inv(Z)
            Yinv = inv(Y)
            Y = 0.5 * (Y + Zinv)
            Z = 0.5 * (Z + Yinv)
        B = Y

    # Taylor series for log(B) about B = I
    BpI_inv = inv(B + I_mat)
    X = (B - I_mat) @ BpI_inv
    X2 = X @ X
    Xpow = X.copy()
    L = X.copy()
    for k in range(3, 17, 2):
        Xpow = Xpow @ X2
        L = L + Xpow / float(k)
    L = 2.0 * L

    # Undo scaling
    L = L * 64.0  # 2**6
    return L


def expm(A):
    """
    Matrix exponential of 3x3 symmetric matrix.

    Scaling-and-squaring with a [6/6] Padé approximant.
    Used for the exponential map: Fp_new = expm(dt*Dp) @ Fp_old.

    Fortran: expm33z(A, E)
    """
    A = np.asarray(A, dtype=np.complex128)
    I_mat = np.eye(3, dtype=np.complex128)

    # Fixed scaling: s = 10  (||A/2^10|| <= 0.01 for ||A|| <= 10)
    s = 10
    B = A / 1024.0  # 2**10

    # Padé [6/6] coefficients
    c1 = 0.5
    c2 = 5.0 / 44.0
    c3 = 1.0 / 66.0
    c4 = 1.0 / 792.0
    c5 = 1.0 / 15840.0
    c6 = 1.0 / 665280.0

    B2 = B @ B
    B3 = B2 @ B
    B4 = B2 @ B2
    B5 = B4 @ B
    B6 = B4 @ B2

    N = (I_mat + c1 * B + c2 * B2 + c3 * B3
         + c4 * B4 + c5 * B5 + c6 * B6)
    D = (I_mat - c1 * B + c2 * B2 - c3 * B3
         + c4 * B4 - c5 * B5 + c6 * B6)
    E = inv(D) @ N

    for _ in range(s):
        E = E @ E
    return E


def polar(F):
    """
    Right polar decomposition: F = R @ U, U = sqrtm(F^T @ F).

    Args:
        F: 3x3 deformation gradient (real or complex)
    Returns:
        R: 3x3 rotation
        U: 3x3 right stretch tensor (symmetric positive definite)

    Fortran: polar33z(F, R, U)
    """
    C = F.T @ F
    U = sqrtm(C)
    R = F @ inv(U)
    return R, U


def cubic_roots(a, b, c, d):
    """
    All roots of a*x**3 + b*x**2 + c*x + d = 0.

    Depressed-cubic + Cardano, all in complex arithmetic (CS-safe).

    Args:
        a, b, c, d: scalar coefficients (real or complex)

    Returns:
        r1, r2, r3: scalar roots (may be complex conjugates)

    Fortran: CALL cs_cubic_roots(a, b, c, d, r1, r2, r3)
    """
    a = np.asarray(a, dtype=np.complex128)
    b = np.asarray(b, dtype=np.complex128)
    c = np.asarray(c, dtype=np.complex128)
    d = np.asarray(d, dtype=np.complex128)

    # Depress the cubic: x = t - b/(3a)
    p = (3.0 * a * c - b ** 2) / (3.0 * a ** 2)
    q = (2.0 * b ** 3 - 9.0 * a * b * c + 27.0 * a ** 2 * d) / (27.0 * a ** 3)

    Delta = (q / 2.0) ** 2 + (p / 3.0) ** 3

    sqrt_Delta = np.sqrt(Delta)
    u = (-q / 2.0 + sqrt_Delta) ** (1.0 / 3.0)
    # v chosen so that u*v = -p/3 (avoids branch-cut ambiguity).
    # Guard u=0 (triple-root / p=0 case) by computing v independently.
    if abs(u) > 1e-30:
        v = -p / (3.0 * u)
    else:
        v = (-q / 2.0 - sqrt_Delta) ** (1.0 / 3.0)

    t1 = u + v
    t2 = -0.5 * (u + v) + 1j * np.sqrt(3.0) / 2.0 * (u - v)
    t3 = -0.5 * (u + v) - 1j * np.sqrt(3.0) / 2.0 * (u - v)

    shift = b / (3.0 * a)
    r1 = t1 - shift
    r2 = t2 - shift
    r3 = t3 - shift

    return r1, r2, r3

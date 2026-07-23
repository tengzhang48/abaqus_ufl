"""
Tensor algebra for constitutive models.

Every function in this module:
  1. Works with both real and complex NumPy arrays (CS-safe)
  2. Has a known Fortran counterpart for code generation

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
  - Generated Fortran currently maps sqrtm/logm/expm to eigendecomposition
    helpers (sqrtm33z/logm33z/expm33z). Those helpers use the same
    trigonometric cubic eigensolver plus diagonal/near-diagonal guards.
  - Iterative Fortran helpers also exist in tensor_ops.for, but the
    generator does not select them by default.

The backend split is intentional for now and is tracked in
the matrix-function templates.  f2py material tests compare generated Fortran
against the Python oracle at representative states.

See the matrix-function templates for details.

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

    Uses the trigonometric solution for the depressed cubic (avoids
    branch cut issues of Cardano's algebraic form). Self-contained —
    no LAPACK dependency.

    Eigenvectors are returned unnormalized (raw cross-product vectors)
    to preserve holomorphic properties for CS differentiation. Use
    inv(V) (not V^T) for reconstruction: f(A) = V @ diag(f(lam)) @ inv(V).

    Args:
        A: 3x3 numpy array (real or complex), assumed symmetric

    Returns:
        lam: length-3 array of eigenvalues (sorted by real part)
        V:   3x3 array where V[:,i] is the eigenvector for lam[i]
             (unnormalized, but linearly independent for distinct eigenvalues)

    Fortran: eig33z(A, lam, V)
    """
    A = np.asarray(A, dtype=np.complex128)

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
        return _eig_near_diagonal(A)

    # For nearly-diagonal matrices the trigonometric formula is unstable
    # when two eigenvalues are nearly equal (|r| ≈ 1).  In that regime
    # tiny rounding errors in det(C) are amplified by the singular
    # derivative dλ/dr and produce spurious splitting.  Fall back to
    # perturbation theory about the diagonal (keeps the off-diagonal CS
    # parts that the old diag/V=I fallback discarded — finding H2).
    p2_abs = abs(p2.real)
    if p2_abs > 0.0 and (
        abs(p1.real) < 1e-14 * p2_abs or abs(p1) < 1e-18
    ):
        return _eig_near_diagonal(A)

    # C = B / p
    C = B / p

    # r = det(C) / 2
    det_C = (C[0, 0] * (C[1, 1] * C[2, 2] - C[1, 2] * C[2, 1])
             - C[0, 1] * (C[1, 0] * C[2, 2] - C[1, 2] * C[2, 0])
             + C[0, 2] * (C[1, 0] * C[2, 1] - C[1, 1] * C[2, 0]))
    r = det_C / 2.0

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


def _eig_near_diagonal(A):
    """Eigendecomposition of a (near-)diagonal complex-symmetric 3x3 by
    quasi-degenerate perturbation theory.

    Used by eig()'s degeneracy guards, where the trigonometric path is
    singular. A = D + E with D = diag(A) and E small (real rounding
    and/or a complex-step perturbation). Diagonal entries are grouped by
    closeness; within each group the dominant real-symmetric part of the
    deviation is diagonalized exactly (numpy.linalg.eigh on a real
    matrix — no recursion into eig), and across groups the eigenvectors
    receive the standard first-order correction

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
        # Deviation from the group mean; rotate by the eigenbasis of its
        # dominant real-symmetric part (exact within the group, keeps
        # the subdominant part to first order on the rotated diagonal).
        dev = sub - (np.trace(sub) / n) * np.eye(n)
        M = dev.imag if np.max(np.abs(dev.imag)) > np.max(np.abs(dev.real)) \
            else dev.real
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

    Returns the RAW (unnormalized) cross-product vector.
    Normalization is not needed because _reconstruct uses inv(V),
    and normalization would be non-holomorphic (breaking CS).

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

    # Degenerate: fall back to coordinate basis
    if abs(best_nsq) < 1e-50:
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

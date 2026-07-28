"""
Tests for matrix functions in tensor.py: eig, sqrtm, logm, expm, polar.

Tests:
  1. Known analytical results (diagonal matrices, identity)
  2. Consistency (expm(logm(A))=A, sqrtm(A)^2=A)
  3. Polar decomposition (R orthogonal, R@U=F, det(R)=1)
  4. CS vs FD tangent of logm(C) at three states
  5. CS vs FD tangent of expm for plasticity use case
"""

import numpy as np

from abaqus_ufl.core.tensor import det, inv, trace, eye, log, sqrt, exp
from abaqus_ufl.core.tensor import eig, sqrtm, logm, expm, normalize, polar


def test_known_results():
    """Test against analytically known results."""
    print("Test 1: Known analytical results")
    all_pass = True

    # logm(diag(a,b,c)) = diag(log(a), log(b), log(c))
    A = np.diag([2.0, 3.0, 5.0])
    L = logm(A)
    L_exact = np.diag([np.log(2.0), np.log(3.0), np.log(5.0)])
    err = np.max(np.abs(L - L_exact))
    ok = err < 1e-12
    all_pass = all_pass and ok
    print(f"  logm(diag): err = {err:.2e} [{'PASS' if ok else 'FAIL'}]")

    # expm(diag(a,b,c)) = diag(exp(a), exp(b), exp(c))
    B = np.diag([0.1, -0.2, 0.3])
    E = expm(B)
    E_exact = np.diag([np.exp(0.1), np.exp(-0.2), np.exp(0.3)])
    err = np.max(np.abs(E - E_exact))
    ok = err < 1e-12
    all_pass = all_pass and ok
    print(f"  expm(diag): err = {err:.2e} [{'PASS' if ok else 'FAIL'}]")

    # sqrtm(I) = I
    S = sqrtm(np.eye(3))
    err = np.max(np.abs(S - np.eye(3)))
    ok = err < 1e-12
    all_pass = all_pass and ok
    print(f"  sqrtm(I):   err = {err:.2e} [{'PASS' if ok else 'FAIL'}]")

    # expm(logm(A)) = A
    A = np.array([[2.0, 0.3, 0.1],
                  [0.3, 3.0, 0.2],
                  [0.1, 0.2, 4.0]])
    EL = expm(logm(A))
    err = np.max(np.abs(EL - A))
    ok = err < 1e-10
    all_pass = all_pass and ok
    print(f"  expm(logm(A))=A: err = {err:.2e} "
          f"[{'PASS' if ok else 'FAIL'}]")

    # sqrtm(A)^2 = A
    S = sqrtm(A)
    err = np.max(np.abs(S @ S - A))
    ok = err < 1e-10
    all_pass = all_pass and ok
    print(f"  sqrtm(A)^2=A:    err = {err:.2e} "
          f"[{'PASS' if ok else 'FAIL'}]")

    assert all_pass


def test_polar():
    """Test polar decomposition."""
    print("\nTest 2: Polar decomposition")
    all_pass = True

    F = np.array([[1.3, 0.1, 0.05],
                  [0.05, 1.2, 0.02],
                  [0.01, 0.03, 1.1]])
    R, U = polar(F)

    err = np.max(np.abs(R.T @ R - np.eye(3)))
    ok = err < 1e-10
    all_pass = all_pass and ok
    print(f"  R^T R = I:  err = {err:.2e} [{'PASS' if ok else 'FAIL'}]")

    err = np.max(np.abs(R @ U - F))
    ok = err < 1e-10
    all_pass = all_pass and ok
    print(f"  R @ U = F:  err = {err:.2e} [{'PASS' if ok else 'FAIL'}]")

    err = np.max(np.abs(U - U.T))
    ok = err < 1e-10
    all_pass = all_pass and ok
    print(f"  U = U^T:    err = {err:.2e} [{'PASS' if ok else 'FAIL'}]")

    detR = np.linalg.det(R.real)
    ok = abs(detR - 1.0) < 1e-10
    all_pass = all_pass and ok
    print(f"  det(R) = {detR:.10f} [{'PASS' if ok else 'FAIL'}]")

    assert all_pass


def test_cs_vs_fd_logm():
    """CS vs FD tangent of logm(F^T @ F) at three states."""
    print("\nTest 3: CS vs FD tangent of logm(F^T @ F)")
    all_pass = True

    h = 1e-20
    eps = 1e-7

    cases = [
        ("F = I (degenerate)", np.eye(3)),
        ("F near I",
         np.eye(3) + 0.001 * np.array([[0.3, 0.1, 0.05],
                                         [0.1, -0.2, 0.02],
                                         [0.05, 0.02, 0.1]])),
        ("F large deformation",
         np.array([[1.5, 0.2, 0.1],
                   [0.05, 0.8, 0.05],
                   [0.02, 0.03, 1.2]])),
    ]

    for name, F in cases:
        dL_cs = np.zeros((3, 3, 3, 3))
        dL_fd = np.zeros((3, 3, 3, 3))
        for k in range(3):
            for l in range(3):
                Fz = F.astype(complex).copy()
                Fz[k, l] += 1j * h
                dL_cs[:, :, k, l] = logm(Fz.T @ Fz).imag / h

                Fp = F.copy(); Fp[k, l] += eps
                Fm = F.copy(); Fm[k, l] -= eps
                dL_fd[:, :, k, l] = (logm(Fp.T @ Fp).real
                                     - logm(Fm.T @ Fm).real) / (2*eps)

        err = np.max(np.abs(dL_cs - dL_fd))
        norm = np.max(np.abs(dL_cs)) + 1e-30
        rel = err / norm
        ok = rel < 1e-5
        all_pass = all_pass and ok
        print(f"  {name}: rel err = {rel:.2e} "
              f"[{'PASS' if ok else 'FAIL'}]")

    assert all_pass


def test_cs_vs_fd_expm():
    """CS vs FD tangent of expm for plasticity."""
    print("\nTest 4: CS vs FD tangent of expm")

    h = 1e-20
    eps = 1e-7

    A = np.array([[0.01, 0.005, 0.0],
                  [0.005, -0.005, 0.0],
                  [0.0, 0.0, -0.005]])

    dE_cs = np.zeros((3, 3, 3, 3))
    dE_fd = np.zeros((3, 3, 3, 3))
    for k in range(3):
        for l in range(3):
            Az = A.astype(complex).copy()
            Az[k, l] += 1j * h; Az[l, k] = Az[k, l]
            dE_cs[:, :, k, l] = expm(Az).imag / h

            Ap = A.copy(); Ap[k, l] += eps; Ap[l, k] = Ap[k, l]
            Am = A.copy(); Am[k, l] -= eps; Am[l, k] = Am[k, l]
            dE_fd[:, :, k, l] = (expm(Ap).real - expm(Am).real) / (2*eps)

    err = np.max(np.abs(dE_cs - dE_fd))
    norm = np.max(np.abs(dE_cs)) + 1e-30
    rel = err / norm
    ok = rel < 1e-5
    print(f"  expm tangent: rel err = {rel:.2e} "
          f"[{'PASS' if ok else 'FAIL'}]")
    assert ok


def test_eig_basic():
    """Basic eigendecomposition checks."""
    print("\nTest 5: Eigendecomposition basics")
    all_pass = True

    A = np.array([[4.0, 1.0, 0.5],
                  [1.0, 3.0, 0.2],
                  [0.5, 0.2, 2.0]])
    lam, V = eig(A)

    # V @ diag(lam) @ inv(V) should reconstruct A
    D = np.diag(lam)
    recon = V @ D @ inv(V)
    err = np.max(np.abs(recon - A))
    ok = err < 1e-12
    all_pass = all_pass and ok
    print(f"  V@D@inv(V) = A: err = {err:.2e} "
          f"[{'PASS' if ok else 'FAIL'}]")

    # V should be invertible (det(V) != 0)
    detV = abs(np.linalg.det(V))
    ok = detV > 1e-10
    all_pass = all_pass and ok
    print(f"  det(V) = {detV:.2e} (invertible) "
          f"[{'PASS' if ok else 'FAIL'}]")

    # Eigenvalues should match numpy
    lam_np = np.sort(np.linalg.eigvalsh(A))
    lam_sorted = np.sort(lam.real)
    err = np.max(np.abs(lam_sorted - lam_np))
    ok = err < 1e-12
    all_pass = all_pass and ok
    print(f"  lam vs numpy: err = {err:.2e} "
          f"[{'PASS' if ok else 'FAIL'}]")

    assert all_pass


def test_normalize_value_and_complex_step_derivative():
    """The public vector normalization uses v^T v, not a Hermitian norm."""
    v = np.array([2.0, -3.0, 6.0])
    direction = np.array([0.3, -0.2, 0.4])
    unit = normalize(v)
    norm = np.sqrt(v @ v)

    np.testing.assert_allclose(unit, v / norm, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(unit @ unit, 1.0, rtol=0.0, atol=1e-14)

    h = 1.0e-30
    got = normalize(v.astype(complex) + 1j * h * direction).imag / h
    expected = (np.eye(3) - np.outer(unit, unit)) @ direction / norm
    np.testing.assert_allclose(got, expected, rtol=2e-14, atol=2e-15)


if __name__ == '__main__':
    print("=" * 60)
    print("tensor.py matrix functions: test suite")
    print("=" * 60)

    r1 = test_known_results()
    r2 = test_polar()
    r3 = test_cs_vs_fd_logm()
    r4 = test_cs_vs_fd_expm()
    r5 = test_eig_basic()

    print()
    all_pass = r1 and r2 and r3 and r4 and r5
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("=" * 60)


def test_eig_cs_at_degenerate_states():
    """Audit finding H2: CS derivatives of spectral functions must survive
    the degeneracy guards (identity / diagonal / repeated-eigenvalue
    states). The old fallbacks returned V=I and discarded off-diagonal
    complex-step perturbations, silently zeroing the tangent at C = I
    (the first increment of every simulation)."""
    h = 1e-20

    def cs_dsqrtm(A0, i, j):
        A = A0.astype(complex)
        A[i, j] += 1j * h
        A[j, i] += 1j * h
        lam, V = eig(A)
        S = V @ np.diag(np.sqrt(lam)) @ np.linalg.inv(V)
        return S[i, j].imag / h

    def dd(a, b):
        # divided difference of sqrt (equals f'(a) when a == b)
        return 1.0 / (np.sqrt(a) + np.sqrt(b))

    cases = [
        (np.eye(3), 0, 1, dd(1.0, 1.0)),               # triple degenerate
        (4.0 * np.eye(3), 0, 2, dd(4.0, 4.0)),         # scaled identity
        (np.diag([1.0, 2.0, 3.0]), 0, 1, dd(1.0, 2.0)),  # distinct diag
        (np.diag([1.0, 2.0, 3.0]), 1, 2, dd(2.0, 3.0)),
        (np.diag([2.0, 2.0, 3.0]), 0, 1, dd(2.0, 2.0)),  # degenerate pair
        (np.diag([2.0, 2.0, 3.0]), 0, 2, dd(2.0, 3.0)),  # cross coupling
        # A near-distinct pair inside the quasi-repeated band is treated in
        # the limiting degenerate basis. The individual basis derivative is
        # intentionally not defined, but the invariant sqrt reconstruction
        # must retain the correct divided-difference derivative.
        (np.diag([2.0, 2.0 + 1e-8, 3.0]), 0, 1,
         dd(2.0, 2.0 + 1e-8)),
        # A repeated pair inside a broader quasi-degenerate cluster must still
        # resolve the CS direction within that pair.
        (np.diag([2.0, 2.0, 2.0 + 1e-8]), 0, 1,
         dd(2.0, 2.0)),
    ]
    for A0, i, j, exact in cases:
        got = cs_dsqrtm(A0, i, j)
        assert abs(got - exact) < 1e-12 * abs(exact), \
            f"CS d(sqrtm)[{i},{j}] at {np.diag(A0)}: {got} != {exact}"


def test_eig_real_path_unchanged_at_guarded_states():
    """The perturbation fallback must reproduce plain eigendecomposition
    for real inputs at the guarded states."""
    for A0 in (np.eye(3), np.diag([3.0, 1.0, 2.0]), np.diag([2.0, 2.0, 5.0])):
        lam, V = eig(A0.astype(complex))
        recon = (V @ np.diag(lam) @ np.linalg.inv(V)).real
        assert np.max(np.abs(recon - A0)) < 1e-12
        assert np.max(np.abs(np.sort(lam.real)
                             - np.sort(np.diag(A0)))) < 1e-12


def test_eig_cs_at_rotated_repeated_and_near_repeated_spectra():
    """The CS fallback must be invariant to the real principal basis."""
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    Q = (
        np.eye(3)
        + np.sin(angle) * K
        + (1.0 - np.cos(angle)) * (K @ K)
    )
    direction = np.array([
        [0.30, -0.20, 0.15],
        [-0.20, -0.10, 0.25],
        [0.15, 0.25, 0.40],
    ])
    h = 1.0e-20

    for gap in (0.0, 1.0e-8):
        diag = np.array([2.0, 2.0 + gap, 3.0])
        A0 = Q @ np.diag(diag) @ Q.T
        Az = A0.astype(complex) + 1j * h * direction
        lam, V = eig(Az)
        got = np.imag(
            V @ np.diag(np.sqrt(lam)) @ np.linalg.inv(V)
        ) / h

        h_principal = Q.T @ direction @ Q
        exact_principal = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                exact_principal[i, j] = (
                    h_principal[i, j]
                    / (np.sqrt(diag[i]) + np.sqrt(diag[j]))
                )
        exact = Q @ exact_principal @ Q.T
        np.testing.assert_allclose(got, exact, rtol=3e-8, atol=3e-10)


def test_eig_value_at_rotated_repeated_and_near_repeated_spectra():
    """Real repeated-spectrum reconstruction must not return singular V."""
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    Q = (
        np.eye(3)
        + np.sin(angle) * K
        + (1.0 - np.cos(angle)) * (K @ K)
    )

    for gap in (0.0, 1.0e-12, 1.0e-8):
        diag = np.array([2.0, 2.0 + gap, 3.0])
        A0 = Q @ np.diag(diag) @ Q.T
        lam, V = eig(A0)
        reconstructed = V @ np.diag(lam) @ np.linalg.inv(V)
        sqrt_reconstructed = (
            V @ np.diag(np.sqrt(lam)) @ np.linalg.inv(V)
        )
        np.testing.assert_allclose(
            reconstructed.real, A0, rtol=3e-12, atol=3e-12)
        np.testing.assert_allclose(
            sqrt_reconstructed.real,
            Q @ np.diag(np.sqrt(diag)) @ Q.T,
            rtol=3e-12,
            atol=3e-12,
        )


def test_eig_value_at_scaled_rotated_distinct_spectra():
    """Guard routing must not depend on the absolute matrix scale.

    The former absolute near-diagonal guard sent sufficiently small rotated
    matrices to a coordinate-basis perturbation formula. That shortcut is only
    valid after rotation to a real principal basis.
    """
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)
    angle = 0.71
    K = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    Q = (
        np.eye(3)
        + np.sin(angle) * K
        + (1.0 - np.cos(angle)) * (K @ K)
    )

    for scale in (
            1.0, 1.0e-9, 1.0e-12, 1.0e-13, 1.0e-20, 1.0e-30, 1.0e-31):
        diag = scale * np.array([0.1, 1.0, 10.0])
        A0 = Q @ np.diag(diag) @ Q.T
        lam, V = eig(A0)
        reconstructed = V @ np.diag(lam) @ np.linalg.inv(V)
        sqrt_reconstructed = (
            V @ np.diag(np.sqrt(lam)) @ np.linalg.inv(V)
        )
        exact_sqrt = Q @ np.diag(np.sqrt(diag)) @ Q.T

        value_rel = (
            np.linalg.norm(reconstructed.real - A0)
            / np.linalg.norm(A0)
        )
        sqrt_rel = (
            np.linalg.norm(sqrt_reconstructed.real - exact_sqrt)
            / np.linalg.norm(exact_sqrt)
        )
        eigenvalue_rel = (
            np.linalg.norm(np.sort(lam.real) - diag)
            / np.linalg.norm(diag)
        )
        assert value_rel < 2.0e-12
        assert sqrt_rel < 2.0e-12
        assert eigenvalue_rel < 2.0e-12

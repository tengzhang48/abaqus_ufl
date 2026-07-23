C======================================================================
C     tensor_ops.for — Complex-arithmetic tensor utilities
C
C     Provides det33, inv33, matmul33, transpose33, outer33 for
C     DOUBLE COMPLEX arguments. Also provides real-only wrappers.
C
C     Design:  1:1 mirror of Python det33/inv33 helpers
C     Purpose: Foundation for complex-step tangent engine
C
C     Convention: All 3x3 matrices stored as (3,3) arrays
C                 Index order: (row, col) = (i, J)
C======================================================================

C----------------------------------------------------------------------
C     det33z: Determinant of 3x3 complex matrix
C----------------------------------------------------------------------
      FUNCTION det33z(A) RESULT(d)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: A(3,3)
      DOUBLE COMPLEX :: d

      d = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
     &  - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
     &  + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      RETURN
      END FUNCTION det33z

C----------------------------------------------------------------------
C     det33d: Determinant of 3x3 real matrix
C----------------------------------------------------------------------
      FUNCTION det33d(A) RESULT(d)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: A(3,3)
      DOUBLE PRECISION :: d

      d = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2))
     &  - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1))
     &  + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))

      RETURN
      END FUNCTION det33d

C----------------------------------------------------------------------
C     inv33z: Inverse of 3x3 complex matrix
C----------------------------------------------------------------------
      SUBROUTINE inv33z(A, Ainv)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: Ainv(3,3)
      DOUBLE COMPLEX :: d, det33z

      d = det33z(A)

      Ainv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / d
      Ainv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) / d
      Ainv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / d
      Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / d
      Ainv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / d
      Ainv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3)) / d
      Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / d
      Ainv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) / d
      Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / d

      RETURN
      END SUBROUTINE inv33z

C----------------------------------------------------------------------
C     inv33d: Inverse of 3x3 real matrix
C----------------------------------------------------------------------
      SUBROUTINE inv33d(A, Ainv)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: A(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: Ainv(3,3)
      DOUBLE PRECISION :: d, det33d

      d = det33d(A)

      Ainv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / d
      Ainv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) / d
      Ainv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / d
      Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / d
      Ainv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / d
      Ainv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3)) / d
      Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / d
      Ainv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) / d
      Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / d

      RETURN
      END SUBROUTINE inv33d

C----------------------------------------------------------------------
C     matmul33z: C = A * B for 3x3 complex matrices
C     (Explicit loop — avoids MATMUL intrinsic for portability)
C----------------------------------------------------------------------
      SUBROUTINE matmul33z(A, B, C)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3), B(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: C(3,3)
      INTEGER :: i, j, k

      DO i = 1, 3
        DO j = 1, 3
          C(i,j) = (0.0d0, 0.0d0)
          DO k = 1, 3
            C(i,j) = C(i,j) + A(i,k) * B(k,j)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE matmul33z

C----------------------------------------------------------------------
C     matmul33d: C = A * B for 3x3 real matrices
C----------------------------------------------------------------------
      SUBROUTINE matmul33d(A, B, C)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: A(3,3), B(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: C(3,3)
      INTEGER :: i, j, k

      DO i = 1, 3
        DO j = 1, 3
          C(i,j) = 0.0d0
          DO k = 1, 3
            C(i,j) = C(i,j) + A(i,k) * B(k,j)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE matmul33d

C----------------------------------------------------------------------
C     transpose33z: B = A^T for 3x3 complex matrix
C----------------------------------------------------------------------
      SUBROUTINE transpose33z(A, AT)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: AT(3,3)
      INTEGER :: i, j

      DO i = 1, 3
        DO j = 1, 3
          AT(i,j) = A(j,i)
        END DO
      END DO

      RETURN
      END SUBROUTINE transpose33z

C----------------------------------------------------------------------
C     transpose33d: B = A^T for 3x3 real matrix
C----------------------------------------------------------------------
      SUBROUTINE transpose33d(A, AT)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: A(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: AT(3,3)
      INTEGER :: i, j

      DO i = 1, 3
        DO j = 1, 3
          AT(i,j) = A(j,i)
        END DO
      END DO

      RETURN
      END SUBROUTINE transpose33d

C----------------------------------------------------------------------
C     matvec33z: y = A * x for 3x3 complex matrix and 3-vector
C----------------------------------------------------------------------
      SUBROUTINE matvec33z(A, x, y)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3), x(3)
      DOUBLE COMPLEX, INTENT(OUT) :: y(3)
      INTEGER :: i, k

      DO i = 1, 3
        y(i) = (0.0d0, 0.0d0)
        DO k = 1, 3
          y(i) = y(i) + A(i,k) * x(k)
        END DO
      END DO

      RETURN
      END SUBROUTINE matvec33z

C----------------------------------------------------------------------
C     matvec33d: y = A * x for 3x3 real matrix and 3-vector
C----------------------------------------------------------------------
      SUBROUTINE matvec33d(A, x, y)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: A(3,3), x(3)
      DOUBLE PRECISION, INTENT(OUT) :: y(3)
      INTEGER :: i, k

      DO i = 1, 3
        y(i) = 0.0d0
        DO k = 1, 3
          y(i) = y(i) + A(i,k) * x(k)
        END DO
      END DO

      RETURN
      END SUBROUTINE matvec33d

C----------------------------------------------------------------------
C     eye33z: A = I (3x3 complex identity)
C----------------------------------------------------------------------
      SUBROUTINE eye33z(A)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(OUT) :: A(3,3)
      INTEGER :: i, j

      DO i = 1, 3
        DO j = 1, 3
          IF (i .EQ. j) THEN
            A(i,j) = (1.0d0, 0.0d0)
          ELSE
            A(i,j) = (0.0d0, 0.0d0)
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE eye33z

C----------------------------------------------------------------------
C     eye33d: A = I (3x3 real identity)
C----------------------------------------------------------------------
      SUBROUTINE eye33d(A)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: A(3,3)
      INTEGER :: i, j

      DO i = 1, 3
        DO j = 1, 3
          IF (i .EQ. j) THEN
            A(i,j) = 1.0d0
          ELSE
            A(i,j) = 0.0d0
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE eye33d

C----------------------------------------------------------------------
C     real2complex33: Convert real 3x3 to complex 3x3
C----------------------------------------------------------------------
      SUBROUTINE real2complex33(A_real, A_cmplx)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: A_real(3,3)
      DOUBLE COMPLEX,   INTENT(OUT) :: A_cmplx(3,3)
      INTEGER :: i, j

      DO i = 1, 3
        DO j = 1, 3
          A_cmplx(i,j) = DCMPLX(A_real(i,j), 0.0d0)
        END DO
      END DO

      RETURN
      END SUBROUTINE real2complex33

C----------------------------------------------------------------------
C     real2complex3: Convert real 3-vector to complex 3-vector
C----------------------------------------------------------------------
      SUBROUTINE real2complex3(v_real, v_cmplx)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: v_real(3)
      DOUBLE COMPLEX,   INTENT(OUT) :: v_cmplx(3)
      INTEGER :: i

      DO i = 1, 3
        v_cmplx(i) = DCMPLX(v_real(i), 0.0d0)
      END DO

      RETURN
      END SUBROUTINE real2complex3

C======================================================================
C     Matrix functions for tensor_ops.for
C     Eigendecomposition-based (Approach A) + iterative fallbacks (B)
C
C     All DOUBLE COMPLEX, CS-safe.
C
C     Subroutines:
C       cross3z        — cross product of complex 3-vectors
C       outer33z       — outer product of two complex 3-vectors
C       sym3z          — symmetric 3x3 tensor from 6 unique entries
C       eig33z         — eigenvalues/vectors of 3x3 symmetric matrix
C       sqrtm33z       — matrix square root via eigendecomposition
C       logm33z        — matrix logarithm via eigendecomposition
C       expm33z        — matrix exponential via eigendecomposition
C       polar33z       — right polar decomposition F = R * U
C       sqrtm33z_iter  — matrix square root (Denman-Beavers)
C       logm33z_iter   — matrix logarithm (inv scaling & squaring)
C       expm33z_iter   — matrix exponential (scaling & squaring + Pade)
C======================================================================

C----------------------------------------------------------------------
C     cross3z: Cross product of two complex 3-vectors
C----------------------------------------------------------------------
      SUBROUTINE cross3z(a, b, c)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: a(3), b(3)
      DOUBLE COMPLEX, INTENT(OUT) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      RETURN
      END SUBROUTINE cross3z

C----------------------------------------------------------------------
C     outer33z: Outer product of two complex 3-vectors, C_ij = a_i b_j
C----------------------------------------------------------------------
      SUBROUTINE outer33z(a, b, C)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: a(3), b(3)
      DOUBLE COMPLEX, INTENT(OUT) :: C(3,3)
      INTEGER :: i, j

      DO i = 1, 3
        DO j = 1, 3
          C(i,j) = a(i) * b(j)
        END DO
      END DO

      RETURN
      END SUBROUTINE outer33z

C----------------------------------------------------------------------
C     sym3z: Symmetric 3x3 tensor from its 6 unique entries
C            T = [[a11, a12, a13],
C                 [a12, a22, a23],
C                 [a13, a23, a33]]
C----------------------------------------------------------------------
      SUBROUTINE sym3z(a11, a22, a33, a12, a13, a23, T)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: a11, a22, a33, a12, a13, a23
      DOUBLE COMPLEX, INTENT(OUT) :: T(3,3)

      T(1,1) = a11
      T(2,2) = a22
      T(3,3) = a33
      T(1,2) = a12
      T(2,1) = a12
      T(1,3) = a13
      T(3,1) = a13
      T(2,3) = a23
      T(3,2) = a23

      RETURN
      END SUBROUTINE sym3z

C----------------------------------------------------------------------
C     eig33z: Eigenvalues and eigenvectors of 3x3 symmetric matrix
C
C     Uses trigonometric solution for the depressed cubic.
C     Eigenvectors are UNNORMALIZED (raw cross products) for CS
C     correctness. Use inv(V) not V^T for reconstruction.
C
C     KNOWN LIMITATION: the degeneracy
C     guards below return lam=q, V=I and DISCARD complex-step
C     perturbations at (near-)diagonal/identity states, silently
C     zeroing CS derivatives of spectral functions there (C = F^T F
C     is exactly I at the first increment). tensor.py eig() was fixed
C     2026-06-12 with a perturbation-theory fallback; this Fortran
C     mirror has NOT been ported (needs a real-symmetric 3x3 Jacobi
C     solver). The default matrix_backend='iterative' does not call
C     eig33z; do not switch materials to backend='eig' until ported.
C----------------------------------------------------------------------
      SUBROUTINE eig33z(A, lam, V)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: lam(3), V(3,3)

C     --- Local variables ---
      DOUBLE COMPLEX :: p1, q, p2, p, r, phi
      DOUBLE COMPLEX :: B(3,3), C(3,3), M(3,3)
      DOUBLE COMPLEX :: zi, zacos, ztmp
      DOUBLE COMPLEX :: det33z
      DOUBLE COMPLEX :: r0(3), r1(3), r2(3)
      DOUBLE COMPLEX :: c1(3), c2(3), c3(3)
      DOUBLE COMPLEX :: nsq1, nsq2, nsq3, best_nsq
      DOUBLE PRECISION, PARAMETER ::
     &  PI = 3.14159265358979323846d0
      DOUBLE PRECISION :: twopi3
      INTEGER :: i, j, k, best_idx

      zi = DCMPLX(0.0d0, 1.0d0)
      twopi3 = 2.0d0 * PI / 3.0d0

C     --- Invariants ---
      p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
      q  = (A(1,1) + A(2,2) + A(3,3)) / 3.0d0

C     B = A - q*I
      DO i = 1, 3
        DO j = 1, 3
          B(i,j) = A(i,j)
        END DO
      END DO
      B(1,1) = B(1,1) - q
      B(2,2) = B(2,2) - q
      B(3,3) = B(3,3) - q

      p2 = B(1,1)**2 + B(2,2)**2 + B(3,3)**2
     &   + 2.0d0 * p1
      p  = SQRT(p2 / 6.0d0)

C     Guard: A is a multiple of identity
      IF (ABS(p) .LT. 1.0d-30) THEN
        lam(1) = q
        lam(2) = q
        lam(3) = q
        CALL eye33z(V)
        RETURN
      END IF

C     Nearly-diagonal guard (matches Python tensor.eig)
      IF (ABS(DBLE(p2)) .GT. 0.0d0 .AND.
     &    (ABS(DBLE(p1)) .LT. 1.0d-14 * ABS(DBLE(p2)) .OR.
     &     ABS(p1) .LT. 1.0d-18)) THEN
        lam(1) = A(1,1)
        lam(2) = A(2,2)
        lam(3) = A(3,3)
        CALL eye33z(V)
        DO i = 1, 2
          DO j = i+1, 3
            IF (DBLE(lam(j)) .LT. DBLE(lam(i))) THEN
              ztmp = lam(i)
              lam(i) = lam(j)
              lam(j) = ztmp
              DO k = 1, 3
                ztmp = V(k,i)
                V(k,i) = V(k,j)
                V(k,j) = ztmp
              END DO
            END IF
          END DO
        END DO
        RETURN
      END IF

C     C = B / p
      DO i = 1, 3
        DO j = 1, 3
          C(i,j) = B(i,j) / p
        END DO
      END DO

C     r = det(C) / 2
      r = det33z(C) / 2.0d0

C     Clamp real part before complex arccos, matching tensor.py eig().
C     Near repeated eigenvalues, roundoff can push DBLE(r) just outside
C     [-1, 1] and poison complex-step tangents with NaNs.
      r = DCMPLX(MAX(-1.0d0, MIN(1.0d0, DBLE(r))), AIMAG(r))

C     Complex arccos: zacos = -i * log(r + i*sqrt(1 - r^2))
      zacos = -zi * LOG(r + zi * SQRT(
     &  DCMPLX(1.0d0, 0.0d0) - r*r))
      phi = zacos / 3.0d0

C     Three eigenvalues
      lam(1) = q + 2.0d0*p*COS(phi)
      lam(2) = q + 2.0d0*p*COS(phi
     &       + DCMPLX(twopi3, 0.0d0))
      lam(3) = q + 2.0d0*p*COS(phi
     &       + DCMPLX(2.0d0*twopi3, 0.0d0))

C     Sort by real part (ascending, 3-element bubble sort)
      DO i = 1, 2
        DO j = i+1, 3
          IF (DBLE(lam(j)) .LT. DBLE(lam(i))) THEN
            ztmp = lam(i)
            lam(i) = lam(j)
            lam(j) = ztmp
          END IF
        END DO
      END DO

C     --- Eigenvectors via cross products of rows of (A - lam*I) ---
      DO k = 1, 3
C       M = A - lam(k) * I
        DO i = 1, 3
          DO j = 1, 3
            M(i,j) = A(i,j)
          END DO
        END DO
        M(1,1) = M(1,1) - lam(k)
        M(2,2) = M(2,2) - lam(k)
        M(3,3) = M(3,3) - lam(k)

C       Extract rows
        DO j = 1, 3
          r0(j) = M(1,j)
          r1(j) = M(2,j)
          r2(j) = M(3,j)
        END DO

C       Three candidate cross products
        CALL cross3z(r0, r1, c1)
        CALL cross3z(r0, r2, c2)
        CALL cross3z(r1, r2, c3)

C       Squared norms (complex)
        nsq1 = c1(1)*c1(1) + c1(2)*c1(2) + c1(3)*c1(3)
        nsq2 = c2(1)*c2(1) + c2(2)*c2(2) + c2(3)*c2(3)
        nsq3 = c3(1)*c3(1) + c3(2)*c3(2) + c3(3)*c3(3)

C       Pick best (largest |nsq|)
        best_idx = 1
        best_nsq = nsq1
        IF (ABS(nsq2) .GT. ABS(best_nsq)) THEN
          best_idx = 2
          best_nsq = nsq2
        END IF
        IF (ABS(nsq3) .GT. ABS(best_nsq)) THEN
          best_idx = 3
          best_nsq = nsq3
        END IF

C       Store eigenvector (unnormalized for CS correctness)
        IF (ABS(best_nsq) .LT. 1.0d-50) THEN
C         Degenerate: coordinate basis fallback
          V(1,k) = DCMPLX(0.0d0, 0.0d0)
          V(2,k) = DCMPLX(0.0d0, 0.0d0)
          V(3,k) = DCMPLX(0.0d0, 0.0d0)
          V(k,k) = DCMPLX(1.0d0, 0.0d0)
        ELSE IF (best_idx .EQ. 1) THEN
          V(1,k) = c1(1)
          V(2,k) = c1(2)
          V(3,k) = c1(3)
        ELSE IF (best_idx .EQ. 2) THEN
          V(1,k) = c2(1)
          V(2,k) = c2(2)
          V(3,k) = c2(3)
        ELSE
          V(1,k) = c3(1)
          V(2,k) = c3(2)
          V(3,k) = c3(3)
        END IF
      END DO

      RETURN
      END SUBROUTINE eig33z

C----------------------------------------------------------------------
C     sqrtm33z: Matrix square root via eigendecomposition
C
C     sqrtm(A) = V * diag(sqrt(lam)) * inv(V)
C     Uses inv(V) (not V^T) for CS correctness.
C----------------------------------------------------------------------
      SUBROUTINE sqrtm33z(A, S)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: S(3,3)

      DOUBLE COMPLEX :: lam(3), V(3,3), Vinv(3,3)
      DOUBLE COMPLEX :: D(3,3), T(3,3)
      INTEGER :: i, j

      CALL eig33z(A, lam, V)
      CALL inv33z(V, Vinv)

C     D = diag(sqrt(lam))
      DO i = 1, 3
        DO j = 1, 3
          D(i,j) = DCMPLX(0.0d0, 0.0d0)
        END DO
      END DO
      D(1,1) = SQRT(lam(1))
      D(2,2) = SQRT(lam(2))
      D(3,3) = SQRT(lam(3))

C     S = V * D * Vinv
      CALL matmul33z(V, D, T)
      CALL matmul33z(T, Vinv, S)

      RETURN
      END SUBROUTINE sqrtm33z

C----------------------------------------------------------------------
C     logm33z: Matrix logarithm via eigendecomposition
C
C     logm(A) = V * diag(log(lam)) * inv(V)
C     For Hencky strain: E = 0.5 * logm(C)
C----------------------------------------------------------------------
      SUBROUTINE logm33z(A, L)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: L(3,3)

      DOUBLE COMPLEX :: lam(3), V(3,3), Vinv(3,3)
      DOUBLE COMPLEX :: D(3,3), T(3,3)
      INTEGER :: i, j

      CALL eig33z(A, lam, V)
      CALL inv33z(V, Vinv)

C     D = diag(log(lam))
      DO i = 1, 3
        DO j = 1, 3
          D(i,j) = DCMPLX(0.0d0, 0.0d0)
        END DO
      END DO
      D(1,1) = LOG(lam(1))
      D(2,2) = LOG(lam(2))
      D(3,3) = LOG(lam(3))

C     L = V * D * Vinv
      CALL matmul33z(V, D, T)
      CALL matmul33z(T, Vinv, L)

      RETURN
      END SUBROUTINE logm33z

C----------------------------------------------------------------------
C     expm33z: Matrix exponential via eigendecomposition
C
C     expm(A) = V * diag(exp(lam)) * inv(V)
C     For exponential map: Fp_new = expm(dt*Dp) * Fp_old
C----------------------------------------------------------------------
      SUBROUTINE expm33z(A, E)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: E(3,3)

      DOUBLE COMPLEX :: lam(3), V(3,3), Vinv(3,3)
      DOUBLE COMPLEX :: D(3,3), T(3,3)
      INTEGER :: i, j

      CALL eig33z(A, lam, V)
      CALL inv33z(V, Vinv)

C     D = diag(exp(lam))
      DO i = 1, 3
        DO j = 1, 3
          D(i,j) = DCMPLX(0.0d0, 0.0d0)
        END DO
      END DO
      D(1,1) = EXP(lam(1))
      D(2,2) = EXP(lam(2))
      D(3,3) = EXP(lam(3))

C     E = V * D * Vinv
      CALL matmul33z(V, D, T)
      CALL matmul33z(T, Vinv, E)

      RETURN
      END SUBROUTINE expm33z

C----------------------------------------------------------------------
C     polar33z: Right polar decomposition F = R * U
C
C     U = sqrtm(F^T * F),  R = F * inv(U)
C----------------------------------------------------------------------
      SUBROUTINE polar33z(F, R, U)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: F(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: R(3,3), U(3,3)

      DOUBLE COMPLEX :: FT(3,3), C(3,3), Uinv(3,3)

C     C = F^T * F
      CALL transpose33z(F, FT)
      CALL matmul33z(FT, F, C)

C     U = sqrtm(C)
      CALL sqrtm33z(C, U)

C     R = F * inv(U)
      CALL inv33z(U, Uinv)
      CALL matmul33z(F, Uinv, R)

      RETURN
      END SUBROUTINE polar33z

C----------------------------------------------------------------------
C     polar33z_iter: Right polar decomposition F = R * U
C
C     Iterative backend: uses sqrtm33z_iter for CS safety.
C----------------------------------------------------------------------
      SUBROUTINE polar33z_iter(F, R, U)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: F(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: R(3,3), U(3,3)

      DOUBLE COMPLEX :: FT(3,3), C(3,3), Uinv(3,3)

C     C = F^T * F
      CALL transpose33z(F, FT)
      CALL matmul33z(FT, F, C)

C     U = sqrtm(C),  Uinv = inv(sqrtm(C))
      CALL sqrtm33z_iter(C, U, Uinv)

C     R = F * inv(U)
      CALL matmul33z(F, Uinv, R)

      RETURN
      END SUBROUTINE polar33z_iter

C----------------------------------------------------------------------
C     sqrtm33z_iter: Matrix square root via Denman-Beavers iteration
C
C     Converges: Y -> sqrtm(A),  Z -> inv(sqrtm(A))
C     Fixed iteration count for CS safety (no branching).
C----------------------------------------------------------------------
      SUBROUTINE sqrtm33z_iter(A, S, Sinv)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: S(3,3), Sinv(3,3)

      DOUBLE COMPLEX :: Y(3,3), Z(3,3)
      DOUBLE COMPLEX :: Yinv(3,3), Zinv(3,3)
      DOUBLE COMPLEX :: half
      INTEGER, PARAMETER :: MAXITER = 20
      INTEGER :: iter, i, j

      half = DCMPLX(0.5d0, 0.0d0)

C     Y = A, Z = I
      DO i = 1, 3
        DO j = 1, 3
          Y(i,j) = A(i,j)
        END DO
      END DO
      CALL eye33z(Z)

      DO iter = 1, MAXITER
        CALL inv33z(Z, Zinv)
        CALL inv33z(Y, Yinv)
        DO i = 1, 3
          DO j = 1, 3
            Y(i,j) = half * (Y(i,j) + Zinv(i,j))
            Z(i,j) = half * (Z(i,j) + Yinv(i,j))
          END DO
        END DO
      END DO

      DO i = 1, 3
        DO j = 1, 3
          S(i,j)    = Y(i,j)
          Sinv(i,j) = Z(i,j)
        END DO
      END DO

      RETURN
      END SUBROUTINE sqrtm33z_iter

C----------------------------------------------------------------------
C     logm33z_iter: Matrix logarithm via inverse scaling & squaring
C
C     Step 1: B = A^{1/2^NSCALE}  (repeated square roots)
C     Step 2: Gregory series for log(B) where B ~ I
C     Step 3: log(A) = 2^NSCALE * log(B)
C     Fixed counts for CS safety.
C----------------------------------------------------------------------
      SUBROUTINE logm33z_iter(A, L)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: L(3,3)

      DOUBLE COMPLEX :: B(3,3), Bsqrt(3,3), Binvsqrt(3,3)
      DOUBLE COMPLEX :: BpI(3,3), BmI(3,3), BpI_inv(3,3)
      DOUBLE COMPLEX :: X(3,3), X2(3,3), Xpow(3,3)
      DOUBLE COMPLEX :: T1(3,3), II(3,3)
      DOUBLE COMPLEX :: scl
      INTEGER, PARAMETER :: NSCALE = 10
      INTEGER, PARAMETER :: NSERIES = 8
      INTEGER :: ks, i, j

C     Step 1: Repeated square roots — B = A^{1/2^NSCALE}
      DO i = 1, 3
        DO j = 1, 3
          B(i,j) = A(i,j)
        END DO
      END DO
      DO ks = 1, NSCALE
        CALL sqrtm33z_iter(B, Bsqrt, Binvsqrt)
        DO i = 1, 3
          DO j = 1, 3
            B(i,j) = Bsqrt(i,j)
          END DO
        END DO
      END DO

C     Step 2: Gregory series for log(B) where B ~ I
C     X = (B - I) * inv(B + I)
      CALL eye33z(II)
      DO i = 1, 3
        DO j = 1, 3
          BpI(i,j) = B(i,j) + II(i,j)
          BmI(i,j) = B(i,j) - II(i,j)
        END DO
      END DO
      CALL inv33z(BpI, BpI_inv)
      CALL matmul33z(BmI, BpI_inv, X)

C     X2 = X * X
      CALL matmul33z(X, X, X2)

C     L = X,  Xpow = X  (start series)
      DO i = 1, 3
        DO j = 1, 3
          L(i,j)    = X(i,j)
          Xpow(i,j) = X(i,j)
        END DO
      END DO

C     L += Xpow * X2 / k  for k = 3, 5, 7, ..., 2*NSERIES+1
      DO ks = 3, 2*NSERIES+1, 2
        CALL matmul33z(Xpow, X2, T1)
        DO i = 1, 3
          DO j = 1, 3
            Xpow(i,j) = T1(i,j)
            L(i,j) = L(i,j)
     &        + Xpow(i,j) / DCMPLX(DBLE(ks), 0.0d0)
          END DO
        END DO
      END DO

C     Step 3: Undo scaling — L = 2 * L * 2^NSCALE
      scl = DCMPLX(DBLE(2**(NSCALE+1)), 0.0d0)
      DO i = 1, 3
        DO j = 1, 3
          L(i,j) = L(i,j) * scl
        END DO
      END DO

      RETURN
      END SUBROUTINE logm33z_iter

C----------------------------------------------------------------------
C     expm33z_iter: Matrix exponential via scaling & squaring + Pade
C
C     Step 1: B = A / 2^NSCALE
C     Step 2: [2/2] Pade approximant for exp(B)
C     Step 3: E = exp(B)^{2^NSCALE}  (repeated squaring)
C     Fixed counts for CS safety.
C----------------------------------------------------------------------
      SUBROUTINE expm33z_iter(A, E)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(OUT) :: E(3,3)

      DOUBLE COMPLEX :: B(3,3), B2(3,3), II(3,3)
      DOUBLE COMPLEX :: N2(3,3), D2(3,3), D2inv(3,3)
      DOUBLE COMPLEX :: T1(3,3)
      DOUBLE COMPLEX :: scl, half, twelfth
      INTEGER, PARAMETER :: NSCALE = 10
      INTEGER :: ks, i, j

      half     = DCMPLX(0.5d0, 0.0d0)
      twelfth  = DCMPLX(1.0d0/12.0d0, 0.0d0)
      scl      = DCMPLX(DBLE(2**NSCALE), 0.0d0)

C     Step 1: B = A / 2^NSCALE
      DO i = 1, 3
        DO j = 1, 3
          B(i,j) = A(i,j) / scl
        END DO
      END DO

C     B2 = B * B
      CALL matmul33z(B, B, B2)

C     [2/2] Pade: exp(B) ~ inv(I - B/2 + B^2/12) * (I + B/2 + B^2/12)
      CALL eye33z(II)
      DO i = 1, 3
        DO j = 1, 3
          N2(i,j) = II(i,j) + half*B(i,j)
     &            + twelfth*B2(i,j)
          D2(i,j) = II(i,j) - half*B(i,j)
     &            + twelfth*B2(i,j)
        END DO
      END DO
      CALL inv33z(D2, D2inv)
      CALL matmul33z(D2inv, N2, E)

C     Step 3: Repeated squaring — E = E^{2^NSCALE}
      DO ks = 1, NSCALE
        CALL matmul33z(E, E, T1)
        DO i = 1, 3
          DO j = 1, 3
            E(i,j) = T1(i,j)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE expm33z_iter

C======================================================================
C     CS-safe scalar trig/hyperbolic functions (DOUBLE COMPLEX)
C     Exact complex formulas — safe for complex-step derivatives
C======================================================================

C----------------------------------------------------------------------
C     cs_sin: complex sine  (analytic continuation of DSIN)
C----------------------------------------------------------------------
      FUNCTION cs_sin(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y
      x = DBLE(z)
      y = AIMAG(z)
      w = DCMPLX(DSIN(x)*DCOSH(y), DCOS(x)*DSINH(y))
      RETURN
      END FUNCTION cs_sin

C----------------------------------------------------------------------
C     cs_cos: complex cosine  (analytic continuation of DCOS)
C----------------------------------------------------------------------
      FUNCTION cs_cos(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y
      x = DBLE(z)
      y = AIMAG(z)
      w = DCMPLX(DCOS(x)*DCOSH(y), -DSIN(x)*DSINH(y))
      RETURN
      END FUNCTION cs_cos

C----------------------------------------------------------------------
C     cs_tan: complex tangent  (analytic continuation of DTAN)
C----------------------------------------------------------------------
      FUNCTION cs_tan(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y, den
      x = DBLE(z)
      y = AIMAG(z)
      den = DCOS(2.0d0*x) + DCOSH(2.0d0*y)
      w = DCMPLX(DSIN(2.0d0*x)/den, DSINH(2.0d0*y)/den)
      RETURN
      END FUNCTION cs_tan

C----------------------------------------------------------------------
C     cs_sinh: complex hyperbolic sine
C----------------------------------------------------------------------
      FUNCTION cs_sinh(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y
      x = DBLE(z)
      y = AIMAG(z)
      w = DCMPLX(DSINH(x)*DCOS(y), DCOSH(x)*DSIN(y))
      RETURN
      END FUNCTION cs_sinh

C----------------------------------------------------------------------
C     cs_cosh: complex hyperbolic cosine
C----------------------------------------------------------------------
      FUNCTION cs_cosh(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y
      x = DBLE(z)
      y = AIMAG(z)
      w = DCMPLX(DCOSH(x)*DCOS(y), DSINH(x)*DSIN(y))
      RETURN
      END FUNCTION cs_cosh

C----------------------------------------------------------------------
C     cs_tanh: complex hyperbolic tangent
C----------------------------------------------------------------------
      FUNCTION cs_tanh(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y, den
      x = DBLE(z)
      y = AIMAG(z)
      den = DCOSH(2.0d0*x) + DCOS(2.0d0*y)
      w = DCMPLX(DSINH(2.0d0*x)/den, DSIN(2.0d0*y)/den)
      RETURN
      END FUNCTION cs_tanh

C----------------------------------------------------------------------
C     cs_erf: complex error function (CS-safe first-order approx)
C----------------------------------------------------------------------
      FUNCTION cs_erf(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE PRECISION :: x, y
      DOUBLE PRECISION, PARAMETER :: two_over_sqrt_pi =
     &  1.1283791670955126d0
      x = DBLE(z)
      y = AIMAG(z)
      w = DCMPLX(ERF(x), two_over_sqrt_pi * EXP(-x*x) * y)
      RETURN
      END FUNCTION cs_erf

C----------------------------------------------------------------------
C     cs_erfc: complex complementary error function
C----------------------------------------------------------------------
      FUNCTION cs_erfc(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      DOUBLE COMPLEX :: cs_erf
      w = DCMPLX(1.0d0, 0.0d0) - cs_erf(z)
      RETURN
      END FUNCTION cs_erfc

C----------------------------------------------------------------------
C     cs_atan: complex arctangent (analytic continuation of DATAN)
C----------------------------------------------------------------------
      FUNCTION cs_atan(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w, zi
      zi = DCMPLX(0.0D0, 1.0D0)
      w = (LOG(1.0D0 + zi*z) - LOG(1.0D0 - zi*z)) / (2.0D0*zi)
      RETURN
      END FUNCTION cs_atan

C----------------------------------------------------------------------
C     cs_acos: complex arccos  (analytic continuation of DACOS)
C----------------------------------------------------------------------
      FUNCTION cs_acos(z) RESULT(w)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN) :: z
      DOUBLE COMPLEX :: w
      w = -DCMPLX(0.0D0, 1.0D0) * LOG(z + DCMPLX(0.0D0, 1.0D0)
     &    * SQRT(1.0D0 - z**2))
      RETURN
      END FUNCTION cs_acos

C----------------------------------------------------------------------
C     cs_cubic_roots: all roots of a*x**3 + b*x**2 + c*x + d = 0
C     Depressed-cubic + Cardano, all in DOUBLE COMPLEX (CS-safe).
C----------------------------------------------------------------------
      SUBROUTINE cs_cubic_roots(a, b, c, d, r1, r2, r3)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: a, b, c, d
      DOUBLE COMPLEX, INTENT(OUT) :: r1, r2, r3
      DOUBLE COMPLEX :: p, q, Delta, sqrt_Delta, u, v, shift
      DOUBLE COMPLEX :: t1, t2, t3
      DOUBLE COMPLEX :: half, third, one, two, three, nine
      DOUBLE COMPLEX :: twentyseven

      half  = DCMPLX(0.5d0, 0.0d0)
      third = DCMPLX(1.0d0/3.0d0, 0.0d0)
      one   = DCMPLX(1.0d0, 0.0d0)
      two   = DCMPLX(2.0d0, 0.0d0)
      three = DCMPLX(3.0d0, 0.0d0)
      nine  = DCMPLX(9.0d0, 0.0d0)
      twentyseven = DCMPLX(27.0d0, 0.0d0)

C     Depress the cubic: x = t - b/(3a)
C     t**3 + p*t + q = 0
      p = (three*a*c - b*b) / (three*a*a)
      q = (two*b**3 - nine*a*b*c + twentyseven*a*a*d)
     &    / (twentyseven*a**3)

C     Discriminant
      Delta = (q/two)**2 + (p/three)**3

C     Cardano terms
      sqrt_Delta = SQRT(Delta)
      u = (-q/two + sqrt_Delta)**third
C     v chosen so that u*v = -p/3 (avoids branch-cut ambiguity).
C     Guard u=0 (triple-root / p=0 case) by computing v independently.
      IF (ABS(u) .GT. 1.0d-30) THEN
        v = -p / (three*u)
      ELSE
        v = (-q/two - sqrt_Delta)**third
      END IF

C     Three roots in t
      t1 = u + v
      t2 = -half*(u+v) + DCMPLX(0.0d0, DSQRT(3.0d0)/2.0d0)*(u-v)
      t3 = -half*(u+v) - DCMPLX(0.0d0, DSQRT(3.0d0)/2.0d0)*(u-v)

C     Shift back
      shift = b / (three*a)
      r1 = t1 - shift
      r2 = t2 - shift
      r3 = t3 - shift

      RETURN
      END SUBROUTINE cs_cubic_roots

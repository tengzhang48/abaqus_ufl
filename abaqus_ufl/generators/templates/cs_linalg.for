C======================================================================
C     cs_linalg.for -- CS-safe small dense linear solve
C
C     Inline LU with partial pivoting for N in {2,3,4,6,9}.
C     No LAPACK dependency -- self-contained for Abaqus UMAT.
C======================================================================

C----------------------------------------------------------------------
C     cs_solve_2: solve A(2,2) * x(2) = b(2)
C----------------------------------------------------------------------
      SUBROUTINE cs_solve_2(A, b, x)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(2,2)
      DOUBLE COMPLEX, INTENT(IN)  :: b(2)
      DOUBLE COMPLEX, INTENT(OUT) :: x(2)
      DOUBLE COMPLEX :: LU(2,2), y(2), temp
      INTEGER :: perm(2)
      INTEGER :: i, j, k, max_row
      DOUBLE PRECISION :: max_val, abs_val

      DO i = 1, 2
        perm(i) = i
        DO j = 1, 2
          LU(i,j) = A(i,j)
        END DO
      END DO

      DO k = 1, 2
        max_val = ABS(LU(k,k))
        max_row = k
        DO i = k+1, 2
          abs_val = ABS(LU(i,k))
          IF (abs_val > max_val) THEN
            max_val = abs_val
            max_row = i
          END IF
        END DO
        IF (max_row /= k) THEN
          DO j = 1, 2
            temp = LU(k,j)
            LU(k,j) = LU(max_row,j)
            LU(max_row,j) = temp
          END DO
          j = perm(k)
          perm(k) = perm(max_row)
          perm(max_row) = j
        END IF
        DO i = k+1, 2
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 2
            LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
          END DO
        END DO
      END DO

      DO i = 1, 2
        y(i) = b(perm(i))
        DO j = 1, i-1
          y(i) = y(i) - LU(i,j) * y(j)
        END DO
      END DO

      DO i = 2, 1, -1
        x(i) = y(i)
        DO j = i+1, 2
          x(i) = x(i) - LU(i,j) * x(j)
        END DO
        x(i) = x(i) / LU(i,i)
      END DO
      RETURN
      END SUBROUTINE cs_solve_2

C----------------------------------------------------------------------
C     cs_solve_3: solve A(3,3) * x(3) = b(3)
C----------------------------------------------------------------------
      SUBROUTINE cs_solve_3(A, b, x)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(3,3)
      DOUBLE COMPLEX, INTENT(IN)  :: b(3)
      DOUBLE COMPLEX, INTENT(OUT) :: x(3)
      DOUBLE COMPLEX :: LU(3,3), y(3), temp
      INTEGER :: perm(3)
      INTEGER :: i, j, k, max_row
      DOUBLE PRECISION :: max_val, abs_val

      DO i = 1, 3
        perm(i) = i
        DO j = 1, 3
          LU(i,j) = A(i,j)
        END DO
      END DO

      DO k = 1, 3
        max_val = ABS(LU(k,k))
        max_row = k
        DO i = k+1, 3
          abs_val = ABS(LU(i,k))
          IF (abs_val > max_val) THEN
            max_val = abs_val
            max_row = i
          END IF
        END DO
        IF (max_row /= k) THEN
          DO j = 1, 3
            temp = LU(k,j)
            LU(k,j) = LU(max_row,j)
            LU(max_row,j) = temp
          END DO
          j = perm(k)
          perm(k) = perm(max_row)
          perm(max_row) = j
        END IF
        DO i = k+1, 3
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 3
            LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
          END DO
        END DO
      END DO

      DO i = 1, 3
        y(i) = b(perm(i))
        DO j = 1, i-1
          y(i) = y(i) - LU(i,j) * y(j)
        END DO
      END DO

      DO i = 3, 1, -1
        x(i) = y(i)
        DO j = i+1, 3
          x(i) = x(i) - LU(i,j) * x(j)
        END DO
        x(i) = x(i) / LU(i,i)
      END DO
      RETURN
      END SUBROUTINE cs_solve_3

C----------------------------------------------------------------------
C     cs_solve_4: solve A(4,4) * x(4) = b(4)
C----------------------------------------------------------------------
      SUBROUTINE cs_solve_4(A, b, x)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(4,4)
      DOUBLE COMPLEX, INTENT(IN)  :: b(4)
      DOUBLE COMPLEX, INTENT(OUT) :: x(4)
      DOUBLE COMPLEX :: LU(4,4), y(4), temp
      INTEGER :: perm(4)
      INTEGER :: i, j, k, max_row
      DOUBLE PRECISION :: max_val, abs_val

      DO i = 1, 4
        perm(i) = i
        DO j = 1, 4
          LU(i,j) = A(i,j)
        END DO
      END DO

      DO k = 1, 4
        max_val = ABS(LU(k,k))
        max_row = k
        DO i = k+1, 4
          abs_val = ABS(LU(i,k))
          IF (abs_val > max_val) THEN
            max_val = abs_val
            max_row = i
          END IF
        END DO
        IF (max_row /= k) THEN
          DO j = 1, 4
            temp = LU(k,j)
            LU(k,j) = LU(max_row,j)
            LU(max_row,j) = temp
          END DO
          j = perm(k)
          perm(k) = perm(max_row)
          perm(max_row) = j
        END IF
        DO i = k+1, 4
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 4
            LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
          END DO
        END DO
      END DO

      DO i = 1, 4
        y(i) = b(perm(i))
        DO j = 1, i-1
          y(i) = y(i) - LU(i,j) * y(j)
        END DO
      END DO

      DO i = 4, 1, -1
        x(i) = y(i)
        DO j = i+1, 4
          x(i) = x(i) - LU(i,j) * x(j)
        END DO
        x(i) = x(i) / LU(i,i)
      END DO
      RETURN
      END SUBROUTINE cs_solve_4

C----------------------------------------------------------------------
C     cs_solve_5: solve A(5,5) * x(5) = b(5)
C----------------------------------------------------------------------
      SUBROUTINE cs_solve_5(A, b, x)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(5,5)
      DOUBLE COMPLEX, INTENT(IN)  :: b(5)
      DOUBLE COMPLEX, INTENT(OUT) :: x(5)
      DOUBLE COMPLEX :: LU(5,5), y(5), temp
      INTEGER :: perm(5)
      INTEGER :: i, j, k, max_row
      DOUBLE PRECISION :: max_val, abs_val

      DO i = 1, 5
        perm(i) = i
        DO j = 1, 5
          LU(i,j) = A(i,j)
        END DO
      END DO

      DO k = 1, 5
        max_val = ABS(LU(k,k))
        max_row = k
        DO i = k+1, 5
          abs_val = ABS(LU(i,k))
          IF (abs_val > max_val) THEN
            max_val = abs_val
            max_row = i
          END IF
        END DO
        IF (max_row /= k) THEN
          DO j = 1, 5
            temp = LU(k,j)
            LU(k,j) = LU(max_row,j)
            LU(max_row,j) = temp
          END DO
          j = perm(k)
          perm(k) = perm(max_row)
          perm(max_row) = j
        END IF
        DO i = k+1, 5
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 5
            LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
          END DO
        END DO
      END DO

      DO i = 1, 5
        y(i) = b(perm(i))
        DO j = 1, i-1
          y(i) = y(i) - LU(i,j) * y(j)
        END DO
      END DO

      DO i = 5, 1, -1
        x(i) = y(i)
        DO j = i+1, 5
          x(i) = x(i) - LU(i,j) * x(j)
        END DO
        x(i) = x(i) / LU(i,i)
      END DO
      RETURN
      END SUBROUTINE cs_solve_5

C----------------------------------------------------------------------
C     cs_solve_6: solve A(6,6) * x(6) = b(6)
C----------------------------------------------------------------------
      SUBROUTINE cs_solve_6(A, b, x)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(6,6)
      DOUBLE COMPLEX, INTENT(IN)  :: b(6)
      DOUBLE COMPLEX, INTENT(OUT) :: x(6)
      DOUBLE COMPLEX :: LU(6,6), y(6), temp
      INTEGER :: perm(6)
      INTEGER :: i, j, k, max_row
      DOUBLE PRECISION :: max_val, abs_val

      DO i = 1, 6
        perm(i) = i
        DO j = 1, 6
          LU(i,j) = A(i,j)
        END DO
      END DO

      DO k = 1, 6
        max_val = ABS(LU(k,k))
        max_row = k
        DO i = k+1, 6
          abs_val = ABS(LU(i,k))
          IF (abs_val > max_val) THEN
            max_val = abs_val
            max_row = i
          END IF
        END DO
        IF (max_row /= k) THEN
          DO j = 1, 6
            temp = LU(k,j)
            LU(k,j) = LU(max_row,j)
            LU(max_row,j) = temp
          END DO
          j = perm(k)
          perm(k) = perm(max_row)
          perm(max_row) = j
        END IF
        DO i = k+1, 6
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 6
            LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
          END DO
        END DO
      END DO

      DO i = 1, 6
        y(i) = b(perm(i))
        DO j = 1, i-1
          y(i) = y(i) - LU(i,j) * y(j)
        END DO
      END DO

      DO i = 6, 1, -1
        x(i) = y(i)
        DO j = i+1, 6
          x(i) = x(i) - LU(i,j) * x(j)
        END DO
        x(i) = x(i) / LU(i,i)
      END DO
      RETURN
      END SUBROUTINE cs_solve_6

C----------------------------------------------------------------------
C     cs_solve_9: solve A(9,9) * x(9) = b(9)
C----------------------------------------------------------------------
      SUBROUTINE cs_solve_9(A, b, x)
      IMPLICIT NONE
      DOUBLE COMPLEX, INTENT(IN)  :: A(9,9)
      DOUBLE COMPLEX, INTENT(IN)  :: b(9)
      DOUBLE COMPLEX, INTENT(OUT) :: x(9)
      DOUBLE COMPLEX :: LU(9,9), y(9), temp
      INTEGER :: perm(9)
      INTEGER :: i, j, k, max_row
      DOUBLE PRECISION :: max_val, abs_val

      DO i = 1, 9
        perm(i) = i
        DO j = 1, 9
          LU(i,j) = A(i,j)
        END DO
      END DO

      DO k = 1, 9
        max_val = ABS(LU(k,k))
        max_row = k
        DO i = k+1, 9
          abs_val = ABS(LU(i,k))
          IF (abs_val > max_val) THEN
            max_val = abs_val
            max_row = i
          END IF
        END DO
        IF (max_row /= k) THEN
          DO j = 1, 9
            temp = LU(k,j)
            LU(k,j) = LU(max_row,j)
            LU(max_row,j) = temp
          END DO
          j = perm(k)
          perm(k) = perm(max_row)
          perm(max_row) = j
        END IF
        DO i = k+1, 9
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 9
            LU(i,j) = LU(i,j) - LU(i,k) * LU(k,j)
          END DO
        END DO
      END DO

      DO i = 1, 9
        y(i) = b(perm(i))
        DO j = 1, i-1
          y(i) = y(i) - LU(i,j) * y(j)
        END DO
      END DO

      DO i = 9, 1, -1
        x(i) = y(i)
        DO j = i+1, 9
          x(i) = x(i) - LU(i,j) * x(j)
        END DO
        x(i) = x(i) / LU(i,i)
      END DO
      RETURN
      END SUBROUTINE cs_solve_9

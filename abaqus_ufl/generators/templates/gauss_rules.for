C======================================================================
C     gauss_rules.for -- Gauss quadrature rules for 2D elements
C
C     Contains:
C       1. gauss_2d_2x2:  4-point (2x2), exact to degree 3
C       2. gauss_2d_3x3:  9-point (3x3), exact to degree 5
C       3. gauss_1d_2pt:  2-point edge rule, exact to degree 3
C       4. gauss_1d_3pt:  3-point edge rule, exact to degree 5
C
C     Recommended pairings:
C       Quad4 volume: 2x2
C       Quad8 volume: 3x3
C       Quad4 edge:   2pt
C       Quad8 edge:   3pt
C======================================================================

C----------------------------------------------------------------------
C     gauss_2d_2x2: 4-point Gauss quadrature (2x2)
C     Exact for polynomials up to degree 3.
C----------------------------------------------------------------------
      SUBROUTINE gauss_2d_2x2(xi_gp, w_gp, ngp)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ngp
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(4,2), w_gp(4)

      DOUBLE PRECISION :: a
      a = 1.0d0 / DSQRT(3.0d0)

      ngp = 4

      xi_gp(1,1) = -a;  xi_gp(1,2) = -a;  w_gp(1) = 1.0d0
      xi_gp(2,1) =  a;  xi_gp(2,2) = -a;  w_gp(2) = 1.0d0
      xi_gp(3,1) =  a;  xi_gp(3,2) =  a;  w_gp(3) = 1.0d0
      xi_gp(4,1) = -a;  xi_gp(4,2) =  a;  w_gp(4) = 1.0d0

      RETURN
      END SUBROUTINE gauss_2d_2x2

C----------------------------------------------------------------------
C     gauss_2d_3x3: 9-point Gauss quadrature (3x3)
C     Exact for polynomials up to degree 5.
C     Recommended for Quad8 elements.
C----------------------------------------------------------------------
      SUBROUTINE gauss_2d_3x3(xi_gp, w_gp, ngp)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ngp
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(9,2), w_gp(9)

      DOUBLE PRECISION :: a, w1, w2, w3
      DOUBLE PRECISION :: pts(3), wts(3)
      INTEGER :: i, j, n

      a  = DSQRT(0.6d0)
      w1 = 5.0d0 / 9.0d0
      w2 = 8.0d0 / 9.0d0
      w3 = 5.0d0 / 9.0d0

      pts(1) = -a;  wts(1) = w1
      pts(2) = 0.0d0; wts(2) = w2
      pts(3) =  a;  wts(3) = w3

      ngp = 9
      n = 0
      DO j = 1, 3
        DO i = 1, 3
          n = n + 1
          xi_gp(n,1) = pts(i)
          xi_gp(n,2) = pts(j)
          w_gp(n) = wts(i) * wts(j)
        END DO
      END DO

      RETURN
      END SUBROUTINE gauss_2d_3x3

C----------------------------------------------------------------------
C     gauss_1d_2pt: 2-point Gauss quadrature for edge integrals
C     Exact for polynomials up to degree 3.
C----------------------------------------------------------------------
      SUBROUTINE gauss_1d_2pt(xi_gp, w_gp, ngp)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ngp
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(2), w_gp(2)

      ngp = 2
      xi_gp(1) = -1.0d0 / DSQRT(3.0d0);  w_gp(1) = 1.0d0
      xi_gp(2) =  1.0d0 / DSQRT(3.0d0);  w_gp(2) = 1.0d0

      RETURN
      END SUBROUTINE gauss_1d_2pt

C----------------------------------------------------------------------
C     gauss_1d_3pt: 3-point Gauss quadrature for edge integrals
C     Exact for polynomials up to degree 5.
C     Recommended for Quad8 edge integration.
C----------------------------------------------------------------------
      SUBROUTINE gauss_1d_3pt(xi_gp, w_gp, ngp)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ngp
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(3), w_gp(3)

      ngp = 3
      xi_gp(1) = -DSQRT(0.6d0);  w_gp(1) = 5.0d0/9.0d0
      xi_gp(2) =  0.0d0;         w_gp(2) = 8.0d0/9.0d0
      xi_gp(3) =  DSQRT(0.6d0);  w_gp(3) = 5.0d0/9.0d0

      RETURN
      END SUBROUTINE gauss_1d_3pt

C======================================================================
C     gauss_hex.for -- Gauss quadrature rules for hexahedral elements
C
C     Provides integration point coordinates and weights for:
C       - 2x2x2 (8-point)  -- Hex8 full integration
C       - 3x3x3 (27-point) -- Hex20 full integration
C       - 1-point           -- Hex8 reduced integration
C       - 2x2x2 (8-point)  -- Hex20 reduced integration
C
C     Generates: gauss_hex8, gauss_hex20, gauss_hex_reduced
C======================================================================

      SUBROUTINE gauss_hex8(xi_gp, w_gp, n_gp)
C     2x2x2 Gauss quadrature for 8-node hexahedral (full integration)
C
C     Output:
C       xi_gp(8,3) -- integration point coordinates (xi,eta,zeta)
C       w_gp(8)    -- integration weights
C       n_gp       -- number of integration points (8)
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(8,3), w_gp(8)
      INTEGER, INTENT(OUT) :: n_gp

      DOUBLE PRECISION :: gp1
      INTEGER :: ii, jj, kk, idx

      n_gp = 8
      gp1 = 1.0d0 / DSQRT(3.0d0)

      idx = 0
      DO kk = 1, 2
        DO jj = 1, 2
          DO ii = 1, 2
            idx = idx + 1
            xi_gp(idx,1) = (-1.0d0)**(ii+1) * gp1
            xi_gp(idx,2) = (-1.0d0)**(jj+1) * gp1
            xi_gp(idx,3) = (-1.0d0)**(kk+1) * gp1
            w_gp(idx) = 1.0d0
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE gauss_hex8


      SUBROUTINE gauss_hex20(xi_gp, w_gp, n_gp)
C     3x3x3 Gauss quadrature for 20-node hexahedral (full integration)
C
C     Output:
C       xi_gp(27,3) -- integration point coordinates
C       w_gp(27)    -- integration weights
C       n_gp        -- number of integration points (27)
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(27,3), w_gp(27)
      INTEGER, INTENT(OUT) :: n_gp

      DOUBLE PRECISION :: gp(3), wt(3)
      INTEGER :: ii, jj, kk, idx

      n_gp = 27

C     3-point Gauss rule: points and weights
      gp(1) = -DSQRT(0.6d0)
      gp(2) =  0.0d0
      gp(3) =  DSQRT(0.6d0)

      wt(1) = 5.0d0 / 9.0d0
      wt(2) = 8.0d0 / 9.0d0
      wt(3) = 5.0d0 / 9.0d0

      idx = 0
      DO kk = 1, 3
        DO jj = 1, 3
          DO ii = 1, 3
            idx = idx + 1
            xi_gp(idx,1) = gp(ii)
            xi_gp(idx,2) = gp(jj)
            xi_gp(idx,3) = gp(kk)
            w_gp(idx) = wt(ii) * wt(jj) * wt(kk)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE gauss_hex20


      SUBROUTINE gauss_hex_reduced(xi_gp, w_gp, n_gp)
C     1-point Gauss quadrature for reduced integration
C
C     Output:
C       xi_gp(1,3) -- integration point at origin
C       w_gp(1)    -- weight = 8.0 (volume of [-1,1]^3)
C       n_gp       -- 1
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(1,3), w_gp(1)
      INTEGER, INTENT(OUT) :: n_gp

      n_gp = 1
      xi_gp(1,1) = 0.0d0
      xi_gp(1,2) = 0.0d0
      xi_gp(1,3) = 0.0d0
      w_gp(1) = 8.0d0

      RETURN
      END SUBROUTINE gauss_hex_reduced

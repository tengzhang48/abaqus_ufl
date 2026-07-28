C======================================================================
C     edge_quad8.for -- Edge (surface) integration for Quad8 elements
C
C     Standalone copy of edge_shape_quad8 for use without shape_quad8.for.
C     Requires: shape_quad8 subroutine (from shape_quad8.for)
C
C     Note: edge_shape_quad8 is also included in shape_quad8.for for
C     convenience. This file exists so the code generator can pull in
C     edge integration independently.
C======================================================================

C----------------------------------------------------------------------
C     edge_shape_quad8: Compute shape functions, edge Jacobian (ds),
C                       and outward unit normal on a specified face.
C
C     face: 1 = bottom (nodes 1-5-2, eta=-1)
C           2 = right  (nodes 2-6-3, xi=+1)
C           3 = top    (nodes 3-7-4, eta=+1)
C           4 = left   (nodes 4-8-1, xi=-1)
C
C     t: parametric coordinate along the edge [-1, +1]
C
C     Returns:
C       sh8(8)    -- Quad8 shape functions evaluated on this face
C       ds        -- edge Jacobian (length element for line integral)
C       normal(2) -- outward unit normal vector (nx, ny)
C
C     The edge is traversed with t from -1 to +1 such that the
C     outward normal is obtained by rotating the tangent (dX/dt, dY/dt)
C     by +90 degrees: normal = (dY/dt, -dX/dt) / ds.
C     For faces 3 and 4, the t-to-(xi,eta) mapping is reversed to
C     ensure the normal points outward.
C----------------------------------------------------------------------
      SUBROUTINE edge_shape_quad8(t, face, coords, sh8, ds, normal)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: t, coords(2,8)
      INTEGER, INTENT(IN) :: face
      DOUBLE PRECISION, INTENT(OUT) :: sh8(8), ds, normal(2)

      DOUBLE PRECISION :: xi, eta, dshxi(8,2)
      DOUBLE PRECISION :: dXdt, dYdt
      INTEGER :: k

C     Map parametric coordinate t to (xi, eta) on the face
      IF (face .EQ. 1) THEN
C       Bottom edge: eta = -1, xi = t
        xi = t;  eta = -1.0d0
      ELSE IF (face .EQ. 2) THEN
C       Right edge: xi = +1, eta = t
        xi = 1.0d0;  eta = t
      ELSE IF (face .EQ. 3) THEN
C       Top edge: eta = +1, xi = -t  (reversed for outward normal)
        xi = -t;  eta = 1.0d0
      ELSE IF (face .EQ. 4) THEN
C       Left edge: xi = -1, eta = -t  (reversed for outward normal)
        xi = -1.0d0;  eta = -t
      END IF

C     Evaluate Quad8 shape functions at (xi, eta)
      CALL shape_quad8(xi, eta, sh8, dshxi)

C     Compute edge tangent: dX/dt = sum_a (dN_a/dxi * dxi/dt + ...) * X_a
      dXdt = 0.0d0
      dYdt = 0.0d0

      IF (face .EQ. 1) THEN
C       dxi/dt = 1, deta/dt = 0
        DO k = 1, 8
          dXdt = dXdt + dshxi(k,1) * coords(1,k)
          dYdt = dYdt + dshxi(k,1) * coords(2,k)
        END DO
      ELSE IF (face .EQ. 2) THEN
C       dxi/dt = 0, deta/dt = 1
        DO k = 1, 8
          dXdt = dXdt + dshxi(k,2) * coords(1,k)
          dYdt = dYdt + dshxi(k,2) * coords(2,k)
        END DO
      ELSE IF (face .EQ. 3) THEN
C       dxi/dt = -1, deta/dt = 0
        DO k = 1, 8
          dXdt = dXdt - dshxi(k,1) * coords(1,k)
          dYdt = dYdt - dshxi(k,1) * coords(2,k)
        END DO
      ELSE IF (face .EQ. 4) THEN
C       dxi/dt = 0, deta/dt = -1
        DO k = 1, 8
          dXdt = dXdt - dshxi(k,2) * coords(1,k)
          dYdt = dYdt - dshxi(k,2) * coords(2,k)
        END DO
      END IF

C     Edge Jacobian
      ds = DSQRT(dXdt*dXdt + dYdt*dYdt)

C     Outward unit normal: rotate tangent by +90 degrees
C     tangent = (dXdt, dYdt), normal = (dYdt, -dXdt) / ds
      IF (ds .GT. 0.0d0) THEN
        normal(1) =  dYdt / ds
        normal(2) = -dXdt / ds
      ELSE
        normal(1) = 0.0d0
        normal(2) = 0.0d0
      END IF

      RETURN
      END SUBROUTINE edge_shape_quad8

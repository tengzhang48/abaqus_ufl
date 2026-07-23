C======================================================================
C     shape_quad8.for — 8-node serendipity shape functions (Quad8)
C
C     Contains:
C       1. shape_quad8:      shape functions and derivatives
C       2. edge_shape_quad8: edge shape functions, Jacobian, normal
C
C     Node numbering:
C
C        4-----7-----3            eta
C        |           |             |
C        8           6             |
C        |           |             +-----xi
C        1-----5-----2
C
C     Corners:  1(-1,-1), 2(+1,-1), 3(+1,+1), 4(-1,+1)
C     Midsides: 5(0,-1),  6(+1,0),  7(0,+1),  8(-1,0)
C
C     For mixed formulation:
C       degree=2 fields (u, mu): use all 8 nodes (Quad8 shapes)
C       degree=1 fields (p):     use corner nodes 1-4 (see shape_quad4)
C======================================================================

C----------------------------------------------------------------------
C     shape_quad8: 8-node serendipity shape functions and derivatives
C
C     sh(8)       = shape function values at (xi, eta)
C     dshxi(8,2)  = derivatives: dshxi(a,1)=dN_a/dxi, dshxi(a,2)=dN_a/deta
C----------------------------------------------------------------------
      SUBROUTINE shape_quad8(xi, eta, sh, dshxi)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: xi, eta
      DOUBLE PRECISION, INTENT(OUT) :: sh(8), dshxi(8,2)

C     Node 1: corner (-1,-1)
      sh(1) = 0.25d0*(1.0d0-xi)*(1.0d0-eta)*(-xi-eta-1.0d0)
      dshxi(1,1) = 0.25d0*(1.0d0-eta)*(2.0d0*xi+eta)
      dshxi(1,2) = 0.25d0*(1.0d0-xi)*(xi+2.0d0*eta)

C     Node 2: corner (+1,-1)
      sh(2) = 0.25d0*(1.0d0+xi)*(1.0d0-eta)*(xi-eta-1.0d0)
      dshxi(2,1) = 0.25d0*(1.0d0-eta)*(2.0d0*xi-eta)
      dshxi(2,2) = 0.25d0*(1.0d0+xi)*(-xi+2.0d0*eta)

C     Node 3: corner (+1,+1)
      sh(3) = 0.25d0*(1.0d0+xi)*(1.0d0+eta)*(xi+eta-1.0d0)
      dshxi(3,1) = 0.25d0*(1.0d0+eta)*(2.0d0*xi+eta)
      dshxi(3,2) = 0.25d0*(1.0d0+xi)*(xi+2.0d0*eta)

C     Node 4: corner (-1,+1)
      sh(4) = 0.25d0*(1.0d0-xi)*(1.0d0+eta)*(-xi+eta-1.0d0)
      dshxi(4,1) = 0.25d0*(1.0d0+eta)*(2.0d0*xi-eta)
      dshxi(4,2) = 0.25d0*(1.0d0-xi)*(-xi+2.0d0*eta)

C     Node 5: midside (0,-1) — between nodes 1 and 2
      sh(5) = 0.5d0*(1.0d0-xi*xi)*(1.0d0-eta)
      dshxi(5,1) = -xi*(1.0d0-eta)
      dshxi(5,2) = -0.5d0*(1.0d0-xi*xi)

C     Node 6: midside (+1,0) — between nodes 2 and 3
      sh(6) = 0.5d0*(1.0d0+xi)*(1.0d0-eta*eta)
      dshxi(6,1) = 0.5d0*(1.0d0-eta*eta)
      dshxi(6,2) = -(1.0d0+xi)*eta

C     Node 7: midside (0,+1) — between nodes 3 and 4
      sh(7) = 0.5d0*(1.0d0-xi*xi)*(1.0d0+eta)
      dshxi(7,1) = -xi*(1.0d0+eta)
      dshxi(7,2) = 0.5d0*(1.0d0-xi*xi)

C     Node 8: midside (-1,0) — between nodes 4 and 1
      sh(8) = 0.5d0*(1.0d0-xi)*(1.0d0-eta*eta)
      dshxi(8,1) = -0.5d0*(1.0d0-eta*eta)
      dshxi(8,2) = -(1.0d0-xi)*eta

      RETURN
      END SUBROUTINE shape_quad8

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
C       sh8(8)    — Quad8 shape functions evaluated on this face
C       ds        — edge Jacobian (length element for line integral)
C       normal(2) — outward unit normal vector (nx, ny)
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

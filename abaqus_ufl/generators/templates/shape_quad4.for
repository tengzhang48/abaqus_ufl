C======================================================================
C     shape_quad4.for — 4-node bilinear shape functions (Quad4)
C
C     Node numbering:
C
C        4-----------3            eta
C        |           |             |
C        |           |             |
C        |           |             +-----xi
C        1-----------2
C
C     Corners: 1(-1,-1), 2(+1,-1), 3(+1,+1), 4(-1,+1)
C
C     Used for:
C       - degree=1 fields (e.g. pressure) in mixed formulation
C       - standalone 4-node elements
C======================================================================

C----------------------------------------------------------------------
C     shape_quad4: 4-node bilinear shape functions and derivatives
C
C     sh(4)       = shape function values at (xi, eta)
C     dshxi(4,2)  = derivatives: dshxi(a,1)=dN_a/dxi, dshxi(a,2)=dN_a/deta
C----------------------------------------------------------------------
      SUBROUTINE shape_quad4(xi, eta, sh, dshxi)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: xi, eta
      DOUBLE PRECISION, INTENT(OUT) :: sh(4), dshxi(4,2)

C     Node 1: (-1,-1)
      sh(1) = 0.25d0*(1.0d0 - xi)*(1.0d0 - eta)
      dshxi(1,1) = -0.25d0*(1.0d0 - eta)
      dshxi(1,2) = -0.25d0*(1.0d0 - xi)

C     Node 2: (+1,-1)
      sh(2) = 0.25d0*(1.0d0 + xi)*(1.0d0 - eta)
      dshxi(2,1) =  0.25d0*(1.0d0 - eta)
      dshxi(2,2) = -0.25d0*(1.0d0 + xi)

C     Node 3: (+1,+1)
      sh(3) = 0.25d0*(1.0d0 + xi)*(1.0d0 + eta)
      dshxi(3,1) = 0.25d0*(1.0d0 + eta)
      dshxi(3,2) = 0.25d0*(1.0d0 + xi)

C     Node 4: (-1,+1)
      sh(4) = 0.25d0*(1.0d0 - xi)*(1.0d0 + eta)
      dshxi(4,1) = -0.25d0*(1.0d0 + eta)
      dshxi(4,2) =  0.25d0*(1.0d0 - xi)

      RETURN
      END SUBROUTINE shape_quad4

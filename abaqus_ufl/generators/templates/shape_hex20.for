C======================================================================
C     shape_hex20.for — Shape functions for 20-node hexahedral element
C
C     Node numbering (Abaqus convention):
C
C     Corners (1-8):
C         8-------7          zeta
C        /|      /|          |  eta
C       5-------6 |          | /
C       | 4-----|-3          |/_____ xi
C       |/      |/
C       1-------2
C
C     Midside nodes (9-20):
C       9:  midpoint of edge 1-2  (xi-dir, eta=-1, zeta=-1)
C      10:  midpoint of edge 2-3  (eta-dir, xi=+1, zeta=-1)
C      11:  midpoint of edge 3-4  (xi-dir, eta=+1, zeta=-1)
C      12:  midpoint of edge 4-1  (eta-dir, xi=-1, zeta=-1)
C      13:  midpoint of edge 5-6  (xi-dir, eta=-1, zeta=+1)
C      14:  midpoint of edge 6-7  (eta-dir, xi=+1, zeta=+1)
C      15:  midpoint of edge 7-8  (xi-dir, eta=+1, zeta=+1)
C      16:  midpoint of edge 8-5  (eta-dir, xi=-1, zeta=+1)
C      17:  midpoint of edge 1-5  (zeta-dir, xi=-1, eta=-1)
C      18:  midpoint of edge 2-6  (zeta-dir, xi=+1, eta=-1)
C      19:  midpoint of edge 3-7  (zeta-dir, xi=+1, eta=+1)
C      20:  midpoint of edge 4-8  (zeta-dir, xi=-1, eta=+1)
C
C     Generates: shape_hex20 (combined N + dN/dxi)
C======================================================================

      SUBROUTINE shape_hex20(xi, eta, zeta, N, dNdxi)
C     20-node hexahedral shape functions and derivatives (serendipity)
C
C     Input:
C       xi, eta, zeta — natural coordinates in [-1,1]
C     Output:
C       N(20)       — shape function values at (xi,eta,zeta)
C       dNdxi(20,3) — dNdxi(a,i) = dN_a/dxi_i
C                      i=1: d/dxi, i=2: d/deta, i=3: d/dzeta
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: xi, eta, zeta
      DOUBLE PRECISION, INTENT(OUT) :: N(20), dNdxi(20,3)

      DOUBLE PRECISION :: xp, xm, ep, em, zp, zm
      DOUBLE PRECISION :: xi2, eta2, zeta2

      xp = 1.0d0 + xi
      xm = 1.0d0 - xi
      ep = 1.0d0 + eta
      em = 1.0d0 - eta
      zp = 1.0d0 + zeta
      zm = 1.0d0 - zeta

      xi2   = 1.0d0 - xi*xi
      eta2  = 1.0d0 - eta*eta
      zeta2 = 1.0d0 - zeta*zeta

C     --- Shape function values ---

C     Corner nodes
      N(1)  = 0.125d0*xm*em*zm*(-xi-eta-zeta-2.0d0)
      N(2)  = 0.125d0*xp*em*zm*( xi-eta-zeta-2.0d0)
      N(3)  = 0.125d0*xp*ep*zm*( xi+eta-zeta-2.0d0)
      N(4)  = 0.125d0*xm*ep*zm*(-xi+eta-zeta-2.0d0)
      N(5)  = 0.125d0*xm*em*zp*(-xi-eta+zeta-2.0d0)
      N(6)  = 0.125d0*xp*em*zp*( xi-eta+zeta-2.0d0)
      N(7)  = 0.125d0*xp*ep*zp*( xi+eta+zeta-2.0d0)
      N(8)  = 0.125d0*xm*ep*zp*(-xi+eta+zeta-2.0d0)

C     Midside nodes on bottom face (zeta=-1)
      N(9)  = 0.25d0*xi2*em*zm
      N(10) = 0.25d0*xp*eta2*zm
      N(11) = 0.25d0*xi2*ep*zm
      N(12) = 0.25d0*xm*eta2*zm

C     Midside nodes on top face (zeta=+1)
      N(13) = 0.25d0*xi2*em*zp
      N(14) = 0.25d0*xp*eta2*zp
      N(15) = 0.25d0*xi2*ep*zp
      N(16) = 0.25d0*xm*eta2*zp

C     Midside nodes on vertical edges
      N(17) = 0.25d0*xm*em*zeta2
      N(18) = 0.25d0*xp*em*zeta2
      N(19) = 0.25d0*xp*ep*zeta2
      N(20) = 0.25d0*xm*ep*zeta2

C     --- Derivatives dN/dxi ---

C     Corner nodes
      dNdxi(1,1)  = 0.125d0*em*zm*( 2.0d0*xi+eta+zeta+1.0d0)
      dNdxi(2,1)  = 0.125d0*em*zm*( 2.0d0*xi-eta-zeta-1.0d0)
      dNdxi(3,1)  = 0.125d0*ep*zm*( 2.0d0*xi+eta-zeta-1.0d0)
      dNdxi(4,1)  = 0.125d0*ep*zm*( 2.0d0*xi-eta+zeta+1.0d0)
      dNdxi(5,1)  = 0.125d0*em*zp*( 2.0d0*xi+eta-zeta+1.0d0)
      dNdxi(6,1)  = 0.125d0*em*zp*( 2.0d0*xi-eta+zeta-1.0d0)
      dNdxi(7,1)  = 0.125d0*ep*zp*( 2.0d0*xi+eta+zeta-1.0d0)
      dNdxi(8,1)  = 0.125d0*ep*zp*( 2.0d0*xi-eta-zeta+1.0d0)

C     Bottom midside
      dNdxi(9,1)  = -0.5d0*xi*em*zm
      dNdxi(10,1) =  0.25d0*eta2*zm
      dNdxi(11,1) = -0.5d0*xi*ep*zm
      dNdxi(12,1) = -0.25d0*eta2*zm

C     Top midside
      dNdxi(13,1) = -0.5d0*xi*em*zp
      dNdxi(14,1) =  0.25d0*eta2*zp
      dNdxi(15,1) = -0.5d0*xi*ep*zp
      dNdxi(16,1) = -0.25d0*eta2*zp

C     Vertical midside
      dNdxi(17,1) = -0.25d0*em*zeta2
      dNdxi(18,1) =  0.25d0*em*zeta2
      dNdxi(19,1) =  0.25d0*ep*zeta2
      dNdxi(20,1) = -0.25d0*ep*zeta2

C     ============================================================
C     dN/deta (column 2)
C     ============================================================

C     Corner nodes
      dNdxi(1,2)  = 0.125d0*xm*zm*( xi+2.0d0*eta+zeta+1.0d0)
      dNdxi(2,2)  = 0.125d0*xp*zm*(-xi+2.0d0*eta+zeta+1.0d0)
      dNdxi(3,2)  = 0.125d0*xp*zm*( xi+2.0d0*eta-zeta-1.0d0)
      dNdxi(4,2)  = 0.125d0*xm*zm*(-xi+2.0d0*eta-zeta-1.0d0)
      dNdxi(5,2)  = 0.125d0*xm*zp*( xi+2.0d0*eta-zeta+1.0d0)
      dNdxi(6,2)  = 0.125d0*xp*zp*(-xi+2.0d0*eta-zeta+1.0d0)
      dNdxi(7,2)  = 0.125d0*xp*zp*( xi+2.0d0*eta+zeta-1.0d0)
      dNdxi(8,2)  = 0.125d0*xm*zp*(-xi+2.0d0*eta+zeta-1.0d0)

C     Bottom midside
      dNdxi(9,2)  = -0.25d0*xi2*zm
      dNdxi(10,2) = -0.5d0*xp*eta*zm
      dNdxi(11,2) =  0.25d0*xi2*zm
      dNdxi(12,2) = -0.5d0*xm*eta*zm

C     Top midside
      dNdxi(13,2) = -0.25d0*xi2*zp
      dNdxi(14,2) = -0.5d0*xp*eta*zp
      dNdxi(15,2) =  0.25d0*xi2*zp
      dNdxi(16,2) = -0.5d0*xm*eta*zp

C     Vertical midside
      dNdxi(17,2) = -0.25d0*xm*zeta2
      dNdxi(18,2) = -0.25d0*xp*zeta2
      dNdxi(19,2) =  0.25d0*xp*zeta2
      dNdxi(20,2) =  0.25d0*xm*zeta2

C     ============================================================
C     dN/dzeta (column 3)
C     ============================================================

C     Corner nodes
      dNdxi(1,3)  = 0.125d0*xm*em*( xi+eta+2.0d0*zeta+1.0d0)
      dNdxi(2,3)  = 0.125d0*xp*em*(-xi+eta+2.0d0*zeta+1.0d0)
      dNdxi(3,3)  = 0.125d0*xp*ep*(-xi-eta+2.0d0*zeta+1.0d0)
      dNdxi(4,3)  = 0.125d0*xm*ep*( xi-eta+2.0d0*zeta+1.0d0)
      dNdxi(5,3)  = 0.125d0*xm*em*(-xi-eta+2.0d0*zeta-1.0d0)
      dNdxi(6,3)  = 0.125d0*xp*em*( xi-eta+2.0d0*zeta-1.0d0)
      dNdxi(7,3)  = 0.125d0*xp*ep*( xi+eta+2.0d0*zeta-1.0d0)
      dNdxi(8,3)  = 0.125d0*xm*ep*(-xi+eta+2.0d0*zeta-1.0d0)

C     Bottom midside
      dNdxi(9,3)  = -0.25d0*xi2*em
      dNdxi(10,3) = -0.25d0*xp*eta2
      dNdxi(11,3) = -0.25d0*xi2*ep
      dNdxi(12,3) = -0.25d0*xm*eta2

C     Top midside
      dNdxi(13,3) =  0.25d0*xi2*em
      dNdxi(14,3) =  0.25d0*xp*eta2
      dNdxi(15,3) =  0.25d0*xi2*ep
      dNdxi(16,3) =  0.25d0*xm*eta2

C     Vertical midside
      dNdxi(17,3) = -0.5d0*xm*em*zeta
      dNdxi(18,3) = -0.5d0*xp*em*zeta
      dNdxi(19,3) = -0.5d0*xp*ep*zeta
      dNdxi(20,3) = -0.5d0*xm*ep*zeta

      RETURN
      END SUBROUTINE shape_hex20

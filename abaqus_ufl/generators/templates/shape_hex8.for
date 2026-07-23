C======================================================================
C     shape_hex8.for — Shape functions for 8-node hexahedral element
C
C     Node numbering (Abaqus convention):
C
C         8-------7          zeta
C        /|      /|          |  eta
C       5-------6 |          | /
C       | 4-----|-3          |/_____ xi
C       |/      |/
C       1-------2
C
C     Corner nodes at (xi,eta,zeta) = (±1,±1,±1):
C       Node 1: (-1,-1,-1)
C       Node 2: (+1,-1,-1)
C       Node 3: (+1,+1,-1)
C       Node 4: (-1,+1,-1)
C       Node 5: (-1,-1,+1)
C       Node 6: (+1,-1,+1)
C       Node 7: (+1,+1,+1)
C       Node 8: (-1,+1,+1)
C
C     Generates: shape_hex8 (combined N + dN/dxi)
C======================================================================

      SUBROUTINE shape_hex8(xi, eta, zeta, N, dNdxi)
C     8-node hexahedral shape functions and derivatives
C
C     Input:
C       xi, eta, zeta — natural coordinates in [-1,1]
C     Output:
C       N(8)      — shape function values at (xi,eta,zeta)
C       dNdxi(8,3) — dNdxi(a,i) = dN_a/dxi_i
C                     i=1: d/dxi, i=2: d/deta, i=3: d/dzeta
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: xi, eta, zeta
      DOUBLE PRECISION, INTENT(OUT) :: N(8), dNdxi(8,3)

      DOUBLE PRECISION :: xp, xm, ep, em, zp, zm

      xp = 1.0d0 + xi
      xm = 1.0d0 - xi
      ep = 1.0d0 + eta
      em = 1.0d0 - eta
      zp = 1.0d0 + zeta
      zm = 1.0d0 - zeta

C     Shape function values
      N(1) = 0.125d0 * xm * em * zm
      N(2) = 0.125d0 * xp * em * zm
      N(3) = 0.125d0 * xp * ep * zm
      N(4) = 0.125d0 * xm * ep * zm
      N(5) = 0.125d0 * xm * em * zp
      N(6) = 0.125d0 * xp * em * zp
      N(7) = 0.125d0 * xp * ep * zp
      N(8) = 0.125d0 * xm * ep * zp

C     dN/dxi
      dNdxi(1,1) = -0.125d0 * em * zm
      dNdxi(2,1) =  0.125d0 * em * zm
      dNdxi(3,1) =  0.125d0 * ep * zm
      dNdxi(4,1) = -0.125d0 * ep * zm
      dNdxi(5,1) = -0.125d0 * em * zp
      dNdxi(6,1) =  0.125d0 * em * zp
      dNdxi(7,1) =  0.125d0 * ep * zp
      dNdxi(8,1) = -0.125d0 * ep * zp

C     dN/deta
      dNdxi(1,2) = -0.125d0 * xm * zm
      dNdxi(2,2) = -0.125d0 * xp * zm
      dNdxi(3,2) =  0.125d0 * xp * zm
      dNdxi(4,2) =  0.125d0 * xm * zm
      dNdxi(5,2) = -0.125d0 * xm * zp
      dNdxi(6,2) = -0.125d0 * xp * zp
      dNdxi(7,2) =  0.125d0 * xp * zp
      dNdxi(8,2) =  0.125d0 * xm * zp

C     dN/dzeta
      dNdxi(1,3) = -0.125d0 * xm * em
      dNdxi(2,3) = -0.125d0 * xp * em
      dNdxi(3,3) = -0.125d0 * xp * ep
      dNdxi(4,3) = -0.125d0 * xm * ep
      dNdxi(5,3) =  0.125d0 * xm * em
      dNdxi(6,3) =  0.125d0 * xp * em
      dNdxi(7,3) =  0.125d0 * xp * ep
      dNdxi(8,3) =  0.125d0 * xm * ep

      RETURN
      END SUBROUTINE shape_hex8

C======================================================================
C     face_hex.for -- Face shape functions for hexahedral elements
C
C     Used for integrating distributed loads (surface tractions)
C     over element faces in the Abaqus UEL NDLOAD interface.
C
C     Abaqus face numbering convention for hex elements:
C
C       Face 1: nodes on eta = -1  (1,2,6,5 for Hex8)
C       Face 2: nodes on eta = +1  (3,4,8,7 for Hex8)
C       Face 3: nodes on xi  = -1  (1,4,8,5 for Hex8)
C       Face 4: nodes on xi  = +1  (2,3,7,6 for Hex8)
C       Face 5: nodes on zeta= -1  (1,2,3,4 for Hex8)
C       Face 6: nodes on zeta= +1  (5,6,7,8 for Hex8)
C
C     Each face is a 2D surface parameterized by two of the three
C     natural coordinates (the third is fixed at +/-1).
C
C     Generates: face_hex8_nodes, face_hex8_shape
C                face_hex20_nodes, face_hex20_shape
C======================================================================

      SUBROUTINE face_hex8_nodes(iface, face_nodes, n_face_nodes)
C     Return local node numbers on a given face of Hex8
C
C     Input:
C       iface -- face number (1-6, Abaqus convention)
C     Output:
C       face_nodes(4) -- local node numbers on this face
C       n_face_nodes  -- 4
C
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: iface
      INTEGER, INTENT(OUT) :: face_nodes(4), n_face_nodes

      n_face_nodes = 4

      IF (iface .EQ. 1) THEN
C       eta = -1: nodes 1,2,6,5
        face_nodes(1) = 1
        face_nodes(2) = 2
        face_nodes(3) = 6
        face_nodes(4) = 5
      ELSE IF (iface .EQ. 2) THEN
C       eta = +1: nodes 3,4,8,7
        face_nodes(1) = 3
        face_nodes(2) = 4
        face_nodes(3) = 8
        face_nodes(4) = 7
      ELSE IF (iface .EQ. 3) THEN
C       xi = -1: nodes 1,4,8,5
        face_nodes(1) = 1
        face_nodes(2) = 4
        face_nodes(3) = 8
        face_nodes(4) = 5
      ELSE IF (iface .EQ. 4) THEN
C       xi = +1: nodes 2,3,7,6
        face_nodes(1) = 2
        face_nodes(2) = 3
        face_nodes(3) = 7
        face_nodes(4) = 6
      ELSE IF (iface .EQ. 5) THEN
C       zeta = -1: nodes 1,2,3,4
        face_nodes(1) = 1
        face_nodes(2) = 2
        face_nodes(3) = 3
        face_nodes(4) = 4
      ELSE IF (iface .EQ. 6) THEN
C       zeta = +1: nodes 5,6,7,8
        face_nodes(1) = 5
        face_nodes(2) = 6
        face_nodes(3) = 7
        face_nodes(4) = 8
      END IF

      RETURN
      END SUBROUTINE face_hex8_nodes


      SUBROUTINE face_hex8_shape(iface, s1, s2,
     &                           N_face, dNds, n_face_nodes)
C     Shape functions on a face of Hex8 (bilinear quad)
C
C     Input:
C       iface -- face number (1-6)
C       s1, s2 -- face natural coordinates in [-1,1]
C     Output:
C       N_face(4) -- shape function values
C       dNds(4,2) -- dN/ds1, dN/ds2
C       n_face_nodes -- 4
C
C     The mapping from face params (s1,s2) to volume params
C     (xi,eta,zeta) depends on which face:
C       Face 1 (eta=-1): s1->xi, s2->zeta
C       Face 2 (eta=+1): s1->xi, s2->zeta
C       Face 3 (xi=-1):  s1->eta, s2->zeta
C       Face 4 (xi=+1):  s1->eta, s2->zeta
C       Face 5 (zeta=-1): s1->xi, s2->eta
C       Face 6 (zeta=+1): s1->xi, s2->eta
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iface
      DOUBLE PRECISION, INTENT(IN)  :: s1, s2
      DOUBLE PRECISION, INTENT(OUT) :: N_face(4), dNds(4,2)
      INTEGER, INTENT(OUT) :: n_face_nodes

      DOUBLE PRECISION :: sp1, sm1, sp2, sm2

      n_face_nodes = 4

      sp1 = 1.0d0 + s1
      sm1 = 1.0d0 - s1
      sp2 = 1.0d0 + s2
      sm2 = 1.0d0 - s2

C     Bilinear quad shape functions on the face
      N_face(1) = 0.25d0 * sm1 * sm2
      N_face(2) = 0.25d0 * sp1 * sm2
      N_face(3) = 0.25d0 * sp1 * sp2
      N_face(4) = 0.25d0 * sm1 * sp2

C     dN/ds1
      dNds(1,1) = -0.25d0 * sm2
      dNds(2,1) =  0.25d0 * sm2
      dNds(3,1) =  0.25d0 * sp2
      dNds(4,1) = -0.25d0 * sp2

C     dN/ds2
      dNds(1,2) = -0.25d0 * sm1
      dNds(2,2) = -0.25d0 * sp1
      dNds(3,2) =  0.25d0 * sp1
      dNds(4,2) =  0.25d0 * sm1

      RETURN
      END SUBROUTINE face_hex8_shape


      SUBROUTINE face_hex20_nodes(iface, face_nodes,
     &                            n_face_nodes)
C     Return local node numbers on a given face of Hex20
C     Each face is an 8-node serendipity quad
C
C     Input:
C       iface -- face number (1-6, Abaqus convention)
C     Output:
C       face_nodes(8) -- local node numbers (4 corners + 4 midside)
C       n_face_nodes  -- 8
C
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: iface
      INTEGER, INTENT(OUT) :: face_nodes(8), n_face_nodes

      n_face_nodes = 8

      IF (iface .EQ. 1) THEN
C       eta = -1: corners 1,2,6,5; midside 9,18,13,17
        face_nodes(1) = 1
        face_nodes(2) = 2
        face_nodes(3) = 6
        face_nodes(4) = 5
        face_nodes(5) = 9
        face_nodes(6) = 18
        face_nodes(7) = 13
        face_nodes(8) = 17
      ELSE IF (iface .EQ. 2) THEN
C       eta = +1: corners 3,4,8,7; midside 11,20,15,19
        face_nodes(1) = 3
        face_nodes(2) = 4
        face_nodes(3) = 8
        face_nodes(4) = 7
        face_nodes(5) = 11
        face_nodes(6) = 20
        face_nodes(7) = 15
        face_nodes(8) = 19
      ELSE IF (iface .EQ. 3) THEN
C       xi = -1: corners 1,4,8,5; midside 12,20,16,17
        face_nodes(1) = 1
        face_nodes(2) = 4
        face_nodes(3) = 8
        face_nodes(4) = 5
        face_nodes(5) = 12
        face_nodes(6) = 20
        face_nodes(7) = 16
        face_nodes(8) = 17
      ELSE IF (iface .EQ. 4) THEN
C       xi = +1: corners 2,3,7,6; midside 10,19,14,18
        face_nodes(1) = 2
        face_nodes(2) = 3
        face_nodes(3) = 7
        face_nodes(4) = 6
        face_nodes(5) = 10
        face_nodes(6) = 19
        face_nodes(7) = 14
        face_nodes(8) = 18
      ELSE IF (iface .EQ. 5) THEN
C       zeta = -1: corners 1,2,3,4; midside 9,10,11,12
        face_nodes(1) = 1
        face_nodes(2) = 2
        face_nodes(3) = 3
        face_nodes(4) = 4
        face_nodes(5) = 9
        face_nodes(6) = 10
        face_nodes(7) = 11
        face_nodes(8) = 12
      ELSE IF (iface .EQ. 6) THEN
C       zeta = +1: corners 5,6,7,8; midside 13,14,15,16
        face_nodes(1) = 5
        face_nodes(2) = 6
        face_nodes(3) = 7
        face_nodes(4) = 8
        face_nodes(5) = 13
        face_nodes(6) = 14
        face_nodes(7) = 15
        face_nodes(8) = 16
      END IF

      RETURN
      END SUBROUTINE face_hex20_nodes


      SUBROUTINE face_hex20_shape(iface, s1, s2,
     &                            N_face, dNds, n_face_nodes)
C     Shape functions on a face of Hex20 (8-node serendipity quad)
C
C     Input:
C       iface -- face number (1-6)
C       s1, s2 -- face natural coordinates in [-1,1]
C     Output:
C       N_face(8) -- shape function values
C       dNds(8,2) -- dN/ds1, dN/ds2
C       n_face_nodes -- 8
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iface
      DOUBLE PRECISION, INTENT(IN)  :: s1, s2
      DOUBLE PRECISION, INTENT(OUT) :: N_face(8), dNds(8,2)
      INTEGER, INTENT(OUT) :: n_face_nodes

      DOUBLE PRECISION :: sp1, sm1, sp2, sm2
      DOUBLE PRECISION :: s1sq, s2sq

      n_face_nodes = 8

      sp1 = 1.0d0 + s1
      sm1 = 1.0d0 - s1
      sp2 = 1.0d0 + s2
      sm2 = 1.0d0 - s2
      s1sq = 1.0d0 - s1*s1
      s2sq = 1.0d0 - s2*s2

C     Corner nodes (serendipity)
      N_face(1) = 0.25d0*sm1*sm2*(-s1-s2-1.0d0)
      N_face(2) = 0.25d0*sp1*sm2*( s1-s2-1.0d0)
      N_face(3) = 0.25d0*sp1*sp2*( s1+s2-1.0d0)
      N_face(4) = 0.25d0*sm1*sp2*(-s1+s2-1.0d0)

C     Midside nodes
      N_face(5) = 0.5d0*s1sq*sm2
      N_face(6) = 0.5d0*sp1*s2sq
      N_face(7) = 0.5d0*s1sq*sp2
      N_face(8) = 0.5d0*sm1*s2sq

C     dN/ds1
      dNds(1,1) = 0.25d0*sm2*( 2.0d0*s1+s2)
      dNds(2,1) = 0.25d0*sm2*( 2.0d0*s1-s2)
      dNds(3,1) = 0.25d0*sp2*( 2.0d0*s1+s2)
      dNds(4,1) = 0.25d0*sp2*( 2.0d0*s1-s2)
      dNds(5,1) = -s1*sm2
      dNds(6,1) =  0.5d0*s2sq
      dNds(7,1) = -s1*sp2
      dNds(8,1) = -0.5d0*s2sq

C     dN/ds2
      dNds(1,2) = 0.25d0*sm1*( s1+2.0d0*s2)
      dNds(2,2) = 0.25d0*sp1*(-s1+2.0d0*s2)
      dNds(3,2) = 0.25d0*sp1*( s1+2.0d0*s2)
      dNds(4,2) = 0.25d0*sm1*(-s1+2.0d0*s2)
      dNds(5,2) = -0.5d0*s1sq
      dNds(6,2) = -sp1*s2
      dNds(7,2) =  0.5d0*s1sq
      dNds(8,2) = -sm1*s2

      RETURN
      END SUBROUTINE face_hex20_shape


      SUBROUTINE face_gauss_2x2(s1_gp, s2_gp, w_gp, n_gp)
C     2x2 Gauss quadrature on a face (for Hex8 faces)
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: s1_gp(4), s2_gp(4), w_gp(4)
      INTEGER, INTENT(OUT) :: n_gp

      DOUBLE PRECISION :: gp1

      n_gp = 4
      gp1 = 1.0d0 / DSQRT(3.0d0)

      s1_gp(1) = -gp1;  s2_gp(1) = -gp1;  w_gp(1) = 1.0d0
      s1_gp(2) =  gp1;  s2_gp(2) = -gp1;  w_gp(2) = 1.0d0
      s1_gp(3) =  gp1;  s2_gp(3) =  gp1;  w_gp(3) = 1.0d0
      s1_gp(4) = -gp1;  s2_gp(4) =  gp1;  w_gp(4) = 1.0d0

      RETURN
      END SUBROUTINE face_gauss_2x2


      SUBROUTINE face_gauss_3x3(s1_gp, s2_gp, w_gp, n_gp)
C     3x3 Gauss quadrature on a face (for Hex20 faces)
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: s1_gp(9), s2_gp(9), w_gp(9)
      INTEGER, INTENT(OUT) :: n_gp

      DOUBLE PRECISION :: gp(3), wt(3)
      INTEGER :: ii, jj, idx

      n_gp = 9

      gp(1) = -DSQRT(0.6d0)
      gp(2) =  0.0d0
      gp(3) =  DSQRT(0.6d0)

      wt(1) = 5.0d0 / 9.0d0
      wt(2) = 8.0d0 / 9.0d0
      wt(3) = 5.0d0 / 9.0d0

      idx = 0
      DO jj = 1, 3
        DO ii = 1, 3
          idx = idx + 1
          s1_gp(idx) = gp(ii)
          s2_gp(idx) = gp(jj)
          w_gp(idx) = wt(ii) * wt(jj)
        END DO
      END DO

      RETURN
      END SUBROUTINE face_gauss_3x3

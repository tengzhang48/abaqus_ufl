C======================================================================
C     shape_tet4.for -- Shape functions + quadrature for 4-node
C     linear tetrahedral element (Abaqus C3D4 convention)
C
C     Reference element: nodes at
C       Node 1: (0,0,0)
C       Node 2: (1,0,0)
C       Node 3: (0,1,0)
C       Node 4: (0,0,1)
C     with natural coordinates (xi,eta,zeta), 0 <= xi+eta+zeta <= 1.
C     Nodes 1-2-3 counterclockwise viewed from node 4 gives detJ > 0
C     (Abaqus C3D4 node ordering).
C
C     N1 = 1 - xi - eta - zeta,  N2 = xi,  N3 = eta,  N4 = zeta
C     All derivatives are constant: the linear tet has an element-wise
C     constant deformation gradient.
C
C     Generates: shape_tet4, gauss_tet4 (4-pt, degree 2),
C                gauss_tet1 (1-pt centroid, degree 1)
C======================================================================

      SUBROUTINE shape_tet4(xi, eta, zeta, N, dNdxi)
C     4-node linear tetrahedron shape functions and derivatives
C
C     Input:
C       xi, eta, zeta -- natural coordinates
C     Output:
C       N(4)       -- shape function values
C       dNdxi(4,3) -- dNdxi(a,i) = dN_a/dxi_i (constant for the tet)
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: xi, eta, zeta
      DOUBLE PRECISION, INTENT(OUT) :: N(4), dNdxi(4,3)

      N(1) = 1.0d0 - xi - eta - zeta
      N(2) = xi
      N(3) = eta
      N(4) = zeta

      dNdxi(1,1) = -1.0d0
      dNdxi(1,2) = -1.0d0
      dNdxi(1,3) = -1.0d0
      dNdxi(2,1) =  1.0d0
      dNdxi(2,2) =  0.0d0
      dNdxi(2,3) =  0.0d0
      dNdxi(3,1) =  0.0d0
      dNdxi(3,2) =  1.0d0
      dNdxi(3,3) =  0.0d0
      dNdxi(4,1) =  0.0d0
      dNdxi(4,2) =  0.0d0
      dNdxi(4,3) =  1.0d0

      RETURN
      END SUBROUTINE shape_tet4


      SUBROUTINE gauss_tet4(xi_gp, w_gp, n_gp)
C     4-point Gauss rule for the reference tetrahedron (degree 2)
C
C     Output:
C       xi_gp(4,3) -- integration point coordinates
C       w_gp(4)    -- weights (sum = 1/6 = reference tet volume)
C       n_gp       -- 4
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(4,3), w_gp(4)
      INTEGER, INTENT(OUT) :: n_gp

      DOUBLE PRECISION :: a, b, w
      INTEGER :: ii

      n_gp = 4
      a = 0.585410196624968454d0
      b = 0.138196601125010515d0
      w = 1.0d0 / 24.0d0

      xi_gp(1,1) = a
      xi_gp(1,2) = b
      xi_gp(1,3) = b
      xi_gp(2,1) = b
      xi_gp(2,2) = a
      xi_gp(2,3) = b
      xi_gp(3,1) = b
      xi_gp(3,2) = b
      xi_gp(3,3) = a
      xi_gp(4,1) = b
      xi_gp(4,2) = b
      xi_gp(4,3) = b

      DO ii = 1, 4
        w_gp(ii) = w
      END DO

      RETURN
      END SUBROUTINE gauss_tet4


      SUBROUTINE gauss_tet1(xi_gp, w_gp, n_gp)
C     1-point (centroid) rule for the reference tetrahedron (degree 1).
C     This is the single-point quadrature highlighted by Scovazzi et al.
C     (2023) for the stabilized linear tet.
C
C     Output:
C       xi_gp(1,3) -- centroid (1/4, 1/4, 1/4)
C       w_gp(1)    -- 1/6 (reference tet volume)
C       n_gp       -- 1
C
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: xi_gp(1,3), w_gp(1)
      INTEGER, INTENT(OUT) :: n_gp

      n_gp = 1
      xi_gp(1,1) = 0.25d0
      xi_gp(1,2) = 0.25d0
      xi_gp(1,3) = 0.25d0
      w_gp(1) = 1.0d0 / 6.0d0

      RETURN
      END SUBROUTINE gauss_tet1

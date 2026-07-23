C======================================================================
C     isoparametric.for — Isoparametric mapping (2D and 3D)
C
C     Contains:
C       map_grad_2d:    Map shape function gradients (2D) + return Jinv
C       apply_jinv:     Map arbitrary dshxi using pre-computed Jinv (2D)
C       map_grad_3d:    Map shape function gradients (3D) + return Jinv
C       apply_jinv_3d:  Map arbitrary dshxi using pre-computed Jinv (3D)
C
C     For mixed-degree elements:
C       1. Call map_grad_Nd with full geometry to get dsh, detJ, Jinv
C       2. Call apply_jinv[_3d] with the SAME Jinv to map lower-degree
C          shape function derivatives
C======================================================================

C----------------------------------------------------------------------
C     map_grad_2d: Map shape function gradients from (xi,eta) to (X,Y)
C
C     Given:
C       dshxi(nNode,2)  — derivatives wrt xi, eta
C       coords(2,nNode) — nodal coordinates (X,Y for each node)
C       nNode           — number of nodes
C
C     Returns:
C       dsh(nNode,2)    — derivatives wrt X, Y
C       detJ            — determinant of the Jacobian
C       Jinv_out(2,2)   — inverse Jacobian (for reuse with apply_jinv)
C       stat            — 0 if detJ <= 0 (degenerate element)
C----------------------------------------------------------------------
      SUBROUTINE map_grad_2d(dshxi, coords, nNode, dsh, detJ,
     &                       Jinv_out, stat)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nNode
      DOUBLE PRECISION, INTENT(IN)  :: dshxi(nNode,2), coords(2,nNode)
      DOUBLE PRECISION, INTENT(OUT) :: dsh(nNode,2), detJ
      DOUBLE PRECISION, INTENT(OUT) :: Jinv_out(2,2)
      INTEGER, INTENT(OUT) :: stat

      DOUBLE PRECISION :: Jac(2,2), detJinv
      INTEGER :: i, k

C     Compute Jacobian: J(i,j) = sum_a dN_a/dxi_i * X_j^a
      Jac(1,1) = 0.0d0
      Jac(1,2) = 0.0d0
      Jac(2,1) = 0.0d0
      Jac(2,2) = 0.0d0
      DO i = 1, 2
        DO k = 1, nNode
          Jac(i,1) = Jac(i,1) + dshxi(k,i) * coords(1,k)
          Jac(i,2) = Jac(i,2) + dshxi(k,i) * coords(2,k)
        END DO
      END DO

C     Determinant
      detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

      IF (detJ .LE. 0.0d0) THEN
        stat = 0
        Jinv_out(1,1) = 0.0d0
        Jinv_out(1,2) = 0.0d0
        Jinv_out(2,1) = 0.0d0
        Jinv_out(2,2) = 0.0d0
        RETURN
      END IF

      stat = 1

C     Inverse of 2x2 Jacobian
      detJinv = 1.0d0 / detJ
      Jinv_out(1,1) =  Jac(2,2) * detJinv
      Jinv_out(1,2) = -Jac(1,2) * detJinv
      Jinv_out(2,1) = -Jac(2,1) * detJinv
      Jinv_out(2,2) =  Jac(1,1) * detJinv

C     Map: dN/dX = Jinv * dN/dxi
      DO k = 1, nNode
        dsh(k,1) = Jinv_out(1,1)*dshxi(k,1)+Jinv_out(1,2)*dshxi(k,2)
        dsh(k,2) = Jinv_out(2,1)*dshxi(k,1)+Jinv_out(2,2)*dshxi(k,2)
      END DO

      RETURN
      END SUBROUTINE map_grad_2d

C----------------------------------------------------------------------
C     apply_jinv: Map shape function derivatives using pre-computed Jinv
C
C     For mixed-degree elements: use the Jacobian inverse from the full
C     geometry mapping to map lower-degree shape function derivatives.
C
C     Given:
C       dshxi(nNode,2)  — derivatives wrt xi, eta
C       Jinv(2,2)       — pre-computed Jacobian inverse
C       nNode           — number of nodes
C
C     Returns:
C       dsh(nNode,2)    — derivatives wrt X, Y
C----------------------------------------------------------------------
      SUBROUTINE apply_jinv(dshxi, Jinv, nNode, dsh)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nNode
      DOUBLE PRECISION, INTENT(IN)  :: dshxi(nNode,2), Jinv(2,2)
      DOUBLE PRECISION, INTENT(OUT) :: dsh(nNode,2)
      INTEGER :: k

      DO k = 1, nNode
        dsh(k,1) = Jinv(1,1)*dshxi(k,1) + Jinv(1,2)*dshxi(k,2)
        dsh(k,2) = Jinv(2,1)*dshxi(k,1) + Jinv(2,2)*dshxi(k,2)
      END DO

      RETURN
      END SUBROUTINE apply_jinv

C----------------------------------------------------------------------
C     map_grad_3d: Map shape function gradients from (xi,eta,zeta)
C                  to (X,Y,Z)
C
C     Given:
C       dshxi(nNode,3)  — derivatives wrt xi, eta, zeta
C       coords(3,nNode) — nodal coordinates (X,Y,Z for each node)
C       nNode           — number of nodes
C
C     Returns:
C       dsh(nNode,3)    — derivatives wrt X, Y, Z
C       detJ            — determinant of the Jacobian
C       Jinv_out(3,3)   — inverse Jacobian (for reuse with apply_jinv_3d)
C       stat            — 0 if detJ <= 0 (degenerate element)
C----------------------------------------------------------------------
      SUBROUTINE map_grad_3d(dshxi, coords, nNode, dsh, detJ,
     &                       Jinv_out, stat)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nNode
      DOUBLE PRECISION, INTENT(IN)  :: dshxi(nNode,3), coords(3,nNode)
      DOUBLE PRECISION, INTENT(OUT) :: dsh(nNode,3), detJ
      DOUBLE PRECISION, INTENT(OUT) :: Jinv_out(3,3)
      INTEGER, INTENT(OUT) :: stat

      DOUBLE PRECISION :: Jac(3,3)
      DOUBLE PRECISION :: det33d
      INTEGER :: i, j, k

C     Compute Jacobian: J(i,j) = sum_a dN_a/dxi_i * X_j^a
      DO i = 1, 3
        DO j = 1, 3
          Jac(i,j) = 0.0d0
          DO k = 1, nNode
            Jac(i,j) = Jac(i,j) + dshxi(k,i) * coords(j,k)
          END DO
        END DO
      END DO

C     Determinant
      detJ = det33d(Jac)

      IF (detJ .LE. 0.0d0) THEN
        stat = 0
        DO i = 1, 3
          DO j = 1, 3
            Jinv_out(i,j) = 0.0d0
          END DO
        END DO
        RETURN
      END IF

      stat = 1

C     Inverse of 3x3 Jacobian
      CALL inv33d(Jac, Jinv_out)

C     Map: dN/dX_j = sum_i Jinv(j,i) * dN/dxi_i
      DO k = 1, nNode
        DO j = 1, 3
          dsh(k,j) = 0.0d0
          DO i = 1, 3
            dsh(k,j) = dsh(k,j) + Jinv_out(j,i) * dshxi(k,i)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE map_grad_3d

C----------------------------------------------------------------------
C     apply_jinv_3d: Map shape function derivatives using
C                    pre-computed 3x3 Jinv
C
C     Given:
C       dshxi(nNode,3)  — derivatives wrt xi, eta, zeta
C       Jinv(3,3)       — pre-computed Jacobian inverse
C       nNode           — number of nodes
C
C     Returns:
C       dsh(nNode,3)    — derivatives wrt X, Y, Z
C----------------------------------------------------------------------
      SUBROUTINE apply_jinv_3d(dshxi, Jinv, nNode, dsh)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nNode
      DOUBLE PRECISION, INTENT(IN)  :: dshxi(nNode,3), Jinv(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: dsh(nNode,3)
      INTEGER :: k, i, j

      DO k = 1, nNode
        DO j = 1, 3
          dsh(k,j) = 0.0d0
          DO i = 1, 3
            dsh(k,j) = dsh(k,j) + Jinv(j,i) * dshxi(k,i)
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE apply_jinv_3d

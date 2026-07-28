C======================================================================
C     tangent_identities.for -- Tensor derivative identity library
C
C     Layer 2 of the symbolic/complex-step tangent system.
C     Pure DOUBLE PRECISION -- no complex arithmetic.
C
C     These are standard closed-form tensor calculus results reused
C     by every model that needs them.
C======================================================================

C-----------------------------------------------------------------------
C     ddetdF33: dJ/dF(i,J) = J * Finv(J,i)
C-----------------------------------------------------------------------
      SUBROUTINE ddetdF33(F, J, dJdF)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3), J
      DOUBLE PRECISION, INTENT(OUT) :: dJdF(3,3)
      DOUBLE PRECISION :: Finv(3,3)
      INTEGER :: ii, JJ1

      CALL inv33d(F, Finv)

      DO ii = 1, 3
        DO JJ1 = 1, 3
          dJdF(ii,JJ1) = J * Finv(JJ1,ii)
        END DO
      END DO

      RETURN
      END SUBROUTINE ddetdF33

C-----------------------------------------------------------------------
C     dlndetdF33: d(ln J)/dF(i,J) = Finv(J,i)
C-----------------------------------------------------------------------
      SUBROUTINE dlndetdF33(F, dlnJdF)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: dlnJdF(3,3)
      DOUBLE PRECISION :: Finv(3,3)
      INTEGER :: ii, JJ1

      CALL inv33d(F, Finv)

      DO ii = 1, 3
        DO JJ1 = 1, 3
          dlnJdF(ii,JJ1) = Finv(JJ1,ii)
        END DO
      END DO

      RETURN
      END SUBROUTINE dlndetdF33

C-----------------------------------------------------------------------
C     dFinvTdF33: d(F^{-T})/dF(i,J,k,L) = -Finv(J,k)*Finv(L,i)
C-----------------------------------------------------------------------
      SUBROUTINE dFinvTdF33(Finv, T4)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: Finv(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: T4(3,3,3,3)
      INTEGER :: ii, JJ1, kk, LL1

      DO ii = 1, 3
        DO JJ1 = 1, 3
          DO kk = 1, 3
            DO LL1 = 1, 3
              T4(ii,JJ1,kk,LL1) = -Finv(JJ1,kk) * Finv(LL1,ii)
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE dFinvTdF33

C-----------------------------------------------------------------------
C     dCinvdC33: d(C^{-1})/dC(I,J,K,L)
C     = -0.5*(Cinv(I,K)*Cinv(J,L) + Cinv(I,L)*Cinv(J,K))
C     for symmetric C.
C-----------------------------------------------------------------------
      SUBROUTINE dCinvdC33(Cinv, T4)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: Cinv(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: T4(3,3,3,3)
      INTEGER :: II1, JJ1, KK, LL1

      DO II1 = 1, 3
        DO JJ1 = 1, 3
          DO KK = 1, 3
            DO LL1 = 1, 3
              T4(II1,JJ1,KK,LL1) = -0.5d0 * (
     &            Cinv(II1,KK) * Cinv(JJ1,LL1)
     &          + Cinv(II1,LL1) * Cinv(JJ1,KK))
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE dCinvdC33

C-----------------------------------------------------------------------
C     dI1dF33: d(tr(F^T F))/dF = 2*F
C-----------------------------------------------------------------------
      SUBROUTINE dI1dF33(F, result)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: result(3,3)
      INTEGER :: ii, JJ1

      DO ii = 1, 3
        DO JJ1 = 1, 3
          result(ii,JJ1) = 2.0d0 * F(ii,JJ1)
        END DO
      END DO

      RETURN
      END SUBROUTINE dI1dF33

C-----------------------------------------------------------------------
C     dI2dF33: d(I2)/dF where I2 = 0.5*(I1^2 - tr(C^2))
C     Result: 2*(I1*F - F*C)
C-----------------------------------------------------------------------
      SUBROUTINE dI2dF33(F, result)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: result(3,3)
      DOUBLE PRECISION :: FT(3,3), C(3,3), FC(3,3), I1
      INTEGER :: ii, JJ1

      CALL transpose33d(F, FT)
      CALL matmul33d(FT, F, C)
      CALL matmul33d(F, C, FC)
      I1 = C(1,1) + C(2,2) + C(3,3)

      DO ii = 1, 3
        DO JJ1 = 1, 3
          result(ii,JJ1) = 2.0d0 * (I1 * F(ii,JJ1) - FC(ii,JJ1))
        END DO
      END DO

      RETURN
      END SUBROUTINE dI2dF33

C-----------------------------------------------------------------------
C     pushforward_PK1_to_spatial:
C     c_spatial(i,j,k,l) = (1/J) * sum(J1,L1)
C                          dPdF(i,J1,k,L1) * F(j,J1) * F(l,L1)
C-----------------------------------------------------------------------
      SUBROUTINE pushforward_PK1_to_spatial(F, J, dPdF, c_spatial)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3), J, dPdF(3,3,3,3)
      DOUBLE PRECISION, INTENT(OUT) :: c_spatial(3,3,3,3)
      DOUBLE PRECISION :: Jinv, AFF
      INTEGER :: ii, jj, kk, ll, JJ1, LL1

      Jinv = 1.0d0 / J

      DO ii = 1, 3
        DO jj = 1, 3
          DO kk = 1, 3
            DO ll = 1, 3
              AFF = 0.0d0
              DO JJ1 = 1, 3
                DO LL1 = 1, 3
                  AFF = AFF
     &              + dPdF(ii,JJ1,kk,LL1) * F(jj,JJ1) * F(ll,LL1)
                END DO
              END DO
              c_spatial(ii,jj,kk,ll) = AFF * Jinv
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE pushforward_PK1_to_spatial

C-----------------------------------------------------------------------
C     jaumann_correction:
C     Add Jaumann rate terms to spatial tangent:
C     c_J(i,j,k,l) = +sigma(i,k) if j==l
C                    -sigma(i,j) if k==l
C-----------------------------------------------------------------------
      SUBROUTINE jaumann_correction(sigma, c_J)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: sigma(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: c_J(3,3,3,3)
      INTEGER :: ii, jj, kk, ll

      DO ii = 1, 3
        DO jj = 1, 3
          DO kk = 1, 3
            DO ll = 1, 3
              c_J(ii,jj,kk,ll) = 0.0d0
              IF (jj .EQ. ll) THEN
                c_J(ii,jj,kk,ll) = c_J(ii,jj,kk,ll) + sigma(ii,kk)
              END IF
              IF (kk .EQ. ll) THEN
                c_J(ii,jj,kk,ll) = c_J(ii,jj,kk,ll) - sigma(ii,jj)
              END IF
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE jaumann_correction

C-----------------------------------------------------------------------
C     voigt_pack66:
C     Pack 4th-order spatial tangent into 6x6 Voigt DDSDDE.
C     Abaqus ordering: 11, 22, 33, 12, 13, 23
C     (kl)-symmetrized, major-symmetrized.
C-----------------------------------------------------------------------
      SUBROUTINE voigt_pack66(c_spatial, D66)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: c_spatial(3,3,3,3)
      DOUBLE PRECISION, INTENT(OUT) :: D66(6,6)
      INTEGER :: i, j
      INTEGER :: voigt_i(6), voigt_j(6)
      DATA voigt_i /1, 2, 3, 1, 1, 2/
      DATA voigt_j /1, 2, 3, 2, 3, 3/

      DO i = 1, 6
        DO j = i, 6
          D66(i,j) = 0.25d0 * (
     &      c_spatial(voigt_i(i),voigt_j(i),voigt_i(j),voigt_j(j))
     &    + c_spatial(voigt_i(i),voigt_j(i),voigt_j(j),voigt_i(j))
     &    + c_spatial(voigt_j(i),voigt_i(i),voigt_i(j),voigt_j(j))
     &    + c_spatial(voigt_j(i),voigt_i(i),voigt_j(j),voigt_i(j)))
          D66(j,i) = D66(i,j)
        END DO
      END DO

      RETURN
      END SUBROUTINE voigt_pack66

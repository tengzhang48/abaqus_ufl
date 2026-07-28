C======================================================================
C     cs_tangent_engine.for -- Complex-step tangent computation
C
C     Computes all tangent blocks for the three-field (u, p, mu)
C     gel formulation by perturbing material function inputs.
C
C     Tangent blocks computed (12 total, 45 perturbations):
C       dP/dF (3,3,3,3), dP/dp (3,3), [dP/dmu skipped: known zero]
C       drp/dF (3,3), drp/dp (scalar), drp/dmu (scalar)
C       djR/dF (3,3,3), djR/dp (3), djR/dmu (3), djR/dgmu (3,3)
C       dcdot/dF (3,3), dcdot/dp (scalar)
C
C     All inputs DOUBLE PRECISION. All outputs DOUBLE PRECISION.
C     Complex arithmetic is internal only.
C
C     Portability: AIMAG (not DIMAG), DBLE, DCMPLX only.
C
C     Safety: Every block resets ALL current-step complex variables
C       unconditionally (Fz, pz, muz, grad_mu_z), even if a given
C       block does not use all of them. This guarantees no cross-block
C       contamination regardless of block ordering or conditional
C       skipping by a code generator.
C======================================================================

      SUBROUTINE cs_gel_tangents(F, p, mu, theta, grad_mu,
     &                           F_old, p_old, props, dt,
     &                           dPdF, dPdp,
     &                           drpdF, drpdp, drpdmu,
     &                           djRdF, djRdp, djRdmu, djRdgmu,
     &                           dcdotdF, dcdotdp)

      IMPLICIT NONE

C     --- Inputs (all real) ---
      DOUBLE PRECISION, INTENT(IN) :: F(3,3), p, mu, theta
      DOUBLE PRECISION, INTENT(IN) :: grad_mu(3)
      DOUBLE PRECISION, INTENT(IN) :: F_old(3,3), p_old
      DOUBLE PRECISION, INTENT(IN) :: props(*)
      DOUBLE PRECISION, INTENT(IN) :: dt

C     --- Outputs (all real tangent values) ---
      DOUBLE PRECISION, INTENT(OUT) :: dPdF(3,3,3,3)
      DOUBLE PRECISION, INTENT(OUT) :: dPdp(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: drpdF(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: drpdp
      DOUBLE PRECISION, INTENT(OUT) :: drpdmu
      DOUBLE PRECISION, INTENT(OUT) :: djRdF(3,3,3)
      DOUBLE PRECISION, INTENT(OUT) :: djRdp(3)
      DOUBLE PRECISION, INTENT(OUT) :: djRdmu(3)
      DOUBLE PRECISION, INTENT(OUT) :: djRdgmu(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: dcdotdF(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: dcdotdp

C     --- Complex-step size ---
      DOUBLE PRECISION, PARAMETER :: CS_H = 1.0d-10

C     --- Complex working variables ---
      DOUBLE COMPLEX :: Fz(3,3), pz, muz, thetaz
      DOUBLE COMPLEX :: grad_mu_z(3)
      DOUBLE COMPLEX :: F_old_z(3,3), p_old_z

C     --- Perturbed outputs ---
      DOUBLE COMPLEX :: Pz(3,3), phiz, rpz, jRz(3), cdotz

C     --- Loop counters ---
      INTEGER :: k, l, i, j

C     --- One-time conversion of time-invariant inputs ---
      CALL real2complex33(F_old, F_old_z)
      p_old_z = DCMPLX(p_old, 0.0d0)
      thetaz  = DCMPLX(theta, 0.0d0)

C======================================================================
C     BLOCK 1: dP/dF (3,3,3,3) -- 9 perturbations of F
C======================================================================
      DO k = 1, 3
        DO l = 1, 3
          CALL real2complex33(F, Fz)
          pz  = DCMPLX(p, 0.0d0)
          muz = DCMPLX(mu, 0.0d0)
          CALL real2complex3(grad_mu, grad_mu_z)
          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)
          CALL gel_stress_PK1(Fz, pz, muz, thetaz, props, Pz, phiz)
          DO i = 1, 3
            DO j = 1, 3
              dPdF(i,j,k,l) = AIMAG(Pz(i,j)) / CS_H
            END DO
          END DO
        END DO
      END DO

C======================================================================
C     BLOCK 2: dP/dp (3,3) -- 1 perturbation of p
C======================================================================
      CALL real2complex33(F, Fz)
      pz  = DCMPLX(p, CS_H)
      muz = DCMPLX(mu, 0.0d0)
      CALL real2complex3(grad_mu, grad_mu_z)
      CALL gel_stress_PK1(Fz, pz, muz, thetaz, props, Pz, phiz)
      DO i = 1, 3
        DO j = 1, 3
          dPdp(i,j) = AIMAG(Pz(i,j)) / CS_H
        END DO
      END DO

C     BLOCK 3: dP/dmu -- SKIPPED (known zero for this model)

C======================================================================
C     BLOCK 4: drp/dF (3,3) -- 9 perturbations of F
C======================================================================
      DO k = 1, 3
        DO l = 1, 3
          CALL real2complex33(F, Fz)
          pz  = DCMPLX(p, 0.0d0)
          muz = DCMPLX(mu, 0.0d0)
          CALL real2complex3(grad_mu, grad_mu_z)
          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)
          CALL gel_pressure_resid(Fz, pz, muz, thetaz, props, rpz)
          drpdF(k,l) = AIMAG(rpz) / CS_H
        END DO
      END DO

C======================================================================
C     BLOCK 5: drp/dp (scalar) -- 1 perturbation of p
C======================================================================
      CALL real2complex33(F, Fz)
      pz  = DCMPLX(p, CS_H)
      muz = DCMPLX(mu, 0.0d0)
      CALL real2complex3(grad_mu, grad_mu_z)
      CALL gel_pressure_resid(Fz, pz, muz, thetaz, props, rpz)
      drpdp = AIMAG(rpz) / CS_H

C======================================================================
C     BLOCK 6: drp/dmu (scalar) -- 1 perturbation of mu
C======================================================================
      CALL real2complex33(F, Fz)
      pz  = DCMPLX(p, 0.0d0)
      muz = DCMPLX(mu, CS_H)
      CALL real2complex3(grad_mu, grad_mu_z)
      CALL gel_pressure_resid(Fz, pz, muz, thetaz, props, rpz)
      drpdmu = AIMAG(rpz) / CS_H

C======================================================================
C     BLOCK 7: djR/dF (3,3,3) -- 9 perturbations of F
C======================================================================
      DO k = 1, 3
        DO l = 1, 3
          CALL real2complex33(F, Fz)
          pz  = DCMPLX(p, 0.0d0)
          muz = DCMPLX(mu, 0.0d0)
          CALL real2complex3(grad_mu, grad_mu_z)
          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)
          CALL gel_flux_ref(Fz, pz, muz, grad_mu_z, thetaz,
     &                      props, jRz)
          DO i = 1, 3
            djRdF(i,k,l) = AIMAG(jRz(i)) / CS_H
          END DO
        END DO
      END DO

C======================================================================
C     BLOCK 8: djR/dp (3) -- 1 perturbation of p
C======================================================================
      CALL real2complex33(F, Fz)
      pz  = DCMPLX(p, CS_H)
      muz = DCMPLX(mu, 0.0d0)
      CALL real2complex3(grad_mu, grad_mu_z)
      CALL gel_flux_ref(Fz, pz, muz, grad_mu_z, thetaz, props, jRz)
      DO i = 1, 3
        djRdp(i) = AIMAG(jRz(i)) / CS_H
      END DO

C======================================================================
C     BLOCK 9: djR/dmu (3) -- 1 perturbation of mu
C======================================================================
      CALL real2complex33(F, Fz)
      pz  = DCMPLX(p, 0.0d0)
      muz = DCMPLX(mu, CS_H)
      CALL real2complex3(grad_mu, grad_mu_z)
      CALL gel_flux_ref(Fz, pz, muz, grad_mu_z, thetaz, props, jRz)
      DO i = 1, 3
        djRdmu(i) = AIMAG(jRz(i)) / CS_H
      END DO

C======================================================================
C     BLOCK 10: djR/d(grad_mu) (3,3) -- 3 perturbations of grad_mu
C======================================================================
      DO k = 1, 3
        CALL real2complex33(F, Fz)
        pz  = DCMPLX(p, 0.0d0)
        muz = DCMPLX(mu, 0.0d0)
        CALL real2complex3(grad_mu, grad_mu_z)
        grad_mu_z(k) = grad_mu_z(k) + DCMPLX(0.0d0, CS_H)
        CALL gel_flux_ref(Fz, pz, muz, grad_mu_z, thetaz, props, jRz)
        DO i = 1, 3
          djRdgmu(i,k) = AIMAG(jRz(i)) / CS_H
        END DO
      END DO

C======================================================================
C     BLOCK 11: dcdot/dF (3,3) -- 9 perturbations of F
C======================================================================
      DO k = 1, 3
        DO l = 1, 3
          CALL real2complex33(F, Fz)
          pz  = DCMPLX(p, 0.0d0)
          muz = DCMPLX(mu, 0.0d0)
          CALL real2complex3(grad_mu, grad_mu_z)
          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)
          CALL gel_conc_rate(Fz, F_old_z, pz, p_old_z, props, dt,
     &                       cdotz)
          dcdotdF(k,l) = AIMAG(cdotz) / CS_H
        END DO
      END DO

C======================================================================
C     BLOCK 12: dcdot/dp (scalar) -- 1 perturbation of p
C======================================================================
      CALL real2complex33(F, Fz)
      pz  = DCMPLX(p, CS_H)
      muz = DCMPLX(mu, 0.0d0)
      CALL real2complex3(grad_mu, grad_mu_z)
      CALL gel_conc_rate(Fz, F_old_z, pz, p_old_z, props, dt, cdotz)
      dcdotdp = AIMAG(cdotz) / CS_H

      RETURN
      END SUBROUTINE cs_gel_tangents

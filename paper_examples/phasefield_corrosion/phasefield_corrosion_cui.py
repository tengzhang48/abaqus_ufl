"""
Cui-style uniform-Q8R phase-field corrosion UEL.

This is the first abaqus_ufl step toward reproducing the authors'
shared ``PhaseFieldSCC.f`` reference UEL.  It matches the original
UEL's interpolation and quadrature choice:

  - Q8 displacement, phase field, and concentration
  - 32 element DOFs
  - 2x2 reduced integration (``Quad8R``)

The material includes the small-strain J2 update and the fatigue/repassivation
cycle state variables used by ``PhaseFieldSCC.f``.  The mech-phase coupling is
explicit through ``xL_old`` in the phase equation, matching the forward-Euler
implementation option described in Cui's documentation.
"""

import os
import shutil
import sys
from pathlib import Path

project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)

import numpy as np

import abaqus_ufl as au
from abaqus_ufl.core.small_strain_plasticity import flow_direction, q_mises
from abaqus_ufl.core.tensor import exp, eye, trace
from abaqus_ufl.generators.uel_gen import generate_uel


class CuiJ2CorrosionMaterial(au.Material):

    # Residual damage stiffness kappa_r: a declaration constant
    # inlined into the generated Fortran (partial-drop default).
    xkap_residual = 1.0e-3
    """Cui-style corrosion model with small-strain J2 plasticity.

    The J2 update follows Cui's power-law hardening relation:
    ``sigma_f = sigma_y * (1 + E * ep / sigma_y)**N``.
    """

    props = dict(
        E=190000.0,
        nu=0.3,
        sigma_y=520.0,
        hardening_n=0.067,
        D=8.5e-4,
        L0=1.0e-3,
        kappa=5.1e-5,
        omega=35.3,
        Achem=53.5,
        k_repassivation=5.0e-4,
        eps_f=3.0e-3,
        t0=10.0,
        c_solid=1.0,
        c_liquid=0.036,
    )
    state_vars = dict(
        ep=0.0,
        epsp=np.zeros((3, 3)),
        deqpl=0.0,
        hydro=0.0,
        xL=1.0e-3,
        ti=0.0,
        ei=0.0,
    )
    state_var_props = dict(xL="L0")

    def stress_PK1(self, F, phi, ep_old, epsp_old, deqpl_old, hydro_old,
                   xL_old, ti_old, ei_old, dt):
        """Small-strain stress returned in the PK1 slot.

        This implements the UEL-side analogue of Cui's `kumat`: a small-strain
        radial-return J2 update with isotropic hardening.  The imposed strains
        in the benchmark are small, so the Cauchy stress and PK1 stress differ
        only by higher-order terms.
        """
        eps = 0.5 * (F + F.T) - eye(3)
        mu = self.E / (2.0 * (1.0 + self.nu))
        lam = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))

        eps_e_trial = eps - epsp_old
        sigma_trial = lam * trace(eps_e_trial) * eye(3) + 2.0 * mu * eps_e_trial
        seq = q_mises(sigma_trial)
        yield_old = self.sigma_y * (
            1.0 + self.E * ep_old / self.sigma_y
        ) ** self.hardening_n
        f_trial = seq - yield_old

        if f_trial.real > 0.0:
            deqpl = 0.0 * f_trial
            yield_new = yield_old
            tangent_hard = self.E * self.hardening_n * (
                1.0 + self.E * ep_old / self.sigma_y
            ) ** (self.hardening_n - 1.0)
            for iteration in range(20):
                resid = seq - 3.0 * mu * deqpl - yield_new
                deqpl = deqpl + resid / (3.0 * mu + tangent_hard)
                yield_new = self.sigma_y * (
                    1.0 + self.E * (ep_old + deqpl) / self.sigma_y
                ) ** self.hardening_n
                tangent_hard = self.E * self.hardening_n * (
                    1.0 + self.E * (ep_old + deqpl) / self.sigma_y
                ) ** (self.hardening_n - 1.0)
            n = flow_direction(sigma_trial, seq)
            sigma0 = sigma_trial - 2.0 * mu * deqpl * n
            ep_new = ep_old + deqpl
            epsp_new = epsp_old + deqpl * n
        else:
            sigma0 = sigma_trial
            ep_new = ep_old
            epsp_new = epsp_old
            deqpl = 0.0 * f_trial

        degr = phi * phi * (3.0 - 2.0 * phi)
        xkap = self.xkap_residual
        sigma_damaged = (degr + xkap) * sigma0
        hydro_new = trace(sigma_damaged) / 3.0
        gas_const = 8314.0
        temp_abs = 300.0
        mech_factor = exp(hydro_new * 7.12e3 / (gas_const * temp_abs)) * (
            1.0 + ep_new / (self.sigma_y / self.E)
        )

        ei_trial = ei_old + deqpl
        if ei_trial.real > self.eps_f.real:
            ti_cycle = 0.0 * ti_old
            ei_cycle = 0.0 * ei_old
        else:
            ti_cycle = ti_old + dt
            ei_cycle = ei_trial

        if ti_cycle.real < self.t0.real:
            repassivation = 1.0 + 0.0 * ti_cycle
        else:
            repassivation = exp(
                -self.k_repassivation * (ti_cycle - self.t0))
        xL_new = self.L0 * mech_factor * repassivation
        return sigma_damaged, {
            'ep': ep_new,
            'epsp': epsp_new,
            'deqpl': deqpl,
            'hydro': hydro_new,
            'xL': xL_new,
            'ti': ti_cycle,
            'ei': ei_cycle,
        }

    def phase_storage(self, F, phi, c, phi_old, xL_old, dt):
        h = phi * phi * (3.0 - 2.0 * phi)
        dh = 6.0 * phi * (1.0 - phi)
        dc_eq = self.c_solid - self.c_liquid
        c_mix = self.c_liquid + h * (self.c_solid - self.c_liquid)
        chem_arg = c - c_mix
        dpsi_chem = -2.0 * self.Achem * chem_arg * dc_eq * dh
        dpsi_dw = (
            2.0 * self.omega * phi * (1.0 - phi) * (1.0 - 2.0 * phi)
        )
        phidot = (phi - phi_old) / dt
        return -phidot / xL_old - dpsi_chem - dpsi_dw

    def phase_flux(self, F, phi, grad_phi):
        return self.kappa * grad_phi

    def species_storage(self, F, c, c_old, dt):
        return (c - c_old) / dt

    def species_flux(self, F, phi, grad_phi, grad_c):
        dh = 6.0 * phi * (1.0 - phi)
        dc_eq = self.c_solid - self.c_liquid
        return -self.D * (grad_c - dc_eq * dh * grad_phi)


class CuiJ2Corrosion(au.WeakForm):
    """Uniform-Q8R ``u, phi, c`` corrosion problem."""

    material = CuiJ2CorrosionMaterial
    ndim = 2

    # Staggered stiffness (residual stays fully coupled).
    # Drop only the mechanics-damage block K_uphi = d(R_u)/d(phi) = g'(phi)*sigma0
    # from AMATRX. In the corrosion front this term is large and one-sided
    # (there is no symmetric d(R_phi)/d(u): the mechanics->phase coupling enters
    # only through the lagged mobility xL_old), so it makes the monolithic tangent
    # unsymmetric/indefinite and inflates the displacement correction. Removing it
    # from the Jacobian only changes the iteration path, not the converged fields
    # (sigma = g(phi)*sigma0 is untouched in the residual). To recover the full
    # monolithic (consistent) tangent for comparison, comment this out.
    drop_tangent_coupling = [('momentum_equation', 'phi')]

    def define_fields(self):
        self.u = au.VectorField("u", degree=2)
        self.phi = au.ScalarField("phi", degree=2, test="eta")
        self.c = au.ScalarField("c", degree=2, test="zeta")

    def momentum_equation(self, v, F, phi):
        return self.material.stress_PK1(F, phi)

    def phase_equation(self, eta, F, phi, c, grad_phi, phi_old, xL_old, dt):
        return (
            self.material.phase_storage(F, phi, c, phi_old, xL_old, dt),
            self.material.phase_flux(F, phi, grad_phi),
        )

    def species_transport_equation(self, zeta, F, phi, c, grad_phi,
                                   grad_c, c_old, dt):
        return (
            self.material.species_storage(F, c, c_old, dt),
            self.material.species_flux(F, phi, grad_phi, grad_c),
        )


def verification_state():
    F = np.array([[1.0015, 0.0002, 0.0],
                  [0.0001, 1.0005, 0.0],
                  [0.0, 0.0, 1.0]])
    return dict(
        F=F,
        phi=0.72,
        phi_old=0.725,
        c=0.68,
        c_old=0.69,
        grad_phi=np.array([0.12, -0.05, 0.0]),
        grad_c=np.array([-0.03, 0.02, 0.0]),
        ep_old=0.0,
        epsp_old=np.zeros((3, 3)),
        deqpl_old=0.0,
        hydro_old=0.0,
        xL_old=1.0e-3,
        ti_old=0.2,
        ei_old=0.0,
        dt=0.01,
    )


def _strip_trailing_whitespace(path):
    with open(path, "r") as f:
        lines = f.readlines()
    with open(path, "w") as f:
        f.writelines(line.rstrip() + "\n" for line in lines)


def _visualization_prelude():
    """Abaqus visualization bridge shared state for the Cui duplicate mesh."""
    return """      MODULE kvisual
      IMPLICIT NONE
      DOUBLE PRECISION :: UserVar(4,25,70000), PitD = 0.0D0
      INTEGER :: nelem = 0
      SAVE
      END MODULE kvisual

      LOGICAL FUNCTION finite_d(x)
      IMPLICIT NONE
      DOUBLE PRECISION x
      finite_d = (x .EQ. x) .AND. (DABS(x) .LT. 1.0D300)
      RETURN
      END

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
      USE kvisual
      IMPLICIT NONE
      INTEGER LOP, LRESTART, KSTEP, KINC
      DOUBLE PRECISION TIME(2), DTIME
      IF (LOP .EQ. 0) THEN
        CALL MutexInit(1)
        UserVar = 0.0D0
        PitD = 0.0D0
        nelem = 0
      END IF
      RETURN
      END

"""


def _uvarm_bridge():
    """Expose UEL integration-point values on Cui's native visualization mesh."""
    return """
C======================================================================
C     UVARM bridge for Abaqus CPE8R visualization elements in the Cui deck.
C     The UEL writes integration-point data to UserVar and UVARM exposes it
C     on the duplicate native-element mesh. The physics remains in the UEL.
C======================================================================
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     &NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     &JMAC,JMATYP,MATLAYO,LACCFLA)
      USE kvisual
      IMPLICIT NONE
      CHARACTER*80 CMNAME, ORNAME, PART_NAME
      CHARACTER*3 FLGRAY(15)
      INTEGER NUVARM, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      INTEGER NDI, NSHR, MATLAYO, LACCFLA
      INTEGER JMAC(*), JMATYP(*)
      DOUBLE PRECISION UVAR(NUVARM), DIRECT(3,3), T(3,3), TIME(2)
      DOUBLE PRECISION DTIME, COORD(*), ARRAY(15)
      INTEGER JARRAY(15)
      DOUBLE PRECISION UVAL, PITD_VAL
      INTEGER I, NOFFSET, UEL_NPT, LOCAL_NOEL, JRCD
      LOGICAL finite_d

C     In an assembly Abaqus passes its internal element number to UVARM.
C     Recover the part-level duplicate label before removing the 10000
C     visualization offset. Flat smoke decks use the label directly.
      PART_NAME = ' '
      LOCAL_NOEL = NOEL
      JRCD = 1
      CALL GETPARTINFO(NOEL,1,PART_NAME,LOCAL_NOEL,JRCD)
      IF (JRCD .EQ. 0) THEN
        NOFFSET = LOCAL_NOEL - 10000
      ELSE
        NOFFSET = NOEL - 10000
      END IF
C     The generated Quad8 UEL reverses native CPE8R points 3 and 4.
      UEL_NPT = NPT
      IF (NPT .EQ. 3) UEL_NPT = 4
      IF (NPT .EQ. 4) UEL_NPT = 3
      DO I = 1, NUVARM
        UVAR(I) = 0.0D0
      END DO
      IF (NOFFSET .GE. 1 .AND. NOFFSET .LE. 70000) THEN
        DO I = 1, MIN(NUVARM, 24)
          UVAL = UserVar(UEL_NPT,I,NOFFSET)
          IF (finite_d(UVAL)) UVAR(I) = UVAL
        END DO
        PITD_VAL = DBLE(PitD)
        IF (NUVARM .GE. 25 .AND. finite_d(PITD_VAL)) UVAR(25) =
     &      PITD_VAL
      END IF
      RETURN
      END
"""


def _apply_visualization_bridge(path):
    """Patch generated Cui UEL Fortran with the reproducible UVARM bridge.

    The bridge is intentionally example-local.  It depends on the Cui deck's
    assembly-to-part element mapping and on the specific UVARM slot layout
    used by the comparison scripts.
    """
    path = Path(path)
    text = path.read_text()
    if "SUBROUTINE UVARM(" in text and "MODULE kvisual" in text:
        return

    marker = "      SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL,"
    if marker not in text:
        raise RuntimeError("Could not find UEL entry point for bridge patch")
    text = text.replace(marker, _visualization_prelude() + marker, 1)

    implicit_marker = "\n      IMPLICIT NONE\n\nC     --- Abaqus UEL interface ---"
    if implicit_marker not in text:
        raise RuntimeError("Could not find UEL IMPLICIT NONE for bridge patch")
    text = text.replace(
        implicit_marker,
        "\n      USE kvisual\n      IMPLICIT NONE\n\n"
        "C     --- Abaqus UEL interface ---",
        1,
    )

    decl_marker = (
        "      DOUBLE PRECISION :: species_flux_real(3)\n"
        "      DOUBLE PRECISION :: species_storage_real\n"
    )
    decl_insert = (
        "      DOUBLE PRECISION :: species_flux_real(3)\n"
        "      DOUBLE PRECISION :: species_storage_real\n"
        "      DOUBLE PRECISION :: eps_vis(4), epsp_vis(4), ee_vis(4)\n"
        "      DOUBLE PRECISION :: y_gp, pitd_candidate, pitd_value\n"
    )
    if decl_marker not in text:
        raise RuntimeError("Could not find material-output declarations")
    text = text.replace(decl_marker, decl_insert, 1)

    integer_marker = "      INTEGER :: idx, ii_v, jj_v, i, j, k, l, kk, row, col\n"
    if integer_marker not in text:
        raise RuntimeError("Could not find integer declaration block")
    text = text.replace(integer_marker, integer_marker + "      LOGICAL finite_d\n", 1)

    dtime_marker = (
        "C     DTIME guard (Abaqus calls with DTIME=0 for initial eval)\n"
        "      dt_safe = DTIME\n"
        "      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14\n\n"
    )
    dtime_insert = dtime_marker + (
        "C     Determine UEL element count for duplicate visualization elements.\n"
        "      IF (DTIME .EQ. 0.0D0) THEN\n"
        "        IF (JELEM .GT. nelem) THEN\n"
        "          CALL MutexLock(1)\n"
        "          IF (JELEM .GT. nelem) nelem = JELEM\n"
        "          CALL MutexUnlock(1)\n"
        "        END IF\n"
        "      END IF\n\n"
    )
    if dtime_marker not in text:
        raise RuntimeError("Could not find DTIME guard block")
    text = text.replace(dtime_marker, dtime_insert, 1)

    state_marker = (
        "        SVARS(LOC + 13) = DBLE(xL_new_z)\n"
        "        SVARS(LOC + 14) = DBLE(ti_new_z)\n"
        "        SVARS(LOC + 15) = DBLE(ei_new_z)\n\n"
    )
    uservar_insert = state_marker + (
        "C       Store fields for duplicate CPE8R visualization elements.\n"
        "        eps_vis(1) = F(1,1) - 1.0D0\n"
        "        eps_vis(2) = F(2,2) - 1.0D0\n"
        "        eps_vis(3) = F(3,3) - 1.0D0\n"
        "        eps_vis(4) = 0.5D0 * (F(1,2) + F(2,1))\n"
        "        epsp_vis(1) = DBLE(epsp_new_z(1,1))\n"
        "        epsp_vis(2) = DBLE(epsp_new_z(2,2))\n"
        "        epsp_vis(3) = DBLE(epsp_new_z(3,3))\n"
        "        epsp_vis(4) = DBLE(epsp_new_z(1,2))\n"
        "        DO i = 1, 4\n"
        "          ee_vis(i) = eps_vis(i) - epsp_vis(i)\n"
        "        END DO\n"
        "        UserVar(kk,1,JELEM) = P_real(1,1)\n"
        "        UserVar(kk,2,JELEM) = P_real(2,2)\n"
        "        UserVar(kk,3,JELEM) = P_real(3,3)\n"
        "        UserVar(kk,4,JELEM) = P_real(1,2)\n"
        "        UserVar(kk,5,JELEM) = eps_vis(1)\n"
        "        UserVar(kk,6,JELEM) = eps_vis(2)\n"
        "        UserVar(kk,7,JELEM) = eps_vis(3)\n"
        "        UserVar(kk,8,JELEM) = eps_vis(4)\n"
        "        UserVar(kk,9,JELEM) = ee_vis(1)\n"
        "        UserVar(kk,10,JELEM) = ee_vis(2)\n"
        "        UserVar(kk,11,JELEM) = ee_vis(3)\n"
        "        UserVar(kk,12,JELEM) = ee_vis(4)\n"
        "        UserVar(kk,13,JELEM) = epsp_vis(1)\n"
        "        UserVar(kk,14,JELEM) = epsp_vis(2)\n"
        "        UserVar(kk,15,JELEM) = epsp_vis(3)\n"
        "        UserVar(kk,16,JELEM) = epsp_vis(4)\n"
        "        UserVar(kk,17,JELEM) = DBLE(ep_new_z)\n"
        "        IF (finite_d(phi_gp)) UserVar(kk,18,JELEM) = phi_gp\n"
        "        UserVar(kk,19,JELEM) = DBLE(xL_new_z)\n"
        "        UserVar(kk,20,JELEM) = DBLE(hydro_new_z)\n"
        "        UserVar(kk,21,JELEM) = DBLE(deqpl_new_z)\n"
        "        IF (finite_d(c_gp)) UserVar(kk,22,JELEM) = c_gp\n"
        "        UserVar(kk,23,JELEM) = DBLE(ti_new_z)\n"
        "        UserVar(kk,24,JELEM) = DBLE(ei_new_z)\n"
        "        y_gp = 0.0D0\n"
        "        DO ii_v = 1, 8\n"
        "          y_gp = y_gp + sh8(ii_v) * coords_2d(2,ii_v)\n"
        "        END DO\n"
        "        IF (phi_gp .LE. 0.5D0) THEN\n"
        "          pitd_candidate = 0.13D0 - y_gp\n"
        "          IF (pitd_candidate .GT. PitD) THEN\n"
        "            CALL MutexLock(1)\n"
        "            IF (pitd_candidate .GT. PitD) PitD = pitd_candidate\n"
        "            CALL MutexUnlock(1)\n"
        "          END IF\n"
        "        END IF\n"
        "        IF (TIME(2) .LT. 3.01D0) PitD = 0.0D0\n"
        "        pitd_value = DBLE(PitD)\n"
        "        IF (finite_d(pitd_value)) UserVar(kk,25,JELEM) =\n"
        "     &    pitd_value\n\n"
    )
    if state_marker not in text:
        raise RuntimeError("Could not find SVARS write block")
    text = text.replace(state_marker, uservar_insert, 1)

    text = text.rstrip() + "\n" + _uvarm_bridge()
    path.write_text(text)


def generate(output_path=None):
    problem = CuiJ2Corrosion()
    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(__file__), "phasefield_corrosion_cui_uel.for")
    generate_uel(problem, output_path, element="Quad8R",
                 formulation="standard")
    _apply_visualization_bridge(output_path)
    _strip_trailing_whitespace(output_path)
    default_output = os.path.join(
        os.path.dirname(__file__), "phasefield_corrosion_cui_uel.for")
    if os.path.abspath(output_path) == os.path.abspath(default_output):
        abaqus_copy = os.path.join(
            os.path.dirname(__file__), "abaqus_test_from_cui",
            "phasefield_corrosion_cui_uel.for")
        shutil.copyfile(output_path, abaqus_copy)
    return output_path, problem


class CuiJ2CorrosionDiagMaterial(CuiJ2CorrosionMaterial):
    """Diagonal-comparison material: Cui's residual stiffness 1e-7.

    The kappa_r value is part of the DECLARATION and is inlined into the
    generated source; no post-generation edit is involved.
    """

    xkap_residual = 1.0e-7


class CuiJ2CorrosionDiag(CuiJ2Corrosion):
    """Purely block-diagonal Cui-style tangent variant.

    At the Python-class level, the residual is unchanged and every
    off-diagonal Jacobian block is dropped, leaving ``K_uu``, ``K_phiphi``,
    and ``K_cc``.  The material declaration carries Cui's residual
    stiffness ``xkap = 1e-7`` (vs the shared default ``1e-3``), so the
    plotted production and diagonal runs differ in both tangent
    structure and one residual parameter; they are not a controlled tangent-
    only comparison.  The diagonal structure reproduces Cui's hand-written
    UEL.  The phase/species ``d/dF`` blocks are dropped too (they are
    identically zero, since neither equation uses ``F``).
    """

    material = CuiJ2CorrosionDiagMaterial

    drop_tangent_coupling = [
        ('momentum_equation', 'phi'),                # K_uphi
        ('phase_equation', 'c'),                     # K_phic
        ('phase_equation', 'F'),                     # zero (phase has no F dep.)
        ('species_transport_equation', 'phi'),       # K_cphi (value part)
        ('species_transport_equation', 'grad_phi'),  # K_cphi (gradient part)
        ('species_transport_equation', 'F'),         # zero (species no F dep.)
    ]


def generate_diag(output_path=None):
    """Generate the block-diagonal Cui-style UEL.

    Writes ``phasefield_corrosion_cui_diag_uel.for`` and copies it next to
    the Abaqus decks.  Run it with the ``*_diag`` deck (which drops
    ``unsymm``).  Cui's ``xkap = 1e-7`` comes from the
    :class:`CuiJ2CorrosionDiagMaterial` declaration.
    """
    problem = CuiJ2CorrosionDiag()
    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(__file__),
            "phasefield_corrosion_cui_diag_uel.for")
    generate_uel(problem, output_path, element="Quad8R",
                 formulation="standard")
    _apply_visualization_bridge(output_path)
    _strip_trailing_whitespace(output_path)
    abaqus_copy = os.path.join(
        os.path.dirname(__file__), "abaqus_test_from_cui",
        "phasefield_corrosion_cui_diag_uel.for")
    shutil.copyfile(output_path, abaqus_copy)
    return output_path, problem


if __name__ == "__main__":
    problem = CuiJ2Corrosion()
    problem.summary()
    ok = problem.verify(state=verification_state(), tol=1.0e-5,
                        verbose=True)
    if not ok:
        raise SystemExit(1)
    path, _ = generate()
    print(f"Generated {path}")
    diag_path, _ = generate_diag()
    print(f"Generated {diag_path}")

"""Stabilized mixed u-theta UEL after Scovazzi, Zorrilla & Rossi (2023).

Paper: "A kinematically stabilized linear tetrahedral finite element for
compressible and nearly incompressible finite elasticity", CMAME 412 (2023)
116076.  PDF: ../1-s2.0-S0045782523002001-main.pdf

Formulation (total Lagrangian, Sec. 5):
  - fields: displacement u and a nodal scalar theta = J(u) carried as the
    discrepancy field thetat = theta - 1 (Sec. 6: "theta-tilde"), so all
    DOFs start at zero;
  - the constitutive law is evaluated at the theta-enriched strain
    Cbar = theta^(2/3) J^(-2/3) C   (Eqs. 14, 25);
  - Simo-Taylor neo-Hookean (Eq. 44):
        S(Cbar) = kappa/2 (theta^2 - 1) Cbar^-1
                + mu theta^(-2/3) (I - tr(Cbar)/3 Cbar^-1)        (Eq. 46)
  - momentum stabilization (Eq. 41a):
        P = F [ S - tau_theta (theta - J) dS/dtheta ],
        dS/dtheta = Dbar/(3 theta),
        Dbar = dS/dE : Cbar = kappa (2 theta^2 + 1) Cbar^-1
             - 2 mu theta^(-2/3) (I - tr(Cbar)/3 Cbar^-1)         (Eq. 49)
  - theta equation (Eqs. 41b, 43): reaction (1 - tau_theta)(theta - J)
    plus the VMS gradient term
        ( grad q , (tau_u/3) (J/theta) Dbar grad theta )
    which is exactly an h^2-scaled GRADIENT PENALTY on the equal-order
    (unstable) theta field: it suppresses the checkerboard modes that a
    linear theta/pressure field otherwise develops.  This is the
    strain-gradient reading of the VMS term: an internal length
    ell^2 ~ tau_u * kappa / 1 with tau_u = c_tau_u h^2/(2 mu) (Eq. 33b),
    vanishing as h^2 (consistency).
  - tau_theta = c_tau_theta mu/(mu + kappa) (Eq. 40).

Element routes and deviations (documented in README.md):
  - Tet4 and single-point Tet4R reproduce the paper's P1/P1 element routes;
    auxiliary bilinear Quad4 / trilinear Hex8 cases apply the same stabilized
    weak form to non-P1 element families, where the interior-residual
    simplification is approximate;
  - 2D runs use plane strain with the 3D volumetric/deviatoric split
    (F33 = 1), not the d = 2 split of Remark 9;
  - h enters through the material property h_elem (uniform meshes);
    mu_min = mu (hyperelastic, no history);
  - the body-force stabilization terms are omitted (b0 = 0 in all tests).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

_home_here = (Path.home() / "abaqus_ufl_dev" / "abaqus_ufl_lab" / "examples"
              / "scovazzi_2023_stabilized_utheta")
_pwd_repo = Path(os.environ.get("PWD", os.getcwd()))
_logical_here = _pwd_repo / "examples" / "scovazzi_2023_stabilized_utheta"
if _home_here.is_dir():
    HERE = _home_here
elif _logical_here.is_dir():
    HERE = _logical_here
else:
    HERE = Path(__file__).absolute().parent
REPO_ROOT = HERE.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, eye, inv
from abaqus_ufl.generators.uel_gen import generate_uel


class ScovazziUThetaMaterial(au.Material):
    """Simo-Taylor neo-Hookean evaluated at the theta-enriched strain.

    PROPS order:
        1. E            Young modulus
        2. nu           Poisson ratio (kappa, mu from linear-elastic relations)
        3. c_tau_u      momentum fine-scale constant (paper: 2.0)
        4. c_tau_theta  theta fine-scale constant (paper: 0.1)
        5. h_elem       element size h used in tau_u (uniform mesh)
    """

    props = dict(
        E=250.0,
        nu=0.25,
        c_tau_u=2.0,
        c_tau_theta=0.1,
        h_elem=1.0,
    )

    # P = F S(E(u, theta)) is deliberately NOT the derivative of an energy in
    # F at fixed theta (paper Eq. 23 approximates F(u+u', th+th') ~ F(u)), so
    # dP/dF is unsymmetric by construction; the theta equation carries the
    # missing coupling.  Skip the hyperelastic major-symmetry check.
    symmetric_tangent = False

    def _S_enriched(self, F, thetat):
        """Second PK stress (Eq. 46) at the theta-enriched strain."""
        I = eye(3)
        mu = self.E / (2.0 * (1.0 + self.nu))
        kappa = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        theta = 1.0 + thetat
        J = det(F)
        C = F.T @ F
        Cbar = (theta ** (2.0 / 3.0)) * (J ** (-2.0 / 3.0)) * C
        Cbar_inv = inv(Cbar)
        trCbar = Cbar[0, 0] + Cbar[1, 1] + Cbar[2, 2]
        S = (0.5 * kappa * (theta ** 2 - 1.0) * Cbar_inv
             + mu * theta ** (-2.0 / 3.0) * (I - (trCbar / 3.0) * Cbar_inv))
        return S, Cbar_inv, trCbar, theta, J, mu, kappa

    def stress_PK1(self, F, thetat):
        """P = F [S - tau_theta (theta - J) dS/dtheta]  (Eq. 41a)."""
        I = eye(3)
        S, Cbar_inv, trCbar, theta, J, mu, kappa = self._S_enriched(F, thetat)
        # Dbar = dS/dE : Cbar (Eq. 49 at the enriched state)
        Dbar = (kappa * (2.0 * theta ** 2 + 1.0) * Cbar_inv
                - 2.0 * mu * theta ** (-2.0 / 3.0)
                * (I - (trCbar / 3.0) * Cbar_inv))
        dS_dtheta = Dbar / (3.0 * theta)
        tau_theta = self.c_tau_theta * mu / (mu + kappa)
        return F @ (S - tau_theta * (theta - J) * dS_dtheta)

    def phase_storage(self, F, thetat):
        """(1 - tau_theta)(theta - J), tested with q (Eq. 41b)."""
        mu = self.E / (2.0 * (1.0 + self.nu))
        kappa = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        theta = 1.0 + thetat
        J = det(F)
        tau_theta = self.c_tau_theta * mu / (mu + kappa)
        return (1.0 - tau_theta) * (theta - J)

    def phase_flux(self, F, thetat, grad_thetat):
        """VMS gradient penalty on theta (Eqs. 41b, 43).

        The residual carries + (grad q, G) with
            G = (tau_u/3) (J/theta) Dbar grad(theta),
        and the package convention assembles storage*N - flux.Grad(N),
        so this adapter returns -G (same pattern as the other phase UELs).
        grad(thetat) = grad(theta) since theta = 1 + thetat.
        """
        I = eye(3)
        mu = self.E / (2.0 * (1.0 + self.nu))
        kappa = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        theta = 1.0 + thetat
        J = det(F)
        C = F.T @ F
        Cbar = (theta ** (2.0 / 3.0)) * (J ** (-2.0 / 3.0)) * C
        Cbar_inv = inv(Cbar)
        trCbar = Cbar[0, 0] + Cbar[1, 1] + Cbar[2, 2]
        Dbar = (kappa * (2.0 * theta ** 2 + 1.0) * Cbar_inv
                - 2.0 * mu * theta ** (-2.0 / 3.0)
                * (I - (trCbar / 3.0) * Cbar_inv))
        tau_u = self.c_tau_u * self.h_elem ** 2 / (2.0 * mu)
        G = (tau_u / 3.0) * (J / theta) * (Dbar @ grad_thetat)
        return -1.0 * G


class ScovazziUThetaNearIncMaterial(ScovazziUThetaMaterial):
    """Nearly incompressible parameter set (paper Sec. 6.1.2 / 6.2)."""

    props = dict(
        E=250.0,
        nu=0.49995,
        c_tau_u=2.0,
        c_tau_theta=0.1,
        h_elem=1.0,
    )


class PlainSimoTaylorMaterial(au.Material):
    """Irreducible (displacement-only) Simo-Taylor baseline, Eq. 46 at C(u).

    Used for the locking comparison ("u" curves in the paper's figures).
    """

    props = dict(E=250.0, nu=0.25)

    def stress_PK1(self, F):
        I = eye(3)
        mu = self.E / (2.0 * (1.0 + self.nu))
        kappa = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        J = det(F)
        C = F.T @ F
        C_inv = inv(C)
        trC = C[0, 0] + C[1, 1] + C[2, 2]
        S = (0.5 * kappa * (J ** 2 - 1.0) * C_inv
             + mu * J ** (-2.0 / 3.0) * (I - (trC / 3.0) * C_inv))
        return F @ S


class IrreducibleSimoTaylorUEL(au.WeakForm):
    material = PlainSimoTaylorMaterial
    ndim = 3

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)

    def momentum_equation(self, v, F):
        return self.material.stress_PK1(F)


class ScovazziUThetaUEL(au.WeakForm):
    """Equal-order u + thetat mixed element with VMS stabilization."""

    material = ScovazziUThetaMaterial
    ndim = 3

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)
        self.thetat = au.ScalarField("thetat", degree=1, test="q")

    def momentum_equation(self, v, F, thetat):
        return self.material.stress_PK1(F, thetat)

    def phase_equation(self, q, F, thetat, grad_thetat):
        return (
            self.material.phase_storage(F, thetat),
            self.material.phase_flux(F, thetat, grad_thetat),
        )


class ScovazziUThetaNearIncUEL(ScovazziUThetaUEL):
    material = ScovazziUThetaNearIncMaterial


def verification_states():
    """States that exercise the stabilization terms, not one benign point.

    thetat is deliberately NOT equal to J - 1 so that the (theta - J)
    stabilization terms and their tangent blocks are nonzero, and
    grad(thetat) is nonzero so the flux block is exercised.
    """
    F1 = np.array([                       # moderate distortion, J > 1
        [1.05, 0.02, 0.01],
        [0.015, 0.98, 0.005],
        [0.01, 0.0, 1.03],
    ])
    F2 = np.array([                       # compressed, J < 1
        [0.94, -0.03, 0.0],
        [0.02, 0.97, 0.01],
        [0.0, 0.015, 0.96],
    ])
    F3 = np.eye(3)                        # reference state, thetat != 0
    states = []
    for F, tt in ((F1, 0.04), (F2, -0.05), (F3, 0.02)):
        states.append({
            "F": F,
            "thetat": tt,
            "grad_thetat": np.array([0.03, -0.02, 0.01]),
            "dt": 1.0,
        })
    return states


def self_check_dS_dtheta(tol=1e-9):
    """Analytic dS/dtheta = Dbar/(3 theta) vs complex-step of S (numpy)."""
    mat = ScovazziUThetaMaterial()
    F = verification_states()[0]["F"].astype(complex)
    thetat = 0.04
    h = 1.0e-25
    S_p, *_ = mat._S_enriched(F, thetat + 1j * h)
    dS_cs = np.asarray(S_p).imag / h

    I = np.eye(3)
    S, Cbar_inv, trCbar, theta, J, mu, kappa = mat._S_enriched(F, thetat)
    Dbar = (kappa * (2.0 * theta ** 2 + 1.0) * Cbar_inv
            - 2.0 * mu * theta ** (-2.0 / 3.0)
            * (I - (trCbar / 3.0) * Cbar_inv))
    dS_an = np.asarray(Dbar / (3.0 * theta), dtype=complex).real
    err = np.max(np.abs(dS_an - dS_cs)) / max(np.max(np.abs(dS_cs)), 1e-30)
    assert err < tol, f"dS/dtheta analytic vs CS mismatch: rel {err:.2e}"
    return err


def self_check_homogeneous(tol=1e-12):
    """At thetat = J - 1 the stabilized P reduces to plain Simo-Taylor."""
    mat = ScovazziUThetaMaterial()
    F = verification_states()[0]["F"]
    J = np.linalg.det(F)
    C = F.T @ F
    Cinv = np.linalg.inv(C)
    mu = mat.E / (2.0 * (1.0 + mat.nu))
    kappa = mat.E / (3.0 * (1.0 - 2.0 * mat.nu))
    S_direct = (0.5 * kappa * (J ** 2 - 1.0) * Cinv
                + mu * J ** (-2.0 / 3.0)
                * (np.eye(3) - (np.trace(C) / 3.0) * Cinv))
    P_direct = F @ S_direct
    P_model = np.asarray(mat.stress_PK1(F.astype(complex), J - 1.0)).real
    err = np.max(np.abs(P_model - P_direct)) / np.max(np.abs(P_direct))
    assert err < tol, f"homogeneous exactness failed: rel {err:.2e}"
    # storage and flux vanish identically at the constrained state
    st = mat.phase_storage(F.astype(complex), J - 1.0)
    assert abs(complex(st)) < 1e-12
    fl = mat.phase_flux(F.astype(complex), J - 1.0, np.zeros(3))
    assert np.max(np.abs(np.asarray(fl))) < 1e-30
    return err


def generate(output_dir=None, element="hex8", near_incompressible=False):
    if output_dir is None:
        output_dir = HERE
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    element_key = element.lower()
    if element_key == "hex8":
        element_name, ndim = "Hex8", 3
    elif element_key == "quad4":
        element_name, ndim = "Quad4", 2
    elif element_key == "tet4":
        element_name, ndim = "tet4", 3        # the paper's element (P1/P1)
    elif element_key == "tet4r":
        element_name, ndim = "tet4r", 3       # single-point quadrature variant
    else:
        raise ValueError("element must be 'hex8', 'quad4', 'tet4' or 'tet4r'")

    cls = ScovazziUThetaNearIncUEL if near_incompressible else ScovazziUThetaUEL
    tag = "nearinc" if near_incompressible else "comp"
    problem = cls(ndim=ndim)
    problem.summary()

    for i, state in enumerate(verification_states()):
        ok = problem.verify(state=state, tol=5e-5, verbose=(i == 0))
        if not ok:
            raise RuntimeError(f"verify() failed at state {i}")
        print(f"verify() passed at state {i}")

    out = output_dir / f"scovazzi_utheta_{tag}_{element_key}.for"
    generate_uel(problem, str(out), element=element_name)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--element", choices=("hex8", "quad4", "tet4", "tet4r"),
                        default="hex8")
    parser.add_argument("--near-incompressible", action="store_true")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--self-check-only", action="store_true")
    args = parser.parse_args()

    e1 = self_check_dS_dtheta()
    print(f"self-check dS/dtheta analytic vs complex-step: rel err {e1:.2e}")
    e2 = self_check_homogeneous()
    print(f"self-check homogeneous exactness:              rel err {e2:.2e}")
    if not args.self_check_only:
        path = generate(args.output_dir, element=args.element,
                        near_incompressible=args.near_incompressible)
        print(f"Generated {path}")

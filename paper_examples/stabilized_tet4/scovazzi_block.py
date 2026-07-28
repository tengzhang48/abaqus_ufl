"""Mixed u-theta Tet4 UEL for the Scovazzi Figure 14 block benchmark.

The block in Section 6.3 uses the compressible neo-Hookean energy

    psi(C, J) = mu/2 * (tr(C) - 3) - mu*ln(J)
                + lam/2 * ln(J)**2,

not the Simo-Taylor energy used by the other examples in the paper.  This
module applies the same stabilized u-theta formulation to that energy and
generates the native Abaqus UEL used by ``abaqus_block_fig14``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]  # paper_examples/<pkg>/ -> repository root
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, eye, inv, log
from abaqus_ufl.generators.uel_gen import generate_uel

from scovazzi_utheta import ScovazziUThetaUEL, verification_states


class ScovazziBlockMaterial(au.Material):
    """Neo-Hookean law and VMS terms for paper Section 6.3.

    PROPS order:
        1. mu           shear modulus [N/mm2]
        2. lam          first Lame parameter [N/mm2]
        3. c_tau_u      displacement fine-scale constant
        4. c_tau_theta  Jacobian fine-scale constant
        5. h_elem       representative element size [mm]
    """

    props = dict(
        mu=80.194,
        lam=400889.806,
        c_tau_u=2.0,
        c_tau_theta=0.1,
        h_elem=0.0625,
    )

    symmetric_tangent = False

    def _S_enriched(self, F, thetat):
        """Second PK stress evaluated at the theta-enriched strain."""
        I = eye(3)
        theta = 1.0 + thetat
        J = det(F)
        C = F.T @ F
        Cbar = theta ** (2.0 / 3.0) * J ** (-2.0 / 3.0) * C
        Cbar_inv = inv(Cbar)
        a = self.lam * log(theta) - self.mu
        S = self.mu * I + a * Cbar_inv
        return S, Cbar_inv, theta, J

    def _Dbar(self, Cbar_inv, theta):
        """Return (dS/dE):Cbar for the benchmark energy."""
        factor = 3.0 * self.lam + 2.0 * self.mu
        factor = factor - 2.0 * self.lam * log(theta)
        return factor * Cbar_inv

    def stress_PK1(self, F, thetat):
        S, Cbar_inv, theta, J = self._S_enriched(F, thetat)
        Dbar = self._Dbar(Cbar_inv, theta)
        dS_dtheta = Dbar / (3.0 * theta)
        bulk = self.lam + 2.0 * self.mu / 3.0
        tau_theta = self.c_tau_theta * self.mu / (self.mu + bulk)
        return F @ (S - tau_theta * (theta - J) * dS_dtheta)

    def phase_storage(self, F, thetat):
        theta = 1.0 + thetat
        J = det(F)
        bulk = self.lam + 2.0 * self.mu / 3.0
        tau_theta = self.c_tau_theta * self.mu / (self.mu + bulk)
        return (1.0 - tau_theta) * (theta - J)

    def phase_flux(self, F, thetat, grad_thetat):
        S_unused, Cbar_inv, theta, J = self._S_enriched(F, thetat)
        Dbar = self._Dbar(Cbar_inv, theta)
        tau_u = self.c_tau_u * self.h_elem ** 2 / (2.0 * self.mu)
        gradient_term = (tau_u / 3.0) * (J / theta)
        gradient_term = gradient_term * (Dbar @ grad_thetat)
        # abaqus_ufl assembles storage*N - flux.grad(N).
        return -1.0 * gradient_term


class ScovazziBlockUEL(ScovazziUThetaUEL):
    material = ScovazziBlockMaterial


def self_check_dS_dtheta(tol=1.0e-9):
    """Check Dbar/(3 theta) against a complex-step derivative."""
    mat = ScovazziBlockMaterial()
    F = verification_states()[0]["F"].astype(complex)
    thetat = 0.04
    step = 1.0e-25
    S_step, _, _, _ = mat._S_enriched(F, thetat + 1j * step)
    derivative_cs = np.asarray(S_step).imag / step

    _, Cbar_inv, theta, _ = mat._S_enriched(F, thetat)
    derivative = np.asarray(mat._Dbar(Cbar_inv, theta) / (3.0 * theta))
    error = np.max(np.abs(derivative - derivative_cs))
    scale = max(np.max(np.abs(derivative_cs)), 1.0e-30)
    relative_error = error / scale
    assert relative_error < tol, relative_error
    return relative_error


def self_check_homogeneous(tol=1.0e-12):
    """At theta=J, recover the original block constitutive law."""
    mat = ScovazziBlockMaterial()
    F = verification_states()[0]["F"].astype(complex)
    J = det(F)
    C = F.T @ F
    C_inv = inv(C)
    S = mat.mu * eye(3) + (mat.lam * log(J) - mat.mu) * C_inv
    expected = F @ S
    actual = mat.stress_PK1(F, J - 1.0)
    error = np.max(np.abs(np.asarray(actual - expected)))
    scale = max(np.max(np.abs(np.asarray(expected))), 1.0e-30)
    relative_error = float(error / scale)
    assert relative_error < tol, relative_error
    return relative_error


def generate(output_dir=None, element="tet4"):
    output_dir = HERE if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    element_key = element.lower()
    if element_key == "tet4":
        element_name = "tet4"
    elif element_key == "tet4r":
        element_name = "tet4r"
    else:
        raise ValueError("element must be 'tet4' or 'tet4r'")

    problem = ScovazziBlockUEL(ndim=3)
    problem.summary()
    for index, state in enumerate(verification_states()):
        if not problem.verify(state=state, tol=5.0e-5, verbose=(index == 0)):
            raise RuntimeError("verify() failed at state {}".format(index))
        print("verify() passed at state {}".format(index))

    output = output_dir / "scovazzi_block_{}.for".format(element_key)
    generate_uel(problem, str(output), element=element_name)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--element", choices=("tet4", "tet4r"), default="tet4")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--self-check-only", action="store_true")
    args = parser.parse_args()

    print("dS/dtheta check: {:.3e}".format(self_check_dS_dtheta()))
    print("homogeneous check: {:.3e}".format(self_check_homogeneous()))
    if not args.self_check_only:
        print("Generated {}".format(generate(args.output_dir, args.element)))


if __name__ == "__main__":
    main()

"""Coupled thermo-mechanical Quad4 UEL (displacement + temperature).

First UEL exemplar of the pipeline. Two nodal fields, a cross-field
coupling (thermal expansion in the stress), and a backward-Euler
transport balance enter through the declared weak form; the generator
emits a self-contained Abaqus UEL with the coupled tangent blocks
evaluated by complex step.

  Momentum:   P = G (F - F^{-T}) + K ln(J) F^{-T} - K alpha T F^{-T}
  Transport:  storage = rho_cp (T - T_old)/dt,   flux = -k grad_T
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log
from abaqus_ufl.generators.uel_gen import generate_uel


class HeatDiffusionMaterial(au.Material):
    """Neo-Hookean mechanics + Fourier heat conduction.

    Properties (Abaqus *UEL PROPERTY order):
        G       shear modulus
        K       logarithmic volumetric coefficient (multiplies ln J)
        alpha   thermal-pressure coupling coefficient (P_th = -K alpha T F^{-T})
        k       thermal conductivity
        rho_cp  volumetric heat capacity
    """

    props = dict(G=1.0, K=10.0, alpha=1e-3, k=0.5, rho_cp=1.0)

    def stress_PK1(self, F, T):
        J = det(F)
        finv_t = inv(F).T
        P_mech = self.G * (F - finv_t) + self.K * log(J) * finv_t
        P_thermal = -self.K * self.alpha * T * finv_t
        return P_mech + P_thermal

    def solvent_storage(self, F, F_old, T, T_old, dt):
        return self.rho_cp * (T - T_old) / dt

    def solvent_flux(self, F, T, grad_T):
        return -self.k * grad_T


class HeatDiffusionProblem(au.WeakForm):
    """2D Quad4 problem with fields u (vector) and T (scalar)."""

    material = HeatDiffusionMaterial
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)
        self.T = au.ScalarField("T", degree=1, test="theta")

    def momentum_equation(self, v, F, T):
        return self.material.stress_PK1(F, T)

    def transport_equation(self, theta, F, T, grad_T, F_old, T_old, dt):
        storage = self.material.solvent_storage(F, F_old, T, T_old, dt)
        flux = self.material.solvent_flux(F, T, grad_T)
        return storage, flux


def verification_state():
    """State for tangent verification with nonzero T and grad_T."""
    F = np.array([
        [1.05, 0.02, 0.0],
        [0.01, 1.03, 0.0],
        [0.0, 0.0, 1.0],
    ])
    return dict(
        F=F,
        F_old=0.98 * F,
        T=10.0,
        T_old=5.0,
        grad_T=np.array([2.0, -1.0, 0.5]),
        dt=0.1,
    )


def generate(output=None):
    """Verify the declared weak form and write the generated UEL."""
    problem = HeatDiffusionProblem()
    if not problem.verify(state=verification_state(), verbose=True):
        raise RuntimeError("WeakForm verification failed")

    if output is None:
        output = Path(__file__).with_name("scalar_diffusion_uel.for")
    output = Path(output)
    generate_uel(problem, str(output), element="Quad4", formulation="standard")
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

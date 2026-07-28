"""Mixed-order thermo-mechanical Quad8 UEL (quadratic u, bilinear T).

Mixed-order exemplar of the pipeline. Displacement is quadratic
serendipity on all eight nodes while the temperature-like scalar is
bilinear on the four corners, so the generated element carries
node-dependent degree-of-freedom maps (corners (u1, u2, T), midsides
(u1, u2); NDOFEL = 20). The scalar flux is pulled back with C^{-1}, so
the transport block is genuinely deformation dependent.

  Momentum:   P = G (F - F^{-T}) + K ln(J) F^{-T} - K alpha T F^{-T}
  Transport:  storage = cT (T - T_old)/dt,  flux = -kappa C^{-1} grad_T
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


class ThermoMechanicalMaterial(au.Material):
    """Neo-Hookean mechanics with thermal pressure and pulled-back flux.

    Properties (Abaqus *UEL PROPERTY order):
        G      shear modulus
        K      logarithmic volumetric coefficient (multiplies ln J)
        alpha  thermal-pressure coupling coefficient (P_th = -K alpha T F^{-T})
        kappa  referential conductivity
        cT     volumetric heat capacity
    """

    props = dict(G=1.0, K=100.0, alpha=1.0e-3, kappa=0.25, cT=1.0)

    def stress_PK1(self, F, T):
        finv_t = inv(F).T
        J = det(F)
        P_mech = self.G * (F - finv_t) + self.K * log(J) * finv_t
        P_thermal = -self.K * self.alpha * T * finv_t
        return P_mech + P_thermal

    def solvent_flux(self, F, T, grad_T):
        Cinv = inv(F.T @ F)
        return -self.kappa * (Cinv @ grad_T)

    def solvent_storage(self, F, F_old, T, T_old, dt):
        return self.cT * (T - T_old) / dt


class ThermoMechanicalQuad8(au.WeakForm):
    """Coupled displacement-temperature problem with mixed-order fields."""

    material = ThermoMechanicalMaterial
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=2)
        self.T = au.ScalarField("T", degree=1, test="theta")

    def momentum_equation(self, v, F, T):
        return self.material.stress_PK1(F, T)

    def transport_equation(self, theta, F, T, grad_T, F_old, T_old, dt):
        storage = self.material.solvent_storage(F, F_old, T, T_old, dt)
        flux = self.material.solvent_flux(F, T, grad_T)
        return storage, flux


def verification_state():
    """State with nonzero temperature, gradient, and deformation."""
    F = np.array([
        [1.08, 0.04, 0.0],
        [0.02, 1.05, 0.0],
        [0.0, 0.0, 1.0],
    ])
    return dict(
        F=F,
        F_old=0.98 * F,
        T=15.0,
        T_old=10.0,
        grad_T=np.array([2.0, -1.0, 0.5]),
        dt=0.1,
    )


def generate(output=None):
    """Verify the declared weak form and write the generated UEL."""
    problem = ThermoMechanicalQuad8()
    if not problem.verify(state=verification_state(), verbose=True):
        raise RuntimeError("WeakForm verification failed")

    if output is None:
        output = Path(__file__).with_name("thermo_mechanics_quad8_uel.for")
    output = Path(output)
    generate_uel(problem, str(output), element="Quad8", formulation="standard")
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

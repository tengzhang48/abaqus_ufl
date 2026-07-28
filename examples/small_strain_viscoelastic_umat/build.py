"""Small-strain Standard Linear Solid (SLS) UMAT.

History exemplar of the pipeline. One Maxwell branch in parallel with an
equilibrium spring; the backward-Euler viscous update is unconditionally
stable in dt/tau, and the deviatoric viscous strain is a full 3x3 TENSOR
state variable threaded through ``STATEV(1..9)`` in column-major order.
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import abaqus_ufl as au
from abaqus_ufl.core.tensor import eye, trace


class SmallStrainViscoelastic(au.SmallStrainMaterial):
    """Compression-positive small-strain SLS (one Maxwell branch).

    Properties (Abaqus *User Material constant order):
        K      equilibrium bulk modulus
        G_inf  equilibrium (long-term) shear modulus
        G_v    Maxwell-branch shear modulus
        tau    Maxwell-branch relaxation time

    State: ``eps_v`` is the deviatoric viscous strain tensor of the
    Maxwell dashpot, stored in ``STATEV(1..9)`` column-major.
    """

    props = dict(K=30.0, G_inf=10.0, G_v=5.0, tau=0.5)
    state_vars = dict(eps_v=np.zeros((3, 3)))

    def stress_update(self, sigma_old, strain_old, dstrain,
                      eps_v_old, dt):
        I = eye(3)
        strain_new = strain_old + dstrain
        ratio = dt / self.tau

        dev_e = strain_new - (1.0 / 3.0) * trace(strain_new) * I
        eps_v_new = (eps_v_old + ratio * dev_e) / (1.0 + ratio)

        sigma_new = (
            self.K * trace(strain_new) * I
            + 2.0 * self.G_inf * dev_e
            + 2.0 * self.G_v * (dev_e - eps_v_new))

        return sigma_new, {'eps_v': eps_v_new}


def generate(output=None):
    """Verify the Python tangent and write the generated UMAT."""
    model = SmallStrainViscoelastic()
    if not model.verify(verbose=True):
        raise RuntimeError("Python tangent verification failed")

    if output is None:
        output = Path(__file__).with_name("small_strain_viscoelastic.for")
    output = Path(output)
    au.generate_small_strain_umat(model, str(output))
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

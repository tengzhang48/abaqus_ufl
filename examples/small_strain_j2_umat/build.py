"""Small-strain J2 plasticity UMAT with linear isotropic hardening.

Stateful exemplar of the pipeline. One scalar state variable, the
equivalent plastic strain, is threaded through ``STATEV(1)``, and the
radial-return update has a genuine elastic/plastic branch.
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import abaqus_ufl as au
from abaqus_ufl.core.small_strain_plasticity import flow_direction, q_mises
from abaqus_ufl.core.tensor import eye, trace


class SmallStrainJ2(au.SmallStrainMaterial):
    """Compression-positive small-strain J2 with linear isotropic hardening.

    The order of ``props`` is the Abaqus *User Material constant order.
    ``STATEV(1)`` holds the equivalent plastic strain.
    """

    props = dict(G=10.0, lam=20.0, sigma_y=0.1, H=5.0)
    state_vars = dict(ep=0.0)

    def stress_update(self, sigma_old, strain_old, dstrain, ep_old, dt):
        sigma_trial = (
            sigma_old + 2.0 * self.G * dstrain
            + self.lam * trace(dstrain) * eye(3))
        seq = q_mises(sigma_trial)
        f = seq - (self.sigma_y + self.H * ep_old)

        if f.real > 0.0:
            dgamma = f / (3.0 * self.G + self.H)
            n = flow_direction(sigma_trial, seq)
            sigma_new = sigma_trial - 2.0 * self.G * dgamma * n
            ep_new = ep_old + dgamma
        else:
            sigma_new = sigma_trial
            ep_new = ep_old

        return sigma_new, {'ep': ep_new}


def generate(output=None):
    """Verify the Python tangent and write the generated UMAT."""
    model = SmallStrainJ2()
    if not model.verify(verbose=True):
        raise RuntimeError("Python tangent verification failed")

    if output is None:
        output = Path(__file__).with_name("small_strain_j2.for")
    output = Path(output)
    au.generate_small_strain_umat(model, str(output))
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

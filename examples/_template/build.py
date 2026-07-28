"""Working template: verify and generate a small-strain elastic UMAT.

The material is deliberately simple.  The reusable part of this directory is
the verification pipeline around it; copied examples must replace both the
model and the model-specific checks.
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import abaqus_ufl as au
from abaqus_ufl.core.tensor import eye, trace


class TemplateElastic(au.SmallStrainMaterial):
    """Isotropic elasticity in the compression-positive material API."""

    # The order is the Abaqus *User Material constant order.
    props = dict(G=10.0, lam=20.0)

    def stress_update(self, sigma_old, strain_old, dstrain):
        I = eye(3)
        return (
            sigma_old
            + 2.0 * self.G * dstrain
            + self.lam * trace(dstrain) * I
        )


def generate(output=None):
    """Verify the Python tangent and write the generated UMAT."""
    model = TemplateElastic()
    if not model.verify(verbose=True):
        raise RuntimeError("Python tangent verification failed")

    if output is None:
        output = Path(__file__).with_name("template_umat.for")
    output = Path(output)
    au.generate_small_strain_umat(model, str(output))
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

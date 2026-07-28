"""Compressible neo-Hookean finite-strain UMAT example.

Stateless finite-strain exemplar of the pipeline. The generated UMAT
receives ``DFGRD1`` and returns Cauchy stress with the Jaumann-rate
tangent that Abaqus/Standard solid elements expect.
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log


class NeoHookean(au.Material):
    """P = G (F - F^{-T}) + K ln(J) F^{-T}.

    The order of ``props`` is the Abaqus *User Material constant order.
    """

    props = dict(G=0.5, K=50.0)

    def stress_PK1(self, F):
        J = det(F)
        FinvT = inv(F).T
        return self.G * (F - FinvT) + self.K * log(J) * FinvT


def generate(output=None):
    """Verify the Python tangent and write the generated UMAT."""
    model = NeoHookean()
    if not model.verify(verbose=True):
        raise RuntimeError("Python tangent verification failed")

    if output is None:
        output = Path(__file__).with_name("neo_hookean_umat.for")
    output = Path(output)
    au.generate_umat(model, str(output))
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

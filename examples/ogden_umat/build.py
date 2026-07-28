"""One-term compressible Ogden UMAT via the explicit spectral path.

Spectral exemplar of the pipeline. The strain energy is expressed in
principal stretches, so the material calls ``eig`` directly and
reconstructs the stress with the eigenspace-invariant form
``V diag(f(lam)) inv(V)``. This is the constitutive route that requires
the scale- and rotation-safe eig fallbacks (repeated principal stretches
occur already in uniaxial stretch).
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, eig, inv, log, sqrt


class OgdenOneTerm(au.Material):
    """W = (2 mu / alpha^2) (lb1^a + lb2^a + lb3^a - 3) + (K/2) ln(J)^2.

    ``lb_i = J^{-1/3} lambda_i`` are isochoric principal stretches. The
    principal Kirchhoff stresses are

        tau_i = (2 mu / alpha) (lb_i^a - (1/3) sum_j lb_j^a) + K ln(J)

    and the second Piola-Kirchhoff stress is the eigenspace-invariant
    reconstruction S = V diag(tau_i / lam_i) inv(V) on eig(C), with
    lam_i = lambda_i^2. Properties are in Abaqus *User Material order.
    """

    props = dict(mu=1.0, alpha=3.5, K=100.0)

    def stress_PK1(self, F):
        J = det(F)
        C = F.T @ F
        lam, V = eig(C)

        scale = J ** (-1.0 / 3.0)
        lb1 = scale * sqrt(lam[0])
        lb2 = scale * sqrt(lam[1])
        lb3 = scale * sqrt(lam[2])

        w1 = lb1 ** self.alpha
        w2 = lb2 ** self.alpha
        w3 = lb3 ** self.alpha
        mean_w = (w1 + w2 + w3) / 3.0

        coef = 2.0 * self.mu / self.alpha
        vol = self.K * log(J)
        t1 = coef * (w1 - mean_w) + vol
        t2 = coef * (w2 - mean_w) + vol
        t3 = coef * (w3 - mean_w) + vol

        D = np.array([
            [t1 / lam[0], 0.0, 0.0],
            [0.0, t2 / lam[1], 0.0],
            [0.0, 0.0, t3 / lam[2]],
        ])
        S = V @ D @ inv(V)
        return F @ S


def generate(output=None):
    """Verify the Python tangent and write the generated UMAT."""
    model = OgdenOneTerm()
    if not model.verify(verbose=True):
        raise RuntimeError("Python tangent verification failed")

    if output is None:
        output = Path(__file__).with_name("ogden_umat.for")
    output = Path(output)
    au.generate_umat(model, str(output))
    return output


if __name__ == "__main__":
    print("Generated {}".format(generate()))

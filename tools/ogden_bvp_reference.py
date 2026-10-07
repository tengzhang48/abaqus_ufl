"""Independent principal-stretch and rectangular-Q1 Ogden references.

These functions do not import the material declaration, generated UEL, or the
package's assembly/shape functions. NumPy's real SVD evaluates the energy's
principal stretches directly. The generated model instead diagonalizes C.
All stresses are first Piola stresses per reference area, tension positive.
"""

import math

import numpy as np


def piola_and_energy(F, *, mu, alpha, K):
    """Differentiate W = 2mu/alpha² Σ(lbᵢ^alpha−1) + K/2 log(J)²."""
    F = np.asarray(F, dtype=float)
    if F.shape[-2:] != (3, 3) or not np.isfinite(F).all():
        raise AssertionError("independent Ogden reference requires finite, positive J")
    J = np.linalg.det(F)
    if np.any(J <= 0):
        raise AssertionError("independent Ogden reference requires finite, positive J")
    left, stretches, right_transpose = np.linalg.svd(F)
    powers = (stretches / J[..., None] ** (1 / 3)) ** alpha
    tau = ((2 * mu / alpha) * (powers - np.mean(powers, axis=-1, keepdims=True))
           + K * np.log(J)[..., None])
    P = (left * (tau / stretches)[..., None, :]) @ right_transpose
    W = (2 * mu / alpha ** 2) * (np.sum(powers, axis=-1) - 3) + K / 2 * np.log(J) ** 2
    return P, float(W) if W.ndim == 0 else W, stretches


def shear_traction(gamma, *, mu, alpha):
    """Closed-form P₁₂ for F = I + gamma e₁⊗e₂, J=1.

    lambda± = (sqrt(gamma²+4) ± gamma)/2, lambda₃=1. Differentiating
    W(gamma) gives 2mu/alpha (lambda+^alpha−lambda−^alpha)/sqrt(gamma²+4).
    """
    gamma = np.asarray(gamma, dtype=float)
    diagonal = np.sqrt(gamma ** 2 + 4)
    plus, minus = (diagonal + gamma) / 2, (diagonal - gamma) / 2
    return 2 * mu / alpha * (plus ** alpha - minus ** alpha) / diagonal


def alpha_two_piola(F, *, mu, K):
    """Eig-free α=2 limit: P=mu J⁻²ᐟ³(F−tr(C)/3 F⁻ᵀ)+K ln(J)F⁻ᵀ."""
    J = float(np.linalg.det(F))
    if not np.isfinite(F).all() or J <= 0:
        raise AssertionError("alpha-two reference requires positive J")
    inverse_transpose = np.linalg.inv(F).T
    return (mu * J ** (-2 / 3) * (F - np.sum(F * F) / 3 * inverse_transpose)
            + K * math.log(J) * inverse_transpose)


def alpha_two_element(coordinates, displacements, *, mu, K):
    """Integrate the α=2 eig-free limit at the same real deformed GP states."""
    displacement = np.asarray(displacements).reshape(4, 2)
    residual = np.zeros((4, 2))
    for gradients, weight in quadrature(coordinates):
        F = np.eye(3)
        F[:2, :2] += displacement.T @ gradients
        P = alpha_two_piola(F, mu=mu, K=K)
        residual += gradients @ P[:2, :2].T * weight
    return residual.ravel()


def quadrature(coordinates):
    """Independent tensor-product, four-point integration on a rectangle."""
    coordinates = np.asarray(coordinates)
    dx = float(coordinates[1, 0] - coordinates[0, 0])
    dy = float(coordinates[3, 1] - coordinates[0, 1])
    expected = coordinates[0] + np.array([[0, 0], [dx, 0], [dx, dy], [0, dy]])
    if dx <= 0 or dy <= 0 or not np.allclose(coordinates, expected, rtol=0, atol=1e-12):
        raise AssertionError("independent quadrature is scoped to rectangular Q1 meshes")
    for xi in (-1 / math.sqrt(3), 1 / math.sqrt(3)):
        for eta in (-1 / math.sqrt(3), 1 / math.sqrt(3)):
            gradients = np.array([
                [-(1 - eta), 1 - eta, 1 + eta, -(1 + eta)],
                [-(1 - xi), -(1 + xi), 1 + xi, 1 - xi],
            ]).T / np.array([2 * dx, 2 * dy])
            yield gradients, dx * dy / 4


def element_reference(coordinates, displacements, *, mu, alpha, K):
    """Integrate P Grad(N), W and stretches without generated assembly."""
    displacements = np.asarray(displacements, dtype=float).reshape(4, 2)
    residual, energy = np.zeros((4, 2)), 0.0
    samples = []
    for gradients, weight in quadrature(coordinates):
        F = np.eye(3)
        F[:2, :2] += displacements.T @ gradients
        P, W, stretches = piola_and_energy(F, mu=mu, alpha=alpha, K=K)
        residual += gradients @ P[:2, :2].T * weight
        energy += W * weight
        samples.append((float(np.linalg.det(F)), float(np.max(stretches)), float(P[0, 1])))
    return residual.ravel(), energy, np.asarray(samples)


def mesh_reference(nodes, elements, solution, *, mu, alpha, K):
    """Assemble the independent nodal residual and retain every GP identity.

    The sample order is element-major, then xi-major/eta-minor Gauss order.
    There are exactly four [J, maximum stretch, P12] samples per element.
    """
    residual = np.zeros((len(nodes), 2))
    displacement = np.asarray(solution).reshape(len(nodes), 2)
    integration = [list(quadrature(coordinates)) for coordinates in nodes[elements]]
    gradients = np.asarray([[point[0] for point in points] for points in integration])
    weights = np.asarray([[point[1] for point in points] for points in integration])
    F = np.broadcast_to(np.eye(3), (len(elements), 4, 3, 3)).copy()
    F[:, :, :2, :2] += np.einsum("eai,eqaJ->eqiJ", displacement[elements], gradients)
    P, W, stretches = piola_and_energy(F, mu=mu, alpha=alpha, K=K)
    local_residual = np.einsum("eqaJ,eqiJ,eq->eai", gradients, P[:, :, :2, :2], weights)
    np.add.at(residual, elements, local_residual)
    samples = np.stack([np.linalg.det(F), np.max(stretches, axis=-1), P[:, :, 0, 1]], axis=-1)
    return residual.ravel(), float(np.sum(W * weights)), samples

"""Small-strain plasticity tensor helpers.

All helpers use the package's small-strain UMAT convention:
compression-positive stress and strain tensors. They are intentionally generic
and small enough to map directly through the Fortran generator.

Full material models should live in examples or higher-level libraries, not in
this core helper module.
"""

from .tensor import eye, sqrt, trace


def mean_stress(sigma):
    """Mean stress p = tr(sigma) / 3, positive in compression."""
    return trace(sigma) / 3.0


def p_eff(sigma):
    """Alias for compression-positive mean effective stress."""
    return mean_stress(sigma)


def dev_stress(sigma):
    """Deviatoric stress s = sigma - p I."""
    return sigma - mean_stress(sigma) * eye(3, dtype=sigma.dtype)


def q_mises(sigma):
    """Mises q = sqrt(3/2 s:s)."""
    s = dev_stress(sigma)
    return sqrt(1.5 * trace(s.T @ s))


def flow_direction(sigma, q, tol=1.0e-30):
    """J2/Drucker-Prager deviatoric direction n = 3s/(2q).

    The small ``tol`` avoids division by zero near hydrostatic states while
    preserving complex-step perturbations away from the singularity.
    """
    return (1.5 / (q + tol)) * dev_stress(sigma)


def bulk_modulus(G, lam):
    """Small-strain isotropic bulk modulus K = lambda + 2G/3."""
    return lam + 2.0 * G / 3.0

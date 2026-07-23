"""Backward-compatible aliases for small-strain plasticity helpers.

New code should import from ``abaqus_ufl.core.small_strain_plasticity``. This
module remains so existing soil examples and user scripts do not break.
"""

from .small_strain_plasticity import (  # noqa: F401
    bulk_modulus,
    dev_stress,
    flow_direction,
    mean_stress,
    p_eff,
    q_mises,
)

__all__ = [
    "bulk_modulus",
    "dev_stress",
    "flow_direction",
    "mean_stress",
    "p_eff",
    "q_mises",
]

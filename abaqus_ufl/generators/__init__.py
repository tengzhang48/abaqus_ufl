"""Code generators: Python Material -> Fortran UMAT/UEL and input scaffolds."""

from .inp_scaffold import (
    UELModelConfig,
    ScaffoldReport,
    generate as generate_inp_scaffold,
    write_job_inp,
)
from .umat_gen import generate_small_strain_umat


def generate_uel(*args, **kwargs):
    """Generate an Abaqus UEL; imported lazily to avoid circular imports."""
    from .uel_gen import generate_uel as _generate_uel
    return _generate_uel(*args, **kwargs)


def generate_umat(*args, **kwargs):
    """Generate an Abaqus UMAT; imported lazily to avoid circular imports."""
    from .umat_gen import generate_umat as _generate_umat
    return _generate_umat(*args, **kwargs)


def generate_uinter(*args, **kwargs):
    """Generate an Abaqus UINTER; imported lazily to avoid circular imports."""
    from .uinter_gen import generate_uinter as _generate_uinter
    return _generate_uinter(*args, **kwargs)


__all__ = [
    "UELModelConfig",
    "ScaffoldReport",
    "generate_inp_scaffold",
    "write_job_inp",
    "generate_small_strain_umat",
    "generate_uel",
    "generate_umat",
    "generate_uinter",
]

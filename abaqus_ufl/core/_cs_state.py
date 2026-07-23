"""Complex-step state preparation utilities."""
import numpy as np


def complexify_state(state):
    """
    Convert state dict values to complex for complex-step evaluation.

    All numerical values except ``dt`` and ``time`` are converted to complex dtype
    (with zero imaginary part). Converting everything to complex prevents
    NumPy from repeatedly upcasting during mixed arithmetic, matching the
    generated Fortran where differentiable arguments are DOUBLE COMPLEX.

    Arrays are always ``.copy()``'d to break aliasing.

    ``dt`` and ``time`` stay real (DOUBLE PRECISION in Fortran, never
    differentiated).
    Non-numeric values are passed through unchanged.
    """
    out = {}
    for key, val in state.items():
        if key in {'dt', 'time'}:
            out[key] = val
        elif isinstance(val, np.ndarray):
            out[key] = np.asarray(val, dtype=complex).copy()
        elif isinstance(val, (int, float)):
            out[key] = complex(val, 0.0)
        else:
            out[key] = val
    return out

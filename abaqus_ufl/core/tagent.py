"""
Complex-step tangent engine.

Given a Material instance and a state (F, scalars, vectors), compute
all tangent blocks by perturbing each input argument by ih and extracting
Im(output) / h.

CS_H = 1e-10: verified safe for Fortran COMPLEX*16, identical to 1e-30
across all 12 tangent blocks of the gel model.

The tangent engine inspects each Material method's signature to determine:
  - Which arguments are DOF variables (need tangent w.r.t.)
  - What shape each argument has (3x3 for F, scalar for p, 3-vector for grad_mu)
  - What shape the output has (3x3 for P, scalar for r_p, 3-vector for j_R)

This information drives the perturbation loops and output array shapes.
"""

import numpy as np

from ._cs_state import complexify_state as _complexify_state

CS_H = 1.0e-10


def _perturb_scalar(x):
    """Return complex copy of scalar with imaginary perturbation."""
    return complex(x, CS_H)


def _perturb_matrix_component(A, k, l):
    """Return complex copy of matrix with A[k,l] perturbed."""
    Ac = A.astype(complex)
    Ac[k, l] += 1j * CS_H
    return Ac


def _perturb_vector_component(v, k):
    """Return complex copy of vector with v[k] perturbed."""
    vc = v.astype(complex)
    vc[k] += 1j * CS_H
    return vc


def _extract_tangent_scalar(output_c):
    """Extract tangent from a complex scalar output."""
    return np.imag(output_c) / CS_H


def _extract_tangent_array(output_c):
    """Extract tangent from a complex array output."""
    return np.imag(output_c) / CS_H


# =====================================================================
# Generic tangent: d(output) / d(input) for any material function
# =====================================================================

def tangent_wrt_matrix(func, wrt_name, state, output_shape):
    """
    Compute d(func)/d(matrix) where matrix is 3x3.

    func: callable(**state) -> array of output_shape
    wrt_name: key in state dict for the 3x3 matrix to differentiate w.r.t.
    state: dict of all arguments {name: value}
    output_shape: shape of func's return value

    Returns: array of shape (*output_shape, 3, 3)
    """
    result = np.zeros(output_shape + (3, 3))
    for k in range(3):
        for l in range(3):
            s = _complexify_state(state)
            s[wrt_name] = _perturb_matrix_component(state[wrt_name], k, l)
            out_c = func(**s)
            if np.isscalar(out_c) or (hasattr(out_c, 'shape') and out_c.shape == ()):
                result[..., k, l] = _extract_tangent_scalar(out_c)
            else:
                result[..., k, l] = _extract_tangent_array(out_c)
    return result


def tangent_wrt_scalar(func, wrt_name, state, output_shape):
    """
    Compute d(func)/d(scalar).

    Returns: array of output_shape (or scalar if output_shape is ())
    """
    s = _complexify_state(state)
    s[wrt_name] = complex(float(state[wrt_name]), CS_H)
    out_c = func(**s)
    if np.isscalar(out_c) or (hasattr(out_c, 'shape') and out_c.shape == ()):
        return _extract_tangent_scalar(out_c)
    else:
        return _extract_tangent_array(out_c)


def tangent_wrt_vector(func, wrt_name, state, output_shape):
    """
    Compute d(func)/d(vector) where vector is n-dimensional.

    Returns: array of shape (*output_shape, n)
    """
    vec = state[wrt_name]
    n = vec.shape[0]
    result = np.zeros(output_shape + (n,))
    for k in range(n):
        s = _complexify_state(state)
        s[wrt_name] = _perturb_vector_component(state[wrt_name], k)
        out_c = func(**s)
        if np.isscalar(out_c) or (hasattr(out_c, 'shape') and out_c.shape == ()):
            result[..., k] = _extract_tangent_scalar(out_c)
        else:
            result[..., k] = _extract_tangent_array(out_c)
    return result


# =====================================================================
# Finite-difference tangent (for cross-verification)
# =====================================================================

FD_EPS = 1e-7


def fd_tangent_wrt_matrix(func, wrt_name, state, output_shape, eps=FD_EPS):
    """FD verification of d(func)/d(matrix)."""
    result = np.zeros(output_shape + (3, 3))
    for k in range(3):
        for l in range(3):
            sp = dict(state)
            sm = dict(state)
            Ap = state[wrt_name].copy()
            Am = state[wrt_name].copy()
            Ap[k, l] += eps
            Am[k, l] -= eps
            sp[wrt_name] = Ap
            sm[wrt_name] = Am
            fp = func(**sp)
            fm = func(**sm)
            result[..., k, l] = (np.asarray(fp) - np.asarray(fm)) / (2 * eps)
    return result


def fd_tangent_wrt_scalar(func, wrt_name, state, output_shape, eps=FD_EPS):
    """FD verification of d(func)/d(scalar)."""
    sp = dict(state)
    sm = dict(state)
    sp[wrt_name] = float(state[wrt_name]) + eps
    sm[wrt_name] = float(state[wrt_name]) - eps
    fp = func(**sp)
    fm = func(**sm)
    return (np.asarray(fp) - np.asarray(fm)) / (2 * eps)


def fd_tangent_wrt_vector(func, wrt_name, state, output_shape, eps=FD_EPS):
    """FD verification of d(func)/d(vector)."""
    vec = state[wrt_name]
    n = vec.shape[0]
    result = np.zeros(output_shape + (n,))
    for k in range(n):
        sp = dict(state)
        sm = dict(state)
        vp = vec.copy(); vp[k] += eps
        vm = vec.copy(); vm[k] -= eps
        sp[wrt_name] = vp
        sm[wrt_name] = vm
        fp = func(**sp)
        fm = func(**sm)
        result[..., k] = (np.asarray(fp) - np.asarray(fm)) / (2 * eps)
    return result

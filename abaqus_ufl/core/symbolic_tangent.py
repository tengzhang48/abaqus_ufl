"""
Symbolic tangent engine with complex-step verification.

Current scope: whole material methods are analyzed with SymPy. If a method is
closed-form and differentiable symbolically, this module can generate pure-real
Fortran tangent expressions and verify them against the existing complex-step
oracle. If symbolic analysis fails, the production generator should fall back
to the existing full complex-step path.

True subexpression-level "hybrid" differentiation, where only an eigensolve,
return map, or other hard node is differentiated by local complex-step and
then chained analytically, is future research work and is not implemented here.
"""

import importlib
import inspect
import numpy as np
import sympy as sp

from . import tensor
from .material import FIELD_ARGS


def _wrap_fortran(line, maxcol=72):
    """Wrap a Fortran fixed-format line at column maxcol."""
    if not line or len(line) <= maxcol:
        return line
    if line and line[0] in ('C', 'c', '!'):
        return line
    if line.strip() == '':
        return line
    result_lines = []
    remaining = line
    while len(remaining) > maxcol:
        breakpt = -1
        for pos in range(maxcol - 1, 6, -1):
            ch = remaining[pos]
            if ch in (' ', ',', '+', '-', '*', '(', ')'):
                breakpt = pos + 1
                break
        if breakpt <= 6:
            breakpt = maxcol
        result_lines.append(remaining[:breakpt])
        remaining = '     & ' + remaining[breakpt:].lstrip()
    if remaining:
        result_lines.append(remaining)
    return '\n'.join(result_lines)


# =====================================================================
# SymPy-compatible tensor function replacements
# =====================================================================

class _SymPyTensor:
    """Drop-in replacements for tensor.py functions that work with SymPy."""

    @staticmethod
    def det(A):
        if isinstance(A, sp.MatrixBase):
            return A.det()
        return _ORIGINAL_TENSOR_FUNCS['det'](A)

    @staticmethod
    def inv(A):
        if isinstance(A, sp.MatrixBase):
            return A.inv()
        return _ORIGINAL_TENSOR_FUNCS['inv'](A)

    @staticmethod
    def log(x):
        if isinstance(x, (sp.Expr, sp.MatrixBase)):
            return sp.log(x)
        return _ORIGINAL_TENSOR_FUNCS['log'](x)

    @staticmethod
    def exp(x):
        if isinstance(x, (sp.Expr, sp.MatrixBase)):
            return sp.exp(x)
        return _ORIGINAL_TENSOR_FUNCS['exp'](x)

    @staticmethod
    def sqrt(x):
        if isinstance(x, (sp.Expr, sp.MatrixBase)):
            return sp.sqrt(x)
        return _ORIGINAL_TENSOR_FUNCS['sqrt'](x)

    @staticmethod
    def trace(A):
        if isinstance(A, sp.MatrixBase):
            return sp.trace(A)
        return _ORIGINAL_TENSOR_FUNCS['trace'](A)

    @staticmethod
    def eye(n, dtype=float):
        if dtype is sp.Expr or dtype == sp.Expr:
            return sp.eye(n)
        return sp.eye(n)

    @staticmethod
    def sym(A):
        if isinstance(A, sp.MatrixBase):
            return (A + A.T) / 2
        return _ORIGINAL_TENSOR_FUNCS['sym'](A)

    @staticmethod
    def dev(A):
        if isinstance(A, sp.MatrixBase):
            tr = sp.trace(A)
            return A - (tr / 3.0) * sp.eye(3)
        return _ORIGINAL_TENSOR_FUNCS['dev'](A)

    @staticmethod
    def I1(A):
        return _SymPyTensor.trace(A)

    @staticmethod
    def I2(A):
        if isinstance(A, sp.MatrixBase):
            trA = sp.trace(A)
            trA2 = sp.trace(A * A)
            return sp.Rational(1, 2) * (trA**2 - trA2)
        return _ORIGINAL_TENSOR_FUNCS['I2'](A)

    @staticmethod
    def dyad(a, b):
        if isinstance(a, sp.MatrixBase) and isinstance(b, sp.MatrixBase):
            return a * b.T
        return _ORIGINAL_TENSOR_FUNCS['dyad'](a, b)

    @staticmethod
    def cross(a, b):
        if isinstance(a, sp.MatrixBase) and isinstance(b, sp.MatrixBase):
            c = sp.zeros(3, 1)
            c[0] = a[1] * b[2] - a[2] * b[1]
            c[1] = a[2] * b[0] - a[0] * b[2]
            c[2] = a[0] * b[1] - a[1] * b[0]
            return c
        return _ORIGINAL_TENSOR_FUNCS['cross'](a, b)


# Names to patch in material module and tensor module
_PATCH_NAMES = [
    'det', 'inv', 'log', 'exp', 'sqrt', 'trace', 'eye',
    'sym', 'dev', 'I1', 'I2', 'dyad', 'cross',
]

# Save original tensor functions for fallback in _SymPyTensor
_ORIGINAL_TENSOR_FUNCS = {}
for _name in _PATCH_NAMES:
    if hasattr(tensor, _name):
        _ORIGINAL_TENSOR_FUNCS[_name] = getattr(tensor, _name)


def _patch_module(mod, sympy_funcs):
    """Patch a module's globals with SymPy-compatible functions."""
    originals = {}
    for name in _PATCH_NAMES:
        if hasattr(mod, name):
            originals[name] = getattr(mod, name)
            setattr(mod, name, getattr(sympy_funcs, name))
    return originals


def _restore_module(mod, originals):
    """Restore a module's patched globals."""
    for name, val in originals.items():
        setattr(mod, name, val)


def _patch_tensor_module(sympy_funcs):
    """Patch abaqus_ufl.core.tensor module."""
    return _patch_module(tensor, sympy_funcs)


# =====================================================================
# SymbolicTangent class
# =====================================================================

class SymbolicTangent:
    """Analyze material functions and generate symbolic tangent code."""

    def __init__(self, material):
        self.material = material
        self.mode = None          # 'symbolic', 'cs', or 'partial'
        self.sym_results = {}     # method_name -> SymPy expression
        self.hard_nodes = []
        self.field_args = getattr(material, '_field_args', FIELD_ARGS)

        # Build symbolic property variables (PROPS_1, PROPS_2, ...)
        # and index mapping for Fortran emission.
        self.prop_symbols = {}
        self.prop_index = {}
        for idx, (name, default) in enumerate(
                material.__class__.props.items(), 1):
            self.prop_symbols[name] = sp.Symbol(f'PROPS_{idx}', real=True)
            self.prop_index[name] = idx

    # -----------------------------------------------------------------
    # Analysis: run with SymPy, see where it breaks
    # -----------------------------------------------------------------

    def analyze(self):
        """Classify each material method as symbolic, cs, or partial."""
        sympy_funcs = _SymPyTensor()

        # Patch tensor module (affects all imports from tensor)
        tensor_originals = _patch_tensor_module(sympy_funcs)

        # Also patch the module where the material class is defined
        # (handles 'from tensor import det' style imports)
        mat_module = importlib.import_module(
            self.material.__class__.__module__)
        mat_originals = _patch_module(mat_module, sympy_funcs)

        # Temporarily replace material properties with symbolic variables
        original_props = {}
        for name, sym in self.prop_symbols.items():
            original_props[name] = getattr(self.material, name)
            setattr(self.material, name, sym)

        try:
            for method_name in self.material.defined_methods():
                method = getattr(self.material, method_name)
                sig = inspect.signature(method)
                params = list(sig.parameters.keys())

                # Build symbolic state for this method
                sym_state = self._build_symbolic_state(params)

                try:
                    result = method(**sym_state)
                    # Check result is SymPy-compatible
                    if isinstance(result, tuple):
                        sym_result = result[0]
                    else:
                        sym_result = result

                    if isinstance(sym_result, (sp.MatrixBase, sp.Expr)):
                        self.sym_results[method_name] = sym_result
                    else:
                        # Returned a numeric type — not symbolic
                        self.mode = 'cs'
                        self.hard_nodes.append(
                            f"{method_name}: returned non-symbolic type")
                except (TypeError, AttributeError, NotImplementedError) as e:
                    self.mode = 'cs'
                    self.hard_nodes.append(f"{method_name}: {e}")
        finally:
            _restore_module(tensor, tensor_originals)
            _restore_module(mat_module, mat_originals)
            # Restore original property values
            for name, val in original_props.items():
                setattr(self.material, name, val)

        # Determine overall mode
        if not self.hard_nodes:
            self.mode = 'symbolic'
        elif self.sym_results:
            self.mode = 'partial'
        else:
            self.mode = 'cs'

        return self.mode

    def _build_symbolic_state(self, params):
        """Build a dict of symbolic arguments for a material method."""
        state = {}
        for p in params:
            if p == 'F':
                state[p] = sp.Matrix(3, 3,
                    lambda i, j: sp.Symbol(f'F_{i+1}{j+1}', real=True))
            elif p in ('p', 'mu'):
                state[p] = sp.Symbol(p, real=True)
            elif p == 'grad_mu':
                state[p] = sp.Matrix(3, 1,
                    lambda i, j: sp.Symbol(f'grad_mu_{i+1}', real=True))
            elif p == 'F_old':
                state[p] = sp.Matrix(3, 3,
                    lambda i, j: sp.Symbol(f'F_old_{i+1}{j+1}', real=True))
            elif p == 'p_old':
                state[p] = sp.Symbol('p_old', real=True)
            elif p == 'dt':
                state[p] = sp.Symbol('dt', real=True)
            elif p in self.field_args:
                shape = self.field_args[p]['shape']
                if shape == ():
                    state[p] = sp.Symbol(p, real=True)
                elif shape == (3,):
                    state[p] = sp.Matrix(3, 1,
                        lambda i, j: sp.Symbol(f'{p}_{i+1}', real=True))
                elif shape == (3, 3):
                    state[p] = sp.Matrix(3, 3,
                        lambda i, j: sp.Symbol(f'{p}_{i+1}{j+1}', real=True))
            elif p in self.material._history_args and p.endswith('_old'):
                base = p[:-4]
                if base not in self.field_args:
                    raise NotImplementedError(
                        f"Symbolic tangent cannot infer the shape of "
                        f"history argument '{p}'")
                shape = self.field_args[base]['shape']
                if shape == ():
                    state[p] = sp.Symbol(p, real=True)
                elif shape == (3,):
                    state[p] = sp.Matrix(3, 1,
                        lambda i, j: sp.Symbol(f'{p}_{i+1}', real=True))
                elif shape == (3, 3):
                    state[p] = sp.Matrix(3, 3,
                        lambda i, j: sp.Symbol(f'{p}_{i+1}{j+1}', real=True))
                else:
                    raise NotImplementedError(
                        f"Unsupported symbolic history shape for '{p}': "
                        f"{shape}")
            elif p in self.material._state_arg_names:
                base = p[:-4]
                init = self.material.state_vars[base]
                if np.isscalar(init):
                    state[p] = sp.Symbol(p, real=True)
                else:
                    arr = np.asarray(init)
                    if arr.shape == (3,):
                        state[p] = sp.Matrix(3, 1,
                            lambda i, j: sp.Symbol(f'{p}_{i+1}', real=True))
                    elif arr.shape == (3, 3):
                        state[p] = sp.Matrix(3, 3,
                            lambda i, j: sp.Symbol(
                                f'{p}_{i+1}{j+1}', real=True))
                    else:
                        raise NotImplementedError(
                            f"Unsupported symbolic state shape for '{p}': "
                            f"{arr.shape}")
            else:
                raise NotImplementedError(
                    f"Symbolic tangent cannot represent argument '{p}'. "
                    "Add it to field/history/state metadata or use the CS "
                    "tangent path.")
        return state

    # -----------------------------------------------------------------
    # Symbolic differentiation
    # -----------------------------------------------------------------

    def symbolic_tangent(self, method_name='stress_PK1', wrt='F'):
        """
        Compute symbolic tangent d(method)/d(wrt).

        Returns a dict mapping (output_idx, input_idx) -> SymPy expression.
        For matrix outputs/inputs, indices are (i,j) and (k,l).
        For vector outputs/inputs, indices are (i,) and (k,).
        """
        if method_name not in self.sym_results:
            raise ValueError(
                f"Method {method_name} is not symbolic. "
                f"Run analyze() first.")

        P_sym = self.sym_results[method_name]
        sig = inspect.signature(getattr(self.material, method_name))
        params = list(sig.parameters.keys())
        sym_state = self._build_symbolic_state(params)

        if wrt not in sym_state:
            raise ValueError(f"Method {method_name} does not take '{wrt}'")

        wrt_sym = sym_state[wrt]

        # Compute shapes
        if isinstance(P_sym, sp.MatrixBase):
            out_shape = P_sym.shape
        else:
            out_shape = ()

        if isinstance(wrt_sym, sp.MatrixBase):
            in_shape = wrt_sym.shape
        else:
            in_shape = ()

        # Build flat dict of all derivatives
        result = {}
        out_indices = list(np.ndindex(out_shape)) if out_shape else [None]
        in_indices = list(np.ndindex(in_shape)) if in_shape else [None]

        for oidx in out_indices:
            for iidx in in_indices:
                if oidx is None and iidx is None:
                    expr = sp.diff(P_sym, wrt_sym)
                elif oidx is None:
                    expr = sp.diff(P_sym, wrt_sym[iidx])
                elif iidx is None:
                    expr = sp.diff(P_sym[oidx], wrt_sym)
                else:
                    expr = sp.diff(P_sym[oidx], wrt_sym[iidx])
                result[(oidx, iidx)] = expr

        return result, out_shape, in_shape

    def symbolic_tangent_cse(self, method_name='stress_PK1', wrt='F'):
        """Differentiate and apply common subexpression elimination."""
        dPdF, raw_out_shape, raw_in_shape = self.symbolic_tangent(method_name, wrt)

        # Normalize shapes: squeeze trailing singleton dimensions
        out_shape = raw_out_shape
        if out_shape and out_shape[-1] == 1:
            out_shape = out_shape[:-1]
        in_shape = raw_in_shape
        if in_shape and in_shape[-1] == 1:
            in_shape = in_shape[:-1]

        # Flatten all expressions for joint CSE
        all_exprs = []
        for key in sorted(dPdF.keys()):
            all_exprs.append(dPdF[key])

        replacements, simplified = sp.cse(all_exprs)
        return replacements, simplified, dPdF, out_shape, in_shape

    # -----------------------------------------------------------------
    # Fortran code generation
    # -----------------------------------------------------------------

    def generate_fortran_tangent(self, method_name='stress_PK1', wrt='F',
                                  result_name='result'):
        """Generate Fortran code string for the tangent."""
        import re
        from sympy.printing import fcode

        replacements, simplified, dPdF, out_shape, in_shape = \
            self.symbolic_tangent_cse(method_name, wrt)

        # Shape assertion: number of components must match out_shape+in_shape
        expected_n = int(np.prod(out_shape + in_shape)) \
            if (out_shape + in_shape) else 1
        actual_n = len(dPdF)
        block_name = f"{method_name}/{wrt}"
        assert actual_n == expected_n, (
            f"Block {block_name}: expected {expected_n} components, "
            f"got {actual_n}. out_shape={out_shape}, in_shape={in_shape}. "
            f"Check squeeze logic in symbolic_tangent_cse().")

        parts = []
        parts.append("C     Generated by abaqus_ufl symbolic tangent")
        parts.append("C     Pure DOUBLE PRECISION - no complex arithmetic")
        parts.append("")

        # Declare CSE temporaries
        if replacements:
            cse_names = [str(var) for var, _ in replacements]
            decl = "      DOUBLE PRECISION :: " + ", ".join(cse_names)
            parts.append(_wrap_fortran(decl))
            parts.append("")

        # Emit temporary variables
        for var, expr in replacements:
            parts.append(fcode(expr, assign_to=str(var),
                               source_format='fixed'))

        if replacements:
            parts.append("")

        # Emit result components using generic shape dispatch
        idx = 0
        for key in sorted(dPdF.keys()):
            oidx, iidx = key
            lhs = self._build_fortran_lhs(
                result_name, out_shape, in_shape, oidx, iidx)
            parts.append(fcode(simplified[idx], assign_to=lhs,
                               source_format='fixed'))
            idx += 1

        code = '\n'.join(parts)

        # Replace symbolic property variables with Fortran PROPS(i) references.
        # MUST run before the indexed-symbol regex below, otherwise PROPS_11
        # gets matched as a 2D index and becomes PROPS(1,1).
        # Sort descending by index so PROPS_11 is replaced before PROPS_1.
        for name, idx in sorted(self.prop_index.items(), key=lambda x: x[1],
                                reverse=True):
            code = code.replace(f'PROPS_{idx}', f'PROPS({idx})')

        # Convert SymPy indexed symbols to Fortran array syntax.
        # e.g. F_11 -> F(1,1), grad_mu_1 -> grad_mu(1)
        # Apply 2D replacement first to avoid partial matches.
        code = re.sub(r'([A-Za-z][A-Za-z0-9_]*)_([1-9])([1-9])',
                      r'\1(\2,\3)', code)
        code = re.sub(r'([A-Za-z][A-Za-z0-9_]*)_([1-9])',
                      r'\1(\2)', code)

        # Normalize continuation characters: SymPy fcode uses '@' in col 6,
        # but the project style uses '&'.
        code = code.replace('\n     @', '\n     &')

        # Re-wrap lines that may have grown beyond 72 chars after regex.
        lines = code.split('\n')
        wrapped = []
        for line in lines:
            wrapped.append(_wrap_fortran(line))
        code = '\n'.join(wrapped)

        return code

    @staticmethod
    def _build_fortran_lhs(result_name, out_shape, in_shape, oidx, iidx):
        """Build Fortran lhs array subscript from shapes and indices.

        Examples:
            out=(3,3), in=(3,3), oidx=(0,1), iidx=(2,0)
                -> result(1,2,3,1)
            out=(3,), in=(3,3), oidx=(0,), iidx=(1,2)
                -> result(1,2,3)
            out=(), in=(), oidx=None, iidx=None
                -> result
        """
        indices = []
        if out_shape:
            for dim, idx in zip(out_shape, oidx):
                indices.append(str(idx + 1))
        if in_shape:
            for dim, idx in zip(in_shape, iidx):
                indices.append(str(idx + 1))
        if indices:
            return f"{result_name}({','.join(indices)})"
        return result_name

    # -----------------------------------------------------------------
    # Numerical verification
    # -----------------------------------------------------------------

    def verify(self, n_tests=5, tol=1e-10):
        """Cross-check symbolic tangent against full complex-step."""
        from .tagent import (
            tangent_wrt_matrix, tangent_wrt_scalar, tangent_wrt_vector
        )

        if self.mode != 'symbolic':
            raise RuntimeError(
                f"verify() only works in symbolic mode, got {self.mode}")

        max_err_P = 0.0
        max_err_dPdF = 0.0
        block_errs = {}

        for _ in range(n_tests):
            # Random deformation gradient
            F_test = _random_deformation_gradient()
            state = self.material._make_state(F=F_test)

            for method_name in self.material.defined_methods():
                func = self.material._call_method_as_func(method_name)
                info = self.material._methods[method_name]

                # Evaluate symbolic expression numerically
                P_sym = self._eval_symbolic(method_name, state)

                # Evaluate via normal (real) function
                P_real = func(**state)

                # Check stress match
                err_P = np.max(np.abs(np.asarray(P_sym) - np.asarray(P_real)))
                max_err_P = max(max_err_P, err_P)

                # Check tangent blocks
                for field in info['field_params']:
                    kind = self.field_args[field]['kind']
                    out_shape = info['output_shape']
                    block_key = f"{method_name}/{field}"

                    if kind == 'matrix':
                        dPdF_cs = tangent_wrt_matrix(
                            func, field, state, out_shape)
                    elif kind == 'scalar':
                        dPdF_cs = tangent_wrt_scalar(
                            func, field, state, out_shape)
                    elif kind == 'vector':
                        dPdF_cs = tangent_wrt_vector(
                            func, field, state, out_shape)
                    else:
                        continue

                    dPdF_sym = self._eval_symbolic_tangent(
                        method_name, field, state)

                    err = np.max(np.abs(np.asarray(dPdF_sym) - np.asarray(dPdF_cs)))
                    max_err_dPdF = max(max_err_dPdF, err)
                    if block_key not in block_errs:
                        block_errs[block_key] = 0.0
                    block_errs[block_key] = max(block_errs[block_key], err)

        if max_err_P > tol:
            raise AssertionError(
                f"Symbolic residual error {max_err_P:.3e} exceeds tol "
                f"{tol:.3e}")

        if max_err_dPdF > tol:
            raise AssertionError(
                f"Symbolic tangent error {max_err_dPdF:.3e} exceeds tol {tol:.3e}")

        print(f"Verified symbolic tangent at {n_tests} states, "
              f"max P err={max_err_P:.3e}, max dPdF err={max_err_dPdF:.3e}")
        for key in sorted(block_errs.keys()):
            print(f"  {key}: max err={block_errs[key]:.3e}")
        return max_err_P, max_err_dPdF, block_errs

    def _eval_symbolic(self, method_name, state):
        """Evaluate a symbolic result numerically at a given state."""
        sym_result = self.sym_results[method_name]
        sig = inspect.signature(getattr(self.material, method_name))
        params = list(sig.parameters.keys())
        sym_state = self._build_symbolic_state(params)

        # Build substitution dict
        subs = {}
        for p in params:
            if p in state:
                val = state[p]
                sym_val = sym_state[p]
                if isinstance(sym_val, sp.MatrixBase):
                    if val.ndim == 1:
                        for i in range(sym_val.shape[0]):
                            for j in range(sym_val.shape[1]):
                                subs[sym_val[i, j]] = float(val[i])
                    else:
                        for i in range(sym_val.shape[0]):
                            for j in range(sym_val.shape[1]):
                                subs[sym_val[i, j]] = float(val[i, j])
                elif isinstance(sym_val, (sp.Expr,)):
                    subs[sym_val] = float(val)

        # Also substitute props
        for name, sym_prop in self.prop_symbols.items():
            subs[sym_prop] = float(getattr(self.material, name))

        num_result = sym_result.subs(subs)
        if isinstance(num_result, sp.MatrixBase):
            arr = np.array(num_result.tolist(), dtype=float)
            if arr.ndim == 2 and arr.shape[1] == 1:
                return arr.flatten()
            return arr
        return float(num_result)

    def _eval_symbolic_tangent(self, method_name, wrt, state):
        """Evaluate symbolic tangent numerically at a given state.

        Returns a numpy array of shape output_shape + input_shape,
        with column vectors normalized to 1D.
        """
        P_sym = self.sym_results[method_name]
        sig = inspect.signature(getattr(self.material, method_name))
        params = list(sig.parameters.keys())
        sym_state = self._build_symbolic_state(params)
        wrt_sym = sym_state[wrt]

        # Build substitution dict
        subs = {}
        for p in params:
            if p in state:
                val = state[p]
                sym_val = sym_state[p]
                if isinstance(sym_val, sp.MatrixBase):
                    if val.ndim == 1:
                        for i in range(sym_val.shape[0]):
                            for j in range(sym_val.shape[1]):
                                subs[sym_val[i, j]] = float(val[i])
                    else:
                        for i in range(sym_val.shape[0]):
                            for j in range(sym_val.shape[1]):
                                subs[sym_val[i, j]] = float(val[i, j])
                elif isinstance(sym_val, (sp.Expr,)):
                    subs[sym_val] = float(val)

        # Also substitute props
        for name, sym_prop in self.prop_symbols.items():
            subs[sym_prop] = float(getattr(self.material, name))

        # Determine shapes (use raw SymPy shapes for indexing)
        if isinstance(P_sym, sp.MatrixBase):
            raw_out_shape = P_sym.shape
        else:
            raw_out_shape = ()

        if isinstance(wrt_sym, sp.MatrixBase):
            raw_in_shape = wrt_sym.shape
        else:
            raw_in_shape = ()

        # Compute all derivatives into a raw-shaped array
        raw_result = np.zeros(raw_out_shape + raw_in_shape)
        out_indices = list(np.ndindex(raw_out_shape)) if raw_out_shape else [None]
        in_indices = list(np.ndindex(raw_in_shape)) if raw_in_shape else [None]

        for oidx in out_indices:
            for iidx in in_indices:
                if oidx is None and iidx is None:
                    expr = sp.diff(P_sym, wrt_sym)
                elif oidx is None:
                    expr = sp.diff(P_sym, wrt_sym[iidx])
                elif iidx is None:
                    expr = sp.diff(P_sym[oidx], wrt_sym)
                else:
                    expr = sp.diff(P_sym[oidx], wrt_sym[iidx])

                val = float(expr.subs(subs))

                if oidx is None and iidx is None:
                    raw_result = np.array(val)
                elif oidx is None:
                    raw_result[iidx] = val
                elif iidx is None:
                    raw_result[oidx] = val
                else:
                    raw_result[oidx + iidx] = val

        # Squeeze trailing singleton dimensions to match numpy conventions
        if isinstance(raw_result, np.ndarray):
            return np.squeeze(raw_result)
        return np.array(raw_result)


# =====================================================================
# Helpers
# =====================================================================

def _random_deformation_gradient():
    """Generate a random but well-conditioned deformation gradient."""
    # Random rotation + stretch
    theta = np.random.uniform(-np.pi, np.pi)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    # Ensure stretches are not too small (avoid phi > 1 in gel model)
    lam = np.diag([np.random.uniform(0.85, 1.3),
                   np.random.uniform(0.85, 1.3),
                   np.random.uniform(0.85, 1.3)])
    F = R @ lam
    F = F + np.random.normal(0, 0.03, (3, 3))
    # Ensure positive determinant
    if np.linalg.det(F) < 0.5:
        F = F @ np.diag([1.1, 1.1, 1.1])
    return F

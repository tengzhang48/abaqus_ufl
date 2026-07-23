"""
UMAT code generator.

Given a Material instance with stress_PK1(self, F), generates a complete
self-contained .for file that works as an Abaqus UMAT subroutine.

The generated file contains:
  1. tensor_ops.for   — det33, inv33, matmul33, transpose33 (real + complex)
  2. stress subroutine — COMPLEX*16, calls tensor_ops
  3. CS tangent engine — computes dP/dF via complex-step perturbation
  4. PK1-to-Cauchy push-forward and Jaumann correction
  5. UMAT wrapper      — standard Abaqus UMAT signature

Internal variables (state_vars):
  Materials may define a state_vars dict mapping variable names to initial
  values (scalar float, 3-vector, or 3x3 tensor). When present, the
  generator adds STATEV read/write to the UMAT wrapper (following
  the standard Abaqus STATEV TIME(2)==0 initialization pattern), extends the
  material subroutine
  signature with state_old/state_new arguments, and passes state_old
  (real, not perturbed) through the CS tangent engine.

  The STATEV layout is determined by iteration order of state_vars.
  Python 3.7+ guarantees dict insertion order, so the layout is stable.
  For explicit control, use collections.OrderedDict.

  dt is passed as DOUBLE PRECISION (never differentiated).

  Small-strain materials may also return a special state key named
  ``"_cutback"`` or ``"request_cutback"``. It is not stored in STATEV; the
  generated UMAT maps values greater than 0.5 to ``PNEWDT = 0.25`` and returns
  before accepting the trial stress/state.

Usage:
    # Stateless material:
    model = NeoHookean(G=0.5, K=50.0)
    model.verify()
    generate_umat(model, 'neo_hookean.for')

    # Material with internal variables:
    model = J2Plasticity(E=200e3, nu=0.3, sigma_y=250, H=1000)
    model.verify()
    generate_umat(model, 'j2_umat.for')
"""

import inspect
import textwrap
import ast
import os
import numpy as np
from collections import OrderedDict

try:
    from ..core.verify import VerificationError
except ImportError:
    VerificationError = AssertionError


# ---------------------------------------------------------------------------
# Field argument names recognized by the framework (from defs.py)
# ---------------------------------------------------------------------------
try:
    from ..core.defs import FIELD_ARGS, PARAM_ARGS
except ImportError:
    # Standalone usage — define the minimal set needed for UMAT
    FIELD_ARGS = {'F', 'p', 'mu', 'grad_mu', 'theta'}
    PARAM_ARGS = {'dt'}


# ===================================================================
# State variable layout helpers
# ===================================================================

def _state_var_info(material):
    """
    Compute state variable layout from material.state_vars.

    Returns OrderedDict:  name -> {'init': value, 'size': int, 'shape': str}
      shape is 'scalar', 'vector', or 'tensor'

    Returns empty OrderedDict if material has no state_vars.
    """
    if not hasattr(material, 'state_vars') or not material.state_vars:
        return OrderedDict()

    info = OrderedDict()
    state_var_props = getattr(material, 'state_var_props', {})
    prop_names = list(getattr(material.__class__, 'props', {}).keys())
    for name, init_val in material.state_vars.items():
        arr = np.asarray(init_val)
        init_prop_index = None
        if name in state_var_props:
            prop_name = state_var_props[name]
            if prop_name not in prop_names:
                raise ValueError(
                    f"state_var_props maps state variable '{name}' to "
                    f"unknown material property '{prop_name}'.")
            init_prop_index = prop_names.index(prop_name) + 1
        if arr.ndim == 0:
            info[name] = {
                'init': float(arr),
                'size': 1,
                'shape': 'scalar',
                'init_prop_index': init_prop_index,
            }
        elif arr.ndim == 1 and arr.shape == (3,):
            info[name] = {
                'init': arr,
                'size': 3,
                'shape': 'vector',
                'init_prop_index': init_prop_index,
            }
        elif arr.ndim == 2 and arr.shape == (3, 3):
            info[name] = {
                'init': arr,
                'size': 9,
                'shape': 'tensor',
                'init_prop_index': init_prop_index,
            }
        else:
            raise ValueError(
                f"state_var '{name}' has unsupported shape {arr.shape}. "
                f"Supported: scalar, (3,) vector, (3,3) tensor.")
    return info


def _nstate_per_gp(sv_info):
    """Total number of STATEV entries per integration point."""
    return sum(v['size'] for v in sv_info.values())


# ===================================================================
# AST helpers
# ===================================================================

def _get_source_body(method):
    """Extract the body of a method as source lines (skip def and docstring)."""
    src = inspect.getsource(method)
    tree = ast.parse(textwrap.dedent(src))
    func_def = tree.body[0]
    # Find where the body starts (skip docstring)
    body_start = 0
    if (func_def.body and isinstance(func_def.body[0], ast.Expr)
            and isinstance(func_def.body[0].value, ast.Constant)
            and isinstance(func_def.body[0].value.value, str)):
        body_start = 1
    body_nodes = func_def.body[body_start:]
    return body_nodes, textwrap.dedent(src)


def _reject_early_returns(body_nodes, method_name='<method>'):
    """Raise if any Return appears inside a nested control-flow block.

    The translator only tracks a single `_return_expr`, so a Return nested
    inside If/For/While would be silently overwritten by the next Return
    it visits. Detect and fail loudly instead of dropping code.

    The single trailing Return in the function body (the last statement
    at top level) is the only allowed form.
    """
    for top in body_nodes:
        for child in ast.walk(top):
            if child is top:
                continue
            if isinstance(child, ast.Return):
                raise NotImplementedError(
                    f"Early `return` inside a control-flow block is not "
                    f"supported in {method_name}() because the translator "
                    f"would silently overwrite earlier return expressions. "
                    f"Assign the result to a local variable through "
                    f"if/elif/else and `return` it once at the end of the "
                    f"function body."
                )


# ===================================================================
# Fortran AST translator
# ===================================================================

class FortranTranslator(ast.NodeVisitor):
    """
    Translate a restricted subset of Python AST to Fortran expressions.

    Every _visit_expr returns (fortran_expr, kind) where kind is:
      'scalar' — a single DOUBLE COMPLEX value
      'tensor' — a 3x3 DOUBLE COMPLEX array name

    Supported operations (the tensor vocabulary):
      det(A)      -> det33z(Az)               [scalar]
      inv(A)      -> inv33z call + temp        [tensor]
      log(x)      -> LOG(xz)                  [scalar]
      exp(x)      -> EXP(xz)                  [scalar]
      sqrt(x)     -> SQRT(xz)                 [scalar]
      trace(A)    -> (A(1,1)+A(2,2)+A(3,3))   [scalar]
      A @ B       -> matmul33z call            [tensor]
      A.T         -> transpose33z call         [tensor]
      x.real      -> DBLE(x)                  [scalar] (for CS-safe branching)
      a * b       -> depends on operand types
      a + b, a - b, a / b, a ** n
      self.G      -> DCMPLX(props(i), 0.0d0)  [scalar]
    """

    def __init__(self, props_names, mat_prefix, state_var_names=None,
                 initial_var_kinds=None, material=None,
                 helper_sources=None, helper_cache=None,
                 helper_method_kinds=None, helper_stack=None,
                 matrix_backend='iterative'):
        if matrix_backend not in ('eig', 'iterative'):
            raise ValueError(
                "matrix_backend must be 'eig' or 'iterative'")
        self.props_names = props_names
        self.mat_prefix = mat_prefix
        self.material = material
        self.matrix_backend = matrix_backend
        self._temp_counter = 0
        self._declarations = []  # (type, name, dims)
        self._statements = []    # Fortran statements
        self._return_expr = None
        self._return_kind = None
        self._return_exprs_kinds = None  # multi-output regular helper
        self._var_kinds = dict(initial_var_kinds or {'F': 'tensor'})
        self._integer_vars = set()  # for-loop indices; emitted as Fortran INTEGER
        self._real_vars = set()    # scalar vars assigned from .real, safe in IF comparisons
        self._state_return_map = None  # set by tuple return handler
        self._state_var_names = state_var_names or []
        self._helper_sources = (
            [] if helper_sources is None else helper_sources)
        self._helper_cache = (
            {} if helper_cache is None else helper_cache)
        self._helper_method_kinds = (
            {} if helper_method_kinds is None else helper_method_kinds)
        self._helper_stack = [] if helper_stack is None else helper_stack

    def _new_temp(self, prefix='t'):
        self._temp_counter += 1
        return f'{prefix}{self._temp_counter}'

    def _trace_expr(self, A):
        return f'({A}(1,1)+{A}(2,2)+{A}(3,3))'

    def _dev_tensor(self, A):
        eye = self._new_temp('EYE')
        self._add_decl('DOUBLE COMPLEX', eye, '(3,3)')
        self._var_kinds[eye] = 'tensor'
        self._add_stmt(f'CALL eye33z({eye})')

        dev = self._new_temp('DEV')
        self._add_decl('DOUBLE COMPLEX', dev, '(3,3)')
        self._var_kinds[dev] = 'tensor'
        self._add_stmt(
            f'{dev} = {A} - ({self._trace_expr(A)} / '
            f'DCMPLX(3.000000000000000d+00, 0.0d0)) * {eye}')
        return dev

    def _add_decl(self, ftype, name, dims=None):
        self._declarations.append((ftype, name, dims))

    def _add_stmt(self, stmt):
        self._statements.append(stmt)

    def _dims_for_kind(self, kind):
        if kind == 'tensor':
            return '(3,3)'
        if kind == 'vector':
            return '(3)'
        if isinstance(kind, str) and kind.startswith('vector'):
            return f'({kind[len("vector"):]})'
        if isinstance(kind, str) and kind.startswith('matrix'):
            n = kind[len("matrix"):]
            return f'({n},{n})'
        return None

    def _declare_temp(self, name, kind):
        dims = self._dims_for_kind(kind)
        self._add_decl('DOUBLE COMPLEX', name, dims)
        self._var_kinds[name] = kind

    def _materialize_call_arg(self, expr, kind):
        """Store a helper-call argument in a local temp before CALL."""
        tmp = self._new_temp('ARG')
        self._declare_temp(tmp, kind)
        self._add_stmt(f'{tmp} = {expr}')
        return tmp

    def _is_self_helper_call(self, node):
        return (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == 'self'
            and node.func.attr.startswith('_')
        )

    def _helper_subroutine_name(self, method_name):
        return f'{self.mat_prefix}_{method_name.lstrip("_")}'

    def _ensure_helper_subroutine(self, method_name, arg_kinds):
        """Translate a self._helper(...) method into a Fortran subroutine."""
        if self.material is None:
            raise NotImplementedError(
                f"Helper method dispatch requires a material instance: "
                f"self.{method_name}()")

        if method_name in self._helper_stack:
            chain = ' -> '.join(self._helper_stack + [method_name])
            raise NotImplementedError(
                f"Recursive helper method dispatch is not supported: {chain}")

        previous_kinds = self._helper_method_kinds.get(method_name)
        if previous_kinds is not None and previous_kinds != tuple(arg_kinds):
            raise NotImplementedError(
                f"Helper self.{method_name}() is called with multiple "
                f"argument kind signatures: {previous_kinds} and "
                f"{tuple(arg_kinds)}")
        self._helper_method_kinds[method_name] = tuple(arg_kinds)

        if method_name in self._helper_cache:
            return self._helper_cache[method_name]

        method = getattr(self.material, method_name, None)
        if method is None or not callable(method):
            raise NotImplementedError(
                f"Unknown helper method self.{method_name}()")

        fh_meta = getattr(method, '__fortran_helper__', None)
        if fh_meta is not None:
            # Hand-written Fortran sidecar — skip body translation
            sub_name = fh_meta.get('subroutine_override')
            if sub_name is None:
                sub_name = self._helper_subroutine_name(method_name)
            outputs = fh_meta.get('outputs', [])
            kind_map = {
                'scalar': 'scalar',
                'vector(3)': 'vector',
                'tensor(3,3)': 'tensor',
            }
            return_kinds = [kind_map.get(o[1], 'scalar') for o in outputs]

            # Arity validation: Python return tuple must match outputs count
            if len(outputs) > 1:
                body_nodes, _ = _get_source_body(method)
                for body_node in body_nodes:
                    for node in ast.walk(body_node):
                        if isinstance(node, ast.Return):
                            if not isinstance(node.value, ast.Tuple):
                                raise NotImplementedError(
                                    f"fortran_helper on self.{method_name}() "
                                    f"declares {len(outputs)} outputs but "
                                    f"returns a non-tuple")
                            if len(node.value.elts) != len(outputs):
                                raise NotImplementedError(
                                    f"fortran_helper on self.{method_name}() "
                                    f"declares {len(outputs)} outputs but "
                                    f"returns {len(node.value.elts)} elements")

            helper_info = {
                'sub_name': sub_name,
                'return_kind': return_kinds[0] if return_kinds else 'scalar',
                'return_kinds': return_kinds,
                '_is_fortran_helper': True,
                '_pass_props': fh_meta.get('pass_props', False),
            }
            self._helper_cache[method_name] = helper_info
            return helper_info

        sig = inspect.signature(method)
        params = list(sig.parameters.keys())
        if len(params) != len(arg_kinds):
            raise NotImplementedError(
                f"Helper self.{method_name}() expected {len(params)} "
                f"arguments, got {len(arg_kinds)}")

        sub_name = self._helper_subroutine_name(method_name)
        helper_info = {'sub_name': sub_name, 'return_kind': None}
        self._helper_cache[method_name] = helper_info

        initial_var_kinds = dict(zip(params, arg_kinds))
        translator = FortranTranslator(
            self.props_names,
            self.mat_prefix,
            self._state_var_names,
            initial_var_kinds=initial_var_kinds,
            material=self.material,
            helper_sources=self._helper_sources,
            helper_cache=self._helper_cache,
            helper_method_kinds=self._helper_method_kinds,
            helper_stack=self._helper_stack + [method_name],
            matrix_backend=self.matrix_backend,
        )
        decls, stmts, return_expr, return_kind = \
            translator.translate_stress(method)
        if translator._state_return_map:
            raise NotImplementedError(
                f"Helper self.{method_name}() returned a state dict. "
                "Helper dispatch supports scalar/tensor return values only.")

        multi_returns = translator._return_exprs_kinds
        if multi_returns:
            # Multi-output regular Python helper: `return a, b, c`.
            return_kinds = [k for (_, k) in multi_returns]
            for rk in return_kinds:
                if rk not in ('scalar', 'tensor', 'vector'):
                    raise NotImplementedError(
                        f"Helper self.{method_name}() returns unsupported "
                        f"kind {rk!r} in a tuple")
            src = self._render_helper_subroutine_multi(
                sub_name, params, arg_kinds, decls, stmts,
                multi_returns, translator._var_kinds)
            self._helper_sources.append(src)
            helper_info['return_kind'] = return_kinds[0]
            helper_info['return_kinds'] = return_kinds
            helper_info['_is_multi_return'] = True
            return helper_info

        if return_kind not in ('scalar', 'tensor', 'vector'):
            raise NotImplementedError(
                f"Helper self.{method_name}() returned unsupported kind "
                f"{return_kind!r}")

        src = self._render_helper_subroutine(
            sub_name, params, arg_kinds, decls, stmts,
            return_expr, return_kind, translator._var_kinds)
        self._helper_sources.append(src)
        helper_info['return_kind'] = return_kind
        return helper_info

    def _render_helper_subroutine(self, sub_name, params, arg_kinds,
                                  decls, stmts, return_expr, return_kind,
                                  var_kinds):
        nprops = len(self.props_names)
        result_dims = self._dims_for_kind(return_kind)
        sig_args = list(params) + ['props', 'result']
        lines = [
            f'      SUBROUTINE {sub_name}(' + ', '.join(sig_args) + ')',
            '      IMPLICIT NONE',
        ]

        for name, kind in zip(params, arg_kinds):
            dims = self._dims_for_kind(kind)
            dim_str = dims if dims else ''
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN)  :: {name}{dim_str}')
        lines.append(
            f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')
        result_dim_str = result_dims if result_dims else ''
        lines.append(
            f'      DOUBLE COMPLEX, INTENT(OUT) :: result{result_dim_str}')
        lines.append('')
        lines.append('      DOUBLE COMPLEX :: det33z')
        lines.append('      DOUBLE COMPLEX :: cs_sin')
        lines.append('      DOUBLE COMPLEX :: cs_cos')
        lines.append('      DOUBLE COMPLEX :: cs_tan')
        lines.append('      DOUBLE COMPLEX :: cs_sinh')
        lines.append('      DOUBLE COMPLEX :: cs_cosh')
        lines.append('      DOUBLE COMPLEX :: cs_tanh')
        lines.append('      DOUBLE COMPLEX :: cs_atan')
        lines.append('      DOUBLE COMPLEX :: cs_acos')
        lines.append('      DOUBLE COMPLEX :: cs_erf')
        lines.append('      DOUBLE COMPLEX :: cs_erfc')

        skip_names = set(params) | {'props', 'result'}
        declared_names = set(skip_names)
        for ftype, name, dims in decls:
            if name in declared_names:
                continue
            dim_str = dims if dims else ''
            lines.append(f'      {ftype} :: {name}{dim_str}')
            declared_names.add(name)

        tensor_vars = {n for n, k in var_kinds.items() if k == 'tensor'}
        for stmt in stmts:
            if '=' in stmt and 'CALL' not in stmt:
                stripped = stmt.strip()
                if stripped.startswith(('IF ', 'ELSE', 'END IF',
                                       'DO ', 'END DO', 'EXIT')):
                    continue
                varname = stmt.split('=')[0].strip()
                if '(' in varname:
                    continue
                if varname not in declared_names:
                    if varname not in tensor_vars:
                        lines.append(f'      DOUBLE COMPLEX :: {varname}')
                    declared_names.add(varname)

        lines.append('      INTEGER :: ii, jj')
        declared_names.update({'ii', 'jj'})
        for stmt in stmts:
            stripped = stmt.strip()
            if stripped.startswith('DO ') and '=' in stripped:
                loop_var = stripped[3:].split('=')[0].strip()
                if loop_var not in declared_names:
                    lines.append(f'      INTEGER :: {loop_var}')
                    declared_names.add(loop_var)
        lines.append('')

        for stmt in stmts:
            lines.append(f'      {stmt}')

        if return_kind == 'tensor':
            lines.append('')
            lines.append('      DO ii = 1, 3')
            lines.append('        DO jj = 1, 3')
            if return_expr in tensor_vars or return_expr in declared_names:
                lines.append(
                    f'          result(ii,jj) = {return_expr}(ii,jj)')
            else:
                expr_ij = _tensorize_expr(
                    return_expr,
                    tensor_vars | {
                        n for n, k in zip(params, arg_kinds)
                        if k == 'tensor'
                    })
                lines.append(f'          result(ii,jj) = {expr_ij}')
            lines.append('        END DO')
            lines.append('      END DO')
        elif return_kind == 'vector':
            lines.append('')
            lines.append('      DO ii = 1, 3')
            lines.append(f'        result(ii) = {return_expr}(ii)')
            lines.append('      END DO')
        else:
            lines.append(f'      result = {return_expr}')

        lines.extend([
            '',
            '      RETURN',
            f'      END SUBROUTINE {sub_name}',
        ])
        return '\n'.join(lines)

    def _render_helper_subroutine_multi(self, sub_name, params, arg_kinds,
                                        decls, stmts, multi_returns,
                                        var_kinds):
        """Render a Python helper that returns a tuple `(a, b, c, ...)`.

        Each return element becomes a trailing INTENT(OUT) argument
        named ``out1, out2, ...``. Sized per the inferred return kind.
        """
        nprops = len(self.props_names)
        return_kinds = [k for (_, k) in multi_returns]
        return_exprs = [e for (e, _) in multi_returns]
        out_names = [f'out{i + 1}' for i in range(len(multi_returns))]

        sig_args = list(params) + ['props'] + out_names
        lines = [
            f'      SUBROUTINE {sub_name}(' + ', '.join(sig_args) + ')',
            '      IMPLICIT NONE',
        ]

        for name, kind in zip(params, arg_kinds):
            dims = self._dims_for_kind(kind)
            dim_str = dims if dims else ''
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN)  :: {name}{dim_str}')
        lines.append(
            f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')
        for oname, okind in zip(out_names, return_kinds):
            odim = self._dims_for_kind(okind) or ''
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(OUT) :: {oname}{odim}')
        lines.append('')
        lines.append('      DOUBLE COMPLEX :: det33z')
        lines.append('      DOUBLE COMPLEX :: cs_sin')
        lines.append('      DOUBLE COMPLEX :: cs_cos')
        lines.append('      DOUBLE COMPLEX :: cs_tan')
        lines.append('      DOUBLE COMPLEX :: cs_sinh')
        lines.append('      DOUBLE COMPLEX :: cs_cosh')
        lines.append('      DOUBLE COMPLEX :: cs_tanh')
        lines.append('      DOUBLE COMPLEX :: cs_atan')
        lines.append('      DOUBLE COMPLEX :: cs_acos')
        lines.append('      DOUBLE COMPLEX :: cs_erf')
        lines.append('      DOUBLE COMPLEX :: cs_erfc')

        skip_names = set(params) | {'props'} | set(out_names)
        declared_names = set(skip_names)
        for ftype, name, dims in decls:
            if name in declared_names:
                continue
            dim_str = dims if dims else ''
            lines.append(f'      {ftype} :: {name}{dim_str}')
            declared_names.add(name)

        tensor_vars = {n for n, k in var_kinds.items() if k == 'tensor'}
        for stmt in stmts:
            if '=' in stmt and 'CALL' not in stmt:
                stripped = stmt.strip()
                if stripped.startswith(('IF ', 'ELSE', 'END IF',
                                       'DO ', 'END DO', 'EXIT')):
                    continue
                varname = stmt.split('=')[0].strip()
                if '(' in varname:
                    continue
                if varname not in declared_names:
                    if varname not in tensor_vars:
                        lines.append(f'      DOUBLE COMPLEX :: {varname}')
                    declared_names.add(varname)

        lines.append('      INTEGER :: ii, jj')
        declared_names.update({'ii', 'jj'})
        for stmt in stmts:
            stripped = stmt.strip()
            if stripped.startswith('DO ') and '=' in stripped:
                loop_var = stripped[3:].split('=')[0].strip()
                if loop_var not in declared_names:
                    lines.append(f'      INTEGER :: {loop_var}')
                    declared_names.add(loop_var)
        lines.append('')

        for stmt in stmts:
            lines.append(f'      {stmt}')

        param_tensors = {
            n for n, k in zip(params, arg_kinds) if k == 'tensor'
        }
        lines.append('')
        for oname, okind, oexpr in zip(out_names, return_kinds, return_exprs):
            if okind == 'tensor':
                lines.append('      DO ii = 1, 3')
                lines.append('        DO jj = 1, 3')
                if oexpr in tensor_vars or oexpr in declared_names:
                    lines.append(
                        f'          {oname}(ii,jj) = {oexpr}(ii,jj)')
                else:
                    expr_ij = _tensorize_expr(
                        oexpr, tensor_vars | param_tensors)
                    lines.append(f'          {oname}(ii,jj) = {expr_ij}')
                lines.append('        END DO')
                lines.append('      END DO')
            elif okind == 'vector':
                lines.append('      DO ii = 1, 3')
                lines.append(f'        {oname}(ii) = {oexpr}(ii)')
                lines.append('      END DO')
            else:
                lines.append(f'      {oname} = {oexpr}')

        lines.extend([
            '',
            '      RETURN',
            f'      END SUBROUTINE {sub_name}',
        ])
        return '\n'.join(lines)

    def _parse_helper_call_args(self, node):
        """Extract materialised arg names and kinds from a self._helper(...) call."""
        if node.keywords:
            raise NotImplementedError(
                "Helper method dispatch does not support keyword arguments")
        arg_names = []
        arg_kinds = []
        for arg in node.args:
            expr, kind = self._visit_expr(arg)
            arg_names.append(self._materialize_call_arg(expr, kind))
            arg_kinds.append(kind)
        return arg_names, arg_kinds

    def _emit_helper_call(self, node):
        method_name = node.func.attr
        arg_names, arg_kinds = self._parse_helper_call_args(node)

        helper_info = self._ensure_helper_subroutine(method_name, arg_kinds)
        return_kind = helper_info['return_kind']
        if return_kind is None:
            raise NotImplementedError(
                f"Recursive helper self.{method_name}() has no known return "
                "kind yet")
        # Single-output helper used in expression context
        result = self._new_temp('H')
        self._declare_temp(result, return_kind)
        if helper_info.get('_is_fortran_helper'):
            # Hand-written sidecar — props is optional
            if helper_info.get('_pass_props'):
                call_args = arg_names + ['props', result]
            else:
                call_args = arg_names + [result]
        else:
            call_args = arg_names + ['props', result]
        self._add_stmt(
            f'CALL {helper_info["sub_name"]}('
            f'{", ".join(call_args)})')
        return result, return_kind

    def translate_stress(self, method):
        """
        Translate stress_PK1(self, F, ...) -> Fortran subroutine body.

        Returns: (declarations, body_statements, return_expr, return_kind)
        """
        nodes, src = _get_source_body(method)
        _reject_early_returns(nodes, method_name=getattr(method, '__name__', '?'))
        for node in nodes:
            self._visit_stmt(node)
        return (self._declarations, self._statements,
                self._return_expr, self._return_kind)

    def _visit_stmt(self, node):
        if isinstance(node, ast.Return):
            # Three return forms:
            #   1. `return P, {'ep': ...}` — stress_update form: tensor + state dict
            #   2. `return a, b, c`        — multi-output regular helper
            #   3. `return expr`           — single value
            if isinstance(node.value, ast.Tuple):
                elts = node.value.elts
                is_stress_state_form = (
                    len(elts) >= 2 and isinstance(elts[1], ast.Dict))
                if is_stress_state_form:
                    # Form 1: stress + state dict
                    expr, kind = self._visit_expr(elts[0])
                    self._return_expr = expr
                    self._return_kind = kind
                    state_node = elts[1]
                    self._state_return_map = {}
                    for key_node, val_node in zip(state_node.keys,
                                                   state_node.values):
                        key = key_node.value if isinstance(
                            key_node, ast.Constant) else key_node.s
                        val_expr, val_kind = self._visit_expr(val_node)
                        self._state_return_map[key] = (val_expr, val_kind)
                else:
                    # Form 2: multi-output regular helper return.
                    self._return_exprs_kinds = [
                        self._visit_expr(e) for e in elts]
                    # Keep _return_expr / _return_kind populated for any
                    # downstream code that only inspects the first slot.
                    if self._return_exprs_kinds:
                        self._return_expr = self._return_exprs_kinds[0][0]
                        self._return_kind = self._return_exprs_kinds[0][1]
            else:
                # Form 3: single value
                expr, kind = self._visit_expr(node.value)
                self._return_expr = expr
                self._return_kind = kind
        elif isinstance(node, ast.Assign):
            target = node.targets[0]
            if isinstance(target, ast.Name):
                varname = target.id
                expr, kind = self._visit_expr(node.value)
                self._var_kinds[varname] = kind
                # Track variables assigned from x.real as real-safe scalars
                if (isinstance(node.value, ast.Attribute)
                        and node.value.attr == 'real'):
                    self._real_vars.add(varname)
                dims = self._dims_for_kind(kind)
                if dims:
                    self._add_decl('DOUBLE COMPLEX', varname, dims)
                    self._add_stmt(f'{varname} = {expr}')
                else:
                    self._add_stmt(f'{varname} = {expr}')
            elif isinstance(target, ast.Tuple):
                # Tuple unpacking: lam, V = eig(C) or R, U = polar(F)
                #                 or r1, r2, r3 = cubic_roots(a,b,c,d)
                #                 or a, b = self._multi_out_helper(...)
                names = [elt.id for elt in target.elts
                         if isinstance(elt, ast.Name)]
                if (isinstance(node.value, ast.Call)
                        and self._is_self_helper_call(node.value)):
                    # Multi-output helper call. Two flavours:
                    #   - @fortran_helper(outputs=[...])  -> hand-written sidecar
                    #   - regular Python helper that `return`s a tuple
                    method_name = node.value.func.attr
                    arg_names, arg_kinds = \
                        self._parse_helper_call_args(node.value)
                    helper_info = self._ensure_helper_subroutine(
                        method_name, arg_kinds)
                    return_kinds = helper_info.get('return_kinds')
                    if not return_kinds:
                        raise NotImplementedError(
                            f"Tuple unpacking from self.{method_name}() "
                            f"requires either @fortran_helper(outputs=[...]) "
                            f"or a Python body that `return`s a tuple of "
                            f"tensors/vectors/scalars."
                        )
                    if len(names) != len(return_kinds):
                        raise NotImplementedError(
                            f"self.{method_name}() returns {len(return_kinds)} "
                            f"value(s) but assignment unpacks {len(names)}")
                    for name, kind in zip(names, return_kinds):
                        self._var_kinds[name] = kind
                        if kind == 'tensor':
                            self._add_decl('DOUBLE COMPLEX', name, '(3,3)')
                        elif kind == 'vector':
                            self._add_decl('DOUBLE COMPLEX', name, '(3)')
                        else:
                            self._add_decl('DOUBLE COMPLEX', name, '')
                    # Regular Python multi-return helpers always take props
                    # (same convention as single-return regular helpers).
                    # @fortran_helper helpers take props only if pass_props=True.
                    is_fortran_helper = helper_info.get('_is_fortran_helper')
                    if is_fortran_helper:
                        pass_props = helper_info.get('_pass_props', False)
                    else:
                        pass_props = True
                    if pass_props:
                        call_args = arg_names + ['props'] + names
                    else:
                        call_args = arg_names + names
                    self._add_stmt(
                        f'CALL {helper_info["sub_name"]}('
                        f'{", ".join(call_args)})')
                elif (isinstance(node.value, ast.Call)
                        and len(names) == 2):
                    func_name = self._get_func_name(node.value)
                    args = [self._visit_expr(a)[0]
                            for a in node.value.args]
                    if func_name == 'eig' and len(args) == 1:
                        # lam, V = eig(A)
                        # -> CALL eig33z(A, lam, V)
                        # lam is vector(3), V is tensor(3,3)
                        self._var_kinds[names[0]] = 'vector'
                        self._var_kinds[names[1]] = 'tensor'
                        self._add_decl('DOUBLE COMPLEX',
                                       names[0], '(3)')
                        self._add_decl('DOUBLE COMPLEX',
                                       names[1], '(3,3)')
                        self._add_stmt(
                            f'CALL eig33z({args[0]}, '
                            f'{names[0]}, {names[1]})')
                    elif func_name == 'polar' and len(args) == 1:
                        # R, U = polar(F)
                        # Both are tensor(3,3)
                        self._var_kinds[names[0]] = 'tensor'
                        self._var_kinds[names[1]] = 'tensor'
                        self._add_decl('DOUBLE COMPLEX',
                                       names[0], '(3,3)')
                        self._add_decl('DOUBLE COMPLEX',
                                       names[1], '(3,3)')
                        if self.matrix_backend == 'iterative':
                            self._add_stmt(
                                f'CALL polar33z_iter({args[0]}, '
                                f'{names[0]}, {names[1]})')
                        else:
                            self._add_stmt(
                                f'CALL polar33z({args[0]}, '
                                f'{names[0]}, {names[1]})')
                    else:
                        raise NotImplementedError(
                            f"Unsupported 2-element tuple-return: "
                            f"{func_name}")
                elif (isinstance(node.value, ast.Call)
                        and len(names) == 3):
                    func_name = self._get_func_name(node.value)
                    args = [self._visit_expr(a)[0]
                            for a in node.value.args]
                    if func_name == 'cubic_roots' and len(args) == 4:
                        # r1, r2, r3 = cubic_roots(a, b, c, d)
                        # -> CALL cs_cubic_roots(a, b, c, d, r1, r2, r3)
                        for n in names:
                            self._var_kinds[n] = 'scalar'
                            self._add_decl('DOUBLE COMPLEX', n, '')
                        self._add_stmt(
                            f'CALL cs_cubic_roots({args[0]}, '
                            f'{args[1]}, {args[2]}, {args[3]}, '
                            f'{names[0]}, {names[1]}, {names[2]})')
                    else:
                        raise NotImplementedError(
                            f"Unsupported 3-element tuple-return: "
                            f"{func_name}")
                else:
                    raise NotImplementedError(
                        f"Unsupported tuple assignment: "
                        f"{ast.dump(target)}")
            else:
                raise NotImplementedError(
                    f"Unsupported assignment target: {ast.dump(target)}")
        elif isinstance(node, ast.If):
            # Support if/else for return mapping (CS-safe branching)
            test_expr, _ = self._visit_expr(node.test)
            self._add_stmt(f'IF ({test_expr}) THEN')
            for child in node.body:
                self._visit_stmt(child)
            if node.orelse:
                self._add_stmt('ELSE')
                for child in node.orelse:
                    self._visit_stmt(child)
            self._add_stmt('END IF')
        elif isinstance(node, ast.For):
            # Support Python range() while preserving Python loop values.
            # Subscript translation already maps Python index i to Fortran
            # (i + 1), so loop variables must remain zero-based.
            # The iterator must be range(N) or range(start, stop)
            target = node.target
            if not isinstance(target, ast.Name):
                raise NotImplementedError(
                    "For loop target must be a simple variable")
            loop_var = target.id
            self._var_kinds[loop_var] = 'scalar'
            # Loop indices are emitted as Fortran INTEGER, not DOUBLE
            # COMPLEX. Tracking them here lets the comparison handler
            # accept `if k < 3` without forcing `k.real < 3`.
            self._integer_vars.add(loop_var)

            iter_node = node.iter
            if (isinstance(iter_node, ast.Call)
                    and self._get_func_name(iter_node) == 'range'):
                args = iter_node.args
                if len(args) == 1:
                    # range(N) -> DO loop_var = 0, N-1
                    arg = args[0]
                    if isinstance(arg, ast.Constant):
                        stop = str(int(arg.value) - 1)
                    else:
                        stop_expr, _ = self._visit_expr(arg)
                        stop = f'INT(DBLE({stop_expr})) - 1'
                    self._add_stmt(f'DO {loop_var} = 0, {stop}')
                elif len(args) == 2:
                    # range(start, stop) -> DO loop_var = start, stop-1
                    start_arg, stop_arg = args
                    if isinstance(start_arg, ast.Constant):
                        start = str(int(start_arg.value))
                    else:
                        se, _ = self._visit_expr(start_arg)
                        start = f'INT(DBLE({se}))'
                    if isinstance(stop_arg, ast.Constant):
                        stop = str(int(stop_arg.value) - 1)
                    else:
                        se, _ = self._visit_expr(stop_arg)
                        stop = f'INT(DBLE({se})) - 1'
                    self._add_stmt(
                        f'DO {loop_var} = {start}, {stop}')
                else:
                    raise NotImplementedError(
                        "range() with 3 arguments not supported")
            else:
                raise NotImplementedError(
                    "For loop must iterate over range()")

            for child in node.body:
                self._visit_stmt(child)
            self._add_stmt('END DO')
        elif isinstance(node, ast.Break):
            # break -> EXIT (Fortran loop exit)
            self._add_stmt('EXIT')
        elif isinstance(node, ast.While):
            raise NotImplementedError(
                "while loops are deliberately rejected in the Fortran "
                "translator because an unbounded loop can hang Abaqus "
                "forever without error. Use the bounded substitute: "
                "for k in range(N): ... if converged.real: break. "
                "See the design notes.")
        elif isinstance(node, ast.AugAssign):
            # x += expr  ->  x = x + expr
            target = node.target
            if not isinstance(target, ast.Name):
                raise NotImplementedError(
                    "Augmented assignment target must be a variable")
            varname = target.id
            kind = self._var_kinds.get(varname, 'scalar')
            expr, expr_kind = self._visit_expr(node.value)

            op_str = {
                ast.Add: '+', ast.Sub: '-',
                ast.Mult: '*', ast.Div: '/',
            }.get(type(node.op))
            if op_str is None:
                raise NotImplementedError(
                    f"Unsupported augmented assignment op: "
                    f"{type(node.op)}")

            self._add_stmt(f'{varname} = {varname} {op_str} {expr}')
        else:
            raise NotImplementedError(
                f"Unsupported statement: {ast.dump(node)}")

    def _visit_expr(self, node):
        """Returns (fortran_expr_str, kind) where kind is 'scalar' or 'tensor'."""

        if isinstance(node, ast.Constant):
            v = node.value
            if isinstance(v, (int, float)):
                return f'DCMPLX({_fortran_float(v)}, 0.0d0)', 'scalar'
            return str(v), 'scalar'

        elif isinstance(node, ast.Name):
            name = node.id
            kind = self._var_kinds.get(name, 'scalar')
            return name, kind

        elif isinstance(node, ast.Subscript):
            value, vk = self._visit_expr(node.value)

            def _fortran_index(idx_node):
                """Convert an index node to a Fortran 1-based integer index."""
                if isinstance(idx_node, ast.Constant):
                    return str(int(idx_node.value) + 1)
                idx_expr, _ = self._visit_expr(idx_node)
                return f'({idx_expr} + 1)'

            # Handle single index: arr[i] -> arr(i+1)
            # Handle multi-index: arr[i,j] -> arr(i+1,j+1)
            if isinstance(node.slice, ast.Tuple):
                if (vk == 'tensor' and len(node.slice.elts) == 2
                        and isinstance(node.slice.elts[0], ast.Slice)
                        and node.slice.elts[0].lower is None
                        and node.slice.elts[0].upper is None
                        and node.slice.elts[0].step is None):
                    col_idx = _fortran_index(node.slice.elts[1])
                    tmp = self._new_temp('COL')
                    self._add_decl('DOUBLE COMPLEX', tmp, '(3)')
                    self._var_kinds[tmp] = 'vector'
                    self._add_stmt(f'DO ii = 1, 3')
                    self._add_stmt(f'  {tmp}(ii) = {value}(ii,{col_idx})')
                    self._add_stmt(f'END DO')
                    return tmp, 'vector'
                indices = [_fortran_index(elt) for elt in node.slice.elts]
                idx_str = ','.join(indices)
                if ((vk == 'tensor' or str(vk).startswith('matrix'))
                        and len(indices) == 2):
                    result_kind = 'scalar'
                elif vk == 'tensor' and len(indices) == 1:
                    result_kind = 'vector'
                elif ((vk == 'vector' or str(vk).startswith('vector'))
                        and len(indices) == 1):
                    result_kind = 'scalar'
                else:
                    result_kind = 'scalar'
                return f'{value}({idx_str})', result_kind
            else:
                fortran_idx = _fortran_index(node.slice)
                result_kind = 'scalar'
                return f'{value}({fortran_idx})', result_kind

        elif isinstance(node, ast.Attribute):
            # np.pi / numpy.pi
            if (isinstance(node.value, ast.Name)
                    and node.value.id in ('np', 'numpy')
                    and node.attr == 'pi'):
                return 'DCMPLX(3.141592653589793d+00, 0.0d0)', 'scalar'
            # self.G -> props(i)
            elif isinstance(node.value, ast.Name) and node.value.id == 'self':
                prop_name = node.attr
                if prop_name in self.props_names:
                    idx = self.props_names.index(prop_name) + 1
                    return f'DCMPLX(props({idx}), 0.0d0)', 'scalar'
                else:
                    raise ValueError(
                        f"Unknown property: self.{prop_name}")
            # A.T -> transpose
            elif node.attr == 'T':
                inner, inner_kind = self._visit_expr(node.value)
                assert inner_kind == 'tensor', \
                    f".T applied to non-tensor: {inner}"
                tmp = self._new_temp('AT')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                self._add_stmt(f'CALL transpose33z({inner}, {tmp})')
                return tmp, 'tensor'
            # x.real -> DBLE(x)  (for CS-safe branching)
            elif node.attr == 'real':
                inner, _ = self._visit_expr(node.value)
                return f'DBLE({inner})', 'scalar'
            else:
                raise NotImplementedError(
                    f"Unsupported attribute: {ast.dump(node)}")

        elif isinstance(node, ast.BinOp):
            left, lk = self._visit_expr(node.left)
            right, rk = self._visit_expr(node.right)

            if isinstance(node.op, ast.MatMult):
                tmp = self._new_temp('MM')
                if rk == 'vector' or (rk == 'scalar' and lk == 'tensor'):
                    # Matrix-vector multiply: A @ x -> y (3-vector)
                    self._add_decl('DOUBLE COMPLEX', tmp, '(3)')
                    self._var_kinds[tmp] = 'vector'
                    self._add_stmt(f'CALL matvec33z({left}, {right}, {tmp})')
                    return tmp, 'vector'
                else:
                    self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                    self._var_kinds[tmp] = 'tensor'
                    self._add_stmt(f'CALL matmul33z({left}, {right}, {tmp})')
                    return tmp, 'tensor'

            op_str = {
                ast.Add: '+', ast.Sub: '-',
                ast.Mult: '*', ast.Div: '/',
                ast.Pow: '**',
            }.get(type(node.op))
            if op_str is None:
                raise NotImplementedError(
                    f"Unsupported binop: {type(node.op)}")

            # Determine result kind
            if lk == 'tensor' or rk == 'tensor':
                result_kind = 'tensor'
            elif lk == 'vector' or rk == 'vector':
                result_kind = 'vector'
            else:
                result_kind = 'scalar'

            return f'({left} {op_str} {right})', result_kind

        elif isinstance(node, ast.UnaryOp):
            operand, ok = self._visit_expr(node.operand)
            if isinstance(node.op, ast.USub):
                return f'(-{operand})', ok
            else:
                raise NotImplementedError(
                    f"Unsupported unaryop: {type(node.op)}")

        elif isinstance(node, ast.Call):
            if self._is_self_helper_call(node):
                return self._emit_helper_call(node)

            func_name = self._get_func_name(node)
            if func_name == 'array':
                return self._emit_np_array(node)
            args_with_kinds = [self._visit_expr(a) for a in node.args]
            args = [a[0] for a in args_with_kinds]

            if func_name == 'det':
                return f'det33z({args[0]})', 'scalar'
            elif func_name == 'inv':
                tmp = self._new_temp('INV')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                self._add_stmt(f'CALL inv33z({args[0]}, {tmp})')
                return tmp, 'tensor'
            elif func_name == 'log':
                return f'LOG({args[0]})', 'scalar'
            elif func_name == 'exp':
                return f'EXP({args[0]})', 'scalar'
            elif func_name == 'sqrt':
                return f'SQRT({args[0]})', 'scalar'
            elif func_name in ('sin', 'cos', 'tan',
                               'sinh', 'cosh', 'tanh',
                               'arctan', 'atan', 'arccos'):
                # CS-safe complex trig/hyperbolic intrinsics
                fn_map = {
                    'sin': 'cs_sin', 'cos': 'cs_cos', 'tan': 'cs_tan',
                    'sinh': 'cs_sinh', 'cosh': 'cs_cosh',
                    'tanh': 'cs_tanh',
                    'arctan': 'cs_atan',
                    'atan': 'cs_atan',
                    'arccos': 'cs_acos',
                }
                fname = fn_map[func_name]
                return f'{fname}({args[0]})', 'scalar'
            elif func_name in ('erf', 'erfc'):
                # CS-safe complex error functions
                fname = 'cs_erf' if func_name == 'erf' else 'cs_erfc'
                return f'{fname}({args[0]})', 'scalar'
            elif func_name == 'isnan':
                # NaN is the only value that compares unequal to itself.
                # Works for DOUBLE COMPLEX in Fortran 90+.
                return f'({args[0]} /= {args[0]})', 'scalar'
            elif func_name == 'solve':
                # np.linalg.solve(A, b) -> inline small dense solve
                if len(args) != 2:
                    raise NotImplementedError(
                        "np.linalg.solve requires exactly 2 arguments")
                # Try to infer size of A from variable kind
                A_node = node.args[0]
                A_kind = None
                if isinstance(A_node, ast.Name):
                    A_kind = self._var_kinds.get(A_node.id)
                size_map = {'tensor': 3}
                N = size_map.get(A_kind)
                if N is None and isinstance(A_kind, str):
                    if A_kind.startswith('matrix'):
                        N = int(A_kind[len('matrix'):])
                if N is None:
                    raise NotImplementedError(
                        "np.linalg.solve matrix size cannot be inferred "
                        "statically. Supported sizes: N=2,3,4,6,9. "
                        "Ensure the matrix is a local variable with a "
                        "known shape (e.g. 3x3 tensor).")
                tmp = self._new_temp('SOL')
                solve_kind = 'vector' if N == 3 else f'vector{N}'
                self._add_decl('DOUBLE COMPLEX', tmp, f'({N})')
                self._var_kinds[tmp] = solve_kind
                self._add_stmt(
                    f'CALL cs_solve_{N}({args[0]}, {args[1]}, {tmp})')
                return tmp, solve_kind
            elif func_name == 'trace':
                A = args[0]
                return self._trace_expr(A), 'scalar'
            elif func_name in ('mean_stress', 'p_eff'):
                A = args[0]
                return (
                    f'({self._trace_expr(A)} / '
                    f'DCMPLX(3.000000000000000d+00, 0.0d0))',
                    'scalar')
            elif func_name in ('dev_stress', 'dev'):
                return self._dev_tensor(args[0]), 'tensor'
            elif func_name == 'q_mises':
                S = self._dev_tensor(args[0])
                ST = self._new_temp('AT')
                self._add_decl('DOUBLE COMPLEX', ST, '(3,3)')
                self._var_kinds[ST] = 'tensor'
                self._add_stmt(f'CALL transpose33z({S}, {ST})')
                MM = self._new_temp('MM')
                self._add_decl('DOUBLE COMPLEX', MM, '(3,3)')
                self._var_kinds[MM] = 'tensor'
                self._add_stmt(f'CALL matmul33z({ST}, {S}, {MM})')
                return (
                    f'SQRT(DCMPLX(1.500000000000000d+00, 0.0d0) * '
                    f'{self._trace_expr(MM)})',
                    'scalar')
            elif func_name == 'flow_direction':
                if len(args) < 2:
                    raise NotImplementedError(
                        "flow_direction requires sigma and q arguments")
                tol = args[2] if len(args) > 2 else \
                    'DCMPLX(1.000000000000000d-30, 0.0d0)'
                S = self._dev_tensor(args[0])
                return (
                    f'(DCMPLX(1.500000000000000d+00, 0.0d0) / '
                    f'({args[1]} + {tol})) * {S}',
                    'tensor')
            elif func_name == 'bulk_modulus':
                if len(args) != 2:
                    raise NotImplementedError(
                        "bulk_modulus requires G and lam arguments")
                return (
                    f'({args[1]} + DCMPLX(2.000000000000000d+00, 0.0d0)'
                    f' * {args[0]} / DCMPLX(3.000000000000000d+00, 0.0d0))',
                    'scalar')
            elif func_name == 'eye':
                # eye(3) -> 3x3 identity matrix
                tmp = self._new_temp('EYE')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                self._add_stmt(f'CALL eye33z({tmp})')
                return tmp, 'tensor'
            elif func_name == 'dyad':
                # dyad(a, b) -> outer product, vector(3) x vector(3) -> tensor(3,3)
                if len(args) != 2:
                    raise NotImplementedError(
                        "dyad requires exactly 2 arguments (two 3-vectors)")
                tmp = self._new_temp('DYAD')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                self._add_stmt(
                    f'CALL outer33z({args[0]}, {args[1]}, {tmp})')
                return tmp, 'tensor'
            elif func_name == 'sym3':
                # sym3(a11, a22, a33, a12, a13, a23) -> symmetric tensor(3,3)
                # Constructs an arbitrary symmetric 3x3 tensor from its 6
                # unique entries. Used for constant tensors (base fabrics,
                # anisotropy directors, projection operators) that can't be
                # built from vector literals (not a translator primitive).
                if len(args) != 6:
                    raise NotImplementedError(
                        "sym3 requires exactly 6 arguments: "
                        "(a11, a22, a33, a12, a13, a23)")
                tmp = self._new_temp('SYM3')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                self._add_stmt(
                    f'CALL sym3z({args[0]}, {args[1]}, {args[2]}, '
                    f'{args[3]}, {args[4]}, {args[5]}, {tmp})')
                return tmp, 'tensor'
            elif func_name == 'abs':
                # abs(x) -> ABS(x)
                # CS safety: should only be applied to real values
                # (e.g., abs(f.real), not abs(complex_var))
                return f'ABS({args[0]})', 'scalar'
            elif func_name == 'sqrtm':
                tmp = self._new_temp('SQRTM')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                if self.matrix_backend == 'iterative':
                    tmp_inv = self._new_temp('SQRTMINV')
                    self._add_decl('DOUBLE COMPLEX', tmp_inv, '(3,3)')
                    self._var_kinds[tmp_inv] = 'tensor'
                    self._add_stmt(
                        f'CALL sqrtm33z_iter({args[0]}, {tmp}, {tmp_inv})')
                else:
                    self._add_stmt(
                        f'CALL sqrtm33z({args[0]}, {tmp})')
                return tmp, 'tensor'
            elif func_name == 'logm':
                tmp = self._new_temp('LOGM')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                if self.matrix_backend == 'iterative':
                    self._add_stmt(
                        f'CALL logm33z_iter({args[0]}, {tmp})')
                else:
                    self._add_stmt(
                        f'CALL logm33z({args[0]}, {tmp})')
                return tmp, 'tensor'
            elif func_name == 'expm':
                tmp = self._new_temp('EXPM')
                self._add_decl('DOUBLE COMPLEX', tmp, '(3,3)')
                self._var_kinds[tmp] = 'tensor'
                if self.matrix_backend == 'iterative':
                    self._add_stmt(
                        f'CALL expm33z_iter({args[0]}, {tmp})')
                else:
                    self._add_stmt(
                        f'CALL expm33z({args[0]}, {tmp})')
                return tmp, 'tensor'
            elif func_name == 'int':
                # int(x) -> INT(x) — used in range() loops
                return f'INT({args[0]})', 'scalar'
            elif func_name == 'float':
                # float(x) -> DBLE(x)
                return f'DBLE({args[0]})', 'scalar'
            elif func_name == 'complex':
                # complex(re, im) -> DCMPLX(re, im)
                if len(args) == 2:
                    return f'DCMPLX({args[0]}, {args[1]})', 'scalar'
                else:
                    return f'DCMPLX({args[0]}, 0.0d0)', 'scalar'
            else:
                raise NotImplementedError(
                    f"Unsupported function: {func_name}")

        elif isinstance(node, ast.Compare):
            # Support comparisons for if-statements (CS-safe branching)
            if len(node.ops) > 1:
                raise NotImplementedError(
                    "Chained comparisons are not supported in material "
                    "methods (e.g. '0 < x < 1'). Rewrite as "
                    "'(a < b) and (b < c)' or use explicit intermediate "
                    "variables."
                )
            left, left_kind = self._visit_expr(node.left)
            # In the complex-step tangent engine all variables are
            # DOUBLE COMPLEX. Comparisons on complex variables are
            # invalid in Fortran. The user must explicitly extract the
            # real part (e.g. x.real > 0) which the generator maps to
            # DBLE(x) .GT. 0.0d0.
            def _is_real_safe(n):
                if isinstance(n, ast.Constant):
                    return True
                if isinstance(n, ast.Name):
                    # for-loop indices are emitted as Fortran INTEGER,
                    # so direct comparison is type-promoted and safe.
                    if n.id in self._integer_vars:
                        return True
                    # Variables assigned from x.real are DBLE(x) and real-safe.
                    if n.id in self._real_vars:
                        return True
                if isinstance(n, ast.Attribute):
                    # x.real is real-safe.
                    if n.attr == 'real':
                        return True
                    # self.<prop> is a material property; the generator
                    # emits it as PROPS(K) which is REAL*8 in Fortran,
                    # never a CS-perturbed complex variable.
                    if isinstance(n.value, ast.Name) and n.value.id == 'self':
                        return True
                if isinstance(n, ast.Call):
                    fn = self._get_func_name(n)
                    if fn in ('abs', 'float', 'int'):
                        return True
                if isinstance(n, ast.BinOp):
                    return _is_real_safe(n.left) and _is_real_safe(n.right)
                if isinstance(n, ast.UnaryOp):
                    return _is_real_safe(n.operand)
                return False

            if not _is_real_safe(node.left):
                var_name = getattr(node.left, 'id', ast.dump(node.left))
                raise NotImplementedError(
                    f"Comparison on potentially-complex expression "
                    f"{var_name!r} is not Fortran-safe. "
                    f"Use {var_name}.real or rewrite to avoid branching."
                )
            op = node.ops[0]
            # For the comparator, use a plain real literal if it's
            # a constant (not DCMPLX wrapped) — comparisons are real
            comp_node = node.comparators[0]
            if isinstance(comp_node, ast.Constant):
                right = _fortran_float(comp_node.value)
            elif (isinstance(comp_node, ast.UnaryOp)
                  and isinstance(comp_node.op, ast.USub)
                  and isinstance(comp_node.operand, ast.Constant)):
                right = '-' + _fortran_float(comp_node.operand.value)
            elif (isinstance(comp_node, ast.Attribute)
                  and isinstance(comp_node.value, ast.Name)
                  and comp_node.value.id == 'self'
                  and comp_node.attr in self.props_names):
                idx = self.props_names.index(comp_node.attr) + 1
                right = f'props({idx})'
            else:
                right, _ = self._visit_expr(comp_node)
            op_map = {
                ast.Gt: '.GT.', ast.GtE: '.GE.',
                ast.Lt: '.LT.', ast.LtE: '.LE.',
                ast.Eq: '.EQ.', ast.NotEq: '.NE.',
            }
            op_str = op_map.get(type(op))
            if op_str is None:
                raise NotImplementedError(
                    f"Unsupported comparison: {type(op)}")
            return f'{left} {op_str} {right}', 'scalar'

        else:
            raise NotImplementedError(
                f"Unsupported expression: {ast.dump(node)}")

    def _emit_np_array(self, node):
        """Translate fixed-size ``np.array`` literals used in local Newton code."""
        if len(node.args) != 1:
            raise NotImplementedError("np.array with one positional argument only")
        value = node.args[0]
        if not isinstance(value, (ast.List, ast.Tuple)):
            raise NotImplementedError("np.array argument must be a literal list")
        rows = value.elts
        if not rows:
            raise NotImplementedError("empty np.array is not supported")

        if all(not isinstance(r, (ast.List, ast.Tuple)) for r in rows):
            n = len(rows)
            kind = 'vector' if n == 3 else f'vector{n}'
            tmp = self._new_temp('ARR')
            self._add_decl('DOUBLE COMPLEX', tmp, f'({n})')
            self._var_kinds[tmp] = kind
            for i, elt in enumerate(rows, start=1):
                expr, _ = self._visit_expr(elt)
                self._add_stmt(f'{tmp}({i}) = {expr}')
            return tmp, kind

        if all(isinstance(r, (ast.List, ast.Tuple)) for r in rows):
            nrow = len(rows)
            ncol = len(rows[0].elts)
            if any(len(r.elts) != ncol for r in rows):
                raise NotImplementedError("ragged np.array is not supported")
            if nrow != ncol:
                raise NotImplementedError("only square np.array matrices supported")
            kind = 'tensor' if nrow == 3 else f'matrix{nrow}'
            tmp = self._new_temp('ARR')
            self._add_decl('DOUBLE COMPLEX', tmp, f'({nrow},{ncol})')
            self._var_kinds[tmp] = kind
            for i, row in enumerate(rows, start=1):
                for j, elt in enumerate(row.elts, start=1):
                    expr, _ = self._visit_expr(elt)
                    self._add_stmt(f'{tmp}({i},{j}) = {expr}')
            return tmp, kind

        raise NotImplementedError("mixed-rank np.array is not supported")

    def _get_func_name(self, call_node):
        if isinstance(call_node.func, ast.Name):
            # Bare 'eigh' is treated as alias for 'eig'
            name = call_node.func.id
            if name == 'eigh':
                return 'eig'
            return name
        elif isinstance(call_node.func, ast.Attribute):
            attr = call_node.func
            # Recognise np.linalg.eigh / numpy.linalg.eigh as alias for eig
            if (attr.attr == 'eigh'
                    and isinstance(attr.value, ast.Attribute)
                    and attr.value.attr == 'linalg'
                    and isinstance(attr.value.value, ast.Name)
                    and attr.value.value.id in ('np', 'numpy')):
                return 'eig'
            return attr.attr
        else:
            raise NotImplementedError(
                f"Cannot extract func name: "
                f"{ast.dump(call_node.func)}")


# ===================================================================
# Fortran code helpers
# ===================================================================

def _fortran_float(v):
    """Convert Python float/int to Fortran double literal."""
    if isinstance(v, int):
        return f'{v}.0d0'
    s = f'{v:.15e}'.replace('e', 'd')
    return s


from ._fortran_format import wrap_fortran_line as _wrap_fortran_line


def _wrap_fortran_source(src):
    """Wrap all lines in a Fortran source string."""
    wrapped = []
    for line in src.split('\n'):
        wrapped.extend(
            part.rstrip() for part in _wrap_fortran_line(line).split('\n'))
    return '\n'.join(wrapped)


def _fortran_call(indent, sub_name, args):
    """
    Generate a Fortran CALL statement with proper continuation lines.

    If the call fits in 72 columns, emits a single line.
    Otherwise, breaks after commas with continuation in column 6.

    Args:
        indent: leading whitespace (e.g. '      ' or '          ')
        sub_name: subroutine name
        args: list of argument strings
    """
    # Try single-line first
    single = f'{indent}CALL {sub_name}({", ".join(args)})'
    if len(single) <= 72:
        return single

    # Multi-line: first arg on same line as CALL, rest continued
    lines = [f'{indent}CALL {sub_name}({args[0]},']
    for i, arg in enumerate(args[1:], 1):
        if i < len(args) - 1:
            lines.append(f'     &  {arg},')
        else:
            lines.append(f'     &  {arg})')
    return '\n'.join(lines)


def _tensorize_expr(expr, tensor_vars):
    """
    Replace tensor variable names with (ii,jj) indexed forms.

    Only variables in tensor_vars get indexed. Scalars are left alone.
    E.g., with tensor_vars={'F', 'AT2'}:
      "(G * (F - AT2))" -> "(G * (F(ii,jj) - AT2(ii,jj)))"
    """
    import re
    sorted_names = sorted(tensor_vars, key=len, reverse=True)
    result = expr
    protected = {}

    # Protect matrix-to-scalar helper calls before the component-wise rewrite below.
    # det33z takes a WHOLE 3x3 tensor and returns a scalar, so its argument must NOT
    # be turned into det33z(F(ii,jj)) -- that compiles (Fortran has no explicit
    # interface here) but passes a single component and returns NaNs.  Shield the
    # ENTIRE call by scanning balanced parentheses, so both a bare-name argument
    # det33z(F) and an expression argument det33z(F - G) are covered (the bare-name-
    # only guard missed the latter, generating det33z(F(ii,jj) - G(ii,jj))).
    def _protect_calls(s, func):
        pat = re.compile(r'\b' + re.escape(func) + r'\s*\(')
        out, i = [], 0
        while i < len(s):
            m = pat.match(s, i)
            if not m:
                out.append(s[i]); i += 1
                continue
            depth, j = 0, m.end() - 1            # j at the opening '('
            while j < len(s):
                if s[j] == '(':
                    depth += 1
                elif s[j] == ')':
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            key = f'__ABQ_UFL_PROTECTED_{len(protected)}__'
            protected[key] = s[i:j + 1]
            out.append(key); i = j + 1
        return ''.join(out)

    for func in ('det33z',):
        result = _protect_calls(result, func)
    for name in sorted_names:
        result = re.sub(r'\b' + re.escape(name) + r'\b(?!\()',
                        f'{name}(ii,jj)', result)
    for key, value in protected.items():
        result = result.replace(key, value)
    return result


def _fortran_init_value(name, sv_entry):
    """
    Generate Fortran lines to initialize a state variable to its default.

    Returns list of Fortran source lines (no leading whitespace).
    """
    shape = sv_entry['shape']
    init = sv_entry['init']
    init_prop_index = sv_entry.get('init_prop_index')
    lines = []
    if shape == 'scalar':
        if init_prop_index is not None:
            lines.append(f'{name} = PROPS({init_prop_index})')
        else:
            lines.append(f'{name} = {_fortran_float(init)}')
    elif shape == 'vector':
        for i in range(3):
            lines.append(
                f'{name}({i+1}) = {_fortran_float(float(init[i]))}')
    elif shape == 'tensor':
        # Column-major: Fortran stores (i,j) with i varying fastest
        for j in range(3):
            for i in range(3):
                lines.append(
                    f'{name}({i+1},{j+1}) = '
                    f'{_fortran_float(float(init[i, j]))}')
    return lines


def _fortran_dims(sv_entry):
    """Return Fortran dimension string for a state variable."""
    shape = sv_entry['shape']
    if shape == 'scalar':
        return ''
    elif shape == 'vector':
        return '(3)'
    elif shape == 'tensor':
        return '(3,3)'


# ===================================================================
# STATEV read/write code generation
# ===================================================================

def _generate_statev_read(sv_info):
    """
    Generate STATEV -> local variable read code.

    Follows the standard Abaqus STATEV pattern: if time(2)==0, use initial
    values; else read from STATEV array.

    For UMAT, there's exactly 1 GP per call, so no GP offset.
    """
    if not sv_info:
        return ''

    lines = []
    lines.append('C     --- Read state variables from STATEV ---')
    lines.append('C     (Standard Abaqus STATEV init pattern, TIME(2)==0)')
    lines.append('      IF (TIME(2) .EQ. 0.0d0) THEN')

    # Initial values
    for name, entry in sv_info.items():
        for line in _fortran_init_value(name + '_old', entry):
            lines.append(f'        {line}')

    lines.append('      ELSE')

    # Read from STATEV
    offset = 0
    for name, entry in sv_info.items():
        shape = entry['shape']
        if shape == 'scalar':
            lines.append(
                f'        {name}_old = STATEV({offset + 1})')
            offset += 1
        elif shape == 'vector':
            for i in range(3):
                lines.append(
                    f'        {name}_old({i+1}) = '
                    f'STATEV({offset + i + 1})')
            offset += 3
        elif shape == 'tensor':
            # Column-major storage: offset + (j-1)*3 + i
            for j in range(3):
                for i in range(3):
                    idx = offset + j * 3 + i + 1
                    lines.append(
                        f'        {name}_old({i+1},{j+1}) = '
                        f'STATEV({idx})')
            offset += 9

    lines.append('      END IF')
    return '\n'.join(lines)


def _generate_statev_write(sv_info):
    """Generate local variable -> STATEV write code (real part only)."""
    if not sv_info:
        return ''

    lines = []
    lines.append(
        'C     --- Write updated state to STATEV (real part only) ---')

    offset = 0
    for name, entry in sv_info.items():
        shape = entry['shape']
        if shape == 'scalar':
            lines.append(
                f'      STATEV({offset + 1}) = '
                f'DBLE({name}_new_z)')
            offset += 1
        elif shape == 'vector':
            for i in range(3):
                lines.append(
                    f'      STATEV({offset + i + 1}) = '
                    f'DBLE({name}_new_z({i+1}))')
            offset += 3
        elif shape == 'tensor':
            for j in range(3):
                for i in range(3):
                    idx = offset + j * 3 + i + 1
                    lines.append(
                        f'      STATEV({idx}) = '
                        f'DBLE({name}_new_z({i+1},{j+1}))')
            offset += 9

    return '\n'.join(lines)


# ===================================================================
# Material stress subroutine generation
# ===================================================================

def _generate_stress_subroutine(material, mat_prefix, matrix_backend='iterative'):
    """
    Generate the COMPLEX*16 stress subroutine from the Python method.

    For stateless materials:
      SUBROUTINE prefix_stress_PK1(F, props, P_out)

    For materials with state_vars:
      SUBROUTINE prefix_stress_PK1(F, ep_old, Fp_old, dt, props,
                                    P_out, ep_new, Fp_new)

    Argument ordering:
      1. Field variables (COMPLEX):     F
      2. State variables old (COMPLEX):  ep_old, Fp_old, ...
      3. Real parameters:               dt, props
      4. Output (COMPLEX):              P_out
      5. State variables new (COMPLEX):  ep_new, Fp_new, ...
    """
    sv_info = _state_var_info(material)
    has_state = bool(sv_info)

    method = material._methods['stress_PK1']['callable']

    # Register state variable kinds in translator
    state_var_names = list(sv_info.keys()) if has_state else []
    translator = FortranTranslator(material.props_names, mat_prefix,
                                   state_var_names, material=material,
                                   matrix_backend=matrix_backend)
    # Pre-register state variable kinds so translator knows their types
    for name, entry in sv_info.items():
        if entry['shape'] == 'tensor':
            translator._var_kinds[name + '_old'] = 'tensor'
        else:
            translator._var_kinds[name + '_old'] = 'scalar'

    decls, stmts, return_expr, return_kind = \
        translator.translate_stress(method)

    nprops = len(material.props_names)
    tensor_vars = {n for n, k in translator._var_kinds.items()
                   if k == 'tensor'}

    lines = []

    # --- Subroutine signature ---
    if has_state:
        sig_args = ['F']
        for name in sv_info:
            sig_args.append(f'{name}_old')
        sig_args.extend(['dt', 'props', 'P_out'])
        for name in sv_info:
            sig_args.append(f'{name}_new')
        sig_str = ', '.join(sig_args)
        lines.append(
            f'      SUBROUTINE {mat_prefix}_stress_PK1({sig_str})')
    else:
        lines.append(
            f'      SUBROUTINE {mat_prefix}_stress_PK1'
            f'(F, props, P_out)')

    lines.append('      IMPLICIT NONE')

    # --- Argument declarations ---
    lines.append('      DOUBLE COMPLEX, INTENT(IN)  :: F(3,3)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN)  :: '
                f'{name}_old{dims}')
        lines.append(
            '      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.append(
        f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')
    lines.append(
        '      DOUBLE COMPLEX, INTENT(OUT) :: P_out(3,3)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(OUT) :: '
                f'{name}_new{dims}')

    lines.append('')

    # External function declarations
    lines.append('      DOUBLE COMPLEX :: det33z')
    lines.append('      DOUBLE COMPLEX :: cs_sin')
    lines.append('      DOUBLE COMPLEX :: cs_cos')
    lines.append('      DOUBLE COMPLEX :: cs_tan')
    lines.append('      DOUBLE COMPLEX :: cs_sinh')
    lines.append('      DOUBLE COMPLEX :: cs_cosh')
    lines.append('      DOUBLE COMPLEX :: cs_tanh')
    lines.append('      DOUBLE COMPLEX :: cs_atan')
    lines.append('      DOUBLE COMPLEX :: cs_acos')
    lines.append('      DOUBLE COMPLEX :: cs_erf')
    lines.append('      DOUBLE COMPLEX :: cs_erfc')

    # Track state variable names to avoid double-declaration
    state_new_names = {f'{n}_new' for n in sv_info}
    state_old_names = {f'{n}_old' for n in sv_info}
    skip_names = state_new_names | state_old_names | {'P_out', 'dt'}

    # Declare temporaries from translator
    # Skip state_new names (already declared as INTENT(OUT) args)
    declared_names = set()
    for ftype, name, dims in decls:
        if name in declared_names:
            continue
        if name in skip_names:
            declared_names.add(name)
            continue
        if name in state_new_names:
            declared_names.add(name)
            continue
        dim_str = dims if dims else ''
        lines.append(f'      {ftype} :: {name}{dim_str}')
        declared_names.add(name)

    # Declare user variables (from assignment statements)
    for stmt in stmts:
        if '=' in stmt and 'CALL' not in stmt:
            stripped = stmt.strip()
            if stripped.startswith(('IF ', 'ELSE', 'END IF',
                                   'DO ', 'END DO', 'EXIT')):
                continue
            varname = stmt.split('=')[0].strip()
            if (varname not in declared_names
                    and varname not in skip_names):
                if '(' in varname:
                    continue
                if varname in tensor_vars:
                    pass  # Already declared by translator
                elif varname in translator._real_vars:
                    lines.append(
                        f'      DOUBLE PRECISION :: {varname}')
                    declared_names.add(varname)
                else:
                    lines.append(
                        f'      DOUBLE COMPLEX :: {varname}')
                    declared_names.add(varname)

    lines.append('      INTEGER :: ii, jj')
    declared_names.update({'ii', 'jj'})

    # Declare loop variables from for-loops as INTEGER
    for stmt in stmts:
        stripped = stmt.strip()
        if stripped.startswith('DO ') and '=' in stripped:
            # Extract loop variable: "DO var = ..."
            loop_var = stripped[3:].split('=')[0].strip()
            if loop_var not in declared_names:
                lines.append(f'      INTEGER :: {loop_var}')
                declared_names.add(loop_var)

    lines.append('')

    # Body statements
    for stmt in stmts:
        lines.append(f'      {stmt}')

    # Return: assign result to P_out
    if return_kind == 'tensor':
        if (return_expr in tensor_vars
                or return_expr in declared_names):
            lines.append('')
            lines.append('      DO ii = 1, 3')
            lines.append('        DO jj = 1, 3')
            lines.append(
                f'          P_out(ii,jj) = {return_expr}(ii,jj)')
            lines.append('        END DO')
            lines.append('      END DO')
        else:
            lines.append('')
            lines.append('      DO ii = 1, 3')
            lines.append('        DO jj = 1, 3')
            expr_ij = _tensorize_expr(
                return_expr, tensor_vars | {'F'})
            lines.append(f'          P_out(ii,jj) = {expr_ij}')
            lines.append('        END DO')
            lines.append('      END DO')
    else:
        lines.append(f'      P_out = {return_expr}')

    # Assign state_new outputs from the return dict
    if has_state and translator._state_return_map:
        has_copy = False
        for sv_name, (val_expr, val_kind) in \
                translator._state_return_map.items():
            # Skip if source is already the output variable
            if val_expr == f'{sv_name}_new':
                continue
            if not has_copy:
                lines.append('')
                lines.append('C     --- Copy state outputs ---')
                has_copy = True
            if val_kind == 'tensor':
                if (val_expr in tensor_vars
                        or val_expr in declared_names):
                    lines.append('      DO ii = 1, 3')
                    lines.append('        DO jj = 1, 3')
                    lines.append(
                        f'          {sv_name}_new(ii,jj) = '
                        f'{val_expr}(ii,jj)')
                    lines.append('        END DO')
                    lines.append('      END DO')
                else:
                    lines.append(
                        f'      {sv_name}_new = {val_expr}')
            else:
                lines.append(
                    f'      {sv_name}_new = {val_expr}')

    lines.append('')
    lines.append('      RETURN')
    lines.append(
        f'      END SUBROUTINE {mat_prefix}_stress_PK1')

    if translator._helper_sources:
        lines.append('')
        lines.append('C' + '=' * 70)
        lines.append('C     Material helper methods')
        lines.append('C' + '=' * 70)
        lines.extend(translator._helper_sources)

    return '\n'.join(lines)


def _predef_var_names(material):
    """Scalar PREDEF arguments requested by a small-strain material."""
    predef_vars = getattr(material, 'predef_vars', ())
    if isinstance(predef_vars, str):
        return [predef_vars]
    return list(predef_vars)


def _small_strain_update_method_info(material, update_method):
    """Return method metadata for the small-strain update method."""
    import inspect

    if update_method == 'stress_update':
        if 'stress_update' not in material._methods:
            raise ValueError(
                "Small-strain UMAT material must define stress_update")
        return material._methods['stress_update']

    method = getattr(material, update_method, None)
    if method is None or not callable(method):
        raise ValueError(
            f"Small-strain UMAT update_method='{update_method}' was not "
            f"found on {material.__class__.__name__}")

    sig = inspect.signature(method)
    params = list(sig.parameters.keys())
    return {
        'callable': method,
        'all_params': params,
        'field_params': [p for p in params if p in material._field_args],
        'history_params': [p for p in params if p in material._history_args],
        'param_params': [p for p in params if p in material._param_args],
        'state_params': [p for p in params if p in material._state_arg_names],
        'output_shape': (3, 3),
        'description': 'small-strain stress update',
    }


def _generate_small_strain_update_subroutine(material, mat_prefix,
                                             predef_names=None,
                                             update_method='stress_update'):
    """Generate COMPLEX*16 small-strain stress_update subroutine."""
    sv_info = _state_var_info(material)
    has_state = bool(sv_info)
    predef_names = list(predef_names or [])
    info = _small_strain_update_method_info(material, update_method)
    method = info['callable']
    input_tensors = {
        'sigma_old': 'tensor',
        'strain_old': 'tensor',
        'dstrain': 'tensor',
    }
    initial_var_kinds = dict(input_tensors)
    for name in predef_names:
        initial_var_kinds[name] = 'scalar'
    state_var_names = list(sv_info.keys()) if has_state else []
    translator = FortranTranslator(
        material.props_names, mat_prefix, state_var_names,
        initial_var_kinds=initial_var_kinds,
        material=material)
    for name, entry in sv_info.items():
        if entry['shape'] == 'tensor':
            translator._var_kinds[name + '_old'] = 'tensor'
        else:
            translator._var_kinds[name + '_old'] = 'scalar'

    decls, stmts, return_expr, return_kind = translator.translate_stress(method)
    nprops = len(material.props_names)
    tensor_vars = {n for n, k in translator._var_kinds.items()
                   if k == 'tensor'}

    sig_args = ['sigma_old', 'strain_old', 'dstrain']
    if has_state:
        for name in sv_info:
            sig_args.append(f'{name}_old')
    for name in predef_names:
        sig_args.append(name)
    if has_state:
        sig_args.append('dt')
    sig_args.extend(['cutback_flag', 'props', 'sigma_new'])
    if has_state:
        for name in sv_info:
            sig_args.append(f'{name}_new')

    lines = [
        f'      SUBROUTINE {mat_prefix}_stress_update('
        + ', '.join(sig_args) + ')',
        '      IMPLICIT NONE',
        '      DOUBLE COMPLEX, INTENT(IN)  :: sigma_old(3,3)',
        '      DOUBLE COMPLEX, INTENT(IN)  :: strain_old(3,3)',
        '      DOUBLE COMPLEX, INTENT(IN)  :: dstrain(3,3)',
    ]
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN)  :: {name}_old{dims}')
    for name in predef_names:
        lines.append(f'      DOUBLE COMPLEX, INTENT(IN)  :: {name}')
    if has_state:
        lines.append('      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.extend([
        '      DOUBLE COMPLEX, INTENT(OUT) :: cutback_flag',
        f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})',
        '      DOUBLE COMPLEX, INTENT(OUT) :: sigma_new(3,3)',
    ])
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(OUT) :: {name}_new{dims}')

    lines.append('')
    lines.append('      DOUBLE COMPLEX :: det33z')
    lines.append('      DOUBLE COMPLEX :: cs_sin')
    lines.append('      DOUBLE COMPLEX :: cs_cos')
    lines.append('      DOUBLE COMPLEX :: cs_tan')
    lines.append('      DOUBLE COMPLEX :: cs_sinh')
    lines.append('      DOUBLE COMPLEX :: cs_cosh')
    lines.append('      DOUBLE COMPLEX :: cs_tanh')
    lines.append('      DOUBLE COMPLEX :: cs_atan')
    lines.append('      DOUBLE COMPLEX :: cs_acos')
    lines.append('      DOUBLE COMPLEX :: cs_erf')
    lines.append('      DOUBLE COMPLEX :: cs_erfc')

    state_new_names = {f'{n}_new' for n in sv_info}
    state_old_names = {f'{n}_old' for n in sv_info}
    input_names = set(initial_var_kinds)
    skip_names = state_new_names | state_old_names | input_names | {
        'sigma_new', 'cutback_flag', 'dt'}

    declared_names = set()
    for ftype, name, dims in decls:
        if name in declared_names:
            continue
        if name in skip_names:
            declared_names.add(name)
            continue
        if name in state_new_names:
            declared_names.add(name)
            continue
        dim_str = dims if dims else ''
        lines.append(f'      {ftype} :: {name}{dim_str}')
        declared_names.add(name)

    for stmt in stmts:
        if '=' in stmt and 'CALL' not in stmt:
            stripped = stmt.strip()
            if stripped.startswith(('IF ', 'ELSE', 'END IF',
                                   'DO ', 'END DO', 'EXIT')):
                continue
            varname = stmt.split('=')[0].strip()
            if '(' in varname:
                continue
            if varname not in declared_names and varname not in skip_names:
                if varname not in tensor_vars:
                    lines.append(f'      DOUBLE COMPLEX :: {varname}')
                    declared_names.add(varname)

    lines.append('      INTEGER :: ii, jj')
    for stmt in stmts:
        stripped = stmt.strip()
        if stripped.startswith('DO ') and '=' in stripped:
            loop_var = stripped[3:].split('=')[0].strip()
            if loop_var not in declared_names:
                lines.append(f'      INTEGER :: {loop_var}')
                declared_names.add(loop_var)

    lines.append('')
    lines.append('      cutback_flag = DCMPLX(0.0d0, 0.0d0)')
    lines.append('')
    for stmt in stmts:
        lines.append(f'      {stmt}')

    if return_kind == 'tensor':
        lines.append('')
        lines.append('      DO ii = 1, 3')
        lines.append('        DO jj = 1, 3')
        if return_expr in tensor_vars or return_expr in declared_names:
            lines.append(
                f'          sigma_new(ii,jj) = {return_expr}(ii,jj)')
        else:
            expr_ij = _tensorize_expr(return_expr, tensor_vars | input_names)
            lines.append(f'          sigma_new(ii,jj) = {expr_ij}')
        lines.append('        END DO')
        lines.append('      END DO')
    else:
        lines.append(f'      sigma_new = {return_expr}')

    if has_state and translator._state_return_map:
        has_copy = False
        for sv_name, (val_expr, val_kind) in \
                translator._state_return_map.items():
            if sv_name in ('_cutback', 'request_cutback'):
                lines.append(f'      cutback_flag = {val_expr}')
                continue
            if sv_name not in sv_info:
                raise ValueError(
                    f"stress_update returned unknown state key '{sv_name}'. "
                    "Use a declared state variable name, or '_cutback' for "
                    "Abaqus PNEWDT cutback requests.")
            if val_expr == f'{sv_name}_new':
                continue
            if not has_copy:
                lines.append('')
                lines.append('C     --- Copy state outputs ---')
                has_copy = True
            if val_kind == 'tensor':
                lines.append('      DO ii = 1, 3')
                lines.append('        DO jj = 1, 3')
                lines.append(
                    f'          {sv_name}_new(ii,jj) = {val_expr}(ii,jj)')
                lines.append('        END DO')
                lines.append('      END DO')
            else:
                lines.append(f'      {sv_name}_new = {val_expr}')

    lines.extend([
        '',
        '      RETURN',
        f'      END SUBROUTINE {mat_prefix}_stress_update',
    ])
    if translator._helper_sources:
        lines.append('')
        lines.append('C' + '=' * 70)
        lines.append('C     Material helper methods')
        lines.append('C' + '=' * 70)
        lines.extend(translator._helper_sources)
    return '\n'.join(lines)


# ===================================================================
# Complex-step tangent engine generation
# ===================================================================

def _generate_cs_dPdF(mat_prefix, nprops, sv_info=None):
    """
    Generate the complex-step dP/dF tangent engine.

    For stateless materials: same as before.
    For materials with state:
      - State_old read from DOUBLE PRECISION args
      - Converted to DOUBLE COMPLEX with zero imaginary part
      - Reset before each perturbation (lesson 4.2)
      - State_new outputs from CS calls are discarded
    """
    if sv_info is None:
        sv_info = OrderedDict()
    has_state = bool(sv_info)

    # Build argument list
    args = ['F']
    if has_state:
        for name in sv_info:
            args.append(f'{name}_old')
        args.append('dt')
    args.extend(['props', 'dPdF'])

    lines = []
    sig = f'      SUBROUTINE {mat_prefix}_cs_dPdF('
    sig += ', '.join(args) + ')'
    lines.append(sig)
    lines.append('      IMPLICIT NONE')
    lines.append(
        '      DOUBLE PRECISION, INTENT(IN)  :: F(3,3)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE PRECISION, INTENT(IN)  :: '
                f'{name}_old{dims}')
        lines.append(
            '      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.append(
        f'      DOUBLE PRECISION, INTENT(IN)  :: props({nprops})')
    lines.append(
        '      DOUBLE PRECISION, INTENT(OUT) :: dPdF(3,3,3,3)')
    lines.append('')
    lines.append(
        '      DOUBLE PRECISION, PARAMETER :: CS_H = 1.0d-10')
    lines.append('      DOUBLE COMPLEX :: Fz(3,3), Pz(3,3)')

    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
        lines.append('      DOUBLE PRECISION :: dt_safe')

    lines.append('      INTEGER :: k, l, i, j')
    lines.append('')

    if has_state:
        lines.append(
            'C     --- Protect against DTIME = 0 ---')
        lines.append('      dt_safe = dt')
        lines.append(
            '      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14')
        lines.append('')

    lines.append('      DO k = 1, 3')
    lines.append('        DO l = 1, 3')
    lines.append('          CALL real2complex33(F, Fz)')

    if has_state:
        # Reset state_old to real before each perturbation (lesson 4.2)
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'          {name}_old_z = '
                    f'DCMPLX({name}_old, 0.0d0)')
            elif shape == 'vector':
                lines.append(
                    f'          CALL real2complex3('
                    f'{name}_old, {name}_old_z)')
            elif shape == 'tensor':
                lines.append(
                    f'          CALL real2complex33('
                    f'{name}_old, {name}_old_z)')

    lines.append(
        '          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)')

    # Call material
    if has_state:
        call_args = ['Fz']
        for name in sv_info:
            call_args.append(f'{name}_old_z')
        call_args.append('dt_safe')
        call_args.append('props')
        call_args.append('Pz')
        for name in sv_info:
            call_args.append(f'{name}_new_z')
        lines.append(_fortran_call(
            '          ', f'{mat_prefix}_stress_PK1', call_args))
    else:
        lines.append(
            f'          CALL {mat_prefix}_stress_PK1'
            f'(Fz, props, Pz)')

    lines.append('          DO i = 1, 3')
    lines.append('            DO j = 1, 3')
    lines.append(
        '              dPdF(i,j,k,l) = AIMAG(Pz(i,j)) / CS_H')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')
    lines.append('      RETURN')
    lines.append(
        f'      END SUBROUTINE {mat_prefix}_cs_dPdF')

    return '\n'.join(lines)


# ===================================================================
# PK1 -> Cauchy + Jaumann (unchanged)
# ===================================================================

def _generate_pk1_to_cauchy_jaumann():
    """
    Generate PK1 -> Cauchy stress + Abaqus Jaumann DDSDDE conversion.

    For Abaqus/Standard solid continuum UMATs, DDSDDE is the spatial
    tangent associated with the Jaumann rate of Kirchhoff stress
    (tau = J * sigma):

        tau^J_ij = J * DDSDDE_ijkl * D_kl

    The chain rule on P_iJ = F_iI * S_IJ produces the identity

        AFF_ijkl / J = c_Truesdell_ijkl + delta_ik * sigma_jl    (*)

    so AFF/J is NOT the Truesdell tangent; it differs by a geometric
    delta-sigma term that is present for every constitutive law
    generated through this push-forward path. The "-delta_ik*sigma_jl"
    in the formula below is the bookkeeping correction required by (*).
    Without it, DDSDDE is wrong by a stress-linear term that is invisible
    to F=I, stress-free-rotation, and uniaxial-bar tests but is ~40%
    relative on finite multi-axial deformation.

    See the Jaumann-rate resolution notes for the derivation and a 50-line
    NumPy reference that catches the entire class of "tangent looks
    right at F=I but is silently wrong at finite stress" bugs.

    References:
        Bonet, Gil & Wood (2016), Nonlinear Solid Mechanics for FEA.
        Tensor Toolbox (adtzlr/ttb): C4 = piola(F, C^e)/J + (S cdya I)
            + (I cdya S).
        Nguyen & Waas (2016), Z. Angew. Math. Phys. 67:35.
    """
    return """
C----------------------------------------------------------------------
C     pk1_to_cauchy_jaumann: Convert PK1 stress and tangent to
C     Cauchy stress (Voigt) and Jaumann DDSDDE (Voigt 6x6)
C
C     For Abaqus/Standard solid continuum UMATs, DDSDDE is the spatial
C     tangent in tau^J_ij = J * DDSDDE_ijkl * D_kl, where tau = J*sigma
C     and ^J denotes the Jaumann rate of Kirchhoff stress.
C
C     DDSDDE_ijkl = (1/J)*AFF_ijkl
C                  - delta_ik*sigma_jl                 (bookkeeping)
C                  + 0.5*(sigma_ik*delta_jl + sigma_jk*delta_il
C                       + sigma_il*delta_jk + sigma_jl*delta_ik)
C
C     where AFF_ijkl = sum_{J,L} dPdF(i,J,k,L) * F(j,J) * F(l,L).
C
C     The bookkeeping line is non-optional. The chain rule on
C     P_iJ = F_iI * S_IJ gives AFF/J = c_Truesdell + delta_ik*sigma_jl;
C     subtracting that delta-sigma piece recovers c_Truesdell, to which
C     the symmetric (sigma cdya I + I cdya sigma) block is then added.
C
C     Abaqus Voigt order: 11, 22, 33, 12, 13, 23
C----------------------------------------------------------------------
      SUBROUTINE pk1_to_cauchy_jaumann(F, P, dPdF, STRESS, DDSDDE)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: F(3,3), P(3,3), dPdF(3,3,3,3)
      DOUBLE PRECISION, INTENT(OUT) :: STRESS(6), DDSDDE(6,6)

            DOUBLE PRECISION :: Jdet, det33d, JdetInv
      DOUBLE PRECISION :: sigma(3,3)
      DOUBLE PRECISION :: c(3,3,3,3)
      DOUBLE PRECISION :: AFF
      DOUBLE PRECISION :: delta(3,3)
      INTEGER :: i, j, k, l, J1, L1

C     Voigt index map: (1,1)->1, (2,2)->2, (3,3)->3,
C                      (1,2)->4, (1,3)->5, (2,3)->6
      INTEGER :: voigt_i(6), voigt_j(6)
      DATA voigt_i /1, 2, 3, 1, 1, 2/
      DATA voigt_j /1, 2, 3, 2, 3, 3/

C     Determinant
      Jdet = det33d(F)
      JdetInv = 1.0d0 / Jdet

C     Cauchy stress: sigma = (1/J) P F^T
      DO i = 1, 3
        DO j = 1, 3
          sigma(i,j) = 0.0d0
          DO J1 = 1, 3
            sigma(i,j) = sigma(i,j) + P(i,J1) * F(j,J1)
          END DO
                    sigma(i,j) = sigma(i,j) * JdetInv
        END DO
      END DO

C     DDSDDE: (1/J)*AFF
C                  - delta_ik*sigma_jl
C                  + 0.5*(sigma_ik*delta_jl + sigma_jk*delta_il
C                       + sigma_il*delta_jk + sigma_jl*delta_ik)
C     where AFF_ijkl = sum_{J,L} dPdF(i,J,k,L) * F(j,J) * F(l,L)
C---------------------------------------------------------------------
      DO i = 1, 3
        DO j = 1, 3
          IF (i .EQ. j) THEN
            delta(i,j) = 1.0d0
          ELSE
            delta(i,j) = 0.0d0
          END IF
        END DO
      END DO

      DO i = 1, 3
        DO j = 1, 3
          DO k = 1, 3
            DO l = 1, 3
C             AFF_ijkl = sum_{J,L} dPdF(i,J,k,L) * F(j,J) * F(l,L)
              AFF = 0.0d0
              DO J1 = 1, 3
                DO L1 = 1, 3
                  AFF = AFF + dPdF(i,J1,k,L1)*F(j,J1)*F(l,L1)
                END DO
              END DO
                            c(i,j,k,l) = AFF * JdetInv
     &          - delta(i,k)*sigma(j,l)
     &          + 0.5d0*(sigma(i,k)*delta(j,l) + sigma(j,k)*delta(i,l)
     &                 + sigma(i,l)*delta(j,k) + sigma(j,l)*delta(i,k))
            END DO
          END DO
        END DO
      END DO

C     Pack symmetric tensor into Voigt 6x6.
C     Major symmetry is satisfied by construction.
      DO i = 1, 6
        STRESS(i) = sigma(voigt_i(i), voigt_j(i))
      END DO

      DO i = 1, 6
        DO j = i, 6
          DDSDDE(i,j) = c(voigt_i(i), voigt_j(i), voigt_i(j), voigt_j(j))
          DDSDDE(j,i) = DDSDDE(i,j)
        END DO
      END DO

      RETURN
      END SUBROUTINE pk1_to_cauchy_jaumann
"""


# ===================================================================
# UMAT wrapper generation
# ===================================================================

def _generate_umat_wrapper(mat_prefix, nprops, sv_info=None):
    """
    Generate the standard Abaqus UMAT interface subroutine.

    For stateless materials: identical to original.
    For materials with state: adds STATEV read/write and dt_safe,
    following the standard Abaqus STATEV initialization pattern (TIME(2)==0).
    """
    if sv_info is None:
        sv_info = OrderedDict()
    has_state = bool(sv_info)
    nstatv = _nstate_per_gp(sv_info)

    lines = []
    lines.append('C' + '=' * 70)
    lines.append(
        'C     UMAT: Standard Abaqus user material interface')
    lines.append('C')
    lines.append(
        f'C     Generated by abaqus_ufl for material: {mat_prefix}')
    lines.append(f'C     NPROPS = {nprops}')
    if has_state:
        lines.append(f'C     NSTATV = {nstatv}')
        lines.append('C')
        lines.append('C     State variables:')
        offset = 0
        for name, entry in sv_info.items():
            size = entry['size']
            if size == 1:
                lines.append(
                    f'C       STATEV({offset+1:3d})'
                    f'       = {name}')
            else:
                lines.append(
                    f'C       STATEV({offset+1:3d}:'
                    f'{offset+size:3d})  = {name} '
                    f'({entry["shape"]}, {size} components)')
            offset += size
    lines.append('C' + '=' * 70)

    # UMAT signature (standard Abaqus — never changes)
    lines.append(
        '      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, '
        'SPD, SCD,')
    lines.append(
        '     &  RPL, DDSDDT, DRPLDE, DRPLDT,')
    lines.append(
        '     &  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, '
        'PREDEF, DPRED,')
    lines.append(
        '     &  CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, '
        'NPROPS,')
    lines.append(
        '     &  COORDS, DROT, PNEWDT, CELENT, DFGRD0, '
        'DFGRD1,')
    lines.append(
        '     &  NOEL, NPT, LAYER, KSPT, JSTEP, KINC)')
    lines.append('')
    lines.append('      IMPLICIT NONE')
    lines.append('')
    lines.append('C     --- Abaqus interface variables ---')
    lines.append('      CHARACTER*80 CMNAME')
    lines.append(
        '      INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS')
    lines.append(
        '      INTEGER NOEL, NPT, LAYER, KSPT, JSTEP(4), KINC')
    lines.append('')
    lines.append(
        '      DOUBLE PRECISION STRESS(NTENS), STATEV(NSTATV)')
    lines.append(
        '      DOUBLE PRECISION DDSDDE(NTENS,NTENS)')
    lines.append(
        '      DOUBLE PRECISION SSE, SPD, SCD, RPL, DTIME, '
        'TEMP, DTEMP')
    lines.append(
        '      DOUBLE PRECISION DDSDDT(NTENS), DRPLDE(NTENS), '
        'DRPLDT')
    lines.append(
        '      DOUBLE PRECISION STRAN(NTENS), DSTRAN(NTENS)')
    lines.append(
        '      DOUBLE PRECISION TIME(2), PREDEF(1), DPRED(1)')
    lines.append(
        '      DOUBLE PRECISION PROPS(NPROPS), COORDS(3)')
    lines.append(
        '      DOUBLE PRECISION DROT(3,3), PNEWDT, CELENT')
    lines.append(
        '      DOUBLE PRECISION DFGRD0(3,3), DFGRD1(3,3)')
    lines.append('')

    # --- Local variable declarations ---
    lines.append('C     --- Local variables ---')
    lines.append(
        '      DOUBLE PRECISION :: F(3,3), P_real(3,3)')
    lines.append(
        '      DOUBLE PRECISION :: dPdF(3,3,3,3)')
    lines.append(
        '      DOUBLE COMPLEX :: Fz(3,3), Pz(3,3)')

    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE PRECISION :: {name}_old{dims}')
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
        lines.append('      DOUBLE PRECISION :: dt_safe')

    lines.append('      INTEGER :: i, j')
    lines.append('')

    # --- Get deformation gradient ---
    lines.append('C     --- Get deformation gradient ---')
    lines.append('      DO i = 1, 3')
    lines.append('        DO j = 1, 3')
    lines.append('          F(i,j) = DFGRD1(i,j)')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')

    if has_state:
        # --- STATEV read (standard Abaqus TIME(2)==0 pattern) ---
        lines.append(_generate_statev_read(sv_info))
        lines.append('')

        # --- dt_safe ---
        lines.append('C     --- Protect against DTIME = 0 ---')
        lines.append('      dt_safe = DTIME')
        lines.append(
            '      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14')
        lines.append('')

        # --- Real evaluation: stress + state update ---
        lines.append(
            'C     --- Compute PK1 stress (real evaluation '
            '+ state update) ---')
        lines.append('      CALL real2complex33(F, Fz)')
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'      {name}_old_z = '
                    f'DCMPLX({name}_old, 0.0d0)')
            elif shape == 'vector':
                lines.append(
                    f'      CALL real2complex3({name}_old, '
                    f'{name}_old_z)')
            elif shape == 'tensor':
                lines.append(
                    f'      CALL real2complex33({name}_old, '
                    f'{name}_old_z)')

        call_args = ['Fz']
        for name in sv_info:
            call_args.append(f'{name}_old_z')
        call_args.append('dt_safe')
        call_args.append('PROPS')
        call_args.append('Pz')
        for name in sv_info:
            call_args.append(f'{name}_new_z')
        lines.append(_fortran_call(
            '      ', f'{mat_prefix}_stress_PK1', call_args))

        lines.append('      DO i = 1, 3')
        lines.append('        DO j = 1, 3')
        lines.append('          P_real(i,j) = DBLE(Pz(i,j))')
        lines.append('        END DO')
        lines.append('      END DO')
        lines.append('')

        # --- STATEV write ---
        lines.append(_generate_statev_write(sv_info))
        lines.append('')

        # --- CS tangent ---
        lines.append(
            'C     --- Compute dP/dF via complex step ---')
        cs_args = ['F']
        for name in sv_info:
            cs_args.append(f'{name}_old')
        cs_args.extend(['dt_safe', 'PROPS', 'dPdF'])
        lines.append(_fortran_call(
            '      ', f'{mat_prefix}_cs_dPdF', cs_args))

    else:
        # --- Stateless path (identical to original) ---
        lines.append(
            'C     --- Compute PK1 stress (real evaluation) ---')
        lines.append('      CALL real2complex33(F, Fz)')
        lines.append(
            f'      CALL {mat_prefix}_stress_PK1'
            f'(Fz, PROPS, Pz)')
        lines.append('      DO i = 1, 3')
        lines.append('        DO j = 1, 3')
        lines.append(
            '          P_real(i,j) = DBLE(Pz(i,j))')
        lines.append('        END DO')
        lines.append('      END DO')
        lines.append('')
        lines.append(
            'C     --- Compute dP/dF via complex step ---')
        lines.append(
            f'      CALL {mat_prefix}_cs_dPdF'
            f'(F, PROPS, dPdF)')

    lines.append('')
    lines.append(
        'C     --- Convert PK1 -> Cauchy + Jaumann tangent '
        '(Voigt) ---')
    lines.append(
        '      CALL pk1_to_cauchy_jaumann(F, P_real, dPdF, '
        'STRESS, DDSDDE)')
    lines.append('')
    lines.append('      RETURN')
    lines.append('      END SUBROUTINE UMAT')

    return '\n'.join(lines)


def _voigt_tensor_helpers():
    return """
C----------------------------------------------------------------------
C     Small-strain Voigt/tensor helpers.
C     Abaqus inputs are tensile-positive. Soil update tensors are
C     compression-positive. Engineering shear strains are halved in tensor
C     form; shear stresses are direct tensor components.
C----------------------------------------------------------------------
      SUBROUTINE abq_voigt_to_tensor_cp(v, A, is_strain)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: v(6)
      DOUBLE PRECISION, INTENT(OUT) :: A(3,3)
      LOGICAL, INTENT(IN) :: is_strain
      DOUBLE PRECISION :: sf

      A = 0.0d0
      sf = 1.0d0
      IF (is_strain) sf = 0.5d0
      A(1,1) = -v(1)
      A(2,2) = -v(2)
      A(3,3) = -v(3)
      A(1,2) = -sf*v(4)
      A(2,1) = A(1,2)
      A(1,3) = -sf*v(5)
      A(3,1) = A(1,3)
      A(2,3) = -sf*v(6)
      A(3,2) = A(2,3)
      RETURN
      END SUBROUTINE abq_voigt_to_tensor_cp

      SUBROUTINE tensor_cp_to_abq_voigt(A, v)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: A(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: v(6)
      v(1) = -A(1,1)
      v(2) = -A(2,2)
      v(3) = -A(3,3)
      v(4) = -A(1,2)
      v(5) = -A(1,3)
      v(6) = -A(2,3)
      RETURN
      END SUBROUTINE tensor_cp_to_abq_voigt
"""


def _generate_small_strain_umat_wrapper(
        mat_prefix, nprops, sv_info=None, predef_names=None):
    """Generate Abaqus UMAT wrapper for compression-positive small strain."""
    if sv_info is None:
        sv_info = OrderedDict()
    predef_names = list(predef_names or [])
    has_state = bool(sv_info)
    nstatv = _nstate_per_gp(sv_info)

    lines = []
    lines.append('C' + '=' * 70)
    lines.append('C     UMAT: small-strain compression-positive material')
    lines.append(f'C     Generated by abaqus_ufl for material: {mat_prefix}')
    lines.append(f'C     NPROPS = {nprops}')
    if has_state:
        lines.append(f'C     NSTATV = {nstatv}')
    lines.append('C' + '=' * 70)
    lines.extend([
        '      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,',
        '     &  RPL, DDSDDT, DRPLDE, DRPLDT,',
        '     &  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,',
        '     &  CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,',
        '     &  COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,',
        '     &  NOEL, NPT, LAYER, KSPT, JSTEP, KINC)',
        '',
        '      IMPLICIT NONE',
        '',
        '      CHARACTER*80 CMNAME',
        '      INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS',
        '      INTEGER NOEL, NPT, LAYER, KSPT, JSTEP(4), KINC',
        '      DOUBLE PRECISION STRESS(NTENS), STATEV(NSTATV)',
        '      DOUBLE PRECISION DDSDDE(NTENS,NTENS)',
        '      DOUBLE PRECISION SSE, SPD, SCD, RPL, DTIME, TEMP, DTEMP',
        '      DOUBLE PRECISION DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT',
        '      DOUBLE PRECISION STRAN(NTENS), DSTRAN(NTENS)',
        '      DOUBLE PRECISION TIME(2), PREDEF(*), DPRED(*)',
        '      DOUBLE PRECISION PROPS(NPROPS), COORDS(3)',
        '      DOUBLE PRECISION DROT(3,3), PNEWDT, CELENT',
        '      DOUBLE PRECISION DFGRD0(3,3), DFGRD1(3,3)',
        '',
        '      DOUBLE PRECISION :: sigma_old(3,3), strain_old(3,3)',
        '      DOUBLE PRECISION :: dstrain(3,3), sigma_new(3,3)',
        '      DOUBLE COMPLEX :: sigma_old_z(3,3), strain_old_z(3,3)',
        '      DOUBLE COMPLEX :: dstrain_z(3,3), sigma_new_z(3,3)',
        '      DOUBLE COMPLEX :: cutback_z',
        '      DOUBLE PRECISION, PARAMETER :: CS_H = 1.0d-10',
    ])
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(f'      DOUBLE PRECISION :: {name}_old{dims}')
    for name in predef_names:
        lines.append(f'      DOUBLE PRECISION :: {name}_real')
        lines.append(f'      DOUBLE COMPLEX :: {name}_z')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
            lines.append(f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
        lines.append('      DOUBLE PRECISION :: dt_safe')
    lines.append('      INTEGER :: i, j, col')
    lines.append('')

    lines.extend([
        '      IF (NTENS .NE. 6) THEN',
        '        PNEWDT = 0.5d0',
        '        RETURN',
        '      END IF',
        '',
        'C     --- Convert Abaqus tension-positive Voigt to soil tensors ---',
        '      CALL abq_voigt_to_tensor_cp(STRESS, sigma_old, .FALSE.)',
        '      CALL abq_voigt_to_tensor_cp(STRAN, strain_old, .TRUE.)',
        '      CALL abq_voigt_to_tensor_cp(DSTRAN, dstrain, .TRUE.)',
        '',
    ])

    if has_state:
        lines.append(_generate_statev_read(sv_info))
        lines.append('')
        lines.extend([
            '      dt_safe = DTIME',
            '      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14',
            '',
        ])
    for i, name in enumerate(predef_names, start=1):
        lines.append(f'      {name}_real = PREDEF({i}) + DPRED({i})')
        lines.append(f'      {name}_z = DCMPLX({name}_real, 0.0d0)')
    if predef_names:
        lines.append('')

    lines.extend([
        'C     --- Real stress update ---',
        '      CALL real2complex33(sigma_old, sigma_old_z)',
        '      CALL real2complex33(strain_old, strain_old_z)',
        '      CALL real2complex33(dstrain, dstrain_z)',
    ])
    if has_state:
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'      {name}_old_z = DCMPLX({name}_old, 0.0d0)')
            elif shape == 'vector':
                lines.append(
                    f'      CALL real2complex3({name}_old, {name}_old_z)')
            elif shape == 'tensor':
                lines.append(
                    f'      CALL real2complex33({name}_old, {name}_old_z)')

    call_args = ['sigma_old_z', 'strain_old_z', 'dstrain_z']
    if has_state:
        for name in sv_info:
            call_args.append(f'{name}_old_z')
    for name in predef_names:
        call_args.append(f'{name}_z')
    if has_state:
        call_args.append('dt_safe')
    call_args.extend(['cutback_z', 'PROPS', 'sigma_new_z'])
    if has_state:
        for name in sv_info:
            call_args.append(f'{name}_new_z')
    lines.append(_fortran_call('      ', f'{mat_prefix}_stress_update',
                               call_args))
    lines.extend([
        '      IF (DBLE(cutback_z) .GT. 0.5d0) THEN',
        '        PNEWDT = 0.25d0',
        '        RETURN',
        '      END IF',
    ])
    lines.extend([
        '      DO i = 1, 3',
        '        DO j = 1, 3',
        '          sigma_new(i,j) = DBLE(sigma_new_z(i,j))',
        '        END DO',
        '      END DO',
        '      CALL tensor_cp_to_abq_voigt(sigma_new, STRESS)',
        '',
    ])

    if has_state:
        lines.append(_generate_statev_write(sv_info))
        lines.append('')

    lines.extend([
        'C     --- Complex-step DDSDDE wrt Abaqus DSTRAN ---',
        '      DO col = 1, 6',
        '        CALL real2complex33(sigma_old, sigma_old_z)',
        '        CALL real2complex33(strain_old, strain_old_z)',
        '        CALL real2complex33(dstrain, dstrain_z)',
    ])
    if has_state:
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'        {name}_old_z = DCMPLX({name}_old, 0.0d0)')
            elif shape == 'vector':
                lines.append(
                    f'        CALL real2complex3({name}_old, {name}_old_z)')
            elif shape == 'tensor':
                lines.append(
                    f'        CALL real2complex33({name}_old, {name}_old_z)')
    for name in predef_names:
        lines.append(f'        {name}_z = DCMPLX({name}_real, 0.0d0)')
    lines.extend([
        '        IF (col .EQ. 1) dstrain_z(1,1) = dstrain_z(1,1)',
        '     &    - DCMPLX(0.0d0, CS_H)',
        '        IF (col .EQ. 2) dstrain_z(2,2) = dstrain_z(2,2)',
        '     &    - DCMPLX(0.0d0, CS_H)',
        '        IF (col .EQ. 3) dstrain_z(3,3) = dstrain_z(3,3)',
        '     &    - DCMPLX(0.0d0, CS_H)',
        '        IF (col .EQ. 4) THEN',
        '          dstrain_z(1,2) = dstrain_z(1,2) - DCMPLX(0.0d0,',
        '     &      0.5d0*CS_H)',
        '          dstrain_z(2,1) = dstrain_z(1,2)',
        '        END IF',
        '        IF (col .EQ. 5) THEN',
        '          dstrain_z(1,3) = dstrain_z(1,3) - DCMPLX(0.0d0,',
        '     &      0.5d0*CS_H)',
        '          dstrain_z(3,1) = dstrain_z(1,3)',
        '        END IF',
        '        IF (col .EQ. 6) THEN',
        '          dstrain_z(2,3) = dstrain_z(2,3) - DCMPLX(0.0d0,',
        '     &      0.5d0*CS_H)',
        '          dstrain_z(3,2) = dstrain_z(2,3)',
        '        END IF',
    ])
    lines.append(_fortran_call('        ', f'{mat_prefix}_stress_update',
                               call_args))
    lines.extend([
        '        DDSDDE(1,col) = -AIMAG(sigma_new_z(1,1)) / CS_H',
        '        DDSDDE(2,col) = -AIMAG(sigma_new_z(2,2)) / CS_H',
        '        DDSDDE(3,col) = -AIMAG(sigma_new_z(3,3)) / CS_H',
        '        DDSDDE(4,col) = -AIMAG(sigma_new_z(1,2)) / CS_H',
        '        DDSDDE(5,col) = -AIMAG(sigma_new_z(1,3)) / CS_H',
        '        DDSDDE(6,col) = -AIMAG(sigma_new_z(2,3)) / CS_H',
        '      END DO',
        '',
        '      RETURN',
        '      END SUBROUTINE UMAT',
    ])
    return '\n'.join(lines)


# ===================================================================
# Top-level generator
# ===================================================================

def generate_umat(material, output_path, mat_prefix=None,
                  matrix_backend='iterative'):
    """
    Generate a complete .for file for Abaqus UMAT.

    Args:
        material: a Material instance with stress_PK1 defined
        output_path: path to write the .for file
        mat_prefix: prefix for generated subroutine names
                    (default: lowercase class name)
        matrix_backend: matrix-function backend for generated Fortran.
                    'iterative' (default) uses fixed-count iterative helpers
                        (CS-safe, no eigendecomposition guards).
                    'eig' uses sqrtm33z/logm33z/expm33z.
                        Not CS-safe at diagonal-F states; available for
                        backward compatibility or explicit choice.
    """
    if matrix_backend not in ('eig', 'iterative'):
        raise ValueError(
            "matrix_backend must be 'eig' or 'iterative'")
    if 'stress_PK1' not in material._methods:
        raise ValueError(
            "Material must define stress_PK1 for UMAT generation")

    # Check that stress_PK1 only depends on F (UMAT limitation)
    info = material._methods['stress_PK1']
    field_params = info['field_params']
    if field_params != ['F']:
        raise ValueError(
            f"UMAT stress_PK1 must only depend on F, "
            f"but depends on {field_params}. "
            f"Use UEL target for multi-field problems."
        )

    if mat_prefix is None:
        mat_prefix = material.__class__.__name__.lower()

    nprops = len(material.props_names)
    sv_info = _state_var_info(material)
    has_state = bool(sv_info)

    # Read tensor_ops.for template
    import os
    template_dir = os.path.join(
        os.path.dirname(__file__), 'templates')
    tensor_ops_path = os.path.join(
        template_dir, 'tensor_ops.for')
    if not os.path.exists(tensor_ops_path):
        raise FileNotFoundError('Packaged Fortran template not found: ' + tensor_ops_path)
    with open(tensor_ops_path, 'r') as f:
        tensor_ops_src = f.read()

    # Generate sections
    stress_src = _generate_stress_subroutine(
        material, mat_prefix, matrix_backend=matrix_backend)
    cs_src = _generate_cs_dPdF(mat_prefix, nprops, sv_info)
    convert_src = _generate_pk1_to_cauchy_jaumann()
    umat_src = _generate_umat_wrapper(
        mat_prefix, nprops, sv_info)

    # Header
    nstatv = _nstate_per_gp(sv_info)
    hdr = []
    hdr.append(f'C{"=" * 70}')
    hdr.append(f'C     {os.path.basename(output_path)}')
    hdr.append('C')
    hdr.append('C     Generated by abaqus_ufl v0.1.0')
    hdr.append(
        f'C     Material: {material.__class__.__name__}')
    hdr.append(f'C     Matrix backend: {matrix_backend}')
    hdr.append(
        f'C     Props: '
        f'{", ".join(f"{k}={getattr(material, k)}" for k in material.props_names)}')
    hdr.append(f'C     NPROPS = {nprops}')
    if has_state:
        hdr.append(f'C     NSTATV = {nstatv}')
        hdr.append(
            f'C     State vars: '
            f'{", ".join(sv_info.keys())}')
    hdr.append('C')
    hdr.append(
        'C     This file is self-contained. '
        'Drop it into your Abaqus job:')
    hdr.append(
        f'C       abaqus job=my_job '
        f'user={os.path.basename(output_path)}')
    hdr.append(f'C{"=" * 70}')
    header = '\n'.join(hdr)

    # Assemble
    full_src = '\n'.join([
        header,
        '',
        'C' + '=' * 70,
        'C     Section 1: Abaqus UMAT interface',
        'C' + '=' * 70,
        umat_src,
        '',
        'C' + '=' * 70,
        f'C     Section 2: Material stress: '
        f'{mat_prefix}_stress_PK1',
        'C' + '=' * 70,
        stress_src,
        '',
        'C' + '=' * 70,
        'C     Section 3: Complex-step tangent: d(PK1)/d(F)',
        'C' + '=' * 70,
        cs_src,
        '',
        'C' + '=' * 70,
        'C     Section 4: PK1 -> Cauchy stress + '
        'Jaumann tangent (Voigt)',
        'C' + '=' * 70,
        convert_src,
        '',
        'C' + '=' * 70,
        'C     Section 5: Tensor operations '
        '(det33, inv33, matmul33, etc.)',
        'C' + '=' * 70,
        tensor_ops_src,
    ])

    # Apply column-72 wrapping to all generated sections
    # (tensor_ops template is already compliant)
    full_src = _wrap_fortran_source(full_src)

    with open(output_path, 'w') as f:
        f.write(full_src)

    subs = [f'{mat_prefix}_stress_PK1',
            f'{mat_prefix}_cs_dPdF',
            'pk1_to_cauchy_jaumann', 'UMAT']
    print(f"  Generated: {output_path}")
    print(f"  Material:  {material.__class__.__name__}")
    props_str = ', '.join(material.props_names)
    print(f"  Props:     {nprops} ({props_str})")
    if has_state:
        sv_str = ', '.join(
            f'{n}({v["size"]})' for n, v in sv_info.items())
        print(f"  State:     {nstatv} ({sv_str})")
    print(f"  Subroutines: {', '.join(subs)}")


def generate_small_strain_umat(material, output_path, mat_prefix=None,
                               update_method='stress_update',
                               extra_fortran_files=None):
    """Generate a small-strain Abaqus UMAT from a stress update method.

    By default the generator translates ``stress_update``.  Advanced models can
    pass ``update_method="stress_update_codegen"`` to keep a readable Python
    reference update while translating a separate generator-friendly method.

    Parameters
    ----------
    extra_fortran_files : list[str | os.PathLike] | None
        Additional Fortran source files to append verbatim to the generated
        ``.for`` output.  This is useful for hand-written helper subroutines
        referenced via ``@fortran_helper``.
    """
    predef_names = _predef_var_names(material)
    if predef_names:
        param_args = set(PARAM_ARGS) | set(predef_names)
        material.refresh_metadata(param_args=param_args)

    info = _small_strain_update_method_info(material, update_method)
    expected = ['sigma_old', 'strain_old', 'dstrain']
    field_params = info['field_params']
    if field_params != expected:
        raise ValueError(
            f"{update_method} must use arguments "
            "(sigma_old, strain_old, dstrain, state_old..., dt); "
            f"recognized field arguments were {field_params}")

    if getattr(material, 'stress_convention', None) != \
            'compression_positive':
        raise ValueError(
            "Small-strain UMAT requires stress_convention="
            "'compression_positive'")

    if mat_prefix is None:
        mat_prefix = material.__class__.__name__.lower()

    nprops = len(material.props_names)
    sv_info = _state_var_info(material)
    has_state = bool(sv_info)

    if extra_fortran_files is None:
        extra_fortran_files = []
    extra_fortran_files = [str(p) for p in extra_fortran_files]

    # --- Validate @fortran_helper sidecars -----------------------------------
    # Scan material for any @fortran_helper methods and verify that each
    # required subroutine name appears in at least one sidecar file.
    fh_methods = []
    for name in dir(material):
        if name.startswith('__'):
            continue
        attr = getattr(material, name, None)
        if callable(attr) and hasattr(attr, '__fortran_helper__'):
            fh_methods.append((name, attr.__fortran_helper__))

    if fh_methods:
        sidecar_text = ''
        for path in extra_fortran_files:
            if os.path.exists(path):
                with open(path, 'r') as f:
                    sidecar_text += f.read() + '\n'
        missing = []
        for method_name, meta in fh_methods:
            sub_name = meta.get('subroutine_override')
            if sub_name is None:
                sub_name = f'{mat_prefix}_{method_name.lstrip("_")}'
            # Simple case-insensitive search for "SUBROUTINE <name>"
            pat = f'SUBROUTINE {sub_name}'
            if pat.upper() not in sidecar_text.upper():
                missing.append(sub_name)
        if missing:
            raise VerificationError(
                f"Missing Fortran subroutine(s) for @fortran_helper methods: "
                f"{', '.join(missing)}.  Ensure each is defined in one of the "
                f"extra_fortran_files: {extra_fortran_files}")

    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    tensor_ops_path = os.path.join(template_dir, 'tensor_ops.for')
    if not os.path.exists(tensor_ops_path):
        raise FileNotFoundError('Packaged Fortran template not found: ' + tensor_ops_path)
    with open(tensor_ops_path, 'r') as f:
        tensor_ops_src = f.read()

    cs_linalg_path = os.path.join(template_dir, 'cs_linalg.for')
    if os.path.exists(cs_linalg_path):
        with open(cs_linalg_path, 'r') as f:
            cs_linalg_src = f.read()
    else:
        cs_linalg_src = ''

    update_src = _generate_small_strain_update_subroutine(
        material, mat_prefix, predef_names, update_method)
    umat_src = _generate_small_strain_umat_wrapper(
        mat_prefix, nprops, sv_info, predef_names)
    helper_src = _voigt_tensor_helpers()

    nstatv = _nstate_per_gp(sv_info)
    hdr = []
    hdr.append(f'C{"=" * 70}')
    hdr.append(f'C     {os.path.basename(output_path)}')
    hdr.append('C')
    hdr.append('C     Generated by abaqus_ufl v0.1.0')
    hdr.append('C     Small-strain UMAT, compression-positive user API')
    hdr.append(f'C     Material: {material.__class__.__name__}')
    hdr.append(f'C     Update method: {update_method}')
    hdr.append(
        f'C     Props: '
        f'{", ".join(f"{k}={getattr(material, k)}" for k in material.props_names)}')
    hdr.append(f'C     NPROPS = {nprops}')
    if has_state:
        hdr.append(f'C     NSTATV = {nstatv}')
        hdr.append(f'C     State vars: {", ".join(sv_info.keys())}')
    hdr.append('C')
    hdr.append(f'C       abaqus job=my_job user={os.path.basename(output_path)}')
    hdr.append(f'C{"=" * 70}')
    header = '\n'.join(hdr)

    sections = [
        header,
        '',
        'C' + '=' * 70,
        'C     Section 1: Abaqus UMAT interface',
        'C' + '=' * 70,
        umat_src,
        '',
        'C' + '=' * 70,
        f'C     Section 2: Material stress update: '
        f'{mat_prefix}_stress_update',
        'C' + '=' * 70,
        update_src,
        '',
        'C' + '=' * 70,
        'C     Section 3: Small-strain Voigt/tensor helpers',
        'C' + '=' * 70,
        helper_src,
        '',
        'C' + '=' * 70,
        'C     Section 4: Tensor operations',
        'C' + '=' * 70,
        tensor_ops_src,
        '',
        'C' + '=' * 70,
        'C     Section 5: Small dense linear solve',
        'C' + '=' * 70,
        cs_linalg_src,
    ]

    for idx, path in enumerate(extra_fortran_files, start=6):
        if os.path.exists(path):
            with open(path, 'r') as f:
                sidecar_src = f.read()
            basename = os.path.basename(path)
            sections.extend([
                '',
                'C' + '=' * 70,
                f'C     Section {idx}: User Fortran helper — {basename}',
                'C' + '=' * 70,
                sidecar_src,
            ])

    full_src = '\n'.join(sections)
    full_src = _wrap_fortran_source(full_src)

    with open(output_path, 'w') as f:
        f.write(full_src)

    subs = [f'{mat_prefix}_stress_update', 'UMAT']
    print(f"  Generated: {output_path}")
    print(f"  Material:  {material.__class__.__name__}")
    print(f"  Kinematics: small_strain")
    print(f"  Convention: compression-positive user stress/strain")
    props_str = ', '.join(material.props_names)
    print(f"  Props:     {nprops} ({props_str})")
    if has_state:
        sv_str = ', '.join(
            f'{n}({v["size"]})' for n, v in sv_info.items())
        print(f"  State:     {nstatv} ({sv_str})")
    print(f"  Subroutines: {', '.join(subs)}")

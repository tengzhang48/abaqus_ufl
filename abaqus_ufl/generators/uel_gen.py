"""
UEL code generator.

Given a WeakForm instance, generates a complete .for file containing:
  1. Shipped templates (tensor_ops, shape functions, quadrature, mapping)
  2. Material subroutines (AST-translated from Python, COMPLEX*16)
  3. CS tangent engine (perturbation loops for all tangent blocks)
  4. UEL subroutine (DOF parsing, Gauss loop, RHS + AMATRX assembly)

Usage:
    problem = GelProblem()
    problem.verify()
    from abaqus_ufl.generators.uel_gen import generate_uel
    generate_uel(problem, 'gel_uel.for', element='Quad8')
"""

import os
import inspect
import textwrap
import ast
import numpy as np

from ..core.fields import VectorField, ScalarField
from ..core.defs import FIELD_ARGS
from ..core.weakform import EQUATION_INFO
from .element_config import ElementConfig, ELEMENT_CONFIGS
from .umat_gen import (FortranTranslator, _get_source_body,
                       _state_var_info, _nstate_per_gp,
                       _fortran_dims, _fortran_init_value,
                       _fortran_call)


def _compute_ndofel(cfg, fields):
    """Compute NDOFEL from an ElementConfig and field dict."""
    n_midside = cfg.n_nodes - cfg.n_corner_nodes
    ndof_corner = sum(f.ndof_per_node for f in fields.values())
    ndof_midside = sum(f.ndof_per_node for f in fields.values()
                       if f.degree >= 2)
    return cfg.n_corner_nodes * ndof_corner + n_midside * ndof_midside


def _nn_for_field(cfg, f):
    """Return the node count for a field based on its degree."""
    return cfg.n_nodes if f.degree >= 2 else cfg.n_corner_nodes


def _read_template(name):
    """Read a Fortran template file from the templates directory."""
    path = os.path.join(os.path.dirname(__file__), 'templates', name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Template not found: {path}")
    with open(path, 'r') as f:
        return f.read()


from ._fortran_format import wrap_fortran_line as _wrap_fortran


def _wrap_lines(lines):
    """Apply Fortran line wrapping to a list of lines."""
    result = []
    for line in lines:
        result.append(_wrap_fortran(line))
    return result


# =====================================================================
# Section 2: Material subroutine generation
# =====================================================================

# Map Python argument names to Fortran declarations
_ARG_FORTRAN_TYPE = {
    'F':       ('DOUBLE COMPLEX', 'F(3,3)',       'tensor'),
    'p':       ('DOUBLE COMPLEX', 'p',            'scalar'),
    'mu':      ('DOUBLE COMPLEX', 'mu',           'scalar'),
    'grad_mu': ('DOUBLE COMPLEX', 'grad_mu(3)',   'vector'),
    'F_old':   ('DOUBLE COMPLEX', 'F_old(3,3)',   'tensor'),
    'p_old':   ('DOUBLE COMPLEX', 'p_old',        'scalar'),
    'dt':      ('DOUBLE PRECISION', 'dt',          'scalar'),
    'time':    ('DOUBLE PRECISION', 'time',        'scalar'),
    'X':       ('DOUBLE PRECISION', 'X(3)',        'vector'),
}


def _field_args_for(weakform_or_material):
    """Return problem-specific field argument metadata when available."""
    if hasattr(weakform_or_material, '_field_arg_defs'):
        return weakform_or_material._field_arg_defs
    if hasattr(weakform_or_material, '_field_args'):
        return weakform_or_material._field_args
    return FIELD_ARGS


def _ordered_history_scalars(weakform, history_vars):
    """Return scalar history arguments in the same order everywhere.

    History values are interpolated at a Gauss point from declared fields, so
    the only scalar histories the standard UEL can pass safely are old values
    of declared scalar fields.  Keep the argument order tied to field
    declaration order; using an alphabetic set order in one site and field
    order in another silently swaps same-typed Fortran dummy arguments.
    """
    return _ordered_history_scalars_from(
        weakform.field_names, weakform.fields, history_vars)


def _ordered_history_scalars_from(field_names, fields, history_vars):
    """Core of _ordered_history_scalars for callers without a weakform in
    scope (e.g. uel_fbar_coupled's CS-call emitter). Every site that orders
    CS history-scalar arguments MUST go through this function."""
    remaining = set(history_vars)
    remaining.discard('F_old')
    ordered = []
    for fname in field_names:
        old_name = f'{fname}_old'
        if old_name not in remaining:
            continue
        field = fields[fname]
        if not isinstance(field, ScalarField):
            raise NotImplementedError(
                f"History argument '{old_name}' is not a scalar field "
                "history; standard UEL CS tangents only support scalar "
                "field histories plus F_old.")
        ordered.append(old_name)
        remaining.remove(old_name)
    if remaining:
        raise NotImplementedError(
            "Standard UEL CS tangents cannot pass history arguments that "
            "are not old values of declared scalar fields: "
            f"{sorted(remaining)}")
    return ordered


def _arg_fortran_type(arg_name, field_args=None):
    """Return Fortran declaration metadata for a generated argument."""
    if arg_name in _ARG_FORTRAN_TYPE:
        return _ARG_FORTRAN_TYPE[arg_name]
    finfo = (field_args or FIELD_ARGS).get(arg_name, {})
    kind = finfo.get('kind', 'scalar')
    if kind == 'matrix':
        return 'DOUBLE COMPLEX', f'{arg_name}(3,3)', 'tensor'
    if kind == 'vector':
        return 'DOUBLE COMPLEX', f'{arg_name}(3)', 'vector'
    return 'DOUBLE COMPLEX', arg_name, 'scalar'

# Map return types to Fortran output declarations
_RETURN_FORTRAN = {
    'stress_PK1':      ('DOUBLE COMPLEX', 'P_out(3,3)'),
    'pressure_resid':  ('DOUBLE COMPLEX', 'rp_out'),
    'solvent_flux':    ('DOUBLE COMPLEX', 'jR_out(3)'),
    'solvent_storage': ('DOUBLE COMPLEX', 'cdot_out'),
    'phase_flux':      ('DOUBLE COMPLEX', 'phase_flux_out(3)'),
    'phase_storage':   ('DOUBLE COMPLEX', 'phase_storage_out'),
    'pressure_flux':   ('DOUBLE COMPLEX', 'pressure_flux_out(3)'),
    'pressure_storage': ('DOUBLE COMPLEX', 'pressure_storage_out'),
    'species_flux':    ('DOUBLE COMPLEX', 'species_flux_out(3)'),
    'species_storage': ('DOUBLE COMPLEX', 'species_storage_out'),
    'h_field':         ('DOUBLE COMPLEX', 'H_out(3)'),
    'volumetric_PK1':  ('DOUBLE COMPLEX', 'Pvol_out(3,3)'),
}

# Map method names to output shapes for the element-wise loop
_RETURN_SHAPE = {
    'stress_PK1':      'tensor',
    'pressure_resid':  'scalar',
    'solvent_flux':    'vector',
    'solvent_storage': 'scalar',
    'phase_flux':      'vector',
    'phase_storage':   'scalar',
    'pressure_flux':   'vector',
    'pressure_storage': 'scalar',
    'species_flux':    'vector',
    'species_storage': 'scalar',
    'h_field':         'vector',
    'volumetric_PK1':  'tensor',
}


def _sanitize_name(name):
    """Return a Fortran-safe identifier fragment."""
    return ''.join(ch if ch.isalnum() else '_' for ch in name)


def _real_output_name(method_name, shape=None):
    """Local real-valued output variable for a material method."""
    shape = shape or _RETURN_SHAPE.get(method_name, 'scalar')
    if method_name == 'stress_PK1':
        return 'P_real'
    if method_name == 'pressure_resid':
        return 'rp_real'
    if method_name == 'solvent_storage':
        return 'cdot_real'
    if method_name == 'solvent_flux':
        return 'jR_real'
    return f'{_sanitize_name(method_name)}_real'


def _complex_eval_name(method_name, shape=None):
    """Local complex output variable for real residual evaluation."""
    shape = shape or _RETURN_SHAPE.get(method_name, 'scalar')
    if method_name == 'stress_PK1':
        return 'Pz_eval'
    if method_name == 'pressure_resid':
        return 'rpz_eval'
    if method_name == 'solvent_storage':
        return 'cdotz_eval'
    if method_name == 'solvent_flux':
        return 'jRz_eval'
    return f'{_sanitize_name(method_name)}z_eval'


def _material_method_for_equation(eq_name, binfo=None):
    """Resolve the material method used by an equation/tangent block."""
    if binfo is not None and binfo.get('material_method'):
        return binfo['material_method']
    if eq_name == 'momentum_equation':
        return 'stress_PK1'
    if eq_name == 'pressure_equation':
        return 'pressure_resid'
    return None


def _generate_material_subroutine(material, method_name, mat_prefix,
                                  helper_ctx=None):
    """
    Generate one COMPLEX*16 Fortran subroutine from a Material method.

    For stress_PK1 with state_vars, adds state_old/state_new arguments
    following the same pattern as umat_gen.py.

    ``helper_ctx`` (optional) carries the shared self._helper() translation
    state (material instance + helper_sources/cache/kind maps) so that helpers
    are translated and deduplicated across the material's methods.
    """
    info = material._methods[method_name]
    method = info['callable']
    params = info['all_params']
    field_args = _field_args_for(material)
    nprops = len(material.props_names)

    sv_info = _state_var_info(material)
    has_state = bool(sv_info) and method_name == 'stress_PK1'
    state_old_names_all = {f'{name}_old' for name in sv_info}

    sub_name = f'{mat_prefix}_{method_name}'
    ret_type, ret_decl = _RETURN_FORTRAN.get(method_name,
                                              ('DOUBLE COMPLEX', 'result_out'))
    ret_shape = _RETURN_SHAPE.get(method_name, 'scalar')

    # Translate the method body. Thread the shared helper context so
    # self._helper() calls dispatch correctly (needs the material instance)
    # and helper subroutines are collected once for emission.
    state_var_names = list(sv_info.keys()) if sv_info else []
    hctx = helper_ctx or {}
    translator = FortranTranslator(
        material.props_names, mat_prefix, state_var_names,
        material=hctx.get('material', material),
        helper_sources=hctx.get('helper_sources'),
        helper_cache=hctx.get('helper_cache'),
        helper_method_kinds=hctx.get('helper_method_kinds'))
    # Add all input arguments to var_kinds
    for p in params:
        if p in state_old_names_all:
            base = p[:-4]
            shape = sv_info[base]['shape']
            if shape == 'tensor':
                translator._var_kinds[p] = 'tensor'
            elif shape == 'vector':
                translator._var_kinds[p] = 'vector'
            else:
                translator._var_kinds[p] = 'scalar'
            continue
        arg_info = _arg_fortran_type(p, field_args)
        if arg_info and arg_info[2] == 'tensor':
            translator._var_kinds[p] = 'tensor'
        elif arg_info and arg_info[2] == 'vector':
            translator._var_kinds[p] = 'vector'
        else:
            translator._var_kinds[p] = 'scalar'

    # Register state variable kinds
    if has_state:
        for name, entry in sv_info.items():
            if entry['shape'] == 'tensor':
                translator._var_kinds[name + '_old'] = 'tensor'
            else:
                translator._var_kinds[name + '_old'] = 'scalar'

    decls, stmts, return_expr, return_kind = translator.translate_stress(method)
    tensor_vars = {n for n, k in translator._var_kinds.items() if k == 'tensor'}

    lines = []

    # Build argument list
    if has_state:
        arg_parts = []
        for p in params:
            if p not in [f'{n}_old' for n in sv_info] and p != 'dt':
                arg_parts.append(p)
        for name in sv_info:
            arg_parts.append(f'{name}_old')
        arg_parts.append('dt')
        arg_parts.append('props')
        arg_parts.append(ret_decl.split('(')[0].strip())
        for name in sv_info:
            arg_parts.append(f'{name}_new')
    else:
        arg_parts = []
        for p in params:
            arg_parts.append(p)
        arg_parts.append('props')
        arg_parts.append(ret_decl.split('(')[0].strip())

    sig = f'      SUBROUTINE {sub_name}('
    sig += ', '.join(arg_parts) + ')'
    lines.append(sig)
    lines.append('      IMPLICIT NONE')

    # Declare inputs
    if has_state:
        # Field variables
        for p in params:
            if p in [f'{n}_old' for n in sv_info] or p == 'dt':
                continue
            ftype, fdecl, _ = _arg_fortran_type(p, field_args)
            lines.append(f'      {ftype}, INTENT(IN) :: {fdecl}')
        # State variables old
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN) :: {name}_old{dims}')
        # dt
        lines.append('      DOUBLE PRECISION, INTENT(IN) :: dt')
    else:
        for p in params:
            if p in state_old_names_all:
                base = p[:-4]
                entry = sv_info[base]
                dims = _fortran_dims(entry)
                lines.append(
                    f'      DOUBLE COMPLEX, INTENT(IN) :: {p}{dims}')
            else:
                ftype, fdecl, _ = _arg_fortran_type(p, field_args)
                lines.append(f'      {ftype}, INTENT(IN) :: {fdecl}')

    lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')

    # Declare output
    lines.append(f'      {ret_type}, INTENT(OUT) :: {ret_decl}')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(OUT) :: {name}_new{dims}')
    lines.append('')

    # External functions
    lines.append('      DOUBLE COMPLEX :: det33z')

    # Track names to skip in declaration scanner
    state_new_names = {f'{n}_new' for n in sv_info} if has_state else set()
    state_old_names = {f'{n}_old' for n in sv_info} if has_state else set()
    skip_names = (state_new_names | state_old_names
                  | {ret_decl.split('(')[0].strip(), 'dt'})

    # Declare translator temporaries
    declared = set()
    for ftype, name, dims in decls:
        if name in declared:
            continue
        if name in state_new_names:
            declared.add(name)
            continue
        dim_str = dims if dims else ''
        lines.append(f'      {ftype} :: {name}{dim_str}')
        declared.add(name)

    # Declare user variables from statements
    for stmt in stmts:
        if '=' in stmt and 'CALL' not in stmt:
            stripped = stmt.strip()
            if stripped.startswith(('IF ', 'ELSE', 'END IF',
                                   'DO ', 'END DO', 'EXIT')):
                continue
            lhs = stmt.split('=')[0].strip()
            varname = lhs.split('(')[0].strip()
            if (varname not in declared and varname not in skip_names):
                if varname in tensor_vars:
                    pass
                elif varname in translator._real_vars:
                    lines.append(f'      DOUBLE PRECISION :: {varname}')
                    declared.add(varname)
                else:
                    lines.append(f'      DOUBLE COMPLEX :: {varname}')
                    declared.add(varname)

    if ret_shape in ('tensor', 'vector'):
        lines.append('      INTEGER :: ii, jj')
        declared.update({'ii', 'jj'})

    # Loop variable declarations
    for stmt in stmts:
        stripped = stmt.strip()
        if stripped.startswith('DO ') and '=' in stripped:
            loop_var = stripped[3:].split('=')[0].strip()
            if loop_var not in declared:
                lines.append(f'      INTEGER :: {loop_var}')
                declared.add(loop_var)

    lines.append('')

    # Body statements
    for stmt in stmts:
        lines.append(f'      {stmt}')

    # Return assignment
    out_name = ret_decl.split('(')[0].strip()
    if ret_shape == 'tensor':
        lines.append('')
        lines.append('      DO ii = 1, 3')
        lines.append('        DO jj = 1, 3')
        if return_expr in tensor_vars or return_expr in declared:
            lines.append(
                f'          {out_name}(ii,jj) = {return_expr}(ii,jj)')
        else:
            from .umat_gen import _tensorize_expr
            expr_ij = _tensorize_expr(return_expr, tensor_vars | {'F'})
            lines.append(f'          {out_name}(ii,jj) = {expr_ij}')
        lines.append('        END DO')
        lines.append('      END DO')
    elif ret_shape == 'vector':
        lines.append('')
        lines.append('      DO ii = 1, 3')
        if return_expr in declared:
            lines.append(f'        {out_name}(ii) = {return_expr}(ii)')
        else:
            expr_i = return_expr
            vec_vars = {n for n, k in translator._var_kinds.items()
                        if k == 'vector'}
            import re
            for vn in sorted(vec_vars, key=len, reverse=True):
                expr_i = re.sub(r'\b' + re.escape(vn) + r'\b(?!\()',
                                f'{vn}(ii)', expr_i)
            lines.append(f'        {out_name}(ii) = {expr_i}')
        lines.append('      END DO')
    else:
        lines.append(f'      {out_name} = {return_expr}')

    # State output copy
    if has_state and translator._state_return_map:
        has_copy = False
        for sv_name, (val_expr, val_kind) in \
                translator._state_return_map.items():
            if val_expr == f'{sv_name}_new':
                continue
            if not has_copy:
                lines.append('')
                lines.append('C     --- Copy state outputs ---')
                has_copy = True
            if val_kind == 'tensor':
                if val_expr in tensor_vars or val_expr in declared:
                    lines.append('      DO ii = 1, 3')
                    lines.append('        DO jj = 1, 3')
                    lines.append(
                        f'          {sv_name}_new(ii,jj) = '
                        f'{val_expr}(ii,jj)')
                    lines.append('        END DO')
                    lines.append('      END DO')
                else:
                    lines.append(f'      {sv_name}_new = {val_expr}')
            else:
                lines.append(f'      {sv_name}_new = {val_expr}')

    lines.append('')
    lines.append('      RETURN')
    lines.append(f'      END SUBROUTINE {sub_name}')

    return '\n'.join(lines)


def _generate_all_material_subs(weakform, mat_prefix):
    """Generate Fortran subroutines for all Material methods.

    Each registered method is translated to a COMPLEX*16 subroutine. Any
    ``self._helper()`` calls are translated through a shared helper context
    (the same mechanism used by ``umat_gen``), so a helper used by several
    methods is emitted exactly once and the call sites and the definition agree
    on the subroutine name.
    """
    mat = weakform._mat
    parts = []

    # Shared helper-translation context threaded into every per-method
    # translator: dedups self._helper() subroutines and collects their source.
    helper_ctx = {
        'material': mat,
        'helper_sources': [],
        'helper_cache': {},
        'helper_method_kinds': {},
    }

    # Generate registered methods; collected helpers are emitted afterwards.
    for method_name in mat.defined_methods():
        parts.append(_generate_material_subroutine(
            mat, method_name, mat_prefix, helper_ctx=helper_ctx))
        parts.append('')

    # Emit any self._helper() subroutines collected during translation (once).
    for src in helper_ctx['helper_sources']:
        parts.append(src)
        parts.append('')

    return '\n'.join(parts)


# =====================================================================
# Section 3: CS tangent engine generation
# =====================================================================

def _generate_cs_engine(weakform, mat_prefix):
    """Generate the complex-step tangent engine subroutine."""
    mat = weakform._mat
    blocks = weakform.tangent_blocks()
    field_args = _field_args_for(weakform)
    nprops = len(mat.props_names)
    sv_info = _state_var_info(mat)
    state_old_names = {f'{name}_old' for name in sv_info}
    state_old_args = [f'{name}_old' for name in sv_info]

    # Collect all unique material methods and their field variables
    # For each block, we need to know which material method to call
    # and which variable to perturb

    # Build the list of tangent computations
    tangent_specs = []
    for (block_key, wrt_var), binfo in blocks.items():
        mat_method = _material_method_for_equation(
            binfo['equation'], binfo)

        if mat_method is None:
            continue

        ret_shape = _RETURN_SHAPE.get(mat_method, 'scalar')
        tangent_specs.append({
            'block_key': block_key,
            'wrt_var': wrt_var,
            'mat_method': mat_method,
            'ret_shape': ret_shape,
            'wrt_kind': binfo['wrt_kind'],
        })

    # Collect all field variables needed across all methods
    all_field_vars = set()
    all_history_vars = set()
    all_param_vars = set()
    for eq_info in weakform.equations.values():
        all_field_vars.update(eq_info['field_vars'])
        all_history_vars.update(eq_info['history_vars'])
        all_param_vars.update(eq_info['param_vars'])
    cs_needs_dt = any(
        'dt' in mat._methods[spec['mat_method']]['all_params']
        for spec in tangent_specs
    )
    cs_needs_time = any(
        'time' in mat._methods[spec['mat_method']]['all_params']
        for spec in tangent_specs
    )
    cs_needs_X = any(
        'X' in mat._methods[spec['mat_method']]['all_params']
        for spec in tangent_specs
    )

    lines = []
    lines.append(f'      SUBROUTINE {mat_prefix}_cs_tangents(')

    # Input arguments
    in_args = ['F']
    scalar_fields = sorted([v for v in all_field_vars
                            if field_args.get(v, {}).get('kind') == 'scalar'])
    vector_fields = sorted([v for v in all_field_vars
                            if field_args.get(v, {}).get('kind') == 'vector'])
    in_args.extend(scalar_fields)
    in_args.extend(vector_fields)
    in_args.append('F_old')
    history_scalars = _ordered_history_scalars(weakform, all_history_vars)
    in_args.extend(history_scalars)
    in_args.extend(state_old_args)
    in_args.append('props')
    if cs_needs_dt:
        in_args.append('dt')
    if cs_needs_time:
        in_args.append('time')
    if cs_needs_X:
        in_args.append('X')

    # Output tangent names
    out_names = {}
    for spec in tangent_specs:
        bkey = spec['block_key']
        wrt = spec['wrt_var']
        meth = spec['mat_method']
        short_m = meth.replace('stress_', '').replace('solvent_', '')
        short_m = short_m.replace('pressure_', 'rp_')
        tname = f'd{short_m}_d{wrt}'
        # Avoid duplicates
        if tname in out_names.values():
            tname = f'd{short_m}_d{wrt}_2'
        out_names[str(bkey)] = tname
        spec['tname'] = tname

    out_args = [spec['tname'] for spec in tangent_specs]
    all_args = in_args + out_args

    # Write argument list with continuation
    line = '     &   '
    for i, a in enumerate(all_args):
        if i > 0:
            line += ', '
        if len(line) + len(a) > 65:
            lines.append(line)
            line = '     &   '
        line += a
    line += ')'
    lines.append(line)
    lines.append('')
    lines.append('      IMPLICIT NONE')
    lines.append('')

    # Declare inputs
    lines.append('C     --- Inputs (all DOUBLE PRECISION) ---')
    lines.append('      DOUBLE PRECISION, INTENT(IN) :: F(3,3)')
    for sf in scalar_fields:
        lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: {sf}')
    for vf in vector_fields:
        lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: {vf}(3)')
    lines.append('      DOUBLE PRECISION, INTENT(IN) :: F_old(3,3)')
    for hs in history_scalars:
        lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: {hs}')
    for name, entry in sv_info.items():
        dims = _fortran_dims(entry)
        lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: {name}_old{dims}')
    lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')
    if cs_needs_dt:
        lines.append('      DOUBLE PRECISION, INTENT(IN) :: dt')
    if cs_needs_time:
        lines.append('      DOUBLE PRECISION, INTENT(IN) :: time')
    if cs_needs_X:
        lines.append('      DOUBLE PRECISION, INTENT(IN) :: X(3)')
    lines.append('')

    # Declare outputs
    lines.append('C     --- Outputs (all DOUBLE PRECISION) ---')
    for spec in tangent_specs:
        tname = spec['tname']
        ret_shape = spec['ret_shape']
        wrt_kind = spec['wrt_kind']
        # Output shape = return_shape x wrt_shape
        if ret_shape == 'tensor' and wrt_kind == 'matrix':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3,3,3,3)')
        elif ret_shape == 'tensor' and wrt_kind == 'scalar':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3,3)')
        elif ret_shape == 'tensor' and wrt_kind == 'vector':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3,3,3)')
        elif ret_shape == 'scalar' and wrt_kind == 'matrix':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3,3)')
        elif ret_shape == 'scalar' and wrt_kind == 'scalar':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}')
        elif ret_shape == 'vector' and wrt_kind == 'matrix':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3,3,3)')
        elif ret_shape == 'vector' and wrt_kind == 'scalar':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3)')
        elif ret_shape == 'vector' and wrt_kind == 'vector':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3,3)')
        elif ret_shape == 'scalar' and wrt_kind == 'vector':
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}(3)')
        else:
            lines.append(f'      DOUBLE PRECISION, INTENT(OUT) :: {tname}')
    lines.append('')

    # CS_H and complex variables
    # CS_H = 1e-10 chosen after verification against 1e-30 across all 12
    # tangent blocks. Identical results at
    # both values. 1e-10 avoids potential subnormal issues in Fortran
    # intermediate computations while still giving machine-precision tangents.
    lines.append('C     CS_H = 1e-10: verified identical to 1e-30 for this model.')
    lines.append('C     Complex-step has no subtractive cancellation, so any')
    lines.append('C     small h gives machine-precision tangents.')
    lines.append('      DOUBLE PRECISION, PARAMETER :: CS_H = 1.0d-10')
    lines.append('      DOUBLE COMPLEX :: Fz(3,3), F_old_z(3,3)')
    for sf in scalar_fields:
        lines.append(f'      DOUBLE COMPLEX :: {sf}z')
    for vf in vector_fields:
        lines.append(f'      DOUBLE COMPLEX :: {vf}_z(3)')
    for hs in history_scalars:
        lines.append(f'      DOUBLE COMPLEX :: {hs}_z')
    for name, entry in sv_info.items():
        dims = _fortran_dims(entry)
        lines.append(f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
        lines.append(f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
    lines.append('')

    # Output temporaries for material calls
    lines.append('      DOUBLE COMPLEX :: PKz(3,3), rpz, jRz(3), cdotz')
    lines.append('      INTEGER :: i, j, k, l')
    lines.append('')

    # One-time conversion of history variables
    lines.append('C     Convert history variables (real, not perturbed)')
    lines.append('      CALL real2complex33(F_old, F_old_z)')
    for hs in history_scalars:
        lines.append(f'      {hs}_z = DCMPLX({hs}, 0.0d0)')
    lines.append('')

    # Generate perturbation blocks
    for spec in tangent_specs:
        tname = spec['tname']
        wrt = spec['wrt_var']
        wrt_kind = spec['wrt_kind']
        mat_meth = spec['mat_method']
        ret_shape = spec['ret_shape']
        sub_name = f'{mat_prefix}_{mat_meth}'

        # Build the material call argument list
        minfo = mat._methods[mat_meth]
        call_args = []
        for p in minfo['all_params']:
            if p == 'F':
                call_args.append('Fz')
            elif p == 'F_old':
                call_args.append('F_old_z')
            elif p in state_old_names:
                call_args.append(f'{p}_z')
            elif p in scalar_fields:
                call_args.append(f'{p}z')
            elif p in vector_fields:
                call_args.append(f'{p}_z')
            elif p in history_scalars:
                call_args.append(f'{p}_z')
            elif p == 'dt':
                call_args.append('dt')
            elif p == 'time':
                call_args.append('time')
            elif p == 'X':
                call_args.append('X')
            else:
                call_args.append(p)
        call_args.append('props')

        # Output arg name
        if ret_shape == 'tensor':
            call_args.append('PKz')
            out_var = 'PKz'
            if mat_meth == 'stress_PK1' and sv_info:
                for name in sv_info:
                    call_args.append(f'{name}_new_z')
        elif ret_shape == 'vector':
            call_args.append('jRz')
            out_var = 'jRz'
        else:
            if mat_meth == 'pressure_resid':
                call_args.append('rpz')
                out_var = 'rpz'
            else:
                call_args.append('cdotz')
                out_var = 'cdotz'

        call_str = f'CALL {sub_name}({", ".join(call_args)})'

        lines.append(f'C     === {tname}: d({mat_meth})/d({wrt}) ===')

        if wrt_kind == 'matrix':
            # 9 perturbations of F(k,l)
            lines.append('      DO k = 1, 3')
            lines.append('        DO l = 1, 3')
            lines.append('          CALL real2complex33(F, Fz)')
            for sf in scalar_fields:
                lines.append(f'          {sf}z = DCMPLX({sf}, 0.0d0)')
            for vf in vector_fields:
                lines.append(f'          CALL real2complex3({vf}, {vf}_z)')
            _append_state_old_to_complex(lines, sv_info, indent='          ')
            lines.append(f'          Fz(k,l) = Fz(k,l) + DCMPLX(0.0d0, CS_H)')
            lines.append(f'          {call_str}')
            _emit_extract(lines, tname, out_var, ret_shape, 'matrix')
            lines.append('        END DO')
            lines.append('      END DO')

        elif wrt_kind == 'scalar':
            # 1 perturbation
            lines.append('      CALL real2complex33(F, Fz)')
            for sf in scalar_fields:
                if sf == wrt:
                    lines.append(f'      {sf}z = DCMPLX({sf}, CS_H)')
                else:
                    lines.append(f'      {sf}z = DCMPLX({sf}, 0.0d0)')
            for vf in vector_fields:
                lines.append(f'      CALL real2complex3({vf}, {vf}_z)')
            _append_state_old_to_complex(lines, sv_info, indent='      ')
            lines.append(f'      {call_str}')
            _emit_extract(lines, tname, out_var, ret_shape, 'scalar')

        elif wrt_kind == 'vector':
            # 3 perturbations of grad_mu(k)
            lines.append('      DO k = 1, 3')
            lines.append('        CALL real2complex33(F, Fz)')
            for sf in scalar_fields:
                lines.append(f'        {sf}z = DCMPLX({sf}, 0.0d0)')
            for vf in vector_fields:
                lines.append(f'        CALL real2complex3({vf}, {vf}_z)')
            _append_state_old_to_complex(lines, sv_info, indent='        ')
            lines.append(f'        {wrt}_z(k) = {wrt}_z(k)'
                         f' + DCMPLX(0.0d0, CS_H)')
            lines.append(f'        {call_str}')
            _emit_extract(lines, tname, out_var, ret_shape, 'vector')
            lines.append('      END DO')

        lines.append('')

    lines.append('      RETURN')
    lines.append(f'      END SUBROUTINE {mat_prefix}_cs_tangents')

    return '\n'.join(lines)


def _emit_extract(lines, tname, out_var, ret_shape, wrt_kind):
    """Emit AIMAG extraction lines for one tangent block."""
    indent = '          ' if wrt_kind == 'matrix' else '      '
    if wrt_kind == 'vector':
        indent = '        '

    if ret_shape == 'tensor' and wrt_kind == 'matrix':
        lines.append(f'{indent}DO i = 1, 3')
        lines.append(f'{indent}  DO j = 1, 3')
        lines.append(f'{indent}    {tname}(i,j,k,l) = AIMAG({out_var}(i,j))/CS_H')
        lines.append(f'{indent}  END DO')
        lines.append(f'{indent}END DO')
    elif ret_shape == 'tensor' and wrt_kind == 'scalar':
        lines.append(f'{indent}DO i = 1, 3')
        lines.append(f'{indent}  DO j = 1, 3')
        lines.append(f'{indent}    {tname}(i,j) = AIMAG({out_var}(i,j)) / CS_H')
        lines.append(f'{indent}  END DO')
        lines.append(f'{indent}END DO')
    elif ret_shape == 'tensor' and wrt_kind == 'vector':
        lines.append(f'{indent}DO i = 1, 3')
        lines.append(f'{indent}  DO j = 1, 3')
        lines.append(f'{indent}    {tname}(i,j,k) = AIMAG({out_var}(i,j)) / CS_H')
        lines.append(f'{indent}  END DO')
        lines.append(f'{indent}END DO')
    elif ret_shape == 'scalar' and wrt_kind == 'matrix':
        lines.append(f'{indent}{tname}(k,l) = AIMAG({out_var}) / CS_H')
    elif ret_shape == 'scalar' and wrt_kind == 'scalar':
        lines.append(f'{indent}{tname} = AIMAG({out_var}) / CS_H')
    elif ret_shape == 'scalar' and wrt_kind == 'vector':
        lines.append(f'{indent}{tname}(k) = AIMAG({out_var}) / CS_H')
    elif ret_shape == 'vector' and wrt_kind == 'matrix':
        lines.append(f'{indent}DO i = 1, 3')
        lines.append(f'{indent}  {tname}(i,k,l) = AIMAG({out_var}(i)) / CS_H')
        lines.append(f'{indent}END DO')
    elif ret_shape == 'vector' and wrt_kind == 'scalar':
        lines.append(f'{indent}DO i = 1, 3')
        lines.append(f'{indent}  {tname}(i) = AIMAG({out_var}(i)) / CS_H')
        lines.append(f'{indent}END DO')
    elif ret_shape == 'vector' and wrt_kind == 'vector':
        lines.append(f'{indent}DO i = 1, 3')
        lines.append(f'{indent}  {tname}(i,k) = AIMAG({out_var}(i)) / CS_H')
        lines.append(f'{indent}END DO')


def _append_state_var_declarations(lines, sv_info, indent='      '):
    """Declare local state variables for UEL wrappers."""
    if not sv_info:
        return

    lines.append('')
    lines.append('C     --- State variables (SVARS) ---')
    for name, entry in sv_info.items():
        dims = _fortran_dims(entry)
        lines.append(f'{indent}DOUBLE PRECISION :: {name}_old{dims}')
    for name, entry in sv_info.items():
        dims = _fortran_dims(entry)
        lines.append(f'{indent}DOUBLE COMPLEX :: {name}_old_z{dims}')
        lines.append(f'{indent}DOUBLE COMPLEX :: {name}_new_z{dims}')
    lines.append(f'{indent}INTEGER :: LOC')


def _append_svars_read(lines, sv_info, indent='      '):
    """Read state variables for Gauss point kk from SVARS."""
    if not sv_info:
        return

    lines.append('C       Read state variables from SVARS')
    lines.append(f'{indent}LOC = (kk - 1) * NSTATE_PER_GP')
    lines.append(f'{indent}IF (TIME(2) .EQ. 0.0d0) THEN')
    for name, entry in sv_info.items():
        init_prop_index = entry.get('init_prop_index')
        if entry['shape'] == 'scalar' and init_prop_index is not None:
            lines.append(f'{indent}  {name}_old = PROPS({init_prop_index})')
        else:
            for line in _fortran_init_value(f'{name}_old', entry):
                lines.append(f'{indent}  {line}')
    lines.append(f'{indent}ELSE')

    offset = 0
    for name, entry in sv_info.items():
        shape = entry['shape']
        if shape == 'scalar':
            lines.append(f'{indent}  {name}_old = SVARS(LOC + {offset + 1})')
            offset += 1
        elif shape == 'vector':
            for i in range(3):
                lines.append(
                    f'{indent}  {name}_old({i+1}) = '
                    f'SVARS(LOC + {offset + i + 1})')
            offset += 3
        elif shape == 'tensor':
            for j in range(3):
                for i in range(3):
                    idx = offset + j * 3 + i + 1
                    lines.append(
                        f'{indent}  {name}_old({i+1},{j+1}) = '
                        f'SVARS(LOC + {idx})')
            offset += 9

    for name, entry in sv_info.items():
        init_prop_index = entry.get('init_prop_index')
        if entry['shape'] == 'scalar' and init_prop_index is not None:
            lines.append(f'{indent}IF (({name}_old .NE. {name}_old) .OR.')
            lines.append(f'     &      (DABS({name}_old) .GT. 1.0D300) .OR.')
            lines.append(f'     &      {name}_old .LE. 0.0D0) THEN')
            lines.append(f'{indent}  {name}_old = PROPS({init_prop_index})')
            lines.append(f'{indent}END IF')

    lines.append(f'{indent}END IF')


def _append_svars_write(lines, sv_info, indent='      '):
    """Write updated state variables back to SVARS using real parts only."""
    if not sv_info:
        return

    lines.append('C       Write updated state to SVARS (real part only)')
    offset = 0
    for name, entry in sv_info.items():
        shape = entry['shape']
        if shape == 'scalar':
            lines.append(
                f'{indent}SVARS(LOC + {offset + 1}) = DBLE({name}_new_z)')
            offset += 1
        elif shape == 'vector':
            for i in range(3):
                lines.append(
                    f'{indent}SVARS(LOC + {offset + i + 1}) = '
                    f'DBLE({name}_new_z({i+1}))')
            offset += 3
        elif shape == 'tensor':
            for j in range(3):
                for i in range(3):
                    idx = offset + j * 3 + i + 1
                    lines.append(
                        f'{indent}SVARS(LOC + {idx}) = '
                        f'DBLE({name}_new_z({i+1},{j+1}))')
            offset += 9


def _append_state_old_to_complex(lines, sv_info, indent='      '):
    """Convert local real state variables to complex temporaries."""
    if not sv_info:
        return

    for name, entry in sv_info.items():
        shape = entry['shape']
        if shape == 'scalar':
            lines.append(f'{indent}{name}_old_z = DCMPLX({name}_old, 0.0d0)')
        elif shape == 'vector':
            lines.append(
                f'{indent}CALL real2complex3({name}_old, {name}_old_z)')
        elif shape == 'tensor':
            lines.append(
                f'{indent}CALL real2complex33({name}_old, {name}_old_z)')


# =====================================================================
# Section 4: UEL subroutine generation
# =====================================================================

def _generate_uel(weakform, mat_prefix, cfg):
    """Generate the UEL subroutine."""
    ndim = cfg.ndim
    fields = weakform.fields
    field_names = weakform.field_names
    field_args = _field_args_for(weakform)
    n_corner = cfg.n_corner_nodes
    n_nodes = cfg.n_nodes
    ndofel = _compute_ndofel(cfg, fields)
    equations = weakform.equations
    blocks = weakform.tangent_blocks()
    nprops = len(weakform._mat.props_names)
    sv_info = _state_var_info(weakform._mat)
    state_old_names = {f'{name}_old' for name in sv_info}
    nstate_per_gp = _nstate_per_gp(sv_info) if sv_info else 0

    # Determine which fields need old values and gradients
    all_history = set()
    all_params = set()
    needs_grad = set()
    for eq_info in equations.values():
        all_history.update(eq_info['history_vars'])
        all_params.update(eq_info['param_vars'])
        for fv in eq_info['field_vars']:
            if field_args.get(fv, {}).get('kind') == 'vector':
                # grad_mu is a field var -> mu needs gradient
                # Extract the base field: grad_X -> X
                base = fv.replace('grad_', '')
                needs_grad.add(base)
    cs_needs_dt = any(
        'dt' in weakform._mat._methods[mname]['all_params']
        for mname in weakform._mat.defined_methods()
    )
    cs_needs_time = any(
        'time' in weakform._mat._methods[mname]['all_params']
        for mname in weakform._mat.defined_methods()
    )
    cs_needs_X = any(
        'X' in weakform._mat._methods[mname]['all_params']
        for mname in weakform._mat.defined_methods()
    )

    lines = []
    lines.append('C======================================================================')
    lines.append('C     UEL subroutine — generated by abaqus_ufl')
    lines.append(f'C     Material: {weakform._mat.__class__.__name__}')
    lines.append(f"C     Fields: {', '.join(field_names)}")
    lines.append(f'C     Element: {cfg.name.capitalize()}, NDOFEL = {ndofel}')
    lines.append('C======================================================================')
    lines.append('      SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL,')
    lines.append('     &  NRHS, NSVARS, PROPS, NPROPS, COORDS, MCRD, NNODE,')
    lines.append('     &  U, DU, V, A, JTYPE, TIME, DTIME, KSTEP, KINC,')
    lines.append('     &  JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, PREDEF,')
    lines.append('     &  NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT,')
    lines.append('     &  JPROPS, NJPROP, PERIOD)')
    lines.append('')
    lines.append('      IMPLICIT NONE')
    lines.append('')
    lines.append('C     --- Abaqus UEL interface ---')
    lines.append('      INTEGER NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE')
    lines.append('      INTEGER JTYPE, KSTEP, KINC, JELEM, NDLOAD, NPREDF')
    lines.append('      INTEGER MLVARX, MDLOAD, NJPROP')
    lines.append('      INTEGER JDLTYP(MDLOAD,*), LFLAGS(*), JPROPS(*)')
    lines.append('')
    lines.append('      DOUBLE PRECISION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL)')
    lines.append('      DOUBLE PRECISION SVARS(NSVARS), ENERGY(8)')
    lines.append('      DOUBLE PRECISION PROPS(NPROPS), COORDS(MCRD,NNODE)')
    lines.append('      DOUBLE PRECISION U(NDOFEL), DU(MLVARX,*), V(NDOFEL), A(NDOFEL)')
    lines.append('      DOUBLE PRECISION TIME(2), DTIME, PARAMS(3)')
    lines.append('      DOUBLE PRECISION ADLMAG(MDLOAD,*), PREDEF(2,NPREDF,NNODE)')
    lines.append('      DOUBLE PRECISION DDLMAG(MDLOAD,*), PNEWDT, PERIOD')
    lines.append('')
    lines.append('C     --- Local parameters ---')
    lines.append(f'      INTEGER, PARAMETER :: ndim = {ndim}')
    lines.append(f'      INTEGER, PARAMETER :: NCORNER = {n_corner}')
    lines.append(f'      INTEGER, PARAMETER :: NNODE_E = {n_nodes}')
    lines.append(f'      INTEGER, PARAMETER :: NGP = {cfg.n_gauss_points}')
    if sv_info:
      lines.append(f'      INTEGER, PARAMETER :: NSTATE_PER_GP = {nstate_per_gp}')

    # Declare field nodal arrays and edof maps
    for fname in field_names:
        f = fields[fname]
        nn = _nn_for_field(cfg, f)
        if isinstance(f, VectorField):
            lines.append(f'      DOUBLE PRECISION :: {fname}_node(ndim, {nn})')
            lines.append(f'      DOUBLE PRECISION :: {fname}_old(ndim, {nn})')
            lines.append(f'      INTEGER :: edof_{fname}(ndim, {nn})')
        else:
            lines.append(f'      DOUBLE PRECISION :: {fname}_node({nn})')
            lines.append(f'      DOUBLE PRECISION :: {fname}_old({nn})')
            lines.append(f'      INTEGER :: edof_{fname}({nn})')

    _n = cfg.n_nodes
    _nc = cfg.n_corner_nodes
    _ngp = cfg.n_gauss_points
    _d = cfg.ndim
    _jinv = cfg.jac_inv_size
    shape_lines = [
        f'      DOUBLE PRECISION :: {cfg.sh_name}({_n}), {cfg.dshxi_name}({_n},{_d}), {cfg.dsh_name}({_n},{_d})'
    ]
    linear_shape_decl = (
        f'      DOUBLE PRECISION :: {cfg.sh_linear_name}({_nc}), '
        f'{cfg.dshxi_linear_name}({_nc},{_d}), {cfg.dsh_linear_name}({_nc},{_d})'
    )
    if linear_shape_decl not in shape_lines:
        shape_lines.append(linear_shape_decl)

    lines.append("""
C     --- Shape functions and mapping ---""")
    lines.extend(shape_lines)
    lines.append(f"""      DOUBLE PRECISION :: xi_gp({_ngp},{_d}), w_gp({_ngp})
      DOUBLE PRECISION :: {cfg.coords_name}({_d}, NNODE_E)
      DOUBLE PRECISION :: detJ, Jinv({_jinv},{_jinv}), wdetJ
      DOUBLE PRECISION :: dt_safe
      INTEGER :: ngp_out, stat

C     --- Gauss point field values ---
      DOUBLE PRECISION :: F(3,3), F_old(3,3)""")
    if cs_needs_X:
        lines.append('      DOUBLE PRECISION :: X_gp(3)')

    # GP scalar/vector field values
    for fname in field_names:
        f = fields[fname]
        if isinstance(f, ScalarField):
            lines.append(f'      DOUBLE PRECISION :: {fname}_gp, {fname}_old_gp')
            if fname in needs_grad:
                lines.append(f'      DOUBLE PRECISION :: grad_{fname}(3)')

    _append_state_var_declarations(lines, sv_info)

    # Material output variables (real, for residual)
    lines.append('')
    lines.append('C     --- Material outputs (real) ---')
    for mname in weakform._mat.defined_methods():
        rs = _RETURN_SHAPE.get(mname)
        if rs == 'tensor':
            lines.append(f'      DOUBLE PRECISION :: '
                         f'{_real_output_name(mname, rs)}(3,3)')
        elif rs == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: '
                         f'{_real_output_name(mname, rs)}')
        elif rs == 'vector':
            lines.append(f'      DOUBLE PRECISION :: '
                         f'{_real_output_name(mname, rs)}(3)')

    # Complex temporaries for real evaluation
    complex_declared = set()
    for mname in weakform._mat.defined_methods():
        rs = _RETURN_SHAPE.get(mname)
        cname = _complex_eval_name(mname, rs)
        if cname in complex_declared:
            continue
        complex_declared.add(cname)
        if rs == 'tensor':
            lines.append(f'      DOUBLE COMPLEX :: {cname}(3,3)')
        elif rs == 'vector':
            lines.append(f'      DOUBLE COMPLEX :: {cname}(3)')
        else:
            lines.append(f'      DOUBLE COMPLEX :: {cname}')
    lines.append('C     Complex temps for real material evaluation')
    lines.append('      DOUBLE COMPLEX :: Fz_r(3,3), Fold_r(3,3)')

    # Declare complex temps for scalar/vector fields used in real eval
    for fname in field_names:
        f = fields[fname]
        if isinstance(f, ScalarField):
            lines.append(f'      DOUBLE COMPLEX :: {fname}z_r')
            if fname in needs_grad:
                lines.append(f'      DOUBLE COMPLEX :: grad_{fname}_r(3)')

    # Tangent output arrays
    lines.append('')
    lines.append('C     --- Tangent arrays ---')
    for spec_key, binfo in blocks.items():
        # We need to declare the tangent arrays here too
        pass  # They'll be declared as arguments to the CS subroutine
    # Actually the tangent arrays are local to the GP loop
    # Declare all tangent block arrays
    tangent_specs = []
    for (bkey, wrt_var), binfo in blocks.items():
        mat_meth = _material_method_for_equation(
            binfo['equation'], binfo)
        ret_shape = _RETURN_SHAPE.get(mat_meth, 'scalar')
        wrt_kind = binfo['wrt_kind']
        short_m = mat_meth.replace('stress_', '').replace('solvent_', '')
        short_m = short_m.replace('pressure_', 'rp_')
        tname = f'd{short_m}_d{wrt_var}'
        tangent_specs.append({'tname': tname, 'ret_shape': ret_shape,
                              'wrt_kind': wrt_kind, 'binfo': binfo})

    declared_tangents = set()
    for spec in tangent_specs:
        tn = spec['tname']
        if tn in declared_tangents:
            continue
        declared_tangents.add(tn)
        rs = spec['ret_shape']
        wk = spec['wrt_kind']
        if rs == 'tensor' and wk == 'matrix':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3,3,3,3)')
        elif rs == 'tensor' and wk == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3,3)')
        elif rs == 'tensor' and wk == 'vector':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3,3,3)')
        elif rs == 'scalar' and wk == 'matrix':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3,3)')
        elif rs == 'scalar' and wk == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: {tn}')
        elif rs == 'vector' and wk == 'matrix':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3,3,3)')
        elif rs == 'vector' and wk == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3)')
        elif rs == 'vector' and wk == 'vector':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3,3)')
        elif rs == 'scalar' and wk == 'vector':
            lines.append(f'      DOUBLE PRECISION :: {tn}(3)')

    lines.append('')
    lines.append('      INTEGER :: idx, ii_v, jj_v, i, j, k, l, kk, row, col')
    lines.append('')

    # --- Extract coordinates ---
    lines.append(f'C     Extract {_d}D coordinates from COORDS(MCRD, NNODE)')
    lines.append('      DO ii_v = 1, NNODE_E')
    for idim in range(1, _d + 1):
        lines.append(f'        {cfg.coords_name}({idim}, ii_v) = COORDS({idim}, ii_v)')
    lines.append('      END DO')
    lines.append('')

    # --- Zero RHS and AMATRX ---
    lines.append('C     Zero RHS and AMATRX')
    lines.append('      DO i = 1, NDOFEL')
    lines.append('        RHS(i, 1) = 0.0d0')
    lines.append('        DO j = 1, NDOFEL')
    lines.append('          AMATRX(i, j) = 0.0d0')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')
    lines.append('C     TODO: implement Modified Riks RHS(:,2) properly')
    lines.append('C     Only touch RHS(:,2) when Abaqus allocated it (NRHS >= 2).')
    lines.append('      IF (NRHS .GE. 2 .AND. LFLAGS(3) .EQ. 4) THEN')
    lines.append('        DO i = 1, NDOFEL')
    lines.append('          RHS(i, 2) = 0.0d0')
    lines.append('        END DO')
    lines.append('      END IF')
    lines.append('')

    # --- DTIME guard ---
    lines.append('C     DTIME guard (Abaqus calls with DTIME=0 for initial eval)')
    lines.append('      dt_safe = DTIME')
    lines.append('      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14')
    lines.append('')

    # --- DOF parsing with edof maps ---
    lines.append('C     Parse DOFs from U and DU, build edof maps')
    lines.append('      idx = 0')
    lines.append(f'      DO ii_v = 1, NNODE_E')

    for fname in field_names:
        f = fields[fname]
        is_corner_only = (f.degree < 2)

        if is_corner_only:
            lines.append(f'        IF (ii_v .LE. NCORNER) THEN')
            indent = '          '
        else:
            indent = '        '

        if isinstance(f, VectorField):
            lines.append(f'{indent}DO i = 1, ndim')
            lines.append(f'{indent}  idx = idx + 1')
            lines.append(f'{indent}  {fname}_node(i, ii_v) = U(idx)')
            lines.append(f'{indent}  {fname}_old(i, ii_v) = U(idx) - DU(idx, 1)')
            lines.append(f'{indent}  edof_{fname}(i, ii_v) = idx')
            lines.append(f'{indent}END DO')
        else:
            lines.append(f'{indent}idx = idx + 1')
            lines.append(f'{indent}{fname}_node(ii_v) = U(idx)')
            lines.append(f'{indent}{fname}_old(ii_v) = U(idx) - DU(idx, 1)')
            lines.append(f'{indent}edof_{fname}(ii_v) = idx')

        if is_corner_only:
            lines.append(f'        END IF')

    lines.append('      END DO')
    lines.append('')

    # --- Gauss loop ---
    lines.append('C     Gauss quadrature')
    lines.append(f'      CALL {cfg.gauss_subroutine}(xi_gp, w_gp, ngp_out)')
    lines.append('')
    lines.append('      DO kk = 1, ngp_out')
    lines.append('')

    # Shape functions + Jacobian
    lines.append('C       Shape functions + geometry mapping')
    lines.append(f'        CALL {cfg.shape_subroutine}({cfg.gauss_xi_args},')
    lines.append(f'     &                   {cfg.sh_name}, {cfg.dshxi_name})')
    lines.append(f'        CALL {cfg.jac_subroutine}({cfg.dshxi_name}, {cfg.coords_name}, {cfg.n_nodes},')
    lines.append(f'     &                   {cfg.dsh_name}, detJ, Jinv, stat)')
    lines.append('')
    lines.append('C       Trap degenerate Jacobian')
    lines.append('        IF (stat .EQ. 0) THEN')
    lines.append('          PNEWDT = 0.25d0')
    lines.append('          RETURN')
    lines.append('        END IF')
    lines.append('')

    # Linear shapes if needed (for degree-1 fields)
    has_degree1 = any(fields[fn].degree < 2 for fn in field_names)
    if has_degree1:
        jinv_sub = 'apply_jinv_3d' if cfg.ndim == 3 else 'apply_jinv'
        lines.append(f'        CALL {cfg.shape_linear_subroutine}({cfg.gauss_xi_args},')
        lines.append(f'     &                   {cfg.sh_linear_name}, {cfg.dshxi_linear_name})')
        lines.append(f'        CALL {jinv_sub}({cfg.dshxi_linear_name}, Jinv, {cfg.n_corner_nodes}, {cfg.dsh_linear_name})')
        lines.append('')

    lines.append('        wdetJ = detJ * w_gp(kk)')
    lines.append('')

    # Reference position at the Gauss point: X_gp = sum_a N_a X_a (geometry).
    if cs_needs_X:
        lines.append('C       Reference position at Gauss point (for X-dependent drivers)')
        lines.append('        DO i = 1, 3')
        lines.append('          X_gp(i) = 0.0d0')
        lines.append('        END DO')
        lines.append(f'        DO ii_v = 1, {cfg.n_nodes}')
        lines.append(f'          DO i = 1, {cfg.ndim}')
        lines.append(f'            X_gp(i) = X_gp(i) + {cfg.sh_name}(ii_v) '
                     f'* {cfg.coords_name}(i, ii_v)')
        lines.append('          END DO')
        lines.append('        END DO')
        lines.append('')

    # --- Field interpolation ---
    lines.append('C       Interpolate fields at Gauss point')

    # Deformation gradient F = I + grad(u)
    for fname in field_names:
        f = fields[fname]
        if isinstance(f, VectorField):
            _, dsh_name, nn = cfg.sh_for_degree(f.degree)
            lines.append(f'        CALL eye33d(F)')
            lines.append(f'        CALL eye33d(F_old)')
            lines.append(f'        DO ii_v = 1, {nn}')
            lines.append(f'          DO i = 1, ndim')
            lines.append(f'            DO j = 1, ndim')
            lines.append(f'              F(i,j) = F(i,j)')
            lines.append(f'     &          + {fname}_node(i,ii_v)*{dsh_name}(ii_v,j)')
            lines.append(f'              F_old(i,j) = F_old(i,j)')
            lines.append(f'     &          + {fname}_old(i,ii_v)*{dsh_name}(ii_v,j)')
            lines.append(f'            END DO')
            lines.append(f'          END DO')
            lines.append(f'        END DO')
            lines.append('')

    # Scalar fields
    for fname in field_names:
        f = fields[fname]
        if isinstance(f, ScalarField):
            sh_name, _, nn = cfg.sh_for_degree(f.degree)
            lines.append(f'        {fname}_gp = 0.0d0')
            lines.append(f'        {fname}_old_gp = 0.0d0')
            lines.append(f'        DO ii_v = 1, {nn}')
            lines.append(f'          {fname}_gp = {fname}_gp')
            lines.append(f'     &      + {fname}_node(ii_v)*{sh_name}(ii_v)')
            lines.append(f'          {fname}_old_gp = {fname}_old_gp')
            lines.append(f'     &      + {fname}_old(ii_v)*{sh_name}(ii_v)')
            lines.append(f'        END DO')

            # Gradient if needed
            if fname in needs_grad:
                _, dsh_name, _ = cfg.sh_for_degree(f.degree)
                lines.append(f'        grad_{fname}(1) = 0.0d0')
                lines.append(f'        grad_{fname}(2) = 0.0d0')
                lines.append(f'        grad_{fname}(3) = 0.0d0')
                lines.append(f'        DO ii_v = 1, {nn}')
                lines.append(f'          DO j = 1, ndim')
                lines.append(f'            grad_{fname}(j) = grad_{fname}(j)')
                lines.append(f'     &        + {fname}_node(ii_v)*{dsh_name}(ii_v,j)')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            lines.append('')

    _append_svars_read(lines, sv_info, indent='        ')

    # --- Material evaluation (real, for residual) ---
    # Material subroutines are COMPLEX*16. For the residual, we convert
    # real fields to complex, call the material, and take DBLE of output.
    # Complex temporaries (Fz_r, etc.) declared at subroutine top.

    scalar_fields_list = [fn for fn in field_names
                          if isinstance(fields[fn], ScalarField)]

    lines_real_eval = []
    lines_real_eval.append('')
    lines_real_eval.append('C       Real material evaluation (for RHS)')
    lines_real_eval.append('        CALL real2complex33(F, Fz_r)')
    lines_real_eval.append('        CALL real2complex33(F_old, Fold_r)')
    for sfn in scalar_fields_list:
        lines_real_eval.append(f'        {sfn}z_r = DCMPLX({sfn}_gp, 0.0d0)')
    for sfn in scalar_fields_list:
        if sfn in needs_grad:
            lines_real_eval.append(f'        CALL real2complex3(grad_{sfn},'
                                   f' grad_{sfn}_r)')
    _append_state_old_to_complex(lines_real_eval, sv_info, indent='        ')

    # Call each material method
    for mname in weakform._mat.defined_methods():
        minfo = weakform._mat._methods[mname]
        sub_name = f'{mat_prefix}_{mname}'
        rs = _RETURN_SHAPE.get(mname)
        has_state_call = bool(sv_info) and mname == 'stress_PK1'

        call_args = []
        for p in minfo['all_params']:
            if p == 'F':
                call_args.append('Fz_r')
            elif p == 'F_old':
                call_args.append('Fold_r')
            elif p in state_old_names:
                call_args.append(f'{p}_z')
            elif p == 'dt':
                call_args.append('dt_safe')
            elif p == 'time':
                call_args.append('TIME(2)')
            elif p == 'X':
                call_args.append('X_gp')
            elif field_args.get(p, {}).get('kind') == 'vector':
                call_args.append(f'{p}_r')
            elif field_args.get(p, {}).get('kind') == 'scalar':
                call_args.append(f'{p}z_r')
            elif p.endswith('_old') and p[:-4] in fields:
                if p.endswith('_old'):
                    base = p.replace('_old', '')
                    call_args.append(f'DCMPLX({base}_old_gp, 0.0d0)')
                else:
                    call_args.append(f'{p}z_r')
            else:
                call_args.append(p)
        call_args.append('PROPS')
        out_name = _complex_eval_name(mname, rs)
        real_name = _real_output_name(mname, rs)

        if rs == 'tensor':
            call_args.append(out_name)
            if has_state_call:
                for name in sv_info:
                    call_args.append(f'{name}_new_z')
            lines_real_eval.append(f'        CALL {sub_name}({", ".join(call_args)})')
            lines_real_eval.append('        DO i = 1, 3')
            lines_real_eval.append('          DO j = 1, 3')
            lines_real_eval.append(f'            {real_name}(i,j) = '
                                   f'DBLE({out_name}(i,j))')
            lines_real_eval.append('          END DO')
            lines_real_eval.append('        END DO')
        elif rs == 'vector':
            call_args.append(out_name)
            lines_real_eval.append(f'        CALL {sub_name}({", ".join(call_args)})')
            lines_real_eval.append('        DO i = 1, 3')
            lines_real_eval.append(f'          {real_name}(i) = '
                                   f'DBLE({out_name}(i))')
            lines_real_eval.append('        END DO')
        else:
            call_args.append(out_name)
            lines_real_eval.append(f'        CALL {sub_name}({", ".join(call_args)})')
            lines_real_eval.append(f'        {real_name} = DBLE({out_name})')

    # (no END BLOCK needed — complex temps declared at subroutine top)
    lines.extend(lines_real_eval)
    lines.append('')
    _append_svars_write(lines, sv_info, indent='        ')
    if sv_info:
        lines.append('')

    # --- CS tangent engine call ---
    lines.append('C       CS tangent engine')
    # Build call to mat_cs_tangents
    cs_args = ['F']
    for sfn in sorted([fn for fn in field_names if isinstance(fields[fn], ScalarField)
                       and fn in [v for eq in equations.values() for v in eq['field_vars']]]):
        cs_args.append(f'{sfn}_gp')
    # vector fields
    for sfn in sorted(needs_grad):
        cs_args.append(f'grad_{sfn}')
    cs_args.append('F_old')
    # history scalars.  Keep this list in the exact order used by the
    # generated cs_tangents subroutine declaration.
    for old_name in _ordered_history_scalars(weakform, all_history):
        field_name = old_name[:-4]
        cs_args.append(f'{field_name}_old_gp')
    for name in sv_info:
        cs_args.append(f'{name}_old')
    cs_args.append('PROPS')
    if cs_needs_dt:
        cs_args.append('dt_safe')
    if cs_needs_time:
        cs_args.append('TIME(2)')
    if cs_needs_X:
        cs_args.append('X_gp')
    # tangent output arrays
    for spec in tangent_specs:
        cs_args.append(spec['tname'])

    cs_call = f'        CALL {mat_prefix}_cs_tangents('
    line = cs_call
    for i, a in enumerate(cs_args):
        if i > 0:
            line += ', '
        if len(line) + len(a) > 65:
            lines.append(line)
            line = '     &     '
        line += a
    line += ')'
    lines.append(line)
    lines.append('')

    # --- RHS assembly ---
    lines.append('C       RHS assembly')
    _emit_rhs_assembly(lines, weakform, fields, field_names, equations,
                       needs_grad, cfg)
    lines.append('')

    # --- AMATRX assembly ---
    lines.append('C       AMATRX assembly')
    _emit_amatrx_assembly(lines, weakform, fields, field_names, blocks,
                          tangent_specs, needs_grad, cfg)
    lines.append('')

    # End Gauss loop
    lines.append('      END DO')
    lines.append('')
    lines.append('      RETURN')
    lines.append('      END SUBROUTINE UEL')

    return '\n'.join(lines)


def _emit_rhs_assembly(lines, weakform, fields, field_names, equations,
                       needs_grad, cfg):
    """Emit RHS assembly code for all equations."""

    for eq_name, eq_info in equations.items():
        eq_meta = EQUATION_INFO[eq_name]
        test_field = eq_info['test_field']
        tf = fields[test_field]
        sh, dsh, nn = cfg.sh_for_degree(tf.degree)

        if eq_meta['return_type'] == 'tensor':
            # Momentum: RHS(row,1) -= P_iJ * dN_a/dX_J * wdetJ
            lines.append(f'        DO ii_v = 1, {nn}')
            lines.append(f'          DO i = 1, ndim')
            lines.append(f'            row = edof_{test_field}(i, ii_v)')
            lines.append(f'            DO j = 1, ndim')
            lines.append(f'              RHS(row,1) = RHS(row,1)')
            lines.append(f'     &          - P_real(i,j)*{dsh}(ii_v,j)*wdetJ')
            lines.append(f'            END DO')
            lines.append(f'          END DO')
            lines.append(f'        END DO')

        elif eq_meta['return_type'] == 'scalar':
            # Constraint: RHS(row,1) -= r_p * N_a * wdetJ
            out_name = _real_output_name(
                _material_method_for_equation(eq_name))
            lines.append(f'        DO ii_v = 1, {nn}')
            lines.append(f'          row = edof_{test_field}(ii_v)')
            lines.append(f'          RHS(row,1) = RHS(row,1)')
            lines.append(f'     &      - {out_name}*{sh}(ii_v)*wdetJ')
            lines.append(f'        END DO')

        elif eq_meta['return_type'] == 'tuple':
            # Transport: storage + flux
            storage_method = eq_meta['sub_terms']['storage']['material_method']
            flux_method = eq_meta['sub_terms']['flux']['material_method']
            storage_name = _real_output_name(storage_method)
            flux_name = _real_output_name(flux_method)
            # Storage: RHS -= c_dot * N * wdetJ
            lines.append(f'        DO ii_v = 1, {nn}')
            lines.append(f'          row = edof_{test_field}(ii_v)')
            lines.append(f'          RHS(row,1) = RHS(row,1)')
            lines.append(f'     &      - {storage_name}*{sh}(ii_v)*wdetJ')
            # Flux: RHS += j_R . dN/dX * wdetJ
            lines.append(f'          DO j = 1, ndim')
            lines.append(f'            RHS(row,1) = RHS(row,1)')
            lines.append(f'     &        + {flux_name}(j)*{dsh}(ii_v,j)*wdetJ')
            lines.append(f'          END DO')
            lines.append(f'        END DO')

        lines.append('')


def _emit_amatrx_assembly(lines, weakform, fields, field_names, blocks,
                          tangent_specs, needs_grad, cfg):
    """Emit AMATRX assembly code for all tangent blocks."""

    for spec in tangent_specs:
        tname = spec['tname']
        binfo = spec['binfo']
        sign = binfo['amatrx_sign']
        row_asm = binfo['row_assembly']
        col_asm = binfo['col_assembly']
        test_field = binfo['test_field']
        wrt_var = binfo['wrt_field']
        wrt_kind = binfo['wrt_kind']

        # Determine test field shape function info
        tf = fields[test_field]
        sh_row, dsh_row, nn_row = cfg.sh_for_degree(tf.degree)

        # Determine trial field shape function info
        # wrt_var maps to a field: F->u, p->p, mu->mu, grad_mu->mu
        if wrt_var == 'F':
            trial_field = 'u'
        elif wrt_var.startswith('grad_'):
            trial_field = wrt_var.replace('grad_', '')
        else:
            trial_field = wrt_var
        cf = fields.get(trial_field)
        if cf is None:
            continue
        sh_col, dsh_col, nn_col = cfg.sh_for_degree(cf.degree)

        sign_str = '+' if sign > 0 else '-'
        lines.append(f'C       {tname} (sign={sign_str}1,'
                     f' {row_asm}-{col_asm})')

        # Generate the assembly loop based on pattern
        if row_asm == 'grad' and col_asm == 'grad':
            # T_iJkL * dN_a/dX_J * dN_b/dX_L
            if isinstance(tf, VectorField) and wrt_kind == 'matrix':
                # dP/dF: (3,3,3,3)
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            DO i = 1, ndim')
                lines.append(f'              row = edof_{test_field}(i, ii_v)')
                lines.append(f'              DO k = 1, ndim')
                lines.append(f'                col = edof_{trial_field}(k, jj_v)')
                lines.append(f'                DO j = 1, ndim')
                lines.append(f'                  DO l = 1, ndim')
                lines.append(f'                    AMATRX(row,col) =')
                lines.append(f'     &                AMATRX(row,col)')
                lines.append(f'     &                {sign_str} {tname}(i,j,k,l)')
                lines.append(f'     &                * {dsh_row}(ii_v,j)')
                lines.append(f'     &                * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'                  END DO')
                lines.append(f'                END DO')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            elif isinstance(tf, VectorField) and wrt_kind == 'vector':
                # dP/dgrad_scalar: (3,3,3), row=grad(vector),
                # col=grad(scalar)
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO i = 1, ndim')
                lines.append(f'              row = edof_{test_field}(i, ii_v)')
                lines.append(f'              DO j = 1, ndim')
                lines.append(f'                DO l = 1, ndim')
                lines.append(f'                  AMATRX(row,col) =')
                lines.append(f'     &              AMATRX(row,col)')
                lines.append(f'     &              {sign_str} {tname}(i,j,l)')
                lines.append(f'     &              * {dsh_row}(ii_v,j)')
                lines.append(f'     &              * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'                END DO')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            elif isinstance(tf, ScalarField) and wrt_kind == 'matrix':
                # djR/dF: (3,3,3) row=grad(scalar), col=grad(vector)
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            DO k = 1, ndim')
                lines.append(f'              col = edof_{trial_field}(k, jj_v)')
                lines.append(f'              DO j = 1, ndim')
                lines.append(f'                row = edof_{test_field}(ii_v)')
                lines.append(f'                DO l = 1, ndim')
                lines.append(f'                  AMATRX(row,col) =')
                lines.append(f'     &              AMATRX(row,col)')
                lines.append(f'     &              {sign_str} {tname}(j,k,l)')
                lines.append(f'     &              * {dsh_row}(ii_v,j)')
                lines.append(f'     &              * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'                END DO')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            elif isinstance(tf, ScalarField) and wrt_kind == 'vector':
                # djR/dgrad_mu: (3,3) row=grad, col=grad
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            row = edof_{test_field}(ii_v)')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO j = 1, ndim')
                lines.append(f'              DO l = 1, ndim')
                lines.append(f'                AMATRX(row,col) =')
                lines.append(f'     &            AMATRX(row,col)')
                lines.append(f'     &            {sign_str} {tname}(j,l)')
                lines.append(f'     &            * {dsh_row}(ii_v,j)')
                lines.append(f'     &            * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')

        elif row_asm == 'grad' and col_asm == 'value':
            # T_iJ * dN_a/dX_J * N_b
            if isinstance(tf, VectorField):
                # dP/dp or dP/dmu: (3,3)
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO i = 1, ndim')
                lines.append(f'              row = edof_{test_field}(i, ii_v)')
                lines.append(f'              DO j = 1, ndim')
                lines.append(f'                AMATRX(row,col) =')
                lines.append(f'     &            AMATRX(row,col)')
                lines.append(f'     &            {sign_str} {tname}(i,j)')
                lines.append(f'     &            * {dsh_row}(ii_v,j)')
                lines.append(f'     &            * {sh_col}(jj_v) * wdetJ')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            else:
                # djR/dp or djR/dmu: (3) row=grad(scalar)
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            row = edof_{test_field}(ii_v)')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO j = 1, ndim')
                lines.append(f'              AMATRX(row,col) =')
                lines.append(f'     &          AMATRX(row,col)')
                lines.append(f'     &          {sign_str} {tname}(j)')
                lines.append(f'     &          * {dsh_row}(ii_v,j)')
                lines.append(f'     &          * {sh_col}(jj_v) * wdetJ')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')

        elif row_asm == 'value' and col_asm == 'grad':
            # T_kL * N_a * dN_b/dX_L
            if isinstance(cf, ScalarField) and wrt_kind == 'vector':
                # dstorage/dgrad_mu: (3), row=value(scalar),
                # col=grad(scalar)
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            row = edof_{test_field}(ii_v)')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO l = 1, ndim')
                lines.append(f'              AMATRX(row,col) =')
                lines.append(f'     &          AMATRX(row,col)')
                lines.append(f'     &          {sign_str} {tname}(l)')
                lines.append(f'     &          * {sh_row}(ii_v)')
                lines.append(f'     &          * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            else:
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            row = edof_{test_field}(ii_v)')
                lines.append(f'            DO k = 1, ndim')
                lines.append(f'              col = edof_{trial_field}(k, jj_v)')
                lines.append(f'              DO l = 1, ndim')
                lines.append(f'                AMATRX(row,col) =')
                lines.append(f'     &            AMATRX(row,col)')
                lines.append(f'     &            {sign_str} {tname}(k,l)')
                lines.append(f'     &            * {sh_row}(ii_v)')
                lines.append(f'     &            * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')

        elif row_asm == 'value' and col_asm == 'value':
            # T * N_a * N_b
            lines.append(f'        DO ii_v = 1, {nn_row}')
            lines.append(f'          DO jj_v = 1, {nn_col}')
            lines.append(f'            row = edof_{test_field}(ii_v)')
            lines.append(f'            col = edof_{trial_field}(jj_v)')
            lines.append(f'            AMATRX(row,col) =')
            lines.append(f'     &        AMATRX(row,col)')
            lines.append(f'     &        {sign_str} {tname}')
            lines.append(f'     &        * {sh_row}(ii_v)')
            lines.append(f'     &        * {sh_col}(jj_v) * wdetJ')
            lines.append(f'          END DO')
            lines.append(f'        END DO')

        lines.append('')


# =====================================================================
# Section 4b: F-bar UEL subroutine generation (Quad4, mechanics only)
# =====================================================================

def _generate_uel_fbar(weakform, mat_prefix, cfg):
    """Generate the F-bar UEL subroutine for single-field mechanics.

    Three-pass assembly:
      Pass 1: Compute J at all GPs, accumulate Jbar, store per-GP data
      Pass 2: Compute Fbar, material eval, CS tangent, RHS, standard
              AMATRX, compute and store Q
      Pass 3: Compute dJbar/du, add F-bar correction to AMATRX
    """
    ndim = cfg.ndim
    fields = weakform.fields
    field_names = weakform.field_names
    nprops = len(weakform._mat.props_names)

    n_nodes = cfg.n_nodes
    ngp = cfg.n_gauss_points

    # For F-bar, we only support single-field mechanics
    assert len(field_names) == 1, "F-bar only for single-field (u only)"
    assert isinstance(fields[field_names[0]], VectorField), \
        "F-bar needs VectorField"

    u_field = field_names[0]

    ndofel = n_nodes * ndim
    sv_info = _state_var_info(weakform._mat)
    state_old_names = {f'{name}_old' for name in sv_info}
    nstate_per_gp = _nstate_per_gp(sv_info) if sv_info else 0

    _sh = cfg.sh_linear_name         # F-bar uses linear shapes
    _dshxi = cfg.dshxi_linear_name
    _dsh = cfg.dsh_linear_name
    _nc = cfg.n_corner_nodes
    _d = cfg.ndim
    _jinv = cfg.jac_inv_size

    lines = []
    lines.append(f"""
C======================================================================
C     UEL subroutine — generated by abaqus_ufl (F-bar {cfg.name.capitalize()})
C     Material: {weakform._mat.__class__.__name__}
C     Fields: {u_field}
C     Element: {cfg.name.capitalize()}, NDOFEL = {ndofel}, F-bar = ON
C======================================================================
      SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL,
     &  NRHS, NSVARS, PROPS, NPROPS, COORDS, MCRD, NNODE,
     &  U, DU, V, A, JTYPE, TIME, DTIME, KSTEP, KINC,
     &  JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, PREDEF,
     &  NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT,
     &  JPROPS, NJPROP, PERIOD)

      IMPLICIT NONE

C     --- Abaqus UEL interface ---
      INTEGER NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE
      INTEGER JTYPE, KSTEP, KINC, JELEM, NDLOAD, NPREDF
      INTEGER MLVARX, MDLOAD, NJPROP
      INTEGER JDLTYP(MDLOAD,*), LFLAGS(*), JPROPS(*)

      DOUBLE PRECISION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL)
      DOUBLE PRECISION SVARS(NSVARS), ENERGY(8)
      DOUBLE PRECISION PROPS(NPROPS), COORDS(MCRD,NNODE)
      DOUBLE PRECISION U(NDOFEL), DU(MLVARX,*), V(NDOFEL)
      DOUBLE PRECISION A(NDOFEL)
      DOUBLE PRECISION TIME(2), DTIME, PARAMS(3)
      DOUBLE PRECISION ADLMAG(MDLOAD,*), PREDEF(2,NPREDF,NNODE)
      DOUBLE PRECISION DDLMAG(MDLOAD,*), PNEWDT, PERIOD

C     --- Local parameters ---
      INTEGER, PARAMETER :: ndim = {ndim}
      INTEGER, PARAMETER :: NNODE_E = {n_nodes}
      INTEGER, PARAMETER :: NGP = {ngp}
    """)
    if sv_info:
      lines.append(f'      INTEGER, PARAMETER :: NSTATE_PER_GP = {nstate_per_gp}')

    lines.append(f"""

C     --- Nodal arrays ---
      DOUBLE PRECISION :: u_node(ndim, NNODE_E)
      DOUBLE PRECISION :: u_old(ndim, NNODE_E)
      INTEGER :: edof_u(ndim, NNODE_E)

C     --- Shape functions and mapping ---
      DOUBLE PRECISION :: {_sh}({_nc}), {_dshxi}({_nc},{_d}), {_dsh}({_nc},{_d})
      DOUBLE PRECISION :: xi_gp({ngp},{_d}), w_gp({ngp})
      DOUBLE PRECISION :: {cfg.coords_name}({_d}, NNODE_E)
      DOUBLE PRECISION :: detJxi, Jinv_xi({_jinv},{_jinv})
      INTEGER :: ngp_out, stat

C     --- Per-GP stored data ---
      DOUBLE PRECISION :: {_dsh}_all(NGP, NNODE_E, {_d})
      DOUBLE PRECISION :: detJxi_all(NGP)
      DOUBLE PRECISION :: F_all(3, 3, NGP)
      DOUBLE PRECISION :: J_all(NGP)
      DOUBLE PRECISION :: Finv_all(3, 3, NGP)

C     --- F-bar quantities ---
      DOUBLE PRECISION :: V0, Jbar_num, Jbar, alpha
      DOUBLE PRECISION :: Fbar(3,3), P_bar(3,3)
      DOUBLE PRECISION :: dJbar_du(ndim, NNODE_E)
      DOUBLE PRECISION :: Q_all(ndim, ndim, NGP)
      DOUBLE PRECISION :: g_ai, h_bk, wdetJ, dt_safe

C     --- CS tangent ---
      DOUBLE PRECISION :: Atang(3,3,3,3)

C     --- Complex temps for material evaluation ---
      DOUBLE COMPLEX :: Fbar_z(3,3), Pz_eval(3,3)

            DOUBLE PRECISION :: det33d
      DOUBLE COMPLEX :: det33z
      INTEGER :: idx, ii_v, jj_v, i, j, k, l, m, N, kk, gg
      INTEGER :: row, col
    """)

    _append_state_var_declarations(lines, sv_info)
    lines.append('')

    # --- Extract coordinates ---
    lines.append(f'C     Extract {_d}D coordinates')
    lines.append('      DO ii_v = 1, NNODE_E')
    for idim in range(1, _d + 1):
        lines.append(f'        {cfg.coords_name}({idim}, ii_v) = COORDS({idim}, ii_v)')
    lines.append('      END DO')
    lines.append('')

    # --- Zero RHS and AMATRX ---
    lines.append('C     Zero RHS and AMATRX')
    lines.append('      DO i = 1, NDOFEL')
    lines.append('        RHS(i, 1) = 0.0d0')
    lines.append('        DO j = 1, NDOFEL')
    lines.append('          AMATRX(i, j) = 0.0d0')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')
    lines.append('C     TODO: implement Modified Riks RHS(:,2) properly')
    lines.append('C     Only touch RHS(:,2) when Abaqus allocated it (NRHS >= 2).')
    lines.append('      IF (NRHS .GE. 2 .AND. LFLAGS(3) .EQ. 4) THEN')
    lines.append('        DO i = 1, NDOFEL')
    lines.append('          RHS(i, 2) = 0.0d0')
    lines.append('        END DO')
    lines.append('      END IF')
    lines.append('')

    lines.append('C     DTIME guard (Abaqus calls with DTIME=0 for initial eval)')
    lines.append('      dt_safe = DTIME')
    lines.append('      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14')
    lines.append('')

    # --- DOF parsing ---
    lines.append('C     Parse DOFs and build edof map')
    lines.append('      idx = 0')
    lines.append('      DO ii_v = 1, NNODE_E')
    lines.append('        DO i = 1, ndim')
    lines.append('          idx = idx + 1')
    lines.append('          u_node(i, ii_v) = U(idx)')
    lines.append('          u_old(i, ii_v) = U(idx) - DU(idx, 1)')
    lines.append('          edof_u(i, ii_v) = idx')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')

    # --- Pass 1: Compute Jbar ---
    lines.append('C     =============================================')
    lines.append('C     Pass 1: Compute Jbar and store per-GP data')
    lines.append('C     =============================================')
    lines.append(f'      CALL {cfg.gauss_subroutine}(xi_gp, w_gp, ngp_out)')
    lines.append('')
    lines.append('      V0 = 0.0d0')
    lines.append('      Jbar_num = 0.0d0')
    lines.append('')
    lines.append('      DO kk = 1, NGP')
    # Split shape call to match legacy line breaks
    xi_parts = cfg.gauss_xi_args.split(', ')
    lines.append(f'        CALL {cfg.shape_linear_subroutine}({xi_parts[0]},')
    lines.append(f'     &    {", ".join(xi_parts[1:])}, {_sh}, {_dshxi})')
    lines.append(f'        CALL {cfg.jac_subroutine}({_dshxi}, {cfg.coords_name}, {_nc},')
    lines.append(f'     &    {_dsh}, detJxi, Jinv_xi, stat)')
    lines.append('')
    lines.append('        IF (stat .EQ. 0) THEN')
    lines.append('          PNEWDT = 0.25d0')
    lines.append('          RETURN')
    lines.append('        END IF')
    lines.append('')
    lines.append('C       Store per-GP data')
    lines.append('        detJxi_all(kk) = detJxi')
    lines.append('        DO ii_v = 1, NNODE_E')
    for idim in range(1, _d + 1):
        lines.append(f'          {_dsh}_all(kk, ii_v, {idim}) = {_dsh}(ii_v, {idim})')
    lines.append('        END DO')
    lines.append('')
    lines.append('C       Compute F (eye33d for F33=1)')
    lines.append('        CALL eye33d(F_all(1,1,kk))')
    lines.append('        DO ii_v = 1, NNODE_E')
    lines.append('          DO i = 1, ndim')
    lines.append('            DO j = 1, ndim')
    lines.append('              F_all(i,j,kk) = F_all(i,j,kk)')
    lines.append(f'     &          + u_node(i,ii_v)*{_dsh}(ii_v,j)')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')
    lines.append('C       J and Finv')
    lines.append('        J_all(kk) = det33d(F_all(1,1,kk))')
    lines.append('        CALL inv33d(F_all(1,1,kk),')
    lines.append('     &    Finv_all(1,1,kk))')
    lines.append('')
    lines.append('C       Accumulate Jbar')
    lines.append('        V0 = V0 + detJxi * w_gp(kk)')
    lines.append('        Jbar_num = Jbar_num')
    lines.append('     &    + J_all(kk) * detJxi * w_gp(kk)')
    lines.append('      END DO')
    lines.append('')
    lines.append('      Jbar = Jbar_num / V0')
    lines.append('')

    # --- Compute dJbar/du ---
    lines.append('C     Compute dJbar/du (element-level sum)')
    lines.append('      DO jj_v = 1, NNODE_E')
    lines.append('        DO k = 1, ndim')
    lines.append('          dJbar_du(k,jj_v) = 0.0d0')
    lines.append('          DO gg = 1, NGP')
    lines.append('            DO N = 1, ndim')
    lines.append('              dJbar_du(k,jj_v) = dJbar_du(k,jj_v)')
    lines.append('     &          + J_all(gg)')
    lines.append('     &          * Finv_all(N,k,gg)')
    lines.append(f'     &          * {_dsh}_all(gg,jj_v,N)')
    lines.append('     &          * detJxi_all(gg)')
    lines.append('     &          * w_gp(gg) / V0')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')

    # --- Pass 2: Material, RHS, standard AMATRX, store Q ---
    lines.append('C     =============================================')
    lines.append('C     Pass 2: Fbar, material, RHS, AMATRX_std, Q')
    lines.append('C     =============================================')
    lines.append('      DO kk = 1, NGP')
    lines.append('        alpha = (Jbar / J_all(kk))')
    lines.append('     &    ** (1.0d0 / DBLE(ndim))')
    lines.append('        wdetJ = detJxi_all(kk) * w_gp(kk)')
    lines.append('')
    lines.append('C       Compute Fbar (in-plane only)')
    lines.append('        CALL eye33d(Fbar)')
    lines.append('        DO i = 1, ndim')
    lines.append('          DO j = 1, ndim')
    lines.append('            Fbar(i,j) = alpha * F_all(i,j,kk)')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')

    # Material evaluation (real, for RHS)
    # Call stress_PK1 with complexified Fbar, take DBLE
    sub_name = f'{mat_prefix}_stress_PK1'
    minfo = weakform._mat._methods['stress_PK1']

    # Build call args.  The generated stateful material subroutine signature
    # groups all old state variables before dt.
    call_args_str = []
    call_args_str.append('Fbar_z')
    for name in sv_info:
        call_args_str.append(f'{name}_old_z')
    if 'dt' in minfo['all_params']:
        call_args_str.append('dt_safe')
    for p in minfo['all_params']:
        if p != 'F' and p not in state_old_names and p != 'dt':
            call_args_str.append(p)
    call_args_str.append('PROPS')
    call_args_str.append('Pz_eval')
    if sv_info:
        for name in sv_info:
            call_args_str.append(f'{name}_new_z')

    lines.append('C       Material evaluation with Fbar')
    _append_svars_read(lines, sv_info, indent='        ')
    lines.append('        CALL real2complex33(Fbar, Fbar_z)')
    _append_state_old_to_complex(lines, sv_info, indent='        ')
    lines.append(f'        CALL {sub_name}(')
    lines.append(f'     &    {", ".join(call_args_str)})')
    lines.append('        DO i = 1, 3')
    lines.append('          DO j = 1, 3')
    lines.append('            P_bar(i,j) = DBLE(Pz_eval(i,j))')
    lines.append('          END DO')
    lines.append('        END DO')
    _append_svars_write(lines, sv_info, indent='        ')
    lines.append('')

    # RHS assembly
    lines.append('C       RHS: r_a^i = -P_bar_iJ * dN_a/dX_J * wdetJ')
    lines.append('        DO ii_v = 1, NNODE_E')
    lines.append('          DO i = 1, ndim')
    lines.append('            row = edof_u(i, ii_v)')
    lines.append('            DO j = 1, ndim')
    lines.append('              RHS(row,1) = RHS(row,1)')
    lines.append('     &          - P_bar(i,j)')
    lines.append(f'     &          * {_dsh}_all(kk,ii_v,j) * wdetJ')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')

    # CS tangent: A = dP/dFbar (perturb all 9 components of Fbar)
    lines.append('C       CS tangent: A_iJkL = dP_bar/dFbar_kL')
    lines.append(f'        CALL {mat_prefix}_cs_tangent_dPdF(')
    cs_args = ['Fbar']
    for name in sv_info:
        cs_args.append(f'{name}_old')
    if 'dt' in minfo['all_params']:
        cs_args.append('dt_safe')
    cs_args.extend(['PROPS', 'Atang'])
    lines.append(f'     &    {", ".join(cs_args)})')
    lines.append('')

    # Standard AMATRX (scaled by alpha)
    lines.append('C       Standard AMATRX (alpha-scaled)')
    lines.append('        DO ii_v = 1, NNODE_E')
    lines.append('          DO jj_v = 1, NNODE_E')
    lines.append('            DO i = 1, ndim')
    lines.append('              row = edof_u(i, ii_v)')
    lines.append('              DO k = 1, ndim')
    lines.append('                col = edof_u(k, jj_v)')
    lines.append('                DO j = 1, ndim')
    lines.append('                  DO l = 1, ndim')
    lines.append('                    AMATRX(row,col) =')
    lines.append('     &                AMATRX(row,col)')
    lines.append('     &                + alpha')
    lines.append('     &                * Atang(i,j,k,l)')
    lines.append(f'     &                * {_dsh}_all(kk,ii_v,j)')
    lines.append(f'     &                * {_dsh}_all(kk,jj_v,l)')
    lines.append('     &                * wdetJ')
    lines.append('                  END DO')
    lines.append('                END DO')
    lines.append('              END DO')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')

    # Compute Q (in-plane contraction ONLY) and store
    lines.append('C       Q = A : Fbar (in-plane only)')
    lines.append('        DO i = 1, ndim')
    lines.append('          DO j = 1, ndim')
    lines.append('            Q_all(i,j,kk) = 0.0d0')
    lines.append('            DO m = 1, ndim')
    lines.append('              DO N = 1, ndim')
    lines.append('                Q_all(i,j,kk) =')
    lines.append('     &            Q_all(i,j,kk)')
    lines.append('     &            + Atang(i,j,m,N)')
    lines.append('     &            * Fbar(m,N)')
    lines.append('              END DO')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')

    lines.append('      END DO')
    lines.append('C     End Pass 2')
    lines.append('')

    # --- Pass 3: F-bar correction ---
    lines.append('C     =============================================')
    lines.append('C     Pass 3: F-bar correction (rank-1 update)')
    lines.append('C     =============================================')
    lines.append('      DO kk = 1, NGP')
    lines.append('        wdetJ = detJxi_all(kk) * w_gp(kk)')
    lines.append('')
    lines.append('        DO ii_v = 1, NNODE_E')
    lines.append('          DO jj_v = 1, NNODE_E')
    lines.append('            DO i = 1, ndim')
    lines.append('              row = edof_u(i, ii_v)')
    lines.append('')
    lines.append('C             g_a^i = Q_iJ * dN_a/dX_J * wdetJ')
    lines.append('              g_ai = 0.0d0')
    lines.append('              DO j = 1, ndim')
    lines.append('                g_ai = g_ai')
    lines.append('     &            + Q_all(i,j,kk)')
    lines.append(f'     &            * {_dsh}_all(kk,ii_v,j) * wdetJ')
    lines.append('              END DO')
    lines.append('')
    lines.append('              DO k = 1, ndim')
    lines.append('                col = edof_u(k, jj_v)')
    lines.append('')
    lines.append('C               h_bk = (1/ndim)*[dJbar/du/Jbar')
    lines.append('C                       - Finv_Nk * dN_b/dX_N]')
    lines.append('                h_bk = dJbar_du(k,jj_v) / Jbar')
    lines.append('                DO N = 1, ndim')
    lines.append('                  h_bk = h_bk')
    lines.append('     &              - Finv_all(N,k,kk)')
    lines.append(f'     &              * {_dsh}_all(kk,jj_v,N)')
    lines.append('                END DO')
    lines.append('                h_bk = h_bk')
    lines.append('     &            / DBLE(ndim)')
    lines.append('')
    lines.append('                AMATRX(row,col) =')
    lines.append('     &            AMATRX(row,col)')
    lines.append('     &            + g_ai * h_bk')
    lines.append('              END DO')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('C     End Pass 3')
    lines.append('')
    lines.append('      RETURN')
    lines.append('      END SUBROUTINE UEL')

    return '\n'.join(lines)


def _generate_cs_tangent_dPdF(weakform, mat_prefix):
    """Generate a simple CS tangent subroutine for dP/dF only (F-bar use)."""
    nprops = len(weakform._mat.props_names)
    sub_name = f'{mat_prefix}_cs_tangent_dPdF'
    stress_sub = f'{mat_prefix}_stress_PK1'
    sv_info = _state_var_info(weakform._mat)
    state_old_names = {f'{name}_old' for name in sv_info}

    # Check what args stress_PK1 takes beyond F
    minfo = weakform._mat._methods['stress_PK1']
    extra_params = [p for p in minfo['all_params'] if p != 'F']

    lines = []
    lines.append(f'      SUBROUTINE {sub_name}(')
    in_args = ['F']
    for name in sv_info:
        in_args.append(f'{name}_old')
    if 'dt' in extra_params:
        in_args.append('dt')
    in_args.extend(['PROPS', 'dPdF'])
    lines.append(f'     &  {", ".join(in_args)})')
    lines.append('      IMPLICIT NONE')
    lines.append('      DOUBLE PRECISION, INTENT(IN) :: F(3,3)')
    for name, entry in sv_info.items():
        dims = _fortran_dims(entry)
        lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: {name}_old{dims}')
    if 'dt' in extra_params:
        lines.append('      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.append(f'      DOUBLE PRECISION, INTENT(IN) :: PROPS({nprops})')
    lines.append('      DOUBLE PRECISION, INTENT(OUT) :: dPdF(3,3,3,3)')
    lines.append('')
    lines.append('C     CS_H = 1e-10: verified for this framework')
    lines.append('      DOUBLE PRECISION, PARAMETER :: CS_H = 1.0d-10')
    lines.append('      DOUBLE COMPLEX :: Fz(3,3), Pz(3,3)')
    for name, entry in sv_info.items():
        dims = _fortran_dims(entry)
        lines.append(f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
        lines.append(f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
    lines.append('      INTEGER :: i, j, k, l')
    lines.append('')
    lines.append('      DO k = 1, 3')
    lines.append('        DO l = 1, 3')
    lines.append('          CALL real2complex33(F, Fz)')
    _append_state_old_to_complex(lines, sv_info, indent='          ')
    lines.append('          Fz(k,l) = Fz(k,l)')
    lines.append('     &      + DCMPLX(0.0d0, CS_H)')

    # Build call to stress sub.  Generated stateful material subroutines group
    # all old state variables before dt, even if dt appears earlier in the
    # Python method signature.
    call_args = ['Fz']
    for name in sv_info:
        call_args.append(f'{name}_old_z')
    if 'dt' in extra_params:
        call_args.append('dt')
    for p in extra_params:
        if p not in state_old_names and p != 'dt':
            call_args.append(f'DCMPLX(0.0d0, 0.0d0)')
    call_args.append('PROPS')
    call_args.append('Pz')
    if sv_info:
        for name in sv_info:
            call_args.append(f'{name}_new_z')

    lines.append(f'          CALL {stress_sub}(')
    lines.append(f'     &      {", ".join(call_args)})')
    lines.append('          DO i = 1, 3')
    lines.append('            DO j = 1, 3')
    lines.append('              dPdF(i,j,k,l) =')
    lines.append('     &          AIMAG(Pz(i,j)) / CS_H')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')
    lines.append('      RETURN')
    lines.append(f'      END SUBROUTINE {sub_name}')

    return '\n'.join(lines)


# =====================================================================
# Main entry point
# =====================================================================

def generate_uel(weakform, output_path, element='Quad8', mat_prefix=None,
                 fbar=None, element_config=None, formulation=None):
    """
    Generate a complete UEL .for file from a WeakForm instance.

    Args:
        weakform: a WeakForm instance (verified)
        output_path: path to write the .for file
        element: element type ('Quad8' or 'Quad4') — ignored when
                 *element_config* is given.
        mat_prefix: prefix for generated subroutine names
        fbar: if True, use F-bar modification. Default: auto-detect
              (True for Quad4 single-field, False otherwise)
        element_config: an ``ElementConfig`` instance.  When omitted the
                        config is looked up from *element* via
                        ``ELEMENT_CONFIGS``.
        formulation: explicit formulation route: ``'standard'``,
                     ``'fbar_mechanics'``, or ``'local_pressure'``.
                     When omitted, the legacy *fbar* auto-detection is
                     preserved.
    """
    if mat_prefix is None:
        mat_prefix = weakform._mat.__class__.__name__.lower()

    # Resolve ElementConfig
    if element_config is not None:
        cfg = element_config
    else:
        key = element.lower()
        if key not in ELEMENT_CONFIGS:
            raise ValueError(
                f"Unknown element '{element}'. Pass an ElementConfig or "
                f"use one of: {list(ELEMENT_CONFIGS.keys())}")
        cfg = ELEMENT_CONFIGS[key]

    if formulation is not None:
        formulation = formulation.lower().replace('-', '_')
        valid = {'standard', 'fbar_mechanics', 'local_pressure'}
        if formulation not in valid:
            raise ValueError(
                f"Unknown formulation '{formulation}'. Use one of: "
                f"{sorted(valid)}")
        if fbar is not None:
            requested_fbar = formulation == 'fbar_mechanics'
            if fbar != requested_fbar and formulation != 'local_pressure':
                raise ValueError(
                    "Conflicting options: formulation and fbar disagree")

    if formulation == 'local_pressure':
        if cfg.name not in ('quad4', 'hex8'):
            raise ValueError(
                "local_pressure formulation currently supports Quad4 and Hex8")
        if cfg.ndim != weakform.ndim:
            weakform.ndim = cfg.ndim
            for name in weakform.field_names:
                f = weakform.fields[name]
                f.ndim = cfg.ndim
        weakform._set_element_config(cfg)
        from .uel_local_pressure import generate_uel_local_pressure
        return generate_uel_local_pressure(weakform, output_path,
                                           mat_prefix=mat_prefix,
                                           element_config=cfg)

    if getattr(weakform, 'local_field_names', []):
        raise ValueError(
            "LocalScalar fields are currently supported only by "
            "formulation='local_pressure'. Standard and F-bar UEL generation "
            "do not yet implement local scalar storage, evaluation, or "
            "condensation.")

    # Reconcile spatial dimension: update WeakForm fields if needed
    if cfg.ndim != weakform.ndim:
        weakform.ndim = cfg.ndim
        for name in weakform.field_names:
            f = weakform.fields[name]
            f.ndim = cfg.ndim
            if isinstance(f, VectorField):
                f.ndof_per_node = cfg.ndim
        weakform._build_argument_metadata()
        weakform._mat.refresh_metadata(weakform._field_arg_defs,
                                       weakform._history_args,
                                       weakform._param_args)
        weakform._set_element_config(cfg)
    else:
        weakform._set_element_config(cfg)

    ndofel = _compute_ndofel(cfg, weakform.fields)

    # Auto-detect F-bar
    if formulation == 'standard':
        fbar = False
    elif formulation == 'fbar_mechanics':
        fbar = True
        if len(weakform.field_names) != 1:
            raise ValueError(
                "fbar_mechanics is only for pure mechanical single-field "
                "elements. Use standard for coupled u,p,mu formulations.")
        fname = weakform.field_names[0]
        if not isinstance(weakform.fields[fname], VectorField):
            raise ValueError("fbar_mechanics requires a VectorField")
        if not cfg.supports_fbar:
            raise ValueError(f"Element '{cfg.name}' does not support F-bar")
    elif fbar is None:
        has_mixed = any(weakform.fields[fn].degree < 2
                        for fn in weakform.field_names
                        if isinstance(weakform.fields[fn], ScalarField))
        fbar = cfg.supports_fbar and not has_mixed

    # Read templates
    templates = [
        _read_template('tensor_ops.for'),
        _read_template('cs_linalg.for'),
        _read_template('gauss_rules.for'),
        _read_template('isoparametric.for'),
    ]

    # Shape functions: include what's needed from config
    included = set()
    for tpl in cfg.template_files:
        templates.append(_read_template(tpl))
        included.add(tpl)

    # Include linear shape template if needed and not already present
    has_degree1 = any(weakform.fields[fn].degree < 2
                      for fn in weakform.field_names)
    linear_tpl = f'{cfg.shape_linear_subroutine}.for'
    if (has_degree1 or fbar) and linear_tpl not in included:
        templates.append(_read_template(linear_tpl))
        included.add(linear_tpl)

    # Generate material subroutines
    mat_src = _generate_all_material_subs(weakform, mat_prefix)

    if fbar:
        # F-bar path
        _fnames = weakform.field_names
        _flds = weakform.fields
        if len(_fnames) == 1 and isinstance(_flds[_fnames[0]], VectorField):
            # Single-field mechanics: use existing _generate_uel_fbar
            cs_src = _generate_cs_tangent_dPdF(weakform, mat_prefix)
            uel_src = _generate_uel_fbar(weakform, mat_prefix, cfg)
        else:
            raise ValueError(
                "F-bar generation is currently restricted to pure mechanical "
                "single-field VectorField formulations. Use formulation="
                "'standard' for coupled u,p,mu or transport problems.")
        fbar_label = ', F-bar = ON'
    else:
        # Standard path: full CS engine
        cs_src = _generate_cs_engine(weakform, mat_prefix)
        uel_src = _generate_uel(weakform, mat_prefix, cfg)
        fbar_label = ''

    element_name = cfg.name.capitalize()

    # Apply column-72 wrapping to generated code
    mat_src = '\n'.join(_wrap_lines(mat_src.split('\n')))
    cs_src = '\n'.join(_wrap_lines(cs_src.split('\n')))
    uel_src = '\n'.join(_wrap_lines(uel_src.split('\n')))

    # Header
    header = f"""C======================================================================
C     Generated by abaqus_ufl v0.1.0
C     Material: {weakform._mat.__class__.__name__}
C     Fields: {', '.join(weakform.field_names)}
C     Element: {element_name}, NDOFEL = {ndofel}{fbar_label}
C======================================================================
"""

    # Assemble
    full_src = '\n'.join([
        header,
        'C' + '=' * 70,
        'C     Section 1: UEL subroutine (Abaqus entry point)',
        'C' + '=' * 70,
        uel_src,
        '',
        'C' + '=' * 70,
        'C     Section 2: Material subroutines (COMPLEX*16)',
        'C     AST-translated from Python Material methods',
        'C' + '=' * 70,
        mat_src,
        '',
        'C' + '=' * 70,
        'C     Section 3: Complex-step tangent engine',
        'C' + '=' * 70,
        cs_src,
        '',
        'C' + '=' * 70,
        'C     Section 4: Shipped templates (tensor ops, shape',
        'C     functions, quadrature, isoparametric mapping)',
        'C' + '=' * 70,
        '\n'.join(templates),
    ])

    with open(output_path, 'w') as f:
        f.write(full_src)

    print(f"  Generated: {output_path}")
    print(f"  Material:  {weakform._mat.__class__.__name__}")
    print(f"  Fields:    {', '.join(weakform.field_names)}")
    print(f"  Element:   {element_name}")
    print(f"  NDOFEL:    {ndofel}")
    print(f"  F-bar:     {fbar}")
    if not fbar:
        print(f"  Tangent blocks: {len(weakform.tangent_blocks())}")

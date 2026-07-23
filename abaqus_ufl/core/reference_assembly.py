"""
Python reference assembly for UEL verification.

Mirrors the generated Fortran UEL exactly:
  - Same DOF layout and edof mapping
  - Same shape functions (Quad8, Quad4)
  - Same isoparametric mapping with shared Jinv
  - Same F = I + grad(u), F_old = I + grad(u_old)
  - Same material evaluation
  - Same CS tangent engine
  - Same RHS and AMATRX assembly with correct signs
  - Same 3x3 Gauss quadrature

Usage:
    from reference_assembly import assemble_element
    RHS, AMATRX = assemble_element(problem, coords, U, DU, DTIME, PROPS)
    RHS, AMATRX = assemble_element(problem, coords4, U, DU, DTIME, PROPS,
                                   STATEV=svars, element="quad4")
    RHS, AMATRX = assemble_element(problem, coords4, U, DU, DTIME, PROPS,
                                   TIME=[step_time, total_time])
"""

import numpy as np
from abaqus_ufl.core.fields import VectorField, ScalarField
from abaqus_ufl.core.defs import FIELD_ARGS
from abaqus_ufl.core.weakform import EQUATION_INFO
from abaqus_ufl.generators.element_config import ELEMENT_CONFIGS


# =====================================================================
# Shape functions (exact match to Fortran templates)
# =====================================================================

def shape_quad8(xi, eta):
    """8-node serendipity shape functions and derivatives."""
    sh = np.zeros(8)
    dshxi = np.zeros((8, 2))

    # Corners
    sh[0] = 0.25*(1-xi)*(1-eta)*(-xi-eta-1)
    dshxi[0, 0] = 0.25*(1-eta)*(2*xi+eta)
    dshxi[0, 1] = 0.25*(1-xi)*(xi+2*eta)

    sh[1] = 0.25*(1+xi)*(1-eta)*(xi-eta-1)
    dshxi[1, 0] = 0.25*(1-eta)*(2*xi-eta)
    dshxi[1, 1] = 0.25*(1+xi)*(-xi+2*eta)

    sh[2] = 0.25*(1+xi)*(1+eta)*(xi+eta-1)
    dshxi[2, 0] = 0.25*(1+eta)*(2*xi+eta)
    dshxi[2, 1] = 0.25*(1+xi)*(xi+2*eta)

    sh[3] = 0.25*(1-xi)*(1+eta)*(-xi+eta-1)
    dshxi[3, 0] = 0.25*(1+eta)*(2*xi-eta)
    dshxi[3, 1] = 0.25*(1-xi)*(-xi+2*eta)

    # Midsides
    sh[4] = 0.5*(1-xi*xi)*(1-eta)
    dshxi[4, 0] = -xi*(1-eta)
    dshxi[4, 1] = -0.5*(1-xi*xi)

    sh[5] = 0.5*(1+xi)*(1-eta*eta)
    dshxi[5, 0] = 0.5*(1-eta*eta)
    dshxi[5, 1] = -(1+xi)*eta

    sh[6] = 0.5*(1-xi*xi)*(1+eta)
    dshxi[6, 0] = -xi*(1+eta)
    dshxi[6, 1] = 0.5*(1-xi*xi)

    sh[7] = 0.5*(1-xi)*(1-eta*eta)
    dshxi[7, 0] = -0.5*(1-eta*eta)
    dshxi[7, 1] = -(1-xi)*eta

    return sh, dshxi


def shape_quad4(xi, eta):
    """4-node bilinear shape functions and derivatives."""
    sh = np.zeros(4)
    dshxi = np.zeros((4, 2))

    sh[0] = 0.25*(1-xi)*(1-eta)
    dshxi[0, 0] = -0.25*(1-eta)
    dshxi[0, 1] = -0.25*(1-xi)

    sh[1] = 0.25*(1+xi)*(1-eta)
    dshxi[1, 0] = 0.25*(1-eta)
    dshxi[1, 1] = -0.25*(1+xi)

    sh[2] = 0.25*(1+xi)*(1+eta)
    dshxi[2, 0] = 0.25*(1+eta)
    dshxi[2, 1] = 0.25*(1+xi)

    sh[3] = 0.25*(1-xi)*(1+eta)
    dshxi[3, 0] = -0.25*(1+eta)
    dshxi[3, 1] = 0.25*(1-xi)

    return sh, dshxi


def gauss_2d_3x3():
    """3x3 Gauss quadrature points and weights."""
    g = np.sqrt(3.0/5.0)
    pts_1d = [-g, 0.0, g]
    wts_1d = [5.0/9.0, 8.0/9.0, 5.0/9.0]

    xi_gp = np.zeros((9, 2))
    w_gp = np.zeros(9)
    idx = 0
    for j in range(3):
        for i in range(3):
            xi_gp[idx, 0] = pts_1d[i]
            xi_gp[idx, 1] = pts_1d[j]
            w_gp[idx] = wts_1d[i] * wts_1d[j]
            idx += 1
    return xi_gp, w_gp


def gauss_2d_2x2():
    """2x2 Gauss quadrature points and weights."""
    g = 1.0 / np.sqrt(3.0)
    xi_gp = np.array([
        [-g, -g],
        [ g, -g],
        [ g,  g],
        [-g,  g],
    ])
    w_gp = np.ones(4)
    return xi_gp, w_gp


def map_grad_2d(dshxi, coords, nNode):
    """Map shape function derivatives from (xi,eta) to (X,Y).

    Returns dsh, detJ, Jinv, stat.
    coords: (2, nNode)
    """
    # Jacobian J(i,j) = sum_a dN_a/dxi_i * X_j^a
    Jac = np.zeros((2, 2))
    for i in range(2):
        for k in range(nNode):
            Jac[i, 0] += dshxi[k, i] * coords[0, k]
            Jac[i, 1] += dshxi[k, i] * coords[1, k]

    detJ = Jac[0, 0]*Jac[1, 1] - Jac[0, 1]*Jac[1, 0]
    if detJ <= 0:
        return None, detJ, None, 0

    Jinv = np.array([[ Jac[1, 1], -Jac[0, 1]],
                     [-Jac[1, 0],  Jac[0, 0]]]) / detJ

    dsh = np.zeros((nNode, 2))
    for k in range(nNode):
        dsh[k, 0] = Jinv[0, 0]*dshxi[k, 0] + Jinv[0, 1]*dshxi[k, 1]
        dsh[k, 1] = Jinv[1, 0]*dshxi[k, 0] + Jinv[1, 1]*dshxi[k, 1]

    return dsh, detJ, Jinv, 1


def apply_jinv(dshxi, Jinv, nNode):
    """Map dshxi using pre-computed Jinv."""
    dsh = np.zeros((nNode, 2))
    for k in range(nNode):
        dsh[k, 0] = Jinv[0, 0]*dshxi[k, 0] + Jinv[0, 1]*dshxi[k, 1]
        dsh[k, 1] = Jinv[1, 0]*dshxi[k, 0] + Jinv[1, 1]*dshxi[k, 1]
    return dsh


# =====================================================================
# CS tangent computation (uses the Python tangent engine)
# =====================================================================

CS_H = 1e-10


def _compute_cs_tangent(func, args, wrt_name, wrt_kind, ret_shape):
    """Compute one CS tangent block.

    func: callable material method
    args: dict of argument name -> value
    wrt_name: which argument to perturb
    wrt_kind: 'matrix', 'scalar', or 'vector'
    ret_shape: 'tensor', 'scalar', or 'vector'
    """
    def _complexify(args_dict, perturb_name=None, perturb_idx=None):
        """Convert all args to complex, optionally perturbing one."""
        cargs = {}
        for k, v in args_dict.items():
            if k in {'dt', 'time'}:
                cargs[k] = v
            elif isinstance(v, np.ndarray):
                cargs[k] = v.astype(complex)
            elif isinstance(v, (int, float)):
                cargs[k] = complex(v)
            else:
                cargs[k] = v
        if perturb_name is not None:
            val = cargs[perturb_name]
            if perturb_idx is not None:
                val[perturb_idx] = val[perturb_idx] + 1j * CS_H
            else:
                cargs[perturb_name] = val + 1j * CS_H
        return cargs

    if wrt_kind == 'matrix':
        # F perturbation: 9 components
        if ret_shape == 'tensor':
            tang = np.zeros((3, 3, 3, 3))
        elif ret_shape == 'scalar':
            tang = np.zeros((3, 3))
        elif ret_shape == 'vector':
            tang = np.zeros((3, 3, 3))

        for k in range(3):
            for l in range(3):
                cargs = _complexify(args)
                cargs[wrt_name][k, l] += 1j * CS_H
                result = func(**cargs)
                if isinstance(result, tuple):
                    result = result[0]
                if ret_shape == 'tensor':
                    tang[:, :, k, l] = result.imag / CS_H
                elif ret_shape == 'scalar':
                    tang[k, l] = result.imag / CS_H
                elif ret_shape == 'vector':
                    tang[:, k, l] = result.imag / CS_H
        return tang

    elif wrt_kind == 'scalar':
        cargs = _complexify(args, wrt_name)
        result = func(**cargs)
        if isinstance(result, tuple):
            result = result[0]
        if ret_shape == 'tensor':
            return result.imag / CS_H
        elif ret_shape == 'scalar':
            return result.imag / CS_H
        elif ret_shape == 'vector':
            return result.imag / CS_H

    elif wrt_kind == 'vector':
        if ret_shape == 'tensor':
            tang = np.zeros((3, 3, 3))
        elif ret_shape == 'scalar':
            tang = np.zeros(3)
        elif ret_shape == 'vector':
            tang = np.zeros((3, 3))

        for k in range(3):
            cargs = _complexify(args)
            cargs[wrt_name][k] += 1j * CS_H
            result = func(**cargs)
            if isinstance(result, tuple):
                result = result[0]
            if ret_shape == 'tensor':
                tang[:, :, k] = result.imag / CS_H
            elif ret_shape == 'scalar':
                tang[k] = result.imag / CS_H
            elif ret_shape == 'vector':
                tang[:, k] = result.imag / CS_H
        return tang


# =====================================================================
# Main assembly function
# =====================================================================

# Map material method names to return shapes
_RET_SHAPE = {
    'stress_PK1': 'tensor',
    'pressure_resid': 'scalar',
    'solvent_flux': 'vector',
    'solvent_storage': 'scalar',
    'phase_flux': 'vector',
    'phase_storage': 'scalar',
    'pressure_flux': 'vector',
    'pressure_storage': 'scalar',
    'species_flux': 'vector',
    'species_storage': 'scalar',
}


def _material_method_for_equation(eq_name, binfo=None):
    if binfo is not None and binfo.get('material_method'):
        return binfo['material_method']
    if eq_name == 'momentum_equation':
        return 'stress_PK1'
    if eq_name == 'pressure_equation':
        return 'pressure_resid'
    return None


def _state_var_sizes(material):
    """Return ordered (name, size, init) tuples for material state vars."""
    info = []
    for name, init in getattr(material, 'state_vars', {}).items():
        arr = np.asarray(init)
        size = 1 if np.isscalar(init) or arr.shape == () else arr.size
        info.append((name, size, init))
    return info


def _state_arg_defaults(material):
    """Default state_old arguments from material.state_vars."""
    state = {}
    for name, _, init in _state_var_sizes(material):
        if np.isscalar(init) or np.asarray(init).shape == ():
            state[f'{name}_old'] = float(init)
        else:
            state[f'{name}_old'] = np.asarray(init, dtype=float).copy()
    return state


def _state_args_for_gp(material, STATEV, kk):
    """Read old state-variable values for one Gauss point."""
    defaults = _state_arg_defaults(material)
    if STATEV is None:
        return defaults

    state = {}
    offset = 0
    for name, size, init in _state_var_sizes(material):
        start = kk * sum(s for _, s, _ in _state_var_sizes(material)) + offset
        vals = np.asarray(STATEV[start:start + size], dtype=float)
        if size == 1:
            state[f'{name}_old'] = vals[0]
        else:
            state[f'{name}_old'] = vals.reshape(np.asarray(init).shape)
        offset += size
    return {**defaults, **state}


def assemble_element(weakform, coords, U, DU, DTIME, PROPS, STATEV=None,
                     element=None, TIME=None, return_state=False):
    """
    Assemble element RHS and AMATRX for a single 2D element.

    Mirrors the generated Fortran UEL exactly.

    Args:
        weakform: WeakForm instance (verified)
        coords: (2, 8) nodal coordinates
        U: (NDOFEL,) total DOFs at end of increment
        DU: (NDOFEL,) incremental DOFs
        DTIME: time increment
        PROPS: material properties array
        STATEV: optional flat Gauss-point state array. If omitted, material
            state_vars defaults are used.
        element: optional element name/config, currently "quad4" or "quad8".
        TIME: optional Abaqus-style [step_time, total_time] sequence. If
            omitted, both values default to 0.0.

    Returns:
        RHS: (NDOFEL,) residual vector
        AMATRX: (NDOFEL, NDOFEL) tangent matrix
    """
    if element is not None:
        cfg = ELEMENT_CONFIGS[element.lower()] if isinstance(element, str) else element
        weakform._set_element_config(cfg)
    else:
        cfg = weakform._layout_config

    ndim = weakform.ndim
    ndofel = weakform.ndofel
    fields = weakform.fields
    field_names = weakform.field_names
    equations = weakform.equations
    blocks = weakform.tangent_blocks()
    mat = weakform._mat

    n_corner = weakform._n_corner
    n_nodes = weakform._n_nodes

    field_args = getattr(weakform, '_field_arg_defs', FIELD_ARGS)

    # Determine which fields need gradients
    needs_grad = set()
    for eq_info in equations.values():
        for fv in eq_info['field_vars']:
            if field_args.get(fv, {}).get('kind') == 'vector':
                needs_grad.add(fv.replace('grad_', ''))

    # ---- DOF parsing (current + old) + edof maps ----
    field_nodes = {}     # field_name -> nodal values
    field_old = {}       # field_name -> old nodal values
    edof = {}            # field_name -> edof array

    idx = 0
    for a in range(n_nodes):
        for fname in field_names:
            f = fields[fname]
            is_corner_only = (f.degree < 2)
            nn = len(weakform._field_nodes[fname])

            if is_corner_only and a >= n_corner:
                continue

            if fname not in field_nodes:
                if isinstance(f, VectorField):
                    field_nodes[fname] = np.zeros((ndim, nn))
                    field_old[fname] = np.zeros((ndim, nn))
                    edof[fname] = np.zeros((ndim, nn), dtype=int)
                else:
                    field_nodes[fname] = np.zeros(nn)
                    field_old[fname] = np.zeros(nn)
                    edof[fname] = np.zeros(nn, dtype=int)

            if isinstance(f, VectorField):
                for i in range(ndim):
                    field_nodes[fname][i, a] = U[idx]
                    field_old[fname][i, a] = U[idx] - DU[idx]
                    edof[fname][i, a] = idx
                    idx += 1
            else:
                # For corner-only fields, a is the absolute node index
                # but we store at the field's local index
                if is_corner_only:
                    local_a = a  # corners are 0..3
                else:
                    local_a = a
                field_nodes[fname][local_a] = U[idx]
                field_old[fname][local_a] = U[idx] - DU[idx]
                edof[fname][local_a] = idx
                idx += 1

    assert idx == ndofel, f"DOF count mismatch: {idx} != {ndofel}"

    # ---- Initialize outputs ----
    RHS = np.zeros(ndofel)
    AMATRX = np.zeros((ndofel, ndofel))
    new_state_by_gp = {}   # kk -> state dict from stress_PK1 (real eval), if any

    # DTIME guard: raise if zero dt is passed to a rate-dependent material
    has_dt_param = any(
        'dt' in minfo['all_params']
        for minfo in weakform._mat._methods.values()
    )
    if DTIME == 0.0 and has_dt_param:
        raise ValueError(
            "DTIME=0 passed to assemble_element, but material has "
            "dt-dependent methods. Use a nonzero time increment."
        )
    dt_safe = DTIME
    if TIME is None:
        time_total = 0.0
    else:
        time_arr = np.asarray(TIME, dtype=float)
        if time_arr.size < 2:
            raise ValueError("TIME must contain [step_time, total_time]")
        time_total = float(time_arr[1])

    # ---- Gauss loop ----
    if cfg.gauss_subroutine == 'gauss_2d_2x2':
        xi_gp, w_gp = gauss_2d_2x2()
    elif cfg.gauss_subroutine == 'gauss_2d_3x3':
        xi_gp, w_gp = gauss_2d_3x3()
    else:
        raise NotImplementedError(
            f"reference_assembly only supports 2D quadrature, got "
            f"{cfg.gauss_subroutine}")

    for kk in range(cfg.n_gauss_points):
        xi, eta = xi_gp[kk]

        # Shape functions
        if cfg.shape_subroutine == 'shape_quad8':
            sh_hi, dshxi_hi = shape_quad8(xi, eta)
        elif cfg.shape_subroutine == 'shape_quad4':
            sh_hi, dshxi_hi = shape_quad4(xi, eta)
        else:
            raise NotImplementedError(
                f"reference_assembly only supports Quad4/Quad8, got "
                f"{cfg.shape_subroutine}")
        dsh_hi, detJ, Jinv, stat = map_grad_2d(dshxi_hi, coords, cfg.n_nodes)
        if stat == 0:
            if return_state:    # degenerate element: keep incoming state
                return RHS, AMATRX, _flatten_gp_state(
                    mat, {}, cfg.n_gauss_points, STATEV)
            return RHS, AMATRX  # degenerate element

        # Quad4 shapes (if needed)
        has_deg1 = any(fields[fn].degree < 2 for fn in field_names)
        if has_deg1:
            sh4, dshxi4 = shape_quad4(xi, eta)
            dsh4 = apply_jinv(dshxi4, Jinv, 4)

        wdetJ = detJ * w_gp[kk]

        # ---- Field interpolation ----
        gp_vals = {}  # field variable -> GP value
        gp_old = {}   # old values

        for fname in field_names:
            f = fields[fname]
            nn = len(weakform._field_nodes[fname])
            sh = sh_hi if f.degree >= 2 else sh4
            dsh = dsh_hi if f.degree >= 2 else dsh4

            if isinstance(f, VectorField):
                # Deformation gradient: F = I + grad(u)
                F = np.eye(3)
                F_old = np.eye(3)
                for a in range(nn):
                    for i in range(ndim):
                        for j in range(ndim):
                            F[i, j] += field_nodes[fname][i, a] * dsh[a, j]
                            F_old[i, j] += field_old[fname][i, a] * dsh[a, j]
                gp_vals['F'] = F
                gp_old['F_old'] = F_old

            else:
                # Scalar field: value + optional gradient
                val = 0.0
                val_old = 0.0
                for a in range(nn):
                    val += field_nodes[fname][a] * sh[a]
                    val_old += field_old[fname][a] * sh[a]
                gp_vals[fname] = val
                gp_old[f'{fname}_old'] = val_old

                if fname in needs_grad:
                    grad = np.zeros(3)
                    for a in range(nn):
                        for j in range(ndim):
                            grad[j] += field_nodes[fname][a] * dsh[a, j]
                    gp_vals[f'grad_{fname}'] = grad

        # ---- Material evaluation (real, for RHS) ----
        mat_outputs = {}  # method_name -> output value
        state_args = _state_args_for_gp(mat, STATEV, kk)

        for mname in mat.defined_methods():
            method = mat._methods[mname]['callable']
            mparams = mat._methods[mname]['all_params']

            call_args = {}
            for p in mparams:
                if p in gp_vals:
                    call_args[p] = gp_vals[p]
                elif p in gp_old:
                    call_args[p] = gp_old[p]
                elif p == 'dt':
                    call_args[p] = dt_safe
                elif p == 'time':
                    call_args[p] = time_total
                elif p == 'F':
                    call_args[p] = gp_vals['F']
                elif p == 'F_old':
                    call_args[p] = gp_old['F_old']
                elif p in state_args:
                    call_args[p] = state_args[p]
                elif p + '_gp' in gp_vals:
                    call_args[p] = gp_vals[p + '_gp']

            result = method(**call_args)
            if isinstance(result, tuple):
                if len(result) > 1 and isinstance(result[1], dict):
                    new_state_by_gp[kk] = result[1]   # updated state (real eval)
                result = result[0]
            if isinstance(result, (complex, np.complexfloating)):
                if abs(result.imag) > 1e-12:
                    raise RuntimeError(
                        f"Material method {mname!r} returned non-trivial "
                        f"complex result ({result}) at a real evaluation point. "
                        f"Imaginary part should be zero."
                    )
                mat_outputs[mname] = np.real(result)
            else:
                mat_outputs[mname] = result

        # ---- CS tangent computation ----
        tangents = {}  # block_key -> tangent array

        for (block_key, wrt_var), binfo in blocks.items():
            mat_method_name = _material_method_for_equation(
                binfo['equation'], binfo)
            if mat_method_name is None:
                continue

            method = mat._methods[mat_method_name]['callable']
            mparams = mat._methods[mat_method_name]['all_params']
            ret_shape = _RET_SHAPE.get(mat_method_name, 'scalar')

            # Build args dict for this method
            call_args = {}
            for p in mparams:
                if p in gp_vals:
                    call_args[p] = gp_vals[p].copy() if isinstance(gp_vals[p], np.ndarray) else gp_vals[p]
                elif p in gp_old:
                    call_args[p] = gp_old[p].copy() if isinstance(gp_old[p], np.ndarray) else gp_old[p]
                elif p == 'dt':
                    call_args[p] = dt_safe
                elif p == 'time':
                    call_args[p] = time_total
                elif p in state_args:
                    call_args[p] = (
                        state_args[p].copy()
                        if isinstance(state_args[p], np.ndarray)
                        else state_args[p]
                    )

            tang = _compute_cs_tangent(method, call_args, wrt_var,
                                       binfo['wrt_kind'], ret_shape)
            tangents[(block_key, wrt_var)] = tang

        # ---- RHS assembly ----
        for eq_name, eq_info in equations.items():
            eq_meta = EQUATION_INFO[eq_name]
            test_field = eq_info['test_field']
            tf = fields[test_field]
            nn_row = len(weakform._field_nodes[test_field])
            sh_row = sh_hi if tf.degree >= 2 else sh4
            dsh_row = dsh_hi if tf.degree >= 2 else dsh4

            if eq_meta['return_type'] == 'tensor':
                # Momentum: RHS -= P_iJ * dN_a/dX_J * wdetJ
                P = mat_outputs['stress_PK1']
                for a in range(nn_row):
                    for i in range(ndim):
                        row = edof[test_field][i, a]
                        for j in range(ndim):
                            RHS[row] -= P[i, j] * dsh_row[a, j] * wdetJ

            elif eq_meta['return_type'] == 'scalar':
                # Constraint: RHS -= r_p * N_a * wdetJ
                rp = mat_outputs[_material_method_for_equation(eq_name)]
                for a in range(nn_row):
                    row = edof[test_field][a]
                    RHS[row] -= rp * sh_row[a] * wdetJ

            elif eq_meta['return_type'] == 'tuple':
                # Transport: storage + flux
                storage_method = eq_meta['sub_terms']['storage']['material_method']
                flux_method = eq_meta['sub_terms']['flux']['material_method']
                cdot = mat_outputs[storage_method]
                jR = mat_outputs[flux_method]
                for a in range(nn_row):
                    row = edof[test_field][a]
                    # Storage: RHS -= c_dot * N * wdetJ
                    RHS[row] -= cdot * sh_row[a] * wdetJ
                    # Flux: RHS += j_R . dN/dX * wdetJ
                    for j in range(ndim):
                        RHS[row] += jR[j] * dsh_row[a, j] * wdetJ

        # ---- AMATRX assembly ----
        for (block_key, wrt_var), binfo in blocks.items():
            tang = tangents.get((block_key, wrt_var))
            if tang is None:
                continue

            sign = binfo['amatrx_sign']
            row_asm = binfo['row_assembly']
            col_asm = binfo['col_assembly']
            test_field = binfo['test_field']
            wrt_kind = binfo['wrt_kind']

            # Test field info
            tf = fields[test_field]
            nn_row = len(weakform._field_nodes[test_field])
            sh_row = sh_hi if tf.degree >= 2 else sh4
            dsh_row = dsh_hi if tf.degree >= 2 else dsh4

            # Trial field info
            if wrt_var == 'F':
                trial_field = 'u'
            elif wrt_var.startswith('grad_'):
                trial_field = wrt_var.replace('grad_', '')
            else:
                trial_field = wrt_var
            cf = fields.get(trial_field)
            if cf is None:
                continue
            nn_col = len(weakform._field_nodes[trial_field])
            sh_col = sh_hi if cf.degree >= 2 else sh4
            dsh_col = dsh_hi if cf.degree >= 2 else dsh4

            # Assembly based on pattern
            if row_asm == 'grad' and col_asm == 'grad':
                if isinstance(tf, VectorField) and wrt_kind == 'matrix':
                    # dP/dF: (3,3,3,3)
                    for a in range(nn_row):
                        for b in range(nn_col):
                            for i in range(ndim):
                                row = edof[test_field][i, a]
                                for k in range(ndim):
                                    col = edof[trial_field][k, b]
                                    for j in range(ndim):
                                        for l in range(ndim):
                                            AMATRX[row, col] += (
                                                sign * tang[i, j, k, l]
                                                * dsh_row[a, j]
                                                * dsh_col[b, l] * wdetJ)

                elif isinstance(tf, ScalarField) and wrt_kind == 'matrix':
                    # djR/dF: (3,3,3)
                    for a in range(nn_row):
                        for b in range(nn_col):
                            for k in range(ndim):
                                col = edof[trial_field][k, b]
                                row = edof[test_field][a]
                                for j in range(ndim):
                                    for l in range(ndim):
                                        AMATRX[row, col] += (
                                            sign * tang[j, k, l]
                                            * dsh_row[a, j]
                                            * dsh_col[b, l] * wdetJ)

                elif isinstance(tf, VectorField) and wrt_kind == 'vector':
                    # dP/dgrad_scalar: (3,3,3)
                    for a in range(nn_row):
                        for b in range(nn_col):
                            col = edof[trial_field][b]
                            for i in range(ndim):
                                row = edof[test_field][i, a]
                                for j in range(ndim):
                                    for l in range(ndim):
                                        AMATRX[row, col] += (
                                            sign * tang[i, j, l]
                                            * dsh_row[a, j]
                                            * dsh_col[b, l] * wdetJ)

                elif isinstance(tf, ScalarField) and wrt_kind == 'vector':
                    # djR/dgrad_mu: (3,3)
                    for a in range(nn_row):
                        for b in range(nn_col):
                            row = edof[test_field][a]
                            col = edof[trial_field][b]
                            for j in range(ndim):
                                for l in range(ndim):
                                    AMATRX[row, col] += (
                                        sign * tang[j, l]
                                        * dsh_row[a, j]
                                        * dsh_col[b, l] * wdetJ)

            elif row_asm == 'grad' and col_asm == 'value':
                if isinstance(tf, VectorField):
                    # dP/dp or dP/dmu: (3,3)
                    for a in range(nn_row):
                        for b in range(nn_col):
                            col = edof[trial_field][b]
                            for i in range(ndim):
                                row = edof[test_field][i, a]
                                for j in range(ndim):
                                    AMATRX[row, col] += (
                                        sign * tang[i, j]
                                        * dsh_row[a, j]
                                        * sh_col[b] * wdetJ)
                else:
                    # djR/dp or djR/dmu: (3,)
                    for a in range(nn_row):
                        for b in range(nn_col):
                            row = edof[test_field][a]
                            col = edof[trial_field][b]
                            for j in range(ndim):
                                AMATRX[row, col] += (
                                    sign * tang[j]
                                    * dsh_row[a, j]
                                    * sh_col[b] * wdetJ)

            elif row_asm == 'value' and col_asm == 'grad':
                # drp/dF or dcdot/dF: (3,3)
                for a in range(nn_row):
                    for b in range(nn_col):
                        row = edof[test_field][a]
                        for k in range(ndim):
                            col = edof[trial_field][k, b]
                            for l in range(ndim):
                                AMATRX[row, col] += (
                                    sign * tang[k, l]
                                    * sh_row[a]
                                    * dsh_col[b, l] * wdetJ)

            elif row_asm == 'value' and col_asm == 'value':
                # drp/dp, drp/dmu, dcdot/dp: scalar
                for a in range(nn_row):
                    for b in range(nn_col):
                        row = edof[test_field][a]
                        col = edof[trial_field][b]
                        AMATRX[row, col] += (
                            sign * tang
                            * sh_row[a]
                            * sh_col[b] * wdetJ)

    if return_state:
        return RHS, AMATRX, _flatten_gp_state(mat, new_state_by_gp,
                                              cfg.n_gauss_points, STATEV)
    return RHS, AMATRX


def _flatten_gp_state(material, new_state_by_gp, n_gp, STATEV):
    """Flatten per-GP stress_PK1 state dicts into a STATEV-layout array (GP-major,
    state_vars order, C-flatten — matching `_state_args_for_gp`'s reshape).  GPs
    with no update fall back to the incoming STATEV (or the material defaults)."""
    sizes = _state_var_sizes(material)
    per_gp = sum(s for _, s, _ in sizes)
    if per_gp == 0:
        return np.zeros(0)
    out = np.zeros(n_gp * per_gp)
    defaults = _state_arg_defaults(material)
    for kk in range(n_gp):
        sd = new_state_by_gp.get(kk)
        off = kk * per_gp
        for name, size, init in sizes:
            if sd is not None and name in sd:
                val = np.real(np.asarray(sd[name], dtype=complex)).ravel()
            elif STATEV is not None:
                val = np.asarray(STATEV[off:off + size], dtype=float).ravel()
            else:
                val = np.asarray(defaults[f'{name}_old'], dtype=float).ravel()
            out[off:off + size] = val
            off += size
    return out

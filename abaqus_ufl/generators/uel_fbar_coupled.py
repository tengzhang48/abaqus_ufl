"""
Experimental F-bar UEL generator for coupled multiphysics problems.

This module is kept as a research prototype because it records the coupled
F-bar experiments from the swelling work. It is deliberately not routed through
the public ``generate_uel(..., formulation=...)`` API. Production coupled gel
work should use the standard mixed ``u,p,mu`` route or a local-pressure
formulation; F-bar is public only for pure mechanical single-field elements.
"""

from ..core.fields import VectorField, ScalarField
from ..core.defs import FIELD_ARGS
from ..core.weakform import EQUATION_INFO
from .umat_gen import _state_var_info, _nstate_per_gp, _fortran_dims


def _emit_amatrx_assembly_fbar(lines, weakform, fields, field_names, blocks,
                               tangent_specs, needs_grad, cfg):
    """Emit AMATRX assembly with F-bar alpha scaling for F-derivative blocks."""
    for spec in tangent_specs:
        tname = spec['tname']
        binfo = spec['binfo']
        sign = binfo['amatrx_sign']
        row_asm = binfo['row_assembly']
        col_asm = binfo['col_assembly']
        test_field = binfo['test_field']
        wrt_var = binfo['wrt_field']
        wrt_kind = binfo['wrt_kind']

        tf = fields[test_field]
        sh_row, dsh_row, nn_row = cfg.sh_for_degree(tf.degree)

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

        alpha_prefix = ''
        if wrt_var == 'F' and binfo['equation'] == 'momentum_equation':
            alpha_prefix = 'alpha * '
        sign_str = '+' if sign > 0 else '-'
        lines.append(f'C       {tname} (sign={sign_str}1,'
                     f' {row_asm}-{col_asm})')

        if row_asm == 'grad' and col_asm == 'grad':
            if isinstance(tf, VectorField) and wrt_kind == 'matrix':
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
                lines.append(f'     &                {sign_str} {alpha_prefix}{tname}(i,j,k,l)')
                lines.append(f'     &                * {dsh_row}(ii_v,j)')
                lines.append(f'     &                * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'                  END DO')
                lines.append(f'                END DO')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            elif isinstance(tf, ScalarField) and wrt_kind == 'matrix':
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            DO k = 1, ndim')
                lines.append(f'              col = edof_{trial_field}(k, jj_v)')
                lines.append(f'              DO j = 1, ndim')
                lines.append(f'                row = edof_{test_field}(ii_v)')
                lines.append(f'                DO l = 1, ndim')
                lines.append(f'                  AMATRX(row,col) =')
                lines.append(f'     &              AMATRX(row,col)')
                lines.append(f'     &              {sign_str} {alpha_prefix}{tname}(j,k,l)')
                lines.append(f'     &              * {dsh_row}(ii_v,j)')
                lines.append(f'     &              * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'                END DO')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            elif isinstance(tf, ScalarField) and wrt_kind == 'vector':
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            row = edof_{test_field}(ii_v)')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO j = 1, ndim')
                lines.append(f'              DO l = 1, ndim')
                lines.append(f'                AMATRX(row,col) =')
                lines.append(f'     &            AMATRX(row,col)')
                lines.append(f'     &            {sign_str} {alpha_prefix}{tname}(j,l)')
                lines.append(f'     &            * {dsh_row}(ii_v,j)')
                lines.append(f'     &            * {dsh_col}(jj_v,l) * wdetJ')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')

        elif row_asm == 'grad' and col_asm == 'value':
            if isinstance(tf, VectorField):
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO i = 1, ndim')
                lines.append(f'              row = edof_{test_field}(i, ii_v)')
                lines.append(f'              DO j = 1, ndim')
                lines.append(f'                AMATRX(row,col) =')
                lines.append(f'     &            AMATRX(row,col)')
                lines.append(f'     &            {sign_str} {alpha_prefix}{tname}(i,j)')
                lines.append(f'     &            * {dsh_row}(ii_v,j)')
                lines.append(f'     &            * {sh_col}(jj_v) * wdetJ')
                lines.append(f'              END DO')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')
            else:
                lines.append(f'        DO ii_v = 1, {nn_row}')
                lines.append(f'          DO jj_v = 1, {nn_col}')
                lines.append(f'            row = edof_{test_field}(ii_v)')
                lines.append(f'            col = edof_{trial_field}(jj_v)')
                lines.append(f'            DO j = 1, ndim')
                lines.append(f'              AMATRX(row,col) =')
                lines.append(f'     &          AMATRX(row,col)')
                lines.append(f'     &          {sign_str} {alpha_prefix}{tname}(j)')
                lines.append(f'     &          * {dsh_row}(ii_v,j)')
                lines.append(f'     &          * {sh_col}(jj_v) * wdetJ')
                lines.append(f'            END DO')
                lines.append(f'          END DO')
                lines.append(f'        END DO')

        elif row_asm == 'value' and col_asm == 'grad':
            lines.append(f'        DO ii_v = 1, {nn_row}')
            lines.append(f'          DO jj_v = 1, {nn_col}')
            lines.append(f'            row = edof_{test_field}(ii_v)')
            lines.append(f'            DO k = 1, ndim')
            lines.append(f'              col = edof_{trial_field}(k, jj_v)')
            lines.append(f'              DO l = 1, ndim')
            lines.append(f'                AMATRX(row,col) =')
            lines.append(f'     &            AMATRX(row,col)')
            lines.append(f'     &            {sign_str} {alpha_prefix}{tname}(k,l)')
            lines.append(f'     &            * {sh_row}(ii_v)')
            lines.append(f'     &            * {dsh_col}(jj_v,l) * wdetJ')
            lines.append(f'              END DO')
            lines.append(f'            END DO')
            lines.append(f'          END DO')
            lines.append(f'        END DO')

        elif row_asm == 'value' and col_asm == 'value':
            lines.append(f'        DO ii_v = 1, {nn_row}')
            lines.append(f'          DO jj_v = 1, {nn_col}')
            lines.append(f'            row = edof_{test_field}(ii_v)')
            lines.append(f'            col = edof_{trial_field}(jj_v)')
            lines.append(f'            AMATRX(row,col) =')
            lines.append(f'     &        AMATRX(row,col)')
            lines.append(f'     &        {sign_str} {alpha_prefix}{tname}')
            lines.append(f'     &        * {sh_row}(ii_v)')
            lines.append(f'     &        * {sh_col}(jj_v) * wdetJ')
            lines.append(f'          END DO')
            lines.append(f'        END DO')

        lines.append('')


def _emit_fbar_two_f_cs_call(lines, mat_prefix, field_names, fields,
                              equations, needs_grad, sv_info, cs_needs_dt,
                              tangent_specs, fbar_name, fbar_old_name,
                              mech_tangents):
    """Emit two CS tangent calls: mechanical at Fbar, transport at physical F.

    The second CS call overwrites all tangent arrays.  We save the mechanical
    tangents from the first call into tmp_* arrays and restore them after.
    """
    scalar_fields = [fn for fn in field_names
                     if isinstance(fields[fn], ScalarField)]
    mech_tnames = {s['tname'] for s in mech_tangents}

    for pass_idx, (f_name, f_old_name) in enumerate([
        (fbar_name, fbar_old_name),
        ('F', 'F_old'),
    ]):
        if pass_idx == 0:
            lines.append('C       --- Call 1: mechanical tangents at Fbar ---')
        else:
            lines.append('C       --- Call 2: transport tangents at physical F ---')

        cs_args = [f_name]
        for sfn in sorted([fn for fn in scalar_fields
                           if fn in [v for eq in equations.values()
                                     for v in eq['field_vars']]]):
            cs_args.append(f'{sfn}_gp')
        for sfn in sorted(needs_grad):
            cs_args.append(f'grad_{sfn}')
        cs_args.append(f_old_name)
        # Single source of truth for history-scalar argument order
        # (audit finding H1): must match the cs_tangents dummy list.
        from .uel_gen import _ordered_history_scalars_from
        all_hist = {v for eq in equations.values()
                    for v in eq['history_vars']}
        for old_name in _ordered_history_scalars_from(
                field_names, fields, all_hist):
            cs_args.append(f'{old_name}_gp')
        for name in sv_info:
            cs_args.append(f'{name}_old')
        cs_args.append('PROPS')
        if cs_needs_dt:
            cs_args.append('dt_safe')
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

        if pass_idx == 0:
            # Save mechanical tangents
            lines.append('C       Save mechanical tangents from Fbar call')
            for spec in mech_tangents:
                tn = spec['tname']
                lines.append(f'        tmp_{tn} = {tn}')
        else:
            # Restore mechanical tangents
            lines.append('C       Restore mechanical tangents from Fbar call')
            for spec in mech_tangents:
                tn = spec['tname']
                lines.append(f'        {tn} = tmp_{tn}')
    lines.append('')


def _generate_uel_fbar_coupled(weakform, mat_prefix, cfg):
    """Generate F-bar UEL for coupled multiphysics (VectorField + ScalarField).

    Three-pass assembly:
      Pass 1: Compute centroid F0 and store per-GP data
      Pass 2: Fbar, material eval, RHS, AMATRX_std (alpha-scaled), store Q
      Pass 3: F-bar correction (rank-1 update) for displacement columns
    """
    # Lazy imports to avoid circular dependency at module load time
    from .uel_gen import (
        _RETURN_SHAPE, _append_state_var_declarations,
        _append_svars_read, _append_svars_write,
        _append_state_old_to_complex,
        _emit_rhs_assembly,
    )

    ndim = cfg.ndim
    fields = weakform.fields
    field_names = weakform.field_names
    equations = weakform.equations
    blocks = weakform.tangent_blocks()
    nprops = len(weakform._mat.props_names)
    sv_info = _state_var_info(weakform._mat)
    state_old_names = {f'{name}_old' for name in sv_info}
    nstate_per_gp = _nstate_per_gp(sv_info) if sv_info else 0

    n_nodes = cfg.n_nodes
    ngp = cfg.n_gauss_points

    _sh = cfg.sh_name
    _dsh = cfg.dsh_name
    _d = cfg.ndim
    _jinv = cfg.jac_inv_size

    vector_fields = [fn for fn in field_names if isinstance(fields[fn], VectorField)]
    scalar_fields = [fn for fn in field_names if isinstance(fields[fn], ScalarField)]
    u_field = vector_fields[0] if vector_fields else None

    needs_grad = set()
    for eq_info in equations.values():
        for fv in eq_info['field_vars']:
            if FIELD_ARGS.get(fv, {}).get('kind') == 'vector':
                base = fv.replace('grad_', '')
                needs_grad.add(base)

    cs_needs_dt = any(
        'dt' in weakform._mat._methods[mname]['all_params']
        for mname in weakform._mat.defined_methods()
    )

    ndofel = n_nodes * sum(f.ndof_per_node for f in fields.values())

    lines = []

    # === SUBROUTINE HEADER ===
    lines.append('C======================================================================')
    lines.append('C     UEL subroutine -- generated by abaqus_ufl (F-bar coupled)')
    lines.append(f'C     Material: {weakform._mat.__class__.__name__}')
    lines.append(f"C     Fields: {', '.join(field_names)}")
    lines.append(f'C     Element: {cfg.name.capitalize()}, NDOFEL = {ndofel}, F-bar = ON')
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

    # Abaqus interface
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

    # Local parameters
    lines.append('C     --- Local parameters ---')
    lines.append(f'      INTEGER, PARAMETER :: ndim = {ndim}')
    lines.append(f'      INTEGER, PARAMETER :: NNODE_E = {n_nodes}')
    lines.append(f'      INTEGER, PARAMETER :: NGP = {ngp}')
    if sv_info:
        lines.append(f'      INTEGER, PARAMETER :: NSTATE_PER_GP = {nstate_per_gp}')
    lines.append('')

    # Nodal arrays
    for fname in field_names:
        f = fields[fname]
        if isinstance(f, VectorField):
            lines.append(f'      DOUBLE PRECISION :: {fname}_node(ndim, NNODE_E)')
            lines.append(f'      DOUBLE PRECISION :: {fname}_old(ndim, NNODE_E)')
            lines.append(f'      INTEGER :: edof_{fname}(ndim, NNODE_E)')
        else:
            lines.append(f'      DOUBLE PRECISION :: {fname}_node(NNODE_E)')
            lines.append(f'      DOUBLE PRECISION :: {fname}_old(NNODE_E)')
            lines.append(f'      INTEGER :: edof_{fname}(NNODE_E)')
    lines.append('')

    # Shape functions
    lines.append('C     --- Shape functions and mapping ---')
    lines.append(f'      DOUBLE PRECISION :: {_sh}(NNODE_E), {cfg.dshxi_name}(NNODE_E,{_d}), {_dsh}(NNODE_E,{_d})')
    lines.append(f'      DOUBLE PRECISION :: xi_gp(NGP,{_d}), w_gp(NGP)')
    lines.append(f'      DOUBLE PRECISION :: {cfg.coords_name}({_d}, NNODE_E)')
    lines.append(f'      DOUBLE PRECISION :: detJxi, Jinv({_jinv},{_jinv})')
    lines.append('      INTEGER :: ngp_out, stat')
    lines.append('')

    # Per-GP stored data
    lines.append('C     --- Per-GP stored data ---')
    lines.append(f'      DOUBLE PRECISION :: {_dsh}_all(NGP, NNODE_E, {_d})')
    lines.append(f'      DOUBLE PRECISION :: {_sh}_all(NGP, NNODE_E)')
    lines.append('      DOUBLE PRECISION :: detJxi_all(NGP)')
    lines.append('      DOUBLE PRECISION :: F_all(3, 3, NGP)')
    lines.append('      DOUBLE PRECISION :: F_old_all(3, 3, NGP)')
    lines.append('      DOUBLE PRECISION :: J_all(NGP)')
    lines.append('      DOUBLE PRECISION :: J_old_all(NGP)')
    lines.append('      DOUBLE PRECISION :: Finv_all(3, 3, NGP)')
    lines.append('')

    # F-bar quantities
    lines.append('C     --- F-bar quantities ---')
    lines.append('      DOUBLE PRECISION :: V0, J0, J0_old')
    lines.append('      DOUBLE PRECISION :: alpha, alpha_old')
    lines.append('      DOUBLE PRECISION :: Fbar(3,3), Fbar_old(3,3)')
    lines.append('      DOUBLE PRECISION :: F0(3,3), F0_old(3,3)')
    lines.append('      DOUBLE PRECISION :: dJ0_du(ndim, NNODE_E)')
    lines.append(f'      DOUBLE PRECISION :: dNx_cent(NNODE_E,{_d}), Finv0(3,3)')
    lines.append('      DOUBLE PRECISION :: Q_all(ndim, ndim, NGP)')
    lines.append('      DOUBLE PRECISION :: Q_flux(ndim, NGP)')
    lines.append('      DOUBLE PRECISION :: Q_storage(NGP)')
    lines.append('      DOUBLE PRECISION :: g_ai, h_bk, wdetJ, dt_safe')
    lines.append('      DOUBLE PRECISION :: elem_diag, inc_abs')
    lines.append('      DOUBLE PRECISION :: g_j, g_c')
    lines.append('')
    lines.append('C     --- Deformation gradient ---')
    lines.append('      DOUBLE PRECISION :: F(3,3), F_old(3,3)')
    lines.append('      DOUBLE PRECISION :: det33d')
    lines.append('      DOUBLE COMPLEX :: det33z')
    lines.append('')

    # Scalar field GP values
    for fname in scalar_fields:
        lines.append(f'      DOUBLE PRECISION :: {fname}_gp, {fname}_old_gp')
        if fname in needs_grad:
            lines.append(f'      DOUBLE PRECISION :: grad_{fname}(3)')
    lines.append('')

    # Material outputs
    lines.append('C     --- Material outputs (real) ---')
    for mname in weakform._mat.defined_methods():
        rs = _RETURN_SHAPE.get(mname)
        if rs == 'tensor':
            lines.append('      DOUBLE PRECISION :: P_real(3,3)')
        elif rs == 'scalar':
            sn = 'rp_real' if mname == 'pressure_resid' else 'cdot_real'
            lines.append(f'      DOUBLE PRECISION :: {sn}')
        elif rs == 'vector':
            lines.append('      DOUBLE PRECISION :: jR_real(3)')
    lines.append('')

    # Complex temps
    lines.append('C     --- Complex temporaries ---')
    lines.append('      DOUBLE COMPLEX :: Fbar_z(3,3), Pz_eval(3,3)')
    lines.append('      DOUBLE COMPLEX :: jRz_eval(3), rpz_eval, cdotz_eval')
    lines.append('      DOUBLE COMPLEX :: Fz_fbar_r(3,3), Fold_fbar_r(3,3)')
    lines.append('      DOUBLE COMPLEX :: Fz_phys_r(3,3), Fold_phys_r(3,3)')
    for fname in scalar_fields:
        lines.append(f'      DOUBLE COMPLEX :: {fname}z_r')
        if fname in needs_grad:
            lines.append(f'      DOUBLE COMPLEX :: grad_{fname}_r(3)')
    lines.append('')

    # Tangent arrays
    lines.append('C     --- Tangent arrays ---')
    tangent_specs = []
    for (bkey, wrt_var), binfo in blocks.items():
        mat_meth = binfo.get('material_method')
        if mat_meth is None:
            eq_name = binfo['equation']
            if eq_name == 'momentum_equation':
                mat_meth = 'stress_PK1'
            elif eq_name == 'pressure_equation':
                mat_meth = 'pressure_resid'
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

    # Temp arrays for mechanical tangents (computed at Fbar, saved before
    # second CS call at physical F overwrites them).
    mech_tangents = []
    for spec in tangent_specs:
        eq_name = spec['binfo'].get('equation', '')
        mat_meth = spec['binfo'].get('material_method', '')
        if eq_name == 'momentum_equation' or mat_meth == 'stress_PK1':
            mech_tangents.append(spec)
    for spec in mech_tangents:
        tn = spec['tname']
        rs = spec['ret_shape']
        wk = spec['wrt_kind']
        if rs == 'tensor' and wk == 'matrix':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3,3,3,3)')
        elif rs == 'tensor' and wk == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3,3)')
        elif rs == 'scalar' and wk == 'matrix':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3,3)')
        elif rs == 'scalar' and wk == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}')
        elif rs == 'vector' and wk == 'matrix':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3,3,3)')
        elif rs == 'vector' and wk == 'scalar':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3)')
        elif rs == 'vector' and wk == 'vector':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3,3)')
        elif rs == 'scalar' and wk == 'vector':
            lines.append(f'      DOUBLE PRECISION :: tmp_{tn}(3)')

    _append_state_var_declarations(lines, sv_info)

    # Hourglass stabilization locals (declared here because Fortran requires
    # declarations before first executable statement)
    if n_nodes == 4:
        lines.append('C     --- Hourglass stabilization locals ---')
        lines.append('      DOUBLE PRECISION :: hg1, hg_k')
        lines.append(f'      DOUBLE PRECISION :: hg_stiff_el, hg_p')
        lines.append(f'      DOUBLE PRECISION :: hg_gamma1({n_nodes})')
        lines.append('      DATA hg_gamma1 / 1.0d0, -1.0d0,  1.0d0, -1.0d0 /')
    lines.append('')

    lines.append('      INTEGER :: idx, ii_v, jj_v, mm_v, i, j, k, l, m, N, kk, gg')
    lines.append('      INTEGER :: row, col')
    lines.append('')

    # Extract coordinates
    lines.append(f'C     Extract {_d}D coordinates')
    lines.append('      DO ii_v = 1, NNODE_E')
    for idim in range(1, _d + 1):
        lines.append(f'        {cfg.coords_name}({idim}, ii_v) = COORDS({idim}, ii_v)')
    lines.append('      END DO')
    lines.append('')

    # Zero RHS and AMATRX
    lines.append('C     Zero RHS and AMATRX')
    lines.append('      DO i = 1, NDOFEL')
    lines.append('        RHS(i, 1) = 0.0d0')
    lines.append('        DO j = 1, NDOFEL')
    lines.append('          AMATRX(i, j) = 0.0d0')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')
    lines.append('C     RHS(:,2) is used by procedures that allocate NRHS >= 2')
    lines.append('C     (e.g. the modified Riks incremental-load vector, which this')
    lines.append('C     element does not implement). Zero it whenever allocated so no')
    lines.append('C     garbage is returned. Note LFLAGS(3)=4 is the Abaqus mass-matrix')
    lines.append('C     request, not Riks; mass matrices are likewise not implemented.')
    lines.append('      IF (NRHS .GE. 2) THEN')
    lines.append('        DO i = 1, NDOFEL')
    lines.append('          RHS(i, 2) = 0.0d0')
    lines.append('        END DO')
    lines.append('      END IF')
    lines.append('')
    lines.append('C     Procedure-type dispatch (Abaqus UEL contract): this element')
    lines.append('C     supports the static/quasi-static residual and stiffness')
    lines.append('C     requests LFLAGS(3)=1,2,5. Mass, damping, and initial-')
    lines.append('C     acceleration requests (3,4,6) are outside the supported')
    lines.append('C     scope and return the zeroed arrays (no dynamics claimed).')
    lines.append('      IF (LFLAGS(3) .EQ. 3 .OR. LFLAGS(3) .EQ. 4 .OR.')
    lines.append('     &    LFLAGS(3) .EQ. 6) THEN')
    lines.append('        RETURN')
    lines.append('      END IF')
    lines.append('')

    # DTIME guard
    lines.append('C     DTIME guard')
    lines.append('      dt_safe = DTIME')
    lines.append('      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14')
    lines.append('')

    # DOF parsing
    lines.append('C     Parse DOFs and build edof maps')
    lines.append('      idx = 0')
    lines.append('      DO ii_v = 1, NNODE_E')
    for fname in field_names:
        f = fields[fname]
        if isinstance(f, VectorField):
            lines.append(f'        DO i = 1, ndim')
            lines.append(f'          idx = idx + 1')
            lines.append(f'          {fname}_node(i, ii_v) = U(idx)')
            lines.append(f'          {fname}_old(i, ii_v) = U(idx) - DU(idx, 1)')
            lines.append(f'          edof_{fname}(i, ii_v) = idx')
            lines.append(f'        END DO')
        else:
            lines.append(f'        idx = idx + 1')
            lines.append(f'        {fname}_node(ii_v) = U(idx)')
            lines.append(f'        {fname}_old(ii_v) = U(idx) - DU(idx, 1)')
            lines.append(f'        edof_{fname}(ii_v) = idx')
    lines.append('      END DO')
    lines.append('')

    # Large trial increments usually indicate that Newton has jumped onto a
    # nonphysical distortion mode.  Ask the driver for a cutback before
    # evaluating material routines at such states.
    lines.append('C     Trial increment guards')
    if n_nodes >= 3:
        lines.append('      elem_diag = 0.0d0')
        for idim in range(1, _d + 1):
            lines.append('      elem_diag = elem_diag')
            lines.append(f'     &  + ({cfg.coords_name}({idim},1)')
            lines.append(f'     &  - {cfg.coords_name}({idim},3))**2')
        lines.append('      elem_diag = DSQRT(elem_diag)')
    else:
        lines.append('      elem_diag = 1.0d0')
    if u_field:
        lines.append('      DO ii_v = 1, NNODE_E')
        lines.append('        DO i = 1, ndim')
        lines.append(f'          inc_abs = DABS(DU(edof_{u_field}(i,ii_v),1))')
        lines.append('          IF (inc_abs .GT. 10.0d0*elem_diag) THEN')
        lines.append('            PNEWDT = 0.5d0')
        lines.append('            RETURN')
        lines.append('          END IF')
        lines.append(f'          inc_abs = DABS({u_field}_node(i,ii_v))')
        lines.append('          IF (inc_abs .GT. 5.0d0*elem_diag) THEN')
        lines.append('            PNEWDT = 0.25d0')
        lines.append('            RETURN')
        lines.append('          END IF')
        lines.append('        END DO')
        lines.append('      END DO')
    for fname in scalar_fields:
        lines.append('      DO ii_v = 1, NNODE_E')
        lines.append(f'        inc_abs = DABS(DU(edof_{fname}(ii_v),1))')
        lines.append('        IF (inc_abs .GT. 1.0d6) THEN')
        lines.append('          PNEWDT = 0.5d0')
        lines.append('          RETURN')
        lines.append('        END IF')
        lines.append('      END DO')
    lines.append('')

    # =====================================================================
    # PASS 1: Compute centroid F0 and store per-GP data
    # =====================================================================
    # Centroid F-bar is a standard choice (cf. Datta & Nguyen 2025).
    # F0 is the deformation gradient at the element centroid (xi=eta=0).
    # Fbar_gp = (detF0 / detF_gp)^(1/ndim) * F_gp
    lines.append('C     =============================================')
    lines.append('C     Pass 1: Compute centroid F0, store per-GP data')
    lines.append('C     =============================================')
    xi_parts = cfg.gauss_xi_args.split(', ')

    lines.append(f'      CALL {cfg.gauss_subroutine}(xi_gp, w_gp, ngp_out)')
    lines.append('')
    lines.append('      V0 = 0.0d0')
    lines.append('')
    lines.append('      DO kk = 1, NGP')
    lines.append(f'        CALL {cfg.shape_subroutine}({xi_parts[0]},')
    lines.append(f'     &    {", ".join(xi_parts[1:])}, {_sh}, {cfg.dshxi_name})')
    lines.append(f'        CALL {cfg.jac_subroutine}({cfg.dshxi_name}, {cfg.coords_name}, NNODE_E,')
    lines.append(f'     &    {_dsh}, detJxi, Jinv, stat)')
    lines.append('')
    lines.append('        IF (stat .EQ. 0) THEN')
    lines.append('          PNEWDT = 0.25d0')
    lines.append('          RETURN')
    lines.append('        END IF')
    lines.append('')
    lines.append('C       Store per-GP data')
    lines.append('        detJxi_all(kk) = detJxi')
    lines.append('        DO ii_v = 1, NNODE_E')
    lines.append(f'          {_sh}_all(kk, ii_v) = {_sh}(ii_v)')
    for idim in range(1, _d + 1):
        lines.append(f'          {_dsh}_all(kk, ii_v, {idim}) = {_dsh}(ii_v, {idim})')
    lines.append('        END DO')
    lines.append('')

    # Compute F at each GP
    lines.append('C       Compute F')
    lines.append('        CALL eye33d(F_all(1,1,kk))')
    lines.append('        CALL eye33d(F_old_all(1,1,kk))')
    if u_field:
        lines.append(f'        DO ii_v = 1, NNODE_E')
        lines.append(f'          DO i = 1, ndim')
        lines.append(f'            DO j = 1, ndim')
        lines.append(f'              F_all(i,j,kk) = F_all(i,j,kk)')
        lines.append(f'     &          + {u_field}_node(i,ii_v)*{_dsh}(ii_v,j)')
        lines.append(f'              F_old_all(i,j,kk) = F_old_all(i,j,kk)')
        lines.append(f'     &          + {u_field}_old(i,ii_v)*{_dsh}(ii_v,j)')
        lines.append(f'            END DO')
        lines.append(f'          END DO')
        lines.append(f'        END DO')
    lines.append('')

    # J and Finv
    lines.append('C       J and Finv')
    lines.append('        J_all(kk) = det33d(F_all(1,1,kk))')
    lines.append('        J_old_all(kk) = det33d(F_old_all(1,1,kk))')
    lines.append('        IF (J_all(kk) .LE. 0.0d0) THEN')
    lines.append('          PNEWDT = 0.25d0')
    lines.append('          RETURN')
    lines.append('        END IF')
    lines.append('        IF (J_old_all(kk) .LE. 0.0d0) THEN')
    lines.append('          PNEWDT = 0.25d0')
    lines.append('          RETURN')
    lines.append('        END IF')
    lines.append('        CALL inv33d(F_all(1,1,kk),')
    lines.append('     &    Finv_all(1,1,kk))')
    lines.append('')
    lines.append('        V0 = V0 + detJxi * w_gp(kk)')
    lines.append('      END DO')
    lines.append('')

    # Compute F at element centroid (xi=eta=0).
    # The centroid F-bar reference is a standard choice (cf. Datta & Nguyen 2025).
    lines.append('C     Compute F at element centroid (F-bar reference)')
    lines.append(f'      CALL {cfg.shape_subroutine}(0.0d0, 0.0d0,')
    lines.append(f'     &    {_sh}, {cfg.dshxi_name})')
    lines.append(f'      CALL {cfg.jac_subroutine}({cfg.dshxi_name}, {cfg.coords_name}, NNODE_E,')
    lines.append(f'     &    dNx_cent, detJxi, Jinv, stat)')
    lines.append('      IF (stat .EQ. 0) THEN')
    lines.append('        PNEWDT = 0.25d0')
    lines.append('        RETURN')
    lines.append('      END IF')
    lines.append('      CALL eye33d(F0)')
    lines.append('      CALL eye33d(F0_old)')
    if u_field:
        lines.append('      DO ii_v = 1, NNODE_E')
        lines.append('        DO i = 1, ndim')
        lines.append('          DO j = 1, ndim')
        lines.append('            F0(i,j) = F0(i,j)')
        lines.append(f'     &        + {u_field}_node(i,ii_v)*dNx_cent(ii_v,j)')
        lines.append('            F0_old(i,j) = F0_old(i,j)')
        lines.append(f'     &        + {u_field}_old(i,ii_v)*dNx_cent(ii_v,j)')
        lines.append('          END DO')
        lines.append('        END DO')
        lines.append('      END DO')
    lines.append('')
    lines.append('      J0 = det33d(F0)')
    lines.append('      J0_old = det33d(F0_old)')
    lines.append('      IF (J0 .LE. 0.0d0) THEN')
    lines.append('        PNEWDT = 0.25d0')
    lines.append('        RETURN')
    lines.append('      END IF')
    lines.append('      IF (J0_old .LE. 0.0d0) THEN')
    lines.append('        PNEWDT = 0.25d0')
    lines.append('        RETURN')
    lines.append('      END IF')
    lines.append('      CALL inv33d(F0, Finv0)')
    lines.append('')
    lines.append('C     dJ0/du(k,a) = J0 * Finv0(N,k) * dNa/dX_N (at centroid)')
    lines.append('      DO jj_v = 1, NNODE_E')
    lines.append('        DO k = 1, ndim')
    lines.append('          dJ0_du(k,jj_v) = 0.0d0')
    lines.append('          DO N = 1, ndim')
    lines.append('            dJ0_du(k,jj_v) = dJ0_du(k,jj_v)')
    lines.append('     &        + J0 * Finv0(N,k) * dNx_cent(jj_v,N)')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('')

    # =====================================================================
    # PASS 2: Fbar, material, RHS, standard AMATRX, store Q
    # =====================================================================
    lines.append('C     =============================================')
    lines.append('C     Pass 2: Fbar, material, RHS, AMATRX_std, Q')
    lines.append('C     =============================================')
    lines.append('      DO kk = 1, NGP')
    lines.append('        IF (J_all(kk) .LE. 0.0d0) THEN')
    lines.append('          PNEWDT = 0.25d0')
    lines.append('          RETURN')
    lines.append('        END IF')
    lines.append('        alpha = (J0 / J_all(kk))')
    lines.append('     &    ** (1.0d0 / DBLE(ndim))')
    lines.append('        alpha_old = (J0_old / J_old_all(kk))')
    lines.append('     &    ** (1.0d0 / DBLE(ndim))')
    lines.append('        wdetJ = detJxi_all(kk) * w_gp(kk)')
    lines.append(f'        {_sh}(1:NNODE_E) = {_sh}_all(kk,1:NNODE_E)')
    lines.append(f'        {_dsh}(1:NNODE_E,1:{_d}) = {_dsh}_all(kk,1:NNODE_E,1:{_d})')
    lines.append('')

    # Compute Fbar
    lines.append('C       Compute Fbar (in-plane only)')
    lines.append('        CALL eye33d(Fbar)')
    lines.append('        CALL eye33d(Fbar_old)')
    lines.append('        DO i = 1, ndim')
    lines.append('          DO j = 1, ndim')
    lines.append('            Fbar(i,j) = alpha * F_all(i,j,kk)')
    lines.append('            Fbar_old(i,j) = alpha_old * F_old_all(i,j,kk)')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')

    # Keep physical F and F_old for non-mechanical methods.
    lines.append('C       Keep physical F and F_old for non-mechanical methods')
    lines.append('        F = F_all(:,:,kk)')
    lines.append('        F_old = F_old_all(:,:,kk)')
    lines.append('')

    # Interpolate scalar fields
    for fname in scalar_fields:
        lines.append(f'C       Interpolate {fname}')
        lines.append(f'        {fname}_gp = 0.0d0')
        lines.append(f'        {fname}_old_gp = 0.0d0')
        lines.append(f'        DO ii_v = 1, NNODE_E')
        lines.append(f'          {fname}_gp = {fname}_gp')
        lines.append(f'     &      + {fname}_node(ii_v)*{_sh}(ii_v)')
        lines.append(f'          {fname}_old_gp = {fname}_old_gp')
        lines.append(f'     &      + {fname}_old(ii_v)*{_sh}(ii_v)')
        lines.append(f'        END DO')
        if fname in needs_grad:
            lines.append(f'        grad_{fname}(1) = 0.0d0')
            lines.append(f'        grad_{fname}(2) = 0.0d0')
            lines.append(f'        grad_{fname}(3) = 0.0d0')
            lines.append(f'        DO ii_v = 1, NNODE_E')
            lines.append(f'          DO j = 1, ndim')
            lines.append(f'            grad_{fname}(j) = grad_{fname}(j)')
            lines.append(f'     &        + {fname}_node(ii_v)*{_dsh}_all(kk,ii_v,j)')
            lines.append(f'          END DO')
            lines.append(f'        END DO')
        lines.append('')

    # Read state variables
    _append_svars_read(lines, sv_info, indent='        ')
    lines.append('')

    # Material evaluation (real, for RHS)
    lines.append('C       Real material evaluation (for RHS)')
    lines.append('C       Fbar kinematics for mechanical stress')
    lines.append('        CALL real2complex33(Fbar, Fz_fbar_r)')
    lines.append('        CALL real2complex33(Fbar_old, Fold_fbar_r)')
    lines.append('C       Physical kinematics for flux, storage, pressure')
    lines.append('        CALL real2complex33(F, Fz_phys_r)')
    lines.append('        CALL real2complex33(F_old, Fold_phys_r)')
    for fname in scalar_fields:
        lines.append(f'        {fname}z_r = DCMPLX({fname}_gp, 0.0d0)')
        if fname in needs_grad:
            lines.append(f'        CALL real2complex3(grad_{fname}, grad_{fname}_r)')
    _append_state_old_to_complex(lines, sv_info, indent='        ')
    lines.append('')

    # Call each material method
    for mname in weakform._mat.defined_methods():
        minfo = weakform._mat._methods[mname]
        sub_name = f'{mat_prefix}_{mname}'
        rs = _RETURN_SHAPE.get(mname)
        has_state_call = bool(sv_info) and mname == 'stress_PK1'

        call_args = []
        is_mech = (mname == 'stress_PK1')
        for p in minfo['all_params']:
            if p == 'F':
                call_args.append('Fz_fbar_r' if is_mech else 'Fz_phys_r')
            elif p == 'F_old':
                call_args.append('Fold_fbar_r' if is_mech else 'Fold_phys_r')
            elif p in state_old_names:
                call_args.append(f'{p}_z')
            elif p == 'dt':
                call_args.append('dt_safe')
            elif FIELD_ARGS.get(p, {}).get('kind') == 'vector':
                call_args.append(f'{p}_r')
            elif p in ('p', 'mu', 'p_old'):
                if p.endswith('_old'):
                    base = p.replace('_old', '')
                    call_args.append(f'DCMPLX({base}_old_gp, 0.0d0)')
                else:
                    call_args.append(f'{p}z_r')
            else:
                call_args.append(p)
        call_args.append('PROPS')

        if rs == 'tensor':
            call_args.append('Pz_eval')
            if has_state_call:
                for name in sv_info:
                    call_args.append(f'{name}_new_z')
            lines.append(f'        CALL {sub_name}(')
            lines.append(f'     &    {", ".join(call_args)})')
            lines.append('        DO i = 1, 3')
            lines.append('          DO j = 1, 3')
            lines.append('            P_real(i,j) = DBLE(Pz_eval(i,j))')
            lines.append('          END DO')
            lines.append('        END DO')
        elif rs == 'vector':
            call_args.append('jRz_eval')
            lines.append(f'        CALL {sub_name}(')
            lines.append(f'     &    {", ".join(call_args)})')
            lines.append('        DO i = 1, 3')
            lines.append('          jR_real(i) = DBLE(jRz_eval(i))')
            lines.append('        END DO')
        else:
            if mname == 'pressure_resid':
                call_args.append('rpz_eval')
                lines.append(f'        CALL {sub_name}(')
                lines.append(f'     &    {", ".join(call_args)})')
                lines.append('        rp_real = DBLE(rpz_eval)')
            else:
                call_args.append('cdotz_eval')
                lines.append(f'        CALL {sub_name}(')
                lines.append(f'     &    {", ".join(call_args)})')
                lines.append('        cdot_real = DBLE(cdotz_eval)')
        lines.append('')

    # Write state variables
    _append_svars_write(lines, sv_info, indent='        ')
    if sv_info:
        lines.append('')

    # CS tangent engine: two calls.
    # Mechanical tangents (stress_PK1): computed with Fbar kinematics.
    # Transport tangents (flux, storage): computed with physical F.
    # The second call overwrites the tangents arrays, so we save/restore
    # the mechanical ones.
    lines.append('C       CS tangent engine')
    _emit_fbar_two_f_cs_call(lines, mat_prefix, field_names, fields,
                             equations, needs_grad, sv_info, cs_needs_dt,
                             tangent_specs, 'Fbar', 'Fbar_old', mech_tangents)
    lines.append('')

    # RHS assembly
    lines.append('C       RHS assembly')
    _emit_rhs_assembly(lines, weakform, fields, field_names, equations,
                       needs_grad, cfg)
    lines.append('')

    # AMATRX assembly with alpha scaling
    lines.append('C       AMATRX assembly (standard, alpha-scaled)')
    _emit_amatrx_assembly_fbar(lines, weakform, fields, field_names, blocks,
                               tangent_specs, needs_grad, cfg)
    lines.append('')

    # Store Q for Pass 3
    lines.append('C       Q = A : Fbar (in-plane only)')
    lines.append('        DO i = 1, ndim')
    lines.append('          DO j = 1, ndim')
    lines.append('            Q_all(i,j,kk) = 0.0d0')
    lines.append('            DO m = 1, ndim')
    lines.append('              DO N = 1, ndim')
    lines.append('                Q_all(i,j,kk) =')
    lines.append('     &            Q_all(i,j,kk)')
    dpdF_name = None
    for spec in tangent_specs:
        if spec['tname'].startswith('dPK1_dF') or spec['tname'].startswith('dP_dF'):
            dpdF_name = spec['tname']
            break
    if dpdF_name:
        lines.append(f'     &            + {dpdF_name}(i,j,m,N)')
        lines.append('     &            * Fbar(m,N)')
    lines.append('              END DO')
    lines.append('            END DO')
    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('')

    # Q_flux and Q_storage for K_mu-u F-bar correction.
    # Deliberately disabled: a common coupled-F-bar approach evaluates coupled
    # residuals with F-bar kinematics, but neglects F-bar corrections in the
    # off-diagonal transport-displacement tangents.
    has_flux = any(spec['tname'] == 'dflux_dF' for spec in tangent_specs)
    has_storage = any(spec['tname'] == 'dstorage_dF' for spec in tangent_specs)
    emit_scalar_fbar_correction = False

    if has_flux and emit_scalar_fbar_correction:
        lines.append('C       Q_flux = dflux_dF : Fbar')
        lines.append('        DO i = 1, ndim')
        lines.append('          Q_flux(i,kk) = 0.0d0')
        lines.append('          DO m = 1, ndim')
        lines.append('            DO N = 1, ndim')
        lines.append('              Q_flux(i,kk) =')
        lines.append('     &          Q_flux(i,kk)')
        lines.append('     &          + dflux_dF(i,m,N)')
        lines.append('     &          * Fbar(m,N)')
        lines.append('            END DO')
        lines.append('          END DO')
        lines.append('        END DO')
        lines.append('')

    if has_storage and emit_scalar_fbar_correction:
        lines.append('C       Q_storage = dstorage_dF : Fbar')
        lines.append('        Q_storage(kk) = 0.0d0')
        lines.append('        DO m = 1, ndim')
        lines.append('          DO N = 1, ndim')
        lines.append('            Q_storage(kk) =')
        lines.append('     &        Q_storage(kk)')
        lines.append('     &        + dstorage_dF(m,N)')
        lines.append('     &        * Fbar(m,N)')
        lines.append('          END DO')
        lines.append('        END DO')
        lines.append('')

    lines.append('      END DO')
    lines.append('C     End Pass 2')
    lines.append('')

    # Diagnostic hourglass stabilization for F-bar Quad4 elements.
    # This is not classic reduced-integration hourglass control: the element
    # still uses 2x2 integration.  Penalize only the bilinear checkerboard mode
    # with a small shear-scale element stiffness so the stabilization does not
    # reinsert an artificial bulk penalty.
    hg_factor = 0.05  # fraction of G * V0 / elem_diag**2
    lines.append('C     =============================================')
    lines.append('C     Hourglass stabilization (anti-hourglass)')
    lines.append('C     =============================================')
    if n_nodes == 4:
        lines.append(f'      IF (elem_diag .GT. 0.0d0) THEN')
        lines.append(f'        hg_stiff_el = {hg_factor}d0 * PROPS(1)')
        lines.append(f'     &    * V0 / (elem_diag*elem_diag)')
        lines.append(f'      ELSE')
        lines.append(f'        hg_stiff_el = {hg_factor}d0 * PROPS(1)')
        lines.append(f'     &    * DSQRT(V0)')
        lines.append(f'      END IF')
        lines.append('')
        # RHS: resisting force against existing hourglass deformation
        lines.append('      DO i = 1, ndim')
        lines.append('        hg1 = 0.0d0')
        lines.append('        DO ii_v = 1, NNODE_E')
        lines.append(f'          hg1 = hg1 + u_node(i, ii_v) * hg_gamma1(ii_v)')
        lines.append('        END DO')
        lines.append('        DO ii_v = 1, NNODE_E')
        lines.append(f'          row = edof_u(i, ii_v)')
        lines.append('          hg_p = hg_stiff_el')
        lines.append('     &      *hg1*hg_gamma1(ii_v)')
        lines.append('          RHS(row,1) = RHS(row,1) - hg_p')
        lines.append('        END DO')
        lines.append('      END DO')
        lines.append('')
        # AMATRX: hourglass stiffness (only for displacement DOFs)
        lines.append('      DO ii_v = 1, NNODE_E')
        lines.append('        DO jj_v = 1, NNODE_E')
        lines.append('          hg_p = hg_gamma1(ii_v) * hg_gamma1(jj_v)')
        lines.append('          hg_k = hg_stiff_el * hg_p')
        lines.append('          DO i = 1, ndim')
        lines.append(f'            row = edof_u(i, ii_v)')
        lines.append(f'            col = edof_u(i, jj_v)')
        lines.append('            AMATRX(row,col) = AMATRX(row,col)'
                      ' + hg_k')
        lines.append('          END DO')
        lines.append('        END DO')
        lines.append('      END DO')
    lines.append('')

    # =====================================================================
    # PASS 3: F-bar correction
    # =====================================================================
    lines.append('C     =============================================')
    lines.append('C     Pass 3: F-bar correction (rank-1 update)')
    lines.append('C     =============================================')
    lines.append('      DO kk = 1, NGP')
    lines.append('        wdetJ = detJxi_all(kk) * w_gp(kk)')
    lines.append('')
    lines.append('C       --- K_uu correction ---')
    lines.append('        DO ii_v = 1, NNODE_E')
    lines.append('          DO jj_v = 1, NNODE_E')
    lines.append('            DO i = 1, ndim')
    lines.append(f'              row = edof_{u_field}(i, ii_v)')
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
    lines.append(f'                col = edof_{u_field}(k, jj_v)')
    lines.append('')
    lines.append('C               h_bk = (1/ndim)*[dJ0/du/Jbar')
    lines.append('C                       - Finv_Nk * dN_b/dX_N]')
    lines.append('                h_bk = dJ0_du(k,jj_v) / J0')
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
    lines.append('')
    lines.append('C       --- K_?u correction (separate nest, no ii_v) ---')
    lines.append('        DO jj_v = 1, NNODE_E')
    lines.append('          DO k = 1, ndim')
    lines.append(f'            col = edof_{u_field}(k, jj_v)')
    lines.append('')
    lines.append('C           h_bk = (1/ndim)*[dJ0/du/Jbar')
    lines.append('C                   - Finv_Nk * dN_b/dX_N]')
    lines.append('            h_bk = dJ0_du(k,jj_v) / J0')
    lines.append('            DO N = 1, ndim')
    lines.append('              h_bk = h_bk')
    lines.append('     &          - Finv_all(N,k,kk)')
    lines.append(f'     &          * {_dsh}_all(kk,jj_v,N)')
    lines.append('            END DO')
    lines.append('            h_bk = h_bk / DBLE(ndim)')
    lines.append('')

    if emit_scalar_fbar_correction:
        for sfn in scalar_fields:
            for eq_name, eq_info in equations.items():
                if eq_info.get('test_field') == sfn:
                    if has_flux:
                        lines.append(f'C           F-bar correction: flux term for {sfn}-equation')
                        lines.append('            DO mm_v = 1, NNODE_E')
                        lines.append(f'              row = edof_{sfn}(mm_v)')
                        lines.append('              g_j = 0.0d0')
                        lines.append('              DO j = 1, ndim')
                        lines.append('                g_j = g_j')
                        lines.append('     &            + Q_flux(j,kk)')
                        lines.append(f'     &            * {_dsh}_all(kk,mm_v,j) * wdetJ')
                        lines.append('              END DO')
                        lines.append('              AMATRX(row,col) =')
                        lines.append('     &          AMATRX(row,col)')
                        lines.append('     &          - g_j * h_bk')
                        lines.append('            END DO')
                        lines.append('')

                    if has_storage:
                        lines.append(f'C           F-bar correction: storage term for {sfn}-equation')
                        lines.append('            DO mm_v = 1, NNODE_E')
                        lines.append(f'              row = edof_{sfn}(mm_v)')
                        lines.append('              g_c = Q_storage(kk)')
                        lines.append(f'     &        * {_sh}_all(kk,mm_v) * wdetJ')
                        lines.append('              AMATRX(row,col) =')
                        lines.append('     &          AMATRX(row,col)')
                        lines.append('     &          + g_c * h_bk')
                        lines.append('            END DO')
                        lines.append('')
                    break

    lines.append('          END DO')
    lines.append('        END DO')
    lines.append('      END DO')
    lines.append('C     End Pass 3')
    lines.append('')
    lines.append('      RETURN')
    lines.append('      END SUBROUTINE UEL')
    lines.append('')

    return '\n'.join(lines)

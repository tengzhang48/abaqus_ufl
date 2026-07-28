"""
Tests for umat_gen.py with internal variable support.

Tests:
  1. Stateless material (NeoHookean) — regression, output unchanged
  2. State variable layout computation (_state_var_info)
  3. STATEV read code (standard Abaqus STATEV pattern: time(2)==0 branch)
  4. STATEV write code (real part only)
  5. Material subroutine signature with state args
  6. CS engine with state_old reset per perturbation
  7. UMAT wrapper with full STATEV read/write/CS flow
  8. Backwards compatibility: no state_vars → no STATEV code
"""

import os
import pytest
import numpy as np
from collections import OrderedDict

from abaqus_ufl.generators import umat_gen as ug


# ===================================================================
# Mock Material classes (mimic the real Material interface)
# ===================================================================

class MockNeoHookean:
    """Stateless material — should produce identical output to original."""
    props = dict(G=0.5, K=50.0)
    props_names = ['G', 'K']

    def __init__(self):
        self.G = 0.5
        self.K = 50.0
        self._methods = {
            'stress_PK1': {
                'callable': self.stress_PK1,
                'field_params': ['F'],
            }
        }

    def stress_PK1(self, F):
        J = det(F)
        Finv = inv(F)
        P = self.G * (F - Finv.T) + self.K * log(J) * inv(F).T
        return P


class MockJ2Plasticity:
    """
    Material with internal variables — simplified J2 for testing
    code generation structure (not a full constitutive model).
    """
    props = dict(E=200e3, nu=0.3, sigma_y=250.0, H=1000.0)
    props_names = ['E', 'nu', 'sigma_y', 'H']
    state_vars = OrderedDict([
        ('ep', 0.0),             # equivalent plastic strain (scalar)
        ('Fp', np.eye(3)),       # plastic deformation gradient (tensor)
    ])

    def __init__(self):
        self.E = 200e3
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1000.0
        self._methods = {
            'stress_PK1': {
                'callable': self.stress_PK1,
                'field_params': ['F'],
            }
        }

    def stress_PK1(self, F, ep_old, Fp_old, dt):
        """Simplified J2 — for structure testing only."""
        Fp_old_inv = inv(Fp_old)
        Fe = F @ Fp_old_inv
        J = det(F)
        P = self.E * (Fe - inv(Fe).T)
        ep_new = ep_old
        Fp_new = Fp_old
        return P, {'ep': ep_new, 'Fp': Fp_new}


class MockScalarOnly:
    """Material with only scalar state vars."""
    props = dict(mu=1.0, lam=10.0)
    props_names = ['mu', 'lam']
    state_vars = OrderedDict([
        ('d', 0.0),    # damage (scalar)
        ('ep', 0.0),   # equiv plastic strain (scalar)
    ])

    def __init__(self):
        self.mu = 1.0
        self.lam = 10.0
        self._methods = {
            'stress_PK1': {
                'callable': self.stress_PK1,
                'field_params': ['F'],
            }
        }

    def stress_PK1(self, F, d_old, ep_old, dt):
        J = det(F)
        P = self.mu * F
        d_new = d_old
        ep_new = ep_old
        return P, {'d': d_new, 'ep': ep_new}


class MockSubscriptAccess:
    """Material that reads individual eigenvalue components."""
    props = dict(G=0.5, K=50.0)
    props_names = ['G', 'K']

    def __init__(self):
        self.G = 0.5
        self.K = 50.0
        self._methods = {
            'stress_PK1': {
                'callable': self.stress_PK1,
                'field_params': ['F'],
            }
        }

    def stress_PK1(self, F):
        sigma_prin, V = eig(F)
        s0 = sigma_prin[0]
        s1 = sigma_prin[1]
        s2 = sigma_prin[2]
        v00 = V[0, 0]
        vcol0 = V[:, 0]
        vcol00 = vcol0[0]
        return F * (s0 + s1 + s2 + v00 + vcol00)


class MockRangeIndexing:
    """Material that indexes arrays with a Python range(N) loop variable."""
    props = dict(G=1.0)
    props_names = ['G']

    def __init__(self):
        self.G = 1.0
        self._methods = {
            'stress_PK1': {
                'callable': self.stress_PK1,
                'field_params': ['F'],
            }
        }

    def stress_PK1(self, F):
        P = 0.0 * F
        for a in range(3):
            col = F[:, a]
            P = P + col[a] * F
        return P


class MockNormalizeVector:
    """Material that normalizes a vector expression in the tensor DSL."""
    props = dict(G=1.0)
    props_names = ['G']

    def __init__(self):
        self.G = 1.0
        self._methods = {
            'stress_PK1': {
                'callable': self.stress_PK1,
                'field_params': ['F'],
            }
        }

    def stress_PK1(self, F):
        sigma_prin, V = eig(F)
        raw = V[:, 0]
        direction = normalize(self.G * raw)
        P = dyad(direction, direction)
        return P


# ===================================================================
# Tests: State variable info
# ===================================================================

class TestStateVarInfo:
    def test_no_state_vars(self):
        mat = MockNeoHookean()
        info = ug._state_var_info(mat)
        assert len(info) == 0

    def test_j2_state_vars(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        assert list(info.keys()) == ['ep', 'Fp']
        assert info['ep']['size'] == 1
        assert info['ep']['shape'] == 'scalar'
        assert info['Fp']['size'] == 9
        assert info['Fp']['shape'] == 'tensor'

    def test_nstate_per_gp(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        assert ug._nstate_per_gp(info) == 10  # 1 + 9

    def test_scalar_only(self):
        mat = MockScalarOnly()
        info = ug._state_var_info(mat)
        assert ug._nstate_per_gp(info) == 2  # 1 + 1


# ===================================================================
# Tests: STATEV read/write code generation
# ===================================================================

class TestStatevReadWrite:
    def test_read_empty(self):
        assert ug._generate_statev_read(OrderedDict()) == ''

    def test_read_j2(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        code = ug._generate_statev_read(info)

        # Init branch: first increment AND all-zero incoming STATEV, so
        # SDVINI / *INITIAL CONDITIONS state is honored, not clobbered.
        assert 'svinit_max = MAX(svinit_max, DABS(STATEV(i)))' in code
        assert ('IF (TIME(2) .EQ. 0.0d0 .AND. '
                'svinit_max .EQ. 0.0d0) THEN') in code
        assert 'ELSE' in code
        assert 'END IF' in code

        # Initial values: ep_old = 0.0, Fp_old = identity
        assert 'ep_old = 0.0d0' in code.replace(
            '0.000000000000000d+00', '0.0d0')

        # Should read from STATEV(1) for ep
        assert 'ep_old = STATEV(1)' in code
        # Should read Fp from STATEV(2:10), column-major
        # Fp_old(1,1) = STATEV(2)  [offset=1, j=0, i=0 -> idx=2]
        assert 'Fp_old(1,1) = STATEV(2)' in code
        # Fp_old(2,1) = STATEV(3)  [offset=1, j=0, i=1 -> idx=3]
        assert 'Fp_old(2,1) = STATEV(3)' in code
        # Fp_old(1,2) = STATEV(5)  [offset=1, j=1, i=0 -> idx=5]
        assert 'Fp_old(1,2) = STATEV(5)' in code
        # Fp_old(3,3) = STATEV(10) [offset=1, j=2, i=2 -> idx=10]
        assert 'Fp_old(3,3) = STATEV(10)' in code

    def test_write_j2(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        code = ug._generate_statev_write(info)

        # Real part only
        assert 'DBLE(ep_new_z)' in code
        assert 'STATEV(1) =' in code
        # Tensor write column-major
        assert 'DBLE(Fp_new_z(1,1))' in code
        assert 'STATEV(2) =' in code
        assert 'STATEV(10) =' in code

    def test_write_empty(self):
        assert ug._generate_statev_write(OrderedDict()) == ''


# ===================================================================
# Tests: Material subroutine signature
# ===================================================================

class TestStressSubroutine:
    def test_stateless_signature(self):
        mat = MockNeoHookean()
        src = ug._generate_stress_subroutine(mat, 'neohookean')
        # Should have original simple signature
        assert 'neohookean_stress_PK1(F, props, P_out)' in src
        # Should NOT have state args or dt
        assert 'ep_old' not in src
        assert 'Fp_old' not in src
        assert 'INTENT(IN) :: dt' not in src

    def test_stateful_signature(self):
        mat = MockJ2Plasticity()
        src = ug._generate_stress_subroutine(mat, 'j2')
        # Should have extended signature with state args
        assert 'j2_stress_PK1(' in src
        assert 'ep_old' in src
        assert 'Fp_old' in src
        assert 'ep_new' in src
        assert 'Fp_new' in src
        assert 'INTENT(IN) :: dt' in src
        # Argument order: F, state_old, dt, props, P_out, state_new
        sig_line = [l for l in src.split('\n')
                    if 'SUBROUTINE j2_stress_PK1' in l][0]
        # ep_old before dt, dt before props, P_out before ep_new
        assert src.index('ep_old') < src.index('INTENT(IN) :: dt')

    def test_stateful_declarations(self):
        mat = MockJ2Plasticity()
        src = ug._generate_stress_subroutine(mat, 'j2')
        # State old args: COMPLEX IN
        assert 'DOUBLE COMPLEX, INTENT(IN)  :: ep_old' in src
        assert 'DOUBLE COMPLEX, INTENT(IN)  :: Fp_old(3,3)' in src
        # dt: DOUBLE PRECISION IN
        assert 'DOUBLE PRECISION, INTENT(IN) :: dt' in src
        # State new: COMPLEX OUT
        assert 'DOUBLE COMPLEX, INTENT(OUT) :: ep_new' in src
        assert 'DOUBLE COMPLEX, INTENT(OUT) :: Fp_new(3,3)' in src


# ===================================================================
# Tests: CS tangent engine
# ===================================================================

class TestCSTangent:
    def test_stateless_cs(self):
        src = ug._generate_cs_dPdF('neohookean', 2)
        assert 'neohookean_cs_dPdF(F, props, dPdF)' in src
        assert 'ep_old' not in src
        assert 'dt' not in src
        # Should call material without state
        assert 'neohookean_stress_PK1(Fz, props, Pz)' in src

    def test_stateful_cs_signature(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_cs_dPdF('j2', 4, info)
        # Signature has state_old + dt
        assert 'j2_cs_dPdF(' in src
        assert 'ep_old' in src
        assert 'Fp_old' in src
        assert 'DOUBLE PRECISION, INTENT(IN) :: dt' in src

    def test_stateful_cs_reset(self):
        """State_old must be reset to real before each perturbation."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_cs_dPdF('j2', 4, info)
        # Inside the DO k/l loop, state must be reset
        assert 'ep_old_z = DCMPLX(ep_old, 0.0d0)' in src
        assert 'CALL real2complex33(Fp_old, Fp_old_z)' in src
        # These resets must appear INSIDE the DO loop
        do_idx = src.index('DO k = 1, 3')
        reset_idx = src.index('ep_old_z = DCMPLX(ep_old, 0.0d0)')
        end_do_idx = src.rindex('END DO')
        assert do_idx < reset_idx < end_do_idx

    def test_stateful_cs_discards_state_new(self):
        """CS engine declares state_new_z but only extracts dP/dF."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_cs_dPdF('j2', 4, info)
        # Should declare state_new_z for material call
        assert 'DOUBLE COMPLEX :: ep_new_z' in src
        assert 'DOUBLE COMPLEX :: Fp_new_z(3,3)' in src
        # But should NOT write them to STATEV or use them
        assert 'STATEV' not in src

    def test_stateful_cs_dt_safe(self):
        """CS engine should protect against DTIME=0."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_cs_dPdF('j2', 4, info)
        assert 'dt_safe = dt' in src
        assert 'IF (dt_safe .LT. 1.0d-14)' in src

    def test_stateful_cs_material_call(self):
        """CS engine calls material with correct arg order."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_cs_dPdF('j2', 4, info)
        # The CALL may span multiple continuation lines.
        # Collect everything from CALL to the closing ')'.
        lines = src.split('\n')
        in_call = False
        call_text = ''
        for line in lines:
            stripped = line.strip()
            if 'CALL j2_stress_PK1' in stripped:
                in_call = True
                call_text = stripped
            elif in_call:
                # Continuation line starts with &
                if stripped.startswith('&'):
                    call_text += ' ' + stripped[1:].strip()
                else:
                    break
            if in_call and ')' in call_text:
                break
        assert call_text, "No CALL j2_stress_PK1 found"
        # Should pass: Fz, ep_old_z, Fp_old_z, dt_safe, props,
        #              Pz, ep_new_z, Fp_new_z
        assert 'Fz' in call_text
        assert 'ep_old_z' in call_text
        assert 'Fp_old_z' in call_text
        assert 'dt_safe' in call_text
        assert 'props' in call_text
        assert 'Pz' in call_text
        assert 'ep_new_z' in call_text
        assert 'Fp_new_z' in call_text


# ===================================================================
# Tests: UMAT wrapper
# ===================================================================

class TestUMATWrapper:
    def test_stateless_wrapper(self):
        src = ug._generate_umat_wrapper('neohookean', 2)
        assert 'SUBROUTINE UMAT(' in src
        # No STATEV read/write
        assert 'TIME(2) .EQ. 0.0d0' not in src
        assert 'dt_safe' not in src
        # Simple material call
        assert 'neohookean_stress_PK1(Fz, PROPS, Pz)' in src
        assert 'neohookean_cs_dPdF(F, PROPS, dPdF)' in src

    def test_optional_abaqus_outputs_are_initialized(self):
        src = ug._generate_umat_wrapper('neohookean', 2)
        for assignment in (
                'RPL = 0.0d0',
                'DDSDDT = 0.0d0',
                'DRPLDE = 0.0d0',
                'DRPLDT = 0.0d0',
                'DDSDDE = 0.0d0'):
            assert assignment in src
        assert 'SSE = 0.0d0' not in src
        assert 'SPD = 0.0d0' not in src
        assert 'SCD = 0.0d0' not in src
        init_idx = src.index('DDSDDE = 0.0d0')
        guard_idx = src.index('IF (NDI .NE. 3')
        assert init_idx < guard_idx

    def test_stateful_wrapper_has_statev_read(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        assert ('IF (TIME(2) .EQ. 0.0d0 .AND. '
                'svinit_max .EQ. 0.0d0) THEN') in src
        assert 'DOUBLE PRECISION :: svinit_max' in src
        assert 'ep_old = STATEV(1)' in src
        assert 'Fp_old(1,1) = STATEV(2)' in src

    def test_stateful_wrapper_has_statev_write(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        assert 'DBLE(ep_new_z)' in src
        assert 'DBLE(Fp_new_z(' in src

    def test_stateful_wrapper_dt_safe(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        assert 'dt_safe = DTIME' in src
        assert 'IF (dt_safe .LT. 1.0d-14)' in src

    def test_stateful_wrapper_declares_state_vars(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        # Real state for STATEV read
        assert 'DOUBLE PRECISION :: ep_old' in src
        assert 'DOUBLE PRECISION :: Fp_old(3,3)' in src
        # Complex state for material call
        assert 'DOUBLE COMPLEX :: ep_old_z' in src
        assert 'DOUBLE COMPLEX :: Fp_old_z(3,3)' in src
        assert 'DOUBLE COMPLEX :: ep_new_z' in src
        assert 'DOUBLE COMPLEX :: Fp_new_z(3,3)' in src

    def test_stateful_wrapper_header_documents_statev(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        assert 'NSTATV = 10' in src
        assert 'STATEV(  1)' in src
        assert 'ep' in src
        assert 'Fp' in src

    def test_stateful_wrapper_cs_call(self):
        """CS engine called with state_old (real) args."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        # Find CS call
        cs_lines = [l.strip() for l in src.split('\n')
                    if 'j2_cs_dPdF' in l]
        assert len(cs_lines) == 1
        cs_call = cs_lines[0]
        # Should pass real state_old (not _z)
        assert 'ep_old,' in cs_call
        assert 'Fp_old,' in cs_call
        assert 'dt_safe' in cs_call
        # Should NOT pass complex versions
        assert 'ep_old_z' not in cs_call
        assert 'Fp_old_z' not in cs_call

    def test_stateful_wrapper_flow_order(self):
        """
        Verify the standard stateful-UMAT flow:
        1. Read F
        2. Read STATEV
        3. dt_safe
        4. Real evaluation (stress + state update)
        5. Write STATEV
        6. CS tangent
        7. Push-forward
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)

        idx_F = src.index('F(i,j) = DFGRD1(i,j)')
        idx_read = src.index('TIME(2) .EQ. 0.0d0')
        idx_dt = src.index('dt_safe = DTIME')
        idx_eval = src.index('j2_stress_PK1(')
        idx_write = src.index('DBLE(ep_new_z)')
        idx_cs = src.index('j2_cs_dPdF(')
        idx_push = src.index('pk1_to_cauchy_jaumann')

        assert idx_F < idx_read < idx_dt < idx_eval
        assert idx_eval < idx_write < idx_cs < idx_push


# ===================================================================
# Tests: Column-major storage order
# ===================================================================

class TestColumnMajor:
    def test_tensor_read_column_major(self):
        """
        Verify tensor STATEV read uses column-major:
        Fp(i,j) at offset + (j-1)*3 + i
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        code = ug._generate_statev_read(info)

        # offset=1 (after ep scalar at position 1)
        # Fp(1,1) -> offset + 0*3 + 0 + 1 = 2
        assert 'Fp_old(1,1) = STATEV(2)' in code
        # Fp(2,1) -> offset + 0*3 + 1 + 1 = 3
        assert 'Fp_old(2,1) = STATEV(3)' in code
        # Fp(3,1) -> offset + 0*3 + 2 + 1 = 4
        assert 'Fp_old(3,1) = STATEV(4)' in code
        # Fp(1,2) -> offset + 1*3 + 0 + 1 = 5
        assert 'Fp_old(1,2) = STATEV(5)' in code
        # Fp(2,2) -> offset + 1*3 + 1 + 1 = 6
        assert 'Fp_old(2,2) = STATEV(6)' in code
        # Fp(3,3) -> offset + 2*3 + 2 + 1 = 10
        assert 'Fp_old(3,3) = STATEV(10)' in code

    def test_tensor_write_column_major(self):
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        code = ug._generate_statev_write(info)

        assert 'STATEV(2) =' in code
        assert 'Fp_new_z(1,1)' in code
        assert 'STATEV(5) =' in code
        assert 'Fp_new_z(1,2)' in code


# ===================================================================
# Tests: Backwards compatibility
# ===================================================================

class TestBackwardsCompat:
    def test_neohookean_no_statev_code(self):
        """Stateless material generates zero STATEV-related code."""
        src = ug._generate_umat_wrapper('neo', 2)
        assert 'dt_safe' not in src
        # Should not have the TIME(2)==0 branching pattern
        assert 'TIME(2) .EQ. 0.0d0' not in src
        # Should not read/write individual STATEV entries
        # (STATEV(NSTATV) in declaration is fine, but no STATEV(1) etc)
        lines = src.split('\n')
        for line in lines:
            stripped = line.strip()
            # Skip declaration lines
            if 'STATEV(NSTATV)' in stripped:
                continue
            if 'STATEV(' in stripped and 'NSTATV' not in stripped:
                assert False, f"Unexpected STATEV access: {stripped}"

    def test_neohookean_cs_no_state(self):
        src = ug._generate_cs_dPdF('neo', 2)
        assert 'ep_old' not in src
        assert 'dt' not in src
        assert 'neo_stress_PK1(Fz, props, Pz)' in src


# ===================================================================
# Tests: Initial values in STATEV read
# ===================================================================

class TestInitialValues:
    def test_identity_tensor_init(self):
        """Fp initialized to identity when time(2)==0."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        code = ug._generate_statev_read(info)

        # In the IF branch (time==0), should set diagonal to 1
        if_block = code.split('ELSE')[0]
        assert 'Fp_old(1,1) = 1.0' in if_block.replace(
            '1.000000000000000d+00', '1.0')
        assert 'Fp_old(2,2) = 1.0' in if_block.replace(
            '1.000000000000000d+00', '1.0')
        assert 'Fp_old(3,3) = 1.0' in if_block.replace(
            '1.000000000000000d+00', '1.0')
        # Off-diagonal should be 0
        assert 'Fp_old(1,2) = 0.0' in if_block.replace(
            '0.000000000000000d+00', '0.0')


# ===================================================================
# Tests: Column 72 compliance
# ===================================================================

def _check_col72(src, section_name=''):
    """Check all code lines <= 72 columns. Returns list of violations."""
    violations = []
    for i, line in enumerate(src.split('\n'), 1):
        if not line:
            continue
        # Comments have no column limit in practice
        if line[0] in ('C', 'c', '!', '*'):
            continue
        if len(line) > 72:
            violations.append(
                f"{section_name} line {i} ({len(line)} chars): "
                f"{line[:80]}...")
    return violations


class TestColumn72:
    def test_wrapper_col72(self):
        """UMAT wrapper must have no lines > 72 columns."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        violations = _check_col72(src, 'wrapper')
        assert violations == [], \
            f"Column 72 violations:\n" + '\n'.join(violations)

    def test_cs_engine_col72(self):
        """CS engine must have no lines > 72 columns."""
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_cs_dPdF('j2', 4, info)
        violations = _check_col72(src, 'cs_engine')
        assert violations == [], \
            f"Column 72 violations:\n" + '\n'.join(violations)

    def test_stateless_wrapper_col72(self):
        """Stateless wrapper also compliant."""
        src = ug._generate_umat_wrapper('neo', 2)
        violations = _check_col72(src, 'stateless_wrapper')
        assert violations == [], \
            f"Column 72 violations:\n" + '\n'.join(violations)

    def test_stateless_cs_col72(self):
        """Stateless CS engine also compliant."""
        src = ug._generate_cs_dPdF('neo', 2)
        violations = _check_col72(src, 'stateless_cs')
        assert violations == [], \
            f"Column 72 violations:\n" + '\n'.join(violations)


class TestSubscriptAccess:
    """Generator supports reading individual array elements."""

    def test_vector_subscript(self):
        """sigma_prin[0] -> sigma_prin(1) in Fortran."""
        mat = MockSubscriptAccess()
        src = ug._generate_stress_subroutine(mat, 'subscr')
        assert 'sigma_prin(1)' in src
        assert 'sigma_prin(2)' in src
        assert 'sigma_prin(3)' in src

    def test_tensor_subscript(self):
        """V[0,0] -> V(1,1) in Fortran."""
        mat = MockSubscriptAccess()
        src = ug._generate_stress_subroutine(mat, 'subscr')
        assert 'V(1,1)' in src

    def test_tensor_column_subscript(self):
        """V[:,0] copies a tensor column into a 3-vector temp."""
        mat = MockSubscriptAccess()
        src = ug._generate_stress_subroutine(mat, 'subscr')
        assert 'DO ii = 1, 3' in src
        assert '= V(ii,1)' in src

    def test_subscript_col72(self):
        """Subscript-using wrapper stays within 72 columns."""
        mat = MockSubscriptAccess()
        src = ug._generate_umat_wrapper('subscr', 2)
        violations = _check_col72(src, 'subscript_wrapper')
        assert violations == [], \
            f"Column 72 violations:\n" + '\n'.join(violations)


class TestRangeLoopIndexing:
    """Python range(N) loop values must stay zero-based in Fortran."""

    def test_range_n_preserves_python_indices(self):
        mat = MockRangeIndexing()
        src = ug._generate_stress_subroutine(mat, 'rangeidx')
        assert 'DO a = 0, 2' in src
        assert '= F(ii,(a + 1))' in src
        assert 'col((a + 1))' in src
        assert 'DO a = 1, 3' not in src
        assert 'F(ii,(a + 2))' not in src


class TestNormalizeVector:
    """CS-safe 3-vector normalization is a first-class DSL operation."""

    def test_normalize_materializes_vector_expression(self):
        mat = MockNormalizeVector()
        src = ug._generate_stress_subroutine(mat, 'normvec')
        assert 'VNORM' in src
        assert 'UNITV' in src
        assert 'SQRT(' in src
        assert 'normalize(' not in src
        # The scaled vector expression must be assigned before component
        # indexing; Fortran cannot subscript a parenthesized array expression.
        assert 'DOUBLE COMPLEX :: ARG' in src


# ===================================================================
# Tests: reference stateful-UMAT flow structure
# ===================================================================

class TestReferenceFlowComparison:
    """
    Compare generated UMAT structure against the reference
    stateful-UMAT flow patterns (TIME(2)==0 init branch).
    """

    def test_reference_flow_flow_statev_init_branch(self):
        """
        The reference stateful flow uses:
          if(time(2).eq.zero) then
            Fp_t = zero; Fp_t(1,1) = one; ...
          else
            Fp_t(1,1) = statev(1); ...
          endif

        We generate the same pattern.
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        code = ug._generate_statev_read(info)
        # Two branches; the init branch additionally requires an all-zero
        # incoming STATEV so SDVINI-provided state is honored.
        assert ('IF (TIME(2) .EQ. 0.0d0 .AND. '
                'svinit_max .EQ. 0.0d0) THEN') in code
        assert 'ELSE' in code
        assert 'END IF' in code
        # Init branch has identity for Fp
        init_block = code.split('ELSE')[0]
        # Verify the 3 diagonal entries are set to 1
        ones = init_block.count('1.000000000000000d+00')
        assert ones == 3, f"Expected 3 identity diagonal entries, got {ones}"
        # ELSE branch reads from STATEV
        else_block = code.split('ELSE')[1]
        assert 'STATEV(1)' in else_block
        assert 'STATEV(10)' in else_block

    def test_reference_flow_flow_call_then_write(self):
        """
        Reference flow: state integration, then STATEV write.
        We: CALL mat_stress_PK1(...) then STATEV(1) = DBLE(ep_new_z).

        The key structural match: material evaluation happens BEFORE
        STATEV write, and only real parts get written.
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        # Material call before STATEV write
        call_idx = src.index('j2_stress_PK1(')
        write_idx = src.index('DBLE(ep_new_z)')
        assert call_idx < write_idx, \
            "Material call must come before STATEV write"

    def test_reference_flow_flow_tangent_after_state(self):
        """
        The reference flow computes the tangent with the update; we compute dP/dF
        via CS AFTER the real evaluation. The CS call must come
        after STATEV write (it uses state_old, not state_new).
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        write_idx = src.index('DBLE(ep_new_z)')
        cs_idx = src.index('j2_cs_dPdF(')
        assert write_idx < cs_idx, \
            "CS tangent must come after STATEV write"

    def test_reference_flow_statev_count_matches(self):
        """
        Reference layout: state scalars and tensors in declaration order.
        Our J2 mock: 10 STATEV (1 ep, 9 Fp).
        Verify the count is computed correctly.
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        assert ug._nstate_per_gp(info) == 10

    def test_reference_flow_voigt_packing_preserved(self):
        """
        The pk1_to_cauchy_jaumann function must still be called
        even with state variables — the push-forward is independent
        of internal state.
        """
        mat = MockJ2Plasticity()
        info = ug._state_var_info(mat)
        src = ug._generate_umat_wrapper('j2', 4, info)
        assert 'pk1_to_cauchy_jaumann(F, P_real, dPdF, ' \
               'STRESS, DDSDDE)' in src


class TestAbaqusTangentContract:
    def test_converter_preserves_every_computed_entry(self):
        src = ug._generate_pk1_to_cauchy_jaumann()
        assert 'DO j = 1, 6' in src
        assert 'DDSDDE(j,i) = DDSDDE(i,j)' not in src

    def test_small_strain_optional_outputs_precede_early_return(self):
        src = ug._generate_small_strain_umat_wrapper('small', 2)
        init_idx = src.index('RPL = 0.0d0')
        early_return_idx = src.index(
            'IF (NDI .NE. 3 .OR. NSHR .NE. 3 .OR. NTENS .NE. 6) THEN')
        assert init_idx < early_return_idx
        for assignment in (
                'DDSDDT = 0.0d0',
                'DRPLDE = 0.0d0',
                'DRPLDT = 0.0d0',
                'DDSDDE = 0.0d0'):
            assert assignment in src
        assert 'SSE = 0.0d0' not in src
        assert 'SPD = 0.0d0' not in src
        assert 'SCD = 0.0d0' not in src


# ===================================================================
# Tests: Fortran line wrapper utility
# ===================================================================

class TestLineWrapper:
    def test_short_line_unchanged(self):
        line = '      x = 1.0d0'
        assert ug._wrap_fortran_line(line) == line

    def test_comment_line_unchanged(self):
        line = 'C' + '=' * 100
        assert ug._wrap_fortran_line(line) == line

    def test_long_line_no_safe_break_raises(self):
        line = '      ' + 'A' * 80  # 86 chars, no break points
        # Should raise rather than emit an ugly/unsafe break
        with pytest.raises(RuntimeError):
            ug._wrap_fortran_line(line)

    def test_long_line_breaks_at_token(self):
        # Line with commas — should break at a comma, not mid-token
        line = '      CALL foo(alpha_long, beta_long, gamma_long, delta_long, epsilon_long, zeta_long)'
        result = ug._wrap_fortran_line(line)
        lines = result.split('\n')
        for l in lines:
            assert len(l) <= 72, f"Line too long ({len(l)}): {l}"
        # Continuation lines start with '     &'
        for l in lines[1:]:
            assert l.startswith('     &'), f"Bad continuation: {l}"

    def test_fortran_call_short(self):
        result = ug._fortran_call('      ', 'foo', ['a', 'b'])
        assert result == '      CALL foo(a, b)'
        assert len(result) <= 72

    def test_fortran_call_long(self):
        args = ['Fz', 'ep_old_z', 'Fp_old_z', 'dt_safe',
                'props', 'Pz', 'ep_new_z', 'Fp_new_z']
        result = ug._fortran_call(
            '          ', 'j2_stress_PK1', args)
        lines = result.split('\n')
        # All lines <= 72 cols
        for line in lines:
            assert len(line) <= 72, \
                f"Line too long ({len(line)}): {line}"
        # First line has CALL and first arg
        assert 'CALL j2_stress_PK1(Fz,' in lines[0]
        # Last line has closing paren
        assert lines[-1].strip().endswith(')')
        # All args present
        full = ' '.join(l.strip() for l in lines)
        for arg in args:
            assert arg in full


if __name__ == '__main__':
    import traceback
    passed = 0
    failed = 0
    for name in sorted(dir()):
        obj = eval(name)
        if isinstance(obj, type) and name.startswith('Test'):
            instance = obj()
            for mname in sorted(dir(instance)):
                if mname.startswith('test_'):
                    method = getattr(instance, mname)
                    try:
                        method()
                        passed += 1
                        print(f'  PASS: {name}.{mname}')
                    except Exception as e:
                        failed += 1
                        print(f'  FAIL: {name}.{mname}: {e}')
                        traceback.print_exc()
    print(f'\n{passed} passed, {failed} failed')

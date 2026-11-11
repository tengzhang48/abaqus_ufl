"""Translator guards: constructs that must fail closed, and emission rules.

Each guard is tested against the construct it exists to stop, so removing a
guard fails a test (audit 2026-10-10: G1, G4, G5, G8).
"""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log
from abaqus_ufl.generators import umat_gen as ug

REPO = Path(__file__).resolve().parents[1]


def stress_source(material):
    return ug._generate_stress_subroutine(material, 'probe')


# ---------------------------------------------------------------------------
# Rejections (G8): constructs the translator cannot represent must raise.
# ---------------------------------------------------------------------------

class EarlyReturn(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        if J.real > 1.0:
            return self.G * F
        return self.G * (F - inv(F).T)


class WhileLoop(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        s = 1.0
        while s.real < 2.0:
            s = s + 0.5
        return self.G * s * (F - inv(F).T)


class ComplexComparison(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        s = 1.0
        if J > 1.0:
            s = 2.0
        return self.G * s * (F - inv(F).T)


class ComplexComparisonOnTheRight(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        s = 1.0
        if 1.0 < J:
            s = 2.0
        return self.G * s * (F - inv(F).T)


class ChainedComparison(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        s = 1.0
        if 0.5 < J.real < 2.0:
            s = 2.0
        return self.G * s * (F - inv(F).T)


@pytest.mark.parametrize('material, message', [
    (EarlyReturn, r'Early `return` inside a control-flow block'),
    (WhileLoop, r'while loops are deliberately rejected'),
    (ComplexComparison, r'Comparison on potentially-complex expression'),
    (ComplexComparisonOnTheRight, r'Comparison on potentially-complex expression'),
    (ChainedComparison, r'Chained comparisons are not supported'),
])
def test_unrepresentable_construct_raises(material, message):
    with pytest.raises(NotImplementedError, match=message):
        stress_source(material())


# ---------------------------------------------------------------------------
# abs() under complex step (G4)
# ---------------------------------------------------------------------------

class AbsOfComplex(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * abs(J - 1.0) * (F - inv(F).T)


class AbsOfRealPart(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * abs(J.real - 1.0) * (F - inv(F).T)


def test_abs_of_possibly_complex_value_raises():
    """ABS(z) of a complex-step variable returns its modulus and drops the
    perturbation, so the generated tangent would be silently wrong."""
    with pytest.raises(NotImplementedError, match=r'abs\(\) of a possibly complex'):
        stress_source(AbsOfComplex())


def test_abs_of_real_part_is_allowed():
    assert 'ABS((DBLE(J) - ' in stress_source(AbsOfRealPart())


class LoopIndexRatio(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        s = 0.0
        for i in range(1, 3):
            for j in range(1, 3):
                s = s + i / j
        return self.G * s * (F - inv(F).T)


def test_division_of_loop_indices_is_real_division():
    """i and j are Fortran INTEGERs; Python's i / j must not truncate."""
    assert '(DBLE(i) / DBLE(j))' in stress_source(LoopIndexRatio())


# ---------------------------------------------------------------------------
# Integer powers (G5)
# ---------------------------------------------------------------------------

class Square(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J ** 2 * (F - inv(F).T)


class SquareFloat(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J ** 2.0 * (F - inv(F).T)


class Cube(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J ** 3 * (F - inv(F).T)


class Reciprocal(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J ** -1 * (F - inv(F).T)


class InverseSquareFloat(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J ** (-2.0) * (F - inv(F).T)


class SquareRoot(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J ** 0.5 * (F - inv(F).T)


@pytest.mark.parametrize('material, emitted', [
    (Square, '(J ** 2)'),
    (SquareFloat, '(J ** 2)'),
    (Cube, '(J ** 3)'),
    (Reciprocal, '(J ** (-1))'),
    (InverseSquareFloat, '(J ** (-2))'),
])
def test_integer_power_uses_an_integer_exponent(material, emitted):
    """A whole-number exponent is emitted as a Fortran integer power, which is
    evaluated by multiplication instead of a complex pow (exp(n log z))."""
    source = stress_source(material())
    assert emitted in source
    assert not re.search(r'\*\* DCMPLX\(', source)


class LoopIndexReciprocal(au.Material):
    props = dict(G=1.0)

    def stress_PK1(self, F):
        s = 0.0
        for k in range(1, 4):
            s = s + k ** (-1)
        return self.G * s * (F - inv(F).T)


def test_negative_power_of_a_loop_index_is_real_division():
    """k is a Fortran INTEGER; k**(-1) must not become integer division."""
    source = stress_source(LoopIndexReciprocal())
    assert '(DBLE(k) ** (-1))' in source
    assert '(k ** (-1))' not in source


def test_fractional_power_keeps_a_complex_exponent():
    assert re.search(r'\(J \*\* DCMPLX\(5\.0+d-01, 0\.0d0\)\)', stress_source(SquareRoot()))


# ---------------------------------------------------------------------------
# det(F) in a direct return expression (G1)
# ---------------------------------------------------------------------------

def component_indexed_det(source):
    """Return det33z calls whose argument contains a component like F(ii,jj).

    det33z takes a whole 3x3 tensor. A component-indexed argument compiles
    (the call has no explicit interface) and returns NaN at run time.
    """
    bad = []
    for match in re.finditer(r'det33z\s*\(', source):
        depth, end = 0, match.end() - 1
        for end in range(match.end() - 1, len(source)):
            depth += {'(': 1, ')': -1}.get(source[end], 0)
            if depth == 0:
                break
        argument = re.sub(r'\s*\n\s*&', ' ', source[match.end():end])
        if re.search(r'[A-Za-z_]\w*\s*\(\s*\w+\s*,\s*\w+\s*\)', argument):
            bad.append(f'det33z({argument})')
    return bad


class DirectReturn(au.Material):
    props = dict(G=1.0, K=10.0)

    def stress_PK1(self, F):
        return self.G * (F - inv(F).T) + self.K * log(det(F)) * inv(F).T


class DirectReturnExpression(au.Material):
    props = dict(G=1.0, K=10.0)

    def stress_PK1(self, F):
        return self.G * (F - inv(F).T) + self.K * log(det(2.0 * F) / 8.0) * inv(F).T


def test_detector_flags_the_historical_bug():
    """Control: the detector must see the generator bug it guards against."""
    assert component_indexed_det('X = LOG(det33z(F(ii,jj)))')
    assert component_indexed_det('X = det33z(F(ii,jj) - G(ii,jj))')
    assert component_indexed_det('X = det33z(F(ii,jj)\n     &  - G(ii,jj))')
    assert not component_indexed_det('X = det33z(F) * det33z((F - G))')


def test_tensorize_shields_det_arguments():
    assert ug._tensorize_expr('(LOG(det33z(F)) * F)', {'F'}) == '(LOG(det33z(F)) * F(ii,jj))'
    assert ug._tensorize_expr('det33z((F - G)) * G', {'F', 'G'}) == 'det33z((F - G)) * G(ii,jj)'


@pytest.mark.parametrize('material', [DirectReturn, DirectReturnExpression])
def test_direct_return_det_passes_the_whole_tensor(material):
    source = stress_source(material())
    assert 'det33z(' in source
    assert component_indexed_det(source) == []


def test_committed_fortran_has_no_component_indexed_det():
    """The documented grep, run over every committed Fortran source."""
    sources = [p for top in ('examples', 'paper_examples', 'abaqus_ufl')
               for p in (REPO / top).rglob('*.for')]
    assert sources
    offenders = {str(p.relative_to(REPO)): bad for p in sources
                 if (bad := component_indexed_det(p.read_text(errors='replace')))}
    assert offenders == {}


# ---------------------------------------------------------------------------
# Comparisons: both operands must be emitted as real values. Fortran has no
# ordering of complex numbers, and constants and properties are DCMPLX(...)
# in arithmetic. Before 2026-11-11 a property or a constant on the left, and
# any expression on the right, were emitted complex and the file did not
# compile (found with `if self.matflag > 2.5:` in a multi-material element).
# ---------------------------------------------------------------------------

class BranchOnProperties(au.Material):
    props = dict(G=1.0, H=0.5, flag=3.0)

    def stress_PK1(self, F):
        J = det(F)
        Jr = J.real
        s = 1.0
        if self.flag > 2.5:
            s = s + 1.0
        if self.G * 2 > 1.0:
            s = s + 1.0
        if 2.0 < J.real:
            s = s + 1.0
        if -self.H < 0.0:
            s = s + 1.0
        if J.real > self.G * 2:
            s = s + 1.0
        if Jr > self.H + 1.0:
            s = s + 1.0
        if J.real > self.G:
            s = s + 1.0
        if abs(J.real) > self.H:
            s = s + 1.0
        if J.real > Jr:
            s = s + 1.0
        return self.G * s * (F - inv(F).T)


EXPECTED_COMPARISONS = [
    'IF (props(3) .GT. 2.500000000000000d+00) THEN',
    'IF (DBLE((DCMPLX(props(1), 0.0d0) * DCMPLX(2.0d0, 0.0d0))) .GT. 1.000000000000000d+00) THEN',
    'IF (2.000000000000000d+00 .LT. DBLE(J)) THEN',
    'IF (DBLE((-DCMPLX(props(2), 0.0d0))) .LT. 0.000000000000000d+00) THEN',
    'IF (DBLE(J) .GT. DBLE((DCMPLX(props(1), 0.0d0) * DCMPLX(2.0d0, 0.0d0)))) THEN',
    'IF (Jr .GT. DBLE((DCMPLX(props(2), 0.0d0) + DCMPLX(1.000000000000000d+00, 0.0d0)))) THEN',
    # operands that were already real are emitted exactly as before
    'IF (DBLE(J) .GT. props(1)) THEN',
    'IF (ABS(DBLE(J)) .GT. props(2)) THEN',
    'IF (DBLE(J) .GT. Jr) THEN',
]


def _joined(source):
    """Fixed-form source with continuation lines folded into one line each."""
    return re.sub(r'\n     &', '', source)


def test_comparison_operands_are_emitted_real():
    source = _joined(stress_source(BranchOnProperties()))
    emitted = [line.strip() for line in source.splitlines() if line.strip().startswith('IF (')]
    assert emitted == EXPECTED_COMPARISONS


def test_no_relational_operator_touches_a_complex_operand():
    """General form of the rule: an operand of .GT. etc. never starts with DCMPLX."""
    source = _joined(stress_source(BranchOnProperties()))
    relational = r'\.(?:GT|GE|LT|LE)\.'
    assert not re.search(r'IF \(\(?-?DCMPLX\(', source)
    assert not re.search(relational + r' \(?-?DCMPLX\(', source)


def test_umat_with_property_comparisons_compiles(tmp_path):
    if shutil.which('gfortran') is None:
        pytest.skip('gfortran not available')
    output = tmp_path / 'umat_branch.for'
    au.generate_umat(BranchOnProperties(), str(output), verify=False)
    result = subprocess.run(['gfortran', '-c', '-std=legacy', '-ffixed-line-length-none', output.name,
                             '-o', 'umat_branch.o'], cwd=str(tmp_path), capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2000:]


# ---------------------------------------------------------------------------
# A vector state variable keeps its kind in stress_PK1. It was registered as
# a scalar there and as a vector in the flux methods, so a helper that takes
# <name>_old could not be shared: "multiple argument kind signatures".
# ---------------------------------------------------------------------------

class VectorStateSharedHelper(au.Material):
    props = dict(G=1.0, K=10.0, c=2.0)
    state_vars = dict(m=np.zeros(3))

    def _relax(self, F, grad_phi, m_old):
        J = det(F)
        return m_old + self.c * J * grad_phi

    def stress_PK1(self, F, grad_phi, m_old, dt):
        J = det(F)
        m_new = self._relax(F, grad_phi, m_old)
        P = self.G * (F - inv(F).T) + self.K * (J - 1.0) * J * inv(F).T
        return P, {'m': m_new}

    def phase_storage(self, F, grad_phi):
        return 0.0 * det(F)

    def phase_flux(self, F, grad_phi, m_old):
        m_new = self._relax(F, grad_phi, m_old)
        return grad_phi - m_new


class VectorStateProblem(au.WeakForm):
    material = VectorStateSharedHelper
    ndim = 3

    def define_fields(self):
        self.u = au.VectorField('u', degree=1)
        self.phi = au.ScalarField('phi', degree=1, test='eta')

    def momentum_equation(self, v, F, grad_phi):
        return self.material.stress_PK1(F, grad_phi)

    def phase_equation(self, eta, F, grad_phi):
        return self.material.phase_storage(F, grad_phi), self.material.phase_flux(F, grad_phi)


def test_helper_with_vector_state_is_shared_by_stress_and_flux(tmp_path):
    output = tmp_path / 'uel_vector_state.for'
    au.generate_uel(VectorStateProblem(), str(output), element='Hex8', formulation='standard', verify=False)
    source = output.read_text()
    assert source.upper().count('SUBROUTINE') > 0
    if shutil.which('gfortran') is None:
        pytest.skip('gfortran not available')
    result = subprocess.run(['gfortran', '-c', '-std=legacy', '-ffixed-line-length-none', output.name,
                             '-o', 'uel_vector_state.o'], cwd=str(tmp_path), capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2000:]

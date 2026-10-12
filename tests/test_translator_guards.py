"""Translator guards: constructs that must fail closed, and emission rules.

Each guard is tested against the construct it exists to stop, so removing a
guard fails a test (audit 2026-10-10: G1, G4, G5, G8).
"""

import re
from pathlib import Path

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


def test_fractional_power_keeps_a_complex_exponent():
    assert '(J ** DCMPLX(5.000000000000000d-01, 0.0d0))' in stress_source(SquareRoot())


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

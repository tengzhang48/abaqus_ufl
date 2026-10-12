"""The complex-step tangent call must pass arguments in the declared order.

In a multi-scalar element, a mismatch between the order in which the tangent
subroutine declares its history arguments and the order in which the caller
passes them still compiles and runs; it shows up only as poor Newton
convergence. Fields here are declared non-alphabetically (phi before c), so an
alphabetical side cannot agree with a declaration-order side by accident
(docs/lessons/coupled-multiphysics.md; audit 2026-10-10: G8).
"""

import re
from pathlib import Path

import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log
from abaqus_ufl.generators.uel_gen import generate_uel

REPO = Path(__file__).resolve().parents[1]


class TwoScalarMaterial(au.Material):
    props = dict(G=1.0, K=10.0, kappa=0.1, D=0.5)

    def stress_PK1(self, F, phi):
        J = det(F)
        g = (1.0 - phi) * (1.0 - phi) + 1.0e-3
        return g * (self.G * (F - inv(F).T) + self.K * log(J) * inv(F).T)

    def phase_storage(self, F, phi, c, phi_old, dt):
        return (phi - phi_old) / dt + phi - 0.5 * c

    def phase_flux(self, F, phi, grad_phi):
        return -self.kappa * grad_phi

    def species_storage(self, F, c, c_old, dt):
        return (c - c_old) / dt

    def species_flux(self, F, phi, grad_phi, grad_c):
        return -self.D * grad_c


class PhiBeforeC(au.WeakForm):
    material = TwoScalarMaterial
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField('u', degree=1)
        self.phi = au.ScalarField('phi', degree=1, test='eta')
        self.c = au.ScalarField('c', degree=1, test='zeta')

    def momentum_equation(self, v, F, phi):
        return self.material.stress_PK1(F, phi)

    def phase_equation(self, eta, F, phi, c, grad_phi, phi_old, dt):
        return (self.material.phase_storage(F, phi, c, phi_old, dt),
                self.material.phase_flux(F, phi, grad_phi))

    def species_transport_equation(self, zeta, F, phi, c, grad_phi, grad_c, c_old, dt):
        return (self.material.species_storage(F, c, c_old, dt),
                self.material.species_flux(F, phi, grad_phi, grad_c))


def argument_list(source, opener):
    """Arguments of the first call or definition starting with ``opener``."""
    start = source.index(opener) + len(opener)
    depth, end = 1, start
    while depth:
        depth += {'(': 1, ')': -1}.get(source[end], 0)
        end += 1
    text = re.sub(r'\n\s*&', ' ', source[start:end - 1])
    return [a.strip() for a in text.split(',')]


def caller_to_dummy(name):
    """Map a call-site argument to the dummy name it should bind to."""
    name = {'PROPS': 'props', 'dt_safe': 'dt'}.get(name, name)
    return re.sub(r'_gp$', '', name)


def tangent_call_pairs(source):
    pairs = []
    for prefix in sorted(set(re.findall(r'SUBROUTINE\s+(\w+)_cs_tangents\s*\(', source))):
        if f'CALL {prefix}_cs_tangents(' not in source:
            continue
        dummies = argument_list(source, f'SUBROUTINE {prefix}_cs_tangents(')
        actuals = argument_list(source, f'CALL {prefix}_cs_tangents(')
        pairs.append((prefix, dummies, actuals))
    return pairs


def test_non_alphabetical_history_order_matches(tmp_path):
    out = tmp_path / 'phi_before_c.for'
    generate_uel(PhiBeforeC(), str(out), element='Quad4', formulation='standard')
    source = out.read_text()
    [(prefix, dummies, actuals)] = tangent_call_pairs(source)
    assert [caller_to_dummy(a) for a in actuals] == dummies
    assert dummies.index('phi_old') < dummies.index('c_old')
    assert actuals.index('phi_old_gp') < actuals.index('c_old_gp')


def test_detector_catches_a_swapped_history_pair(tmp_path):
    """Control: swapping the two history arguments at the call site is caught."""
    out = tmp_path / 'phi_before_c.for'
    generate_uel(PhiBeforeC(), str(out), element='Quad4', formulation='standard')
    swapped = out.read_text().replace('phi_old_gp, c_old_gp', 'c_old_gp, phi_old_gp', 1)
    [(_, dummies, actuals)] = tangent_call_pairs(swapped)
    assert [caller_to_dummy(a) for a in actuals] != dummies


COMMITTED_UELS = sorted(
    p for top in ('examples', 'paper_examples') for p in (REPO / top).rglob('*.for')
    if '_cs_tangents(' in p.read_text(errors='replace'))


@pytest.mark.parametrize('path', COMMITTED_UELS, ids=lambda p: str(p.relative_to(REPO)))
def test_committed_tangent_calls_match_their_signatures(path):
    pairs = tangent_call_pairs(path.read_text(errors='replace'))
    assert pairs
    for prefix, dummies, actuals in pairs:
        assert [caller_to_dummy(a) for a in actuals] == dummies, prefix

"""Operator-level regression for the storage/flux sign-pairing trap.

The generator assembles r = storage*eta - flux.grad(eta), so a scalar
equation with a gradient term produces the single-field element block

    K_ff = (d storage / d f) * M  -  grad(N) . (d flux / d grad_f) . grad(N)

(M is the mass matrix). Once the element resolves the regularization length
l_eff = sqrt(|dflux/dgrad| / |dstorage|), this block must be definite: all
positive, or all negative for a consistently negated equation. An indefinite
block means an anti-diffusive gradient term, typically a flipped flux sign.
It diverges exactly where a phase field localizes, and complex-step versus
finite-difference checks cannot see it.

verify() screens for this and emits OperatorSignWarning; it does not fail,
because mixed formulations can legitimately warn (see verify.py). The tests
below check the screen against an independent assembly of the block, with a
deliberately flipped control (audit 2026-10-10: G2).
"""

import importlib.util
import itertools
import warnings
from pathlib import Path

import numpy as np
import pytest

import abaqus_ufl as au
from abaqus_ufl.core.defs import FIELD_ARGS
from abaqus_ufl.core.verify import OperatorSignWarning, _cs_tangent

REPO = Path(__file__).resolve().parents[1]
FTEST = np.array([[1.02, 0.01, 0.0],
                  [0.005, 1.01, 0.0],
                  [0.0, 0.0, 1.0]])
PHASE_STATE = {'F': FTEST, 'd': 0.1, 'grad_d': np.array([0.05, -0.03, 0.0])}


# ---------------------------------------------------------------------------
# A regularized phase-field equation: (Gc/ell) d - div(Gc ell grad d) = 0
# ---------------------------------------------------------------------------

class PhaseField(au.Material):
    props = dict(Gc=1.0, ell=0.1)

    def phase_storage(self, d):
        return self.Gc / self.ell * d

    def phase_flux(self, grad_d):
        return -self.Gc * self.ell * grad_d


class FlippedPhaseField(PhaseField):
    def phase_flux(self, grad_d):
        return self.Gc * self.ell * grad_d          # the historical bug


class NegatedPhaseField(PhaseField):
    """The whole equation multiplied by -1: still a definite block."""

    def phase_storage(self, d):
        return -self.Gc / self.ell * d

    def phase_flux(self, grad_d):
        return self.Gc * self.ell * grad_d


def phase_material(cls):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')               # unknown-argument notices
        material = cls()
        material.refresh_metadata(field_args={
            **FIELD_ARGS,
            'd': {'shape': (), 'kind': 'scalar'},
            'grad_d': {'shape': (3,), 'kind': 'vector'},
        })
    return material


# ---------------------------------------------------------------------------
# Independent block assembly (does not use verify's screening code)
# ---------------------------------------------------------------------------

def quad4_matrices(h):
    g = 1.0 / np.sqrt(3.0)
    M = np.zeros((4, 4))
    L4 = np.zeros((4, 4, 2, 2))
    for xi, eta in itertools.product((-g, g), (-g, g)):
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)],
                              [(1 - eta), -(1 + xi)],
                              [(1 + eta), (1 + xi)],
                              [-(1 + eta), (1 - xi)]]) * (2.0 / h)
        detJ = (h / 2.0) ** 2
        M += np.outer(N, N) * detJ
        L4 += np.einsum('ai,bj->abij', dN, dN) * detJ
    return M, L4


def single_field_block(material, storage, flux, field, grad, state):
    """Assemble K_ff at h = l_eff/10 from the material's complex-step tangents."""
    s = float(np.real(_cs_tangent(material, storage, field, state)))
    Df = np.real(_cs_tangent(material, flux, grad, state))
    Df2 = 0.5 * (Df + Df.T)[:2, :2]
    l_eff = np.sqrt(np.max(np.abs(np.linalg.eigvalsh(Df2))) / abs(s))
    M, L4 = quad4_matrices(l_eff / 10.0)
    return s * M - np.einsum('abij,ij->ab', L4, Df2)


def is_indefinite(K):
    w = np.linalg.eigvalsh(0.5 * (K + K.T))
    tol = 1e-9 * np.max(np.abs(w))
    return bool((w > tol).any() and (w < -tol).any())


@pytest.mark.parametrize('cls, indefinite', [
    (PhaseField, False), (NegatedPhaseField, False), (FlippedPhaseField, True)])
def test_block_assembly_classifies_the_sign(cls, indefinite):
    K = single_field_block(phase_material(cls), 'phase_storage', 'phase_flux',
                           'd', 'grad_d', PHASE_STATE)
    assert is_indefinite(K) is indefinite


# ---------------------------------------------------------------------------
# verify() screening agrees with the assembly
# ---------------------------------------------------------------------------

def operator_sign_warnings(material, state):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        assert material.verify(state=state, verbose=False)
    return [w for w in caught if issubclass(w.category, OperatorSignWarning)]


def test_verify_warns_on_a_flipped_flux_sign():
    with pytest.warns(OperatorSignWarning, match='INDEFINITE'):
        phase_material(FlippedPhaseField).verify(state=PHASE_STATE, verbose=False)


@pytest.mark.parametrize('cls', [PhaseField, NegatedPhaseField])
def test_verify_is_silent_on_a_definite_pairing(cls):
    assert operator_sign_warnings(phase_material(cls), PHASE_STATE) == []


# ---------------------------------------------------------------------------
# A shipped example
# ---------------------------------------------------------------------------

def test_heat_diffusion_example_block_is_definite():
    spec = importlib.util.spec_from_file_location(
        '_pfop_scalar_diffusion', REPO / 'examples' / 'scalar_diffusion_uel' / 'build.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    material = module.HeatDiffusionProblem()._mat
    state = {'F': FTEST, 'F_old': np.eye(3), 'T': 0.3, 'T_old': 0.2,
             'grad_T': np.array([0.4, -0.1, 0.0]), 'dt': 0.1}
    K = single_field_block(material, 'solvent_storage', 'solvent_flux', 'T', 'grad_T', state)
    assert not is_indefinite(K)
    assert np.all(np.linalg.eigvalsh(0.5 * (K + K.T)) > 0)

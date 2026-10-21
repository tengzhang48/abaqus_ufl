"""verify() as a generation gate, and its failure report.

A deliberately wrong tangent must be rejected by verify() and must stop
generation unless the caller opts out (audit 2026-10-10: G3, G8, G12).
"""

import numpy as np
import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, eye, inv, log, trace
from abaqus_ufl.core.verify import VerificationError
from abaqus_ufl.generators.uel_gen import generate_uel
from abaqus_ufl.generators.umat_gen import generate_small_strain_umat, generate_umat

SHEAR_SKEW = np.array([[1.0, 0.3, 0.0],
                       [0.0, 1.0, 0.0],
                       [0.0, 0.0, 1.0]])


class NeoHookean(au.Material):
    props = dict(G=1.0, K=10.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * (F - inv(F).T) + self.K * log(J) * inv(F).T


class Swelling(au.Material):
    props = dict(G=1.0, beta=0.2, c0=0.7)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.refresh_metadata(field_args={**self._field_args,
                                         'c': {'kind': 'scalar', 'shape': ()}})

    def stress_PK1(self, F, c):
        return self.G * (F - eye(3)) - self.beta * (c - self.c0) * eye(3)


def test_swelling_reference_uses_supplied_concentration_without_mutating_state():
    material = Swelling()
    state = {'F': 1.2 * np.eye(3), 'c': material.c0}
    before = state['F'].copy()
    assert material.verify(state=state, verbose=False)
    assert np.array_equal(state['F'], before)
    assert state['c'] == material.c0


def test_stress_free_failure_reports_its_assumptions_even_when_silent():
    with pytest.raises(VerificationError, match=r'default reference') as failure:
        Swelling().verify(verbose=False)
    assert 'c=0.' in str(failure.value)
    with pytest.raises(VerificationError, match=r'verify\(state\)') as failure:
        Swelling().verify(state={'c': 0.4}, verbose=False)
    assert 'c=0.4' in str(failure.value)


def test_separate_and_declared_reference_states_allow_nonreference_tangent_probes():
    material = Swelling()
    assert material.verify(state={'c': 0.4}, reference_state={'c': material.c0}, verbose=False)
    material.reference_state = {'c': material.c0}
    assert material.verify(state={'c': 0.4}, verbose=False)
    with pytest.raises(VerificationError, match=r'verify\(reference_state\)'):
        material.verify(reference_state={'c': 0.0}, verbose=False)


def test_generation_accepts_a_swelling_reference_state(tmp_path):
    class Problem(au.WeakForm):
        material = Swelling
        ndim = 2

        def define_fields(self):
            self.u = au.VectorField('u', degree=1)
            self.c = au.ScalarField('c', degree=1)

        def momentum_equation(self, v, F, c):
            return self.material.stress_PK1(F, c)

    generate_uel(Problem(), str(tmp_path / 'swelling.for'),
                 element='Quad4', formulation='standard',
                 verify={'reference_state': {'c': 0.7}})


class DroppedTerm(au.Material):
    """J.real removes dJ/dF from the complex-step tangent. It translates and
    compiles, but the generated tangent is wrong."""
    props = dict(G=1.0, K=10.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * (F - inv(F).T) + self.K * log(J.real) * inv(F).T


class SkewDroppedTerm(au.Material):
    """Non-symmetric dP/dF and a wrong tangent block."""
    props = dict(G=1.0)

    def stress_PK1(self, F):
        J = det(F)
        return self.G * J.real * (F - eye(3)) @ SHEAR_SKEW


class SkewDroppedTermOptedOut(SkewDroppedTerm):
    symmetric_tangent = False


class LinearElastic(au.SmallStrainMaterial):
    props = dict(G=10.0, lam=20.0)
    state_vars = dict(a=0.0)

    def stress_update(self, sigma_old, strain_old, dstrain, a_old, dt):
        sigma_new = sigma_old + 2.0 * self.G * dstrain + self.lam * trace(dstrain) * eye(3)
        return sigma_new, {'a': a_old}


class LinearElasticDroppedTerm(au.SmallStrainMaterial):
    props = dict(G=10.0, lam=20.0)
    state_vars = dict(a=0.0)

    def stress_update(self, sigma_old, strain_old, dstrain, a_old, dt):
        sigma_new = sigma_old + 2.0 * self.G * dstrain + self.lam * trace(dstrain).real * eye(3)
        return sigma_new, {'a': a_old}


def displacement_problem(material):
    class Problem(au.WeakForm):
        ndim = 2

        def define_fields(self):
            self.u = au.VectorField('u', degree=1)

        def momentum_equation(self, v, F):
            return self.material.stress_PK1(F)

    Problem.material = material
    return Problem()


# --- verify() itself ------------------------------------------------------

def test_verify_accepts_a_consistent_tangent():
    assert NeoHookean().verify(verbose=False)


def test_verify_rejects_a_wrong_tangent():
    with pytest.raises(VerificationError, match=r'd\(stress_PK1\)/d\(F\)'):
        DroppedTerm().verify(verbose=False)


def test_failure_report_lists_symmetry_when_it_is_checked():
    with pytest.raises(VerificationError) as failure:
        SkewDroppedTerm().verify(verbose=False)
    assert 'major symmetry dP/dF' in str(failure.value)


def test_failure_report_respects_the_symmetric_tangent_opt_out():
    """A material that declares symmetric_tangent = False is not checked for
    major symmetry, so the failure report must not list it (G12)."""
    with pytest.raises(VerificationError) as failure:
        SkewDroppedTermOptedOut().verify(verbose=False)
    report = str(failure.value)
    assert 'd(stress_PK1)/d(F)' in report
    assert 'major symmetry' not in report


# --- generation runs verify() by default (G3) ------------------------------

def test_generate_umat_verifies_first(tmp_path):
    out = tmp_path / 'umat.for'
    with pytest.raises(VerificationError):
        generate_umat(DroppedTerm(), str(out))
    assert not out.exists()
    generate_umat(DroppedTerm(), str(out), verify=False)
    assert out.exists()
    generate_umat(NeoHookean(), str(tmp_path / 'ok.for'))


def test_generate_small_strain_umat_verifies_first(tmp_path):
    out = tmp_path / 'umat.for'
    with pytest.raises(VerificationError):
        generate_small_strain_umat(LinearElasticDroppedTerm(), str(out))
    assert not out.exists()
    generate_small_strain_umat(LinearElasticDroppedTerm(), str(out), verify=False)
    assert out.exists()
    generate_small_strain_umat(LinearElastic(), str(tmp_path / 'ok.for'))


def test_generate_uel_verifies_first(tmp_path):
    out = tmp_path / 'uel.for'
    with pytest.raises(VerificationError):
        generate_uel(displacement_problem(DroppedTerm), str(out), element='Quad4', formulation='standard')
    assert not out.exists()
    generate_uel(displacement_problem(DroppedTerm), str(out), element='Quad4', formulation='standard',
                 verify=False)
    assert out.exists()
    generate_uel(displacement_problem(NeoHookean), str(tmp_path / 'ok.for'), element='Quad4',
                 formulation='standard')


class VerifyCalled(Exception):
    pass


class Unverified:
    """Stands in for a weak form; only records that verify() was called."""

    def verify(self, **kwargs):
        raise VerifyCalled(kwargs)


@pytest.mark.parametrize('formulation', ['standard', 'local_pressure'])
def test_every_uel_route_verifies_before_generating(tmp_path, formulation):
    with pytest.raises(VerifyCalled):
        generate_uel(Unverified(), str(tmp_path / 'uel.for'), element='Quad4', formulation=formulation)


def test_local_pressure_generator_verifies_before_generating(tmp_path):
    from abaqus_ufl.generators.uel_local_pressure import generate_uel_local_pressure
    with pytest.raises(VerifyCalled):
        generate_uel_local_pressure(Unverified(), str(tmp_path / 'uel.for'))


def test_uel_generator_passes_a_model_specific_verification_state(tmp_path):
    """A model whose representative state is not the default (e.g. the
    corrosion example at phi = 0.72) passes its own state and tolerance."""
    state = {'F': np.eye(3)}
    with pytest.raises(VerifyCalled) as called:
        generate_uel(Unverified(), str(tmp_path / 'uel.for'), element='Quad4',
                     formulation='standard', verify=dict(state=state, tol=1e-5))
    kwargs = called.value.args[0]
    assert kwargs['state'] is state and kwargs['tol'] == 1e-5 and kwargs['verbose'] is False


def test_umat_generator_passes_a_model_specific_verification_state(tmp_path):
    material = NeoHookean()
    seen = {}

    def record(**kwargs):
        seen.update(kwargs)
        return True

    material.verify = record
    state = material._make_state()
    generate_umat(material, str(tmp_path / 'umat.for'), verify=dict(state=state, tol=1e-5))
    assert seen['state'] is state and seen['tol'] == 1e-5


@pytest.mark.parametrize('value', [None, 0, 'yes'])
def test_verify_option_must_be_bool_or_dict(tmp_path, value):
    with pytest.raises(TypeError, match='verify must be True, False, or a dict'):
        generate_umat(NeoHookean(), str(tmp_path / 'umat.for'), verify=value)


def test_build_compiled_uel_verifies_and_forwards_the_option(monkeypatch):
    from abaqus_ufl.fe import compiled_element

    with pytest.raises(VerificationError):
        compiled_element.build_compiled_uel(displacement_problem(DroppedTerm))
    problem = displacement_problem(NeoHookean)
    seen = {}

    def record(**kwargs):
        seen.update(kwargs)
        raise VerifyCalled()

    monkeypatch.setattr(problem, 'verify', record)
    state = {'F': np.eye(3)}
    with pytest.raises(VerifyCalled):
        compiled_element.build_compiled_uel(problem, verify=dict(state=state, tol=1e-5))
    assert seen['state'] is state and seen['tol'] == 1e-5

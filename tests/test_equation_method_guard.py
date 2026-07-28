"""Guard regression: unrecognized storage/flux names must fail fast.

The generator emits calls to the CANONICAL sub-term material methods
(``solvent_storage``/``solvent_flux`` and friends). Before this guard, a
material defining e.g. ``heat_storage``/``heat_flux`` passed ``verify()``
and generation, and only the Fortran compile failed (calls to
subroutines that were never generated). The WeakForm constructor now
raises immediately with the expected names in the message.
"""

import numpy as np
import pytest

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log


class _BadNameMaterial(au.Material):
    props = dict(G=1.0, K=10.0, alpha=1e-3, rho_cp=1.0, k=0.5)

    def stress_PK1(self, F, T):
        J = det(F)
        FinvT = inv(F).T
        return (self.G * (F - FinvT) + self.K * log(J) * FinvT
                - self.K * self.alpha * T * FinvT)

    def heat_storage(self, F, F_old, T, T_old, dt):
        return self.rho_cp * (T - T_old) / dt

    def heat_flux(self, F, T, grad_T):
        return -self.k * grad_T


class _GoodNameMaterial(_BadNameMaterial):
    def solvent_storage(self, F, F_old, T, T_old, dt):
        return self.rho_cp * (T - T_old) / dt

    def solvent_flux(self, F, T, grad_T):
        return -self.k * grad_T


def _problem_class(material_cls):
    class Thermo(au.WeakForm):
        material = material_cls
        ndim = 2

        def define_fields(self):
            self.u = au.VectorField("u", degree=1)
            self.T = au.ScalarField("T", degree=1, test="theta")

        def momentum_equation(self, v, F, T):
            return self.material.stress_PK1(F, T)

        def transport_equation(self, theta, F, T, grad_T, F_old, T_old, dt):
            storage = self.material.solvent_storage(F, F_old, T, T_old, dt) \
                if hasattr(self.material, 'solvent_storage') else \
                self.material.heat_storage(F, F_old, T, T_old, dt)
            flux = self.material.solvent_flux(F, T, grad_T) \
                if hasattr(self.material, 'solvent_flux') else \
                self.material.heat_flux(F, T, grad_T)
            return storage, flux

    return Thermo


def test_unrecognized_storage_flux_names_raise():
    with pytest.raises(TypeError) as excinfo:
        _problem_class(_BadNameMaterial)()
    message = str(excinfo.value)
    assert "solvent_storage" in message
    assert "solvent_flux" in message
    assert "transport_equation" in message


def test_recognized_names_still_construct_and_verify():
    problem = _problem_class(_GoodNameMaterial)()
    F = np.array([
        [1.05, 0.02, 0.0],
        [0.01, 1.03, 0.0],
        [0.0, 0.0, 1.0],
    ])
    state = dict(
        F=F, F_old=0.98 * F, T=10.0, T_old=5.0,
        grad_T=np.array([2.0, -1.0, 0.5]), dt=0.1,
    )
    assert problem.verify(state=state, verbose=False)

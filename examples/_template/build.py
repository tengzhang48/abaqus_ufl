"""Template example: define a material model in Python, verify its tangent,
and generate an Abaqus UMAT (or UEL).

Copy this whole folder to `examples/<your_name>/`, rename the class and the
output `.for`, fill in the physics, then run:

    python build.py            # verifies tangents, writes <your_name>.for

This is the AUTHORING half of an example. The VALIDATION half (running the
generated .for in Abaqus and comparing against a reference) lives in ./abaqus/.
"""
from pathlib import Path

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, inv, log  # noqa: F401 (available helpers)


# --------------------------------------------------------------------------- #
# UMAT path: subclass au.Material and return the first Piola-Kirchhoff stress. #
# The framework differentiates stress_PK1 to build the consistent tangent and #
# checks it (complex-step vs finite difference) inside model.verify().         #
# --------------------------------------------------------------------------- #
class TemplateMaterial(au.Material):
    # Property names + defaults. These map, in order, to the values under
    # *User Material in abaqus/job.inp.
    props = dict(G=0.5, K=50.0)

    def stress_PK1(self, F):
        """Return P (1st PK stress) as a function of the deformation gradient F.

        Replace this neo-Hookean placeholder with your model. If your model
        has history/state, add arguments (F_old, state variables, dt, ...) and
        the framework will thread them through — see the tested examples:
          examples/small_strain_j2_umat  (plastic history)
          examples/small_strain_viscoelastic_umat  (rate/time)
        """
        J = det(F)
        FinvT = inv(F).T
        return self.G * (F - FinvT) + self.K * log(J) * FinvT


if __name__ == "__main__":
    model = TemplateMaterial()

    # Tangent verification is the FIRST gate — do not generate code that fails it.
    assert model.verify(verbose=True), "Tangent verification failed"

    out = Path(__file__).parent / "template_umat.for"
    au.generate_umat(model, str(out))
    print(f"Generated {out}")


# --------------------------------------------------------------------------- #
# UEL path (element-level weak form) instead of a UMAT?                        #
# Subclass au.WeakForm, declare fields in define_fields(), write the residual  #
# equations, call model.verify(state=...), then:                              #
#     from abaqus_ufl.generators.uel_gen import generate_uel                   #
#     generate_uel(problem, out, element="Quad4", formulation="standard")      #
# Worked, Abaqus-tested reference: examples/thermo_mechanics_quad8/            #
# --------------------------------------------------------------------------- #

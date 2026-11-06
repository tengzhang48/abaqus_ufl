"""A UEL compiled outside abaqus_ufl: build_external_uel and CompiledAbaqusElement.from_module.

The compiled checks build a generated UEL twice: once by build_compiled_uel,
and once from its Fortran source as if it were external.  Both adapters must
give the same residual, tangent and committed SVARS bit for bit, and the same
staggered solve.  The external objects are compiled at -O3, the optimization
of the meson release build that build_compiled_uel uses.
"""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from abaqus_ufl.fe import (CompiledAbaqusElement, CompiledElementState,
                           build_compiled_uel, build_external_uel,
                           solve_steps_staggered, structured_quad_mesh)
from abaqus_ufl.generators.uel_gen import generate_uel
from fe_stepping_toys import RecordingKernel
from test_fe_batch_cache import AccumulatorProblem

ROOT = Path(__file__).resolve().parents[1]


def _heat_problem():
    spec = importlib.util.spec_from_file_location(
        "_external_uel_heat", str(ROOT / "examples/scalar_diffusion_uel/build.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.HeatDiffusionProblem(G=3.5)


@pytest.fixture(scope="module")
def builds(tmp_path_factory):
    """(problem, generated module, external module) for a stateless example
    (transient heat with mechanics) and a stateful one (history accumulator)."""
    out = {}
    for name, problem in (("heat", _heat_problem()), ("accumulator", AccumulatorProblem())):
        root = tmp_path_factory.mktemp(name)
        generated = build_compiled_uel(problem, workdir=root / "generated")
        source = root / "uel.for"
        generate_uel(problem, str(source), element="Quad4", formulation="standard",
                     verify=False)
        external = build_external_uel(source, workdir=root / "external",
                                      compile_flags=("-O3",))
        out[name] = (problem, generated, external)
    return out


def _adapters(build, n_elem, dt=0.2):
    problem, generated, external = build
    gen = CompiledAbaqusElement(generated, problem, dt=dt, n_elem=n_elem)
    nsvars = 0 if gen.svars is None else gen.svars.shape[1]
    ext = CompiledAbaqusElement.from_module(
        external, dof_per_node=gen.dof_per_node, props=gen.props, n_elem=n_elem,
        nsvars=nsvars, jprops=(0,), dt=dt)
    return gen, ext


@pytest.mark.parametrize("model", ["heat", "accumulator"])
def test_external_build_matches_the_generated_adapter_bit_for_bit(builds, model):
    nodes, elems = structured_quad_mesh(2, 2)
    gen, ext = _adapters(builds[model], len(elems))
    assert sorted(vars(gen)) == sorted(vars(ext))          # one shared setup
    dpn = gen.dof_per_node
    gm = (elems[:, :, None] * dpn + np.arange(dpn)).reshape(len(elems), -1)
    rng = np.random.default_rng(7)
    U_prev = rng.normal(scale=0.03, size=len(nodes) * dpn)
    U = U_prev + rng.normal(scale=0.02, size=U_prev.shape)
    for increment in (1, 3):        # TIME(2) = 0 (initialization), then stored state
        for element in (gen, ext):
            element.begin_increment(increment, time_begin=0.2 * (increment - 1))
        for ei, el in enumerate(elems):
            a = gen.element_rk(ei, nodes[el], U[gm[ei]], (U - U_prev)[gm[ei]])
            b = ext.element_rk(ei, nodes[el], U[gm[ei]], (U - U_prev)[gm[ei]])
            assert all(np.array_equal(x, y) for x, y in zip(a, b))
        a = gen.element_rk_batch(nodes[elems], U[gm], (U - U_prev)[gm])
        b = ext.element_rk_batch(nodes[elems], U[gm], (U - U_prev)[gm])
        assert all(np.array_equal(x, y) for x, y in zip(a, b))
        for element in (gen, ext):
            element.commit_step(nodes, elems, dpn, U, U_prev)
        if model == "accumulator":
            assert np.array_equal(gen.svars, ext.svars) and np.count_nonzero(ext.svars)
        else:
            assert gen.svars is None and ext.svars is None


def test_external_adapter_drives_solve_steps_staggered(builds):
    nodes, elems = structured_quad_mesh(4, 2)
    left, right = np.flatnonzero(nodes[:, 0] == 0.0), np.flatnonzero(nodes[:, 0] == 1.0)

    def dirichlet(frac):          # (ux, uy, T): clamped left edge, heated right edge
        bc = {3 * int(n) + c: 0.0 for n in left for c in (0, 1, 2)}
        bc.update({3 * int(n) + 2: 2.0 * frac for n in right})
        return bc

    results = []
    for element in _adapters(builds["heat"], len(elems)):
        results.append(solve_steps_staggered(
            nodes, elems, 3, element.element_rk, dirichlet, [[0, 1, 2]], 3,
            batch_fn=element.element_rk_batch, commit_fn=element.commit_step,
            state=CompiledElementState([element]), dt=0.2, tol=1e-10,
            return_report=True))
    (U_gen, report_gen), (U_ext, report_ext) = results
    np.testing.assert_array_equal(U_ext, U_gen)
    assert report_ext == report_gen and report_ext.steps_accepted == 3
    assert np.abs(U_ext[2::3]).max() > 1.0 and np.abs(U_ext[0::3]).max() > 0.0


# --------------------------------------------------------------------------
# from_module without compilation
# --------------------------------------------------------------------------

def _capture_module():
    calls = []

    def drive_uel(*args):
        calls.append(args)
        return np.zeros(8), np.eye(8), np.array(args[0], copy=True), 1.0

    return SimpleNamespace(drive_uel=drive_uel, drive_uel_batch=lambda *a: None), calls


def test_from_module_passes_props_jprops_lflags_and_the_beginning_time():
    module, calls = _capture_module()
    element = CompiledAbaqusElement.from_module(
        module, dof_per_node=2, props=[3.5, 1.0], jprops=[7, 8],
        lflags=(1, 0, 1, 0, 0, 0), dt=0.25, time=1.0)
    element.element_rk(0, np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]]), np.zeros(8), np.zeros(8))
    (sv, coords, u, du, props, jprops, time, dtime, pnewdt, lflags, params,
     jtype, kstep, kinc, jelem, period) = calls[0]
    assert props.tolist() == [3.5, 1.0]
    assert jprops.dtype == np.int32 and jprops.tolist() == [7, 8]
    assert lflags.dtype == np.int32 and lflags.tolist() == [1, 0, 1, 0, 0, 0]
    assert time.tolist() == [1.0, 1.0] and dtime == 0.25 and (kstep, kinc) == (1, 1)
    assert element.svars is None and element.svars_size == 1


def test_from_module_stateless_and_stateful_elements():
    stateless = CompiledAbaqusElement.from_module(RecordingKernel(), dof_per_node=1,
                                                  props=[0.0], n_elem=2)
    assert stateless.svars is None
    nodes, elems = structured_quad_mesh(2, 1)
    U = np.linspace(0.0, 0.1, len(nodes))
    stateless.commit_step(nodes, elems, 1, U, np.zeros_like(U))
    assert stateless.svars is None

    stateful = CompiledAbaqusElement.from_module(RecordingKernel(), dof_per_node=1,
                                                 props=[0.0], n_elem=2, nsvars=1,
                                                 svars_initial=[0.25], dt=0.5)
    np.testing.assert_array_equal(stateful.svars, [[0.25], [0.25]])
    stateful.commit_step(nodes, elems, 1, U, np.zeros_like(U))
    np.testing.assert_array_equal(stateful.svars, [[0.75], [0.75]])   # + DTIME per commit
    zeros = CompiledAbaqusElement.from_module(RecordingKernel(), dof_per_node=1,
                                              props=[0.0], n_elem=3, nsvars=2)
    np.testing.assert_array_equal(zeros.svars, np.zeros((3, 2)))


def test_a_subclass_may_override_the_time_convention():
    class EndOfIncrement(CompiledAbaqusElement):
        def begin_increment(self, kinc, time_begin=None, dt=None):
            super().begin_increment(kinc, time_begin=time_begin, dt=dt)
            self.time = self.time + self.dt

    module, calls = _capture_module()
    element = EndOfIncrement.from_module(module, dof_per_node=2, props=[1.0], dt=0.25)
    assert isinstance(element, EndOfIncrement) and element.time == 0.25
    element.begin_increment(3, time_begin=0.5, dt=0.25)
    element.element_rk(0, np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]]), np.zeros(8), np.zeros(8))
    assert calls[-1][6].tolist() == [0.75, 0.75] and calls[-1][13] == 3


_FAKE = SimpleNamespace(drive_uel=lambda *a: None, drive_uel_batch=lambda *a: None)


@pytest.mark.parametrize("changes, message", [
    ({"module": SimpleNamespace(drive_uel=lambda *a: None)}, "drive_uel_batch"),
    ({"dof_per_node": 0}, "dof_per_node"),
    ({"dof_per_node": 1.5}, "dof_per_node"),
    ({"dof_per_node": True}, "dof_per_node"),
    ({"props": [[1.0]]}, "props"),
    ({"props": [1.0, np.nan]}, "props"),
    ({"props": [np.inf]}, "props"),
    ({"nsvars": -1}, "nsvars"),
    ({"nsvars": 2.0}, "nsvars"),
    ({"nsvars": True}, "nsvars"),
    ({"nsvars": 2}, "n_elem"),
    ({"nsvars": 2, "n_elem": -1}, "n_elem"),
    ({"n_elem": 1.5}, "n_elem"),
    ({"nsvars": 2, "n_elem": 1, "svars_initial": [1.0]}, "svars_initial"),
    ({"nsvars": 2, "n_elem": 1, "svars_initial": [1.0, np.nan]}, "svars_initial"),
    ({"jprops": [1.5]}, "jprops"),
    ({"jprops": [[1]]}, "jprops"),
    ({"jprops": [2 ** 31]}, "jprops"),
    ({"lflags": [1, 0, 1, 0, 0]}, "lflags"),
    ({"lflags": [1.0, 0, 1, 0, 0, 0]}, "lflags"),
    ({"dt": 0.0}, "dt"),
    ({"dt": -1.0}, "dt"),
    ({"dt": np.nan}, "dt"),
    ({"dt": np.inf}, "dt"),
    ({"dt": True}, "dt"),
    ({"time": -1.0}, "time_begin"),
    ({"time": np.nan}, "time_begin"),
])
def test_from_module_rejects_invalid_arguments(changes, message):
    arguments = dict(module=_FAKE, dof_per_node=2, props=[1.0, 2.0])
    arguments.update(changes)
    module = arguments.pop("module")
    with pytest.raises(ValueError, match=message):
        CompiledAbaqusElement.from_module(module, **arguments)


def test_build_external_uel_rejects_invalid_inputs(tmp_path):
    with pytest.raises(FileNotFoundError, match="not found"):
        build_external_uel(tmp_path / "missing.f90")
    with pytest.raises(ValueError, match="at least one"):
        build_external_uel([])
    source = tmp_path / "uel.f90"
    source.write_text("subroutine uel()\nend subroutine uel\n")
    with pytest.raises(ValueError, match="compile_flags"):
        build_external_uel(source, compile_flags="-O3")
    with pytest.raises(ValueError, match="module_name"):
        build_external_uel(source, module_name="not a name")


def test_build_external_uel_reports_a_compiler_error(tmp_path):
    source = tmp_path / "broken.f90"
    source.write_text("subroutine uel(\nend\n")
    with pytest.raises(RuntimeError, match="gfortran failed on"):
        build_external_uel(source, workdir=tmp_path / "build")
    assert "broken.f90" in (tmp_path / "build" / "compile.log").read_text()


# --------------------------------------------------------------------------
# multi-source builds, includes and rebuilds
# --------------------------------------------------------------------------

_UEL_USING_A_MODULE_AND_AN_INCLUDE = """\
subroutine uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, props, nprops, &
               coords, mcrd, nnode, u, du, v, a, jtype, time, dtime, kstep, kinc, &
               jelem, params, ndload, jdltyp, adlmag, predef, npredf, lflags, &
               mlvarx, ddlmag, mdload, pnewdt, jprops, njprop, period)
  use uel_constants, only: scale
  implicit none
  include 'factor.inc'
  integer :: ndofel, nrhs, nsvars, nprops, mcrd, nnode, jtype, kstep, kinc, jelem
  integer :: ndload, npredf, mlvarx, mdload, njprop
  integer :: jdltyp(mdload, *), lflags(*), jprops(*)
  double precision :: rhs(mlvarx, *), amatrx(ndofel, ndofel), svars(*), energy(8)
  double precision :: props(*), coords(mcrd, nnode), u(ndofel), du(mlvarx, *)
  double precision :: v(ndofel), a(ndofel), time(2), dtime, params(*)
  double precision :: adlmag(mdload, *), predef(2, npredf, nnode), ddlmag(mdload, *)
  double precision :: pnewdt, period
  integer :: i
  amatrx = 0d0
  do i = 1, ndofel
    rhs(i, 1) = -scale * factor * props(1) * u(i)
    amatrx(i, i) = scale * factor * props(1)
  end do
  svars(1) = svars(1) + dtime
end subroutine uel
"""


def _write_external_uel(folder, scale, factor):
    """Two sources: a Fortran module (compiled first) and the UEL using it;
    the UEL includes inc/factor.inc, found through the relative flag -Iinc."""
    (folder / "inc").mkdir(parents=True, exist_ok=True)
    (folder / "constants.f90").write_text(
        "module uel_constants\n  implicit none\n"
        f"  double precision, parameter :: scale = {scale}d0\nend module uel_constants\n")
    (folder / "inc" / "factor.inc").write_text(
        f"  double precision, parameter :: factor = {factor}d0\n")
    (folder / "uel.f90").write_text(_UEL_USING_A_MODULE_AND_AN_INCLUDE)
    return [folder / "constants.f90", folder / "uel.f90"]


def _gain(module):
    element = CompiledAbaqusElement.from_module(module, dof_per_node=1, props=[3.0],
                                                n_elem=1, nsvars=1, dt=0.5)
    r, k = element.element_rk(0, np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]]),
                              np.ones(4), np.zeros(4))
    assert np.array_equal(r, k @ np.ones(4))
    return r[0] / 3.0


def test_multi_source_build_with_a_module_and_a_relative_include(tmp_path):
    sources = _write_external_uel(tmp_path / "src", scale=2.0, factor=5.0)
    module = build_external_uel(sources, workdir=tmp_path / "build", compile_flags=["-Iinc"])
    assert _gain(module) == 10.0
    log = (tmp_path / "build" / "compile.log").read_text()
    assert log.index("constants.f90") < log.index("uel.f90")


def test_a_rebuild_with_the_same_name_and_workdir_runs_the_new_kernel(tmp_path):
    """Rebuilding must never reuse native code already loaded under that name
    or path, and must not replace an unrelated sys.modules entry."""
    import sys
    sentinel = object()
    sys.modules["external_kernel"] = sentinel
    try:
        build = dict(workdir=tmp_path / "build", compile_flags=["-Iinc"],
                     module_name="external_kernel")
        first = build_external_uel(_write_external_uel(tmp_path / "src", 2.0, 5.0), **build)
        second = build_external_uel(_write_external_uel(tmp_path / "src", 7.0, 5.0), **build)
        assert _gain(second) == 35.0, "the rebuild ran the previously loaded kernel"
        assert _gain(first) == 10.0, "the first module changed under its user"
        assert first is not second and first.__name__ != second.__name__
        assert sys.modules["external_kernel"] is sentinel
    finally:
        sys.modules.pop("external_kernel", None)

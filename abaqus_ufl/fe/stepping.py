"""Serial incremental stepping drivers: monolithic, staggered and multi-material.

Every driver here assembles through a private element loop or batch call with
the standard sign convention ``K = dR/dU`` (the compiled adapter already
converts Abaqus's ``RHS``/``AMATRX``). ``DU = U - U_previous`` is measured from
the accepted state at the start of the increment, so stateful and
rate-dependent elements see their history.

Unlike the strict public :func:`~abaqus_ufl.fe.driver.assemble`, the stepping
assembler accepts repeated node indices in one element (a collapsed Quad4, for
example at the centre of a polar mesh) and scatter-adds their contributions.
Whether such an element is numerically admissible remains the element's
responsibility.

A nonfinite element residual or tangent raises
:class:`~abaqus_ufl.fe.compiled_element.NonfiniteElementOutput`. Line-search
trials treat it as an infinite residual and backtrack; at a Newton iterate it
fails the step (or, in :func:`solve_steps_staggered`, the attempt, which may be
cut back). Any other exception from an element propagates unchanged.

Compiled elements receive the Abaqus beginning-of-increment ``TIME`` through
their ``begin_increment(kinc, time_begin, dt)`` method when a ``state`` object
(see :class:`~abaqus_ufl.fe.bounds.CompiledElementState`) is supplied.
"""

from dataclasses import asdict, dataclass, field
from typing import List, Optional

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from .compiled_element import IncrementTooLarge, NonfiniteElementOutput

__all__ = [
    "StaggeredAttempt", "StaggeredFieldReport", "StaggeredReport",
    "StaggeredSolveFailed", "StaggeredSweepReport", "assemble_groups",
    "global_dof", "solve_steps", "solve_steps_multimaterial",
    "solve_steps_staggered",
]


def _scipy_solve(K, b):
    return spsolve(K, b)


def global_dof(node, comp, dof_per_node):
    """Global index of component ``comp`` at ``node`` (interleaved layout)."""
    return node * dof_per_node + comp


# --------------------------------------------------------------------------
# input validation and the private stepping assembler
# --------------------------------------------------------------------------

def _check_positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _stepping_nodes(nodes):
    nodes = np.asarray(nodes, dtype=float)
    if nodes.ndim != 2 or not len(nodes) or not np.isfinite(nodes).all():
        raise ValueError("nodes must be a finite, nonempty two-dimensional array")
    return nodes


def _stepping_elems(elems, n_node, *, allow_empty=False):
    elems = np.asarray(elems)
    if (elems.ndim != 2 or elems.shape[1] < 1
            or (not allow_empty and not len(elems))
            or not (np.issubdtype(elems.dtype, np.integer) or elems.size == 0)
            or (elems.size and (elems.min() < 0 or elems.max() >= n_node))):
        raise ValueError(
            "elems must be a two-dimensional integer array of valid node "
            "indices, one row per element; repeated (collapsed) nodes are "
            "allowed")
    return elems.astype(int, copy=False)


def _stepping_mesh(nodes, elems, dof_per_node):
    nodes = _stepping_nodes(nodes)
    elems = _stepping_elems(elems, len(nodes))
    _check_positive_int(dof_per_node, "dof_per_node")
    return nodes, elems, len(nodes) * int(dof_per_node)


def _dof_map(elems, comps, dof_per_node):
    """``(nelem, nne*len(comps))`` global DOFs, node-major within an element."""
    comps = np.asarray(comps, dtype=int)
    return (elems[:, :, None] * dof_per_node + comps[None, None, :]).reshape(
        len(elems), elems.shape[1] * len(comps))


def _initial_vector(U0, ndof, name="U0"):
    U = np.zeros(ndof) if U0 is None else np.array(U0, dtype=float, copy=True)
    if U.shape != (ndof,):
        raise ValueError(f"{name} has shape {U.shape}, expected ({ndof},)")
    return U


def _check_element_output(R, K, n, context):
    R = np.asarray(R, dtype=float)
    K = np.asarray(K, dtype=float)
    if R.shape != n[0] or K.shape != n[1]:
        raise ValueError(
            f"{context}: element output shapes {R.shape}, {K.shape} disagree "
            f"with the DOF layout {n[0]}, {n[1]}")
    if not np.isfinite(R).all() or not np.isfinite(K).all():
        raise NonfiniteElementOutput(
            f"{context}: nonfinite element residual or tangent")
    return R, K


def _make_assembler(nodes, elems, dof_per_node, element_fn, batch_fn):
    """Return ``asm(U, DU) -> (R, K)`` over validated mesh arrays.

    With ``batch_fn(coords_all, U_all, DU_all)`` all elements are evaluated in
    one call; otherwise ``element_fn(ei, coords, U_e, DU_e)`` is called per
    element. Both scatter-add, so repeated connectivity accumulates.
    """
    gm = _dof_map(elems, np.arange(dof_per_node), dof_per_node)
    nelem, nd = gm.shape
    ndof = len(nodes) * dof_per_node
    rows_b = np.broadcast_to(gm[:, :, None], (nelem, nd, nd)).ravel()
    cols_b = np.broadcast_to(gm[:, None, :], (nelem, nd, nd)).ravel()

    if batch_fn is not None:
        def asm(U, DU):
            R_all, K_all = batch_fn(nodes[elems], U[gm], DU[gm])
            R_all, K_all = _check_element_output(
                R_all, K_all, ((nelem, nd), (nelem, nd, nd)), "batch assembly")
            R = np.zeros(ndof)
            np.add.at(R, gm.ravel(), R_all.ravel())
            K = sp.coo_matrix((K_all.ravel(), (rows_b, cols_b)),
                              shape=(ndof, ndof)).tocsr()
            return R, K
        return asm

    if element_fn is None:
        raise ValueError("pass element_fn or batch_fn")

    def asm(U, DU):
        R = np.zeros(ndof)
        rows, cols, vals = [], [], []
        for ei, (el, dofs) in enumerate(zip(elems, gm)):
            Re, Ke = _check_element_output(
                *element_fn(ei, nodes[el], U[dofs], DU[dofs]),
                ((nd,), (nd, nd)), f"element {ei}")
            # ``R[dofs] += Re`` would overwrite repeated indices of a
            # collapsed element; ``add.at`` sums them.
            np.add.at(R, dofs, Re)
            gg = np.broadcast_to(dofs, (nd, nd))
            rows.append(gg.T.ravel())
            cols.append(gg.ravel())
            vals.append(Ke.ravel())
        K = sp.coo_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))),
                          shape=(ndof, ndof)).tocsr()
        return R, K
    return asm


def _make_step_assembler(asm, external_fn, step, frac, U_prev):
    """Add an optional external residual/tangent to one load-step assembler.

    ``external_fn(step, frac, U, U_prev)`` returns ``(R_ext, K_ext)`` in the
    same global DOF layout as the internal system.  The equilibrium convention
    is ``R = R_internal - R_ext`` and ``K = K_internal - K_ext``.  ``K_ext`` may
    be ``None`` for a dead load.  Keeping this operation inside the assembler
    makes Newton, line-search trials, and final convergence checks evaluate the
    same residual.
    """
    if external_fn is None:
        return asm

    def step_asm(U, DU):
        R, K = asm(U, DU)
        R_ext, K_ext = external_fn(step, frac, U, U_prev)
        R_ext = np.asarray(R_ext, dtype=float)
        if R_ext.shape != R.shape:
            raise ValueError(
                "external_fn residual has shape "
                f"{R_ext.shape}; expected {R.shape}")
        R = R - R_ext
        if K_ext is not None:
            if sp.issparse(K_ext):
                K_ext_csr = K_ext.tocsr()
            else:
                K_ext_csr = sp.csr_matrix(np.asarray(K_ext, dtype=float))
            if K_ext_csr.shape != K.shape:
                raise ValueError(
                    "external_fn tangent has shape "
                    f"{K_ext_csr.shape}; expected {K.shape}")
            K = K - K_ext_csr
        return R, K

    return step_asm


def _apply_dirichlet_rows(K, R, con):
    """Row-replacement Dirichlet: constrained rows -> identity, R[con]=0, so the
    Newton increment is 0 there (the prescribed value lives in U).  Columns are
    left intact — fine because ``du[con]=0`` kills their contribution."""
    K = K.tolil()
    for c in con:
        K.rows[c] = [int(c)]
        K.data[c] = [1.0]
    R = R.copy()
    R[con] = 0.0
    return K.tocsr(), R


def _assert_finite_system(context, R, K):
    """Fail before the sparse solve if a residual/tangent entry is NaN/Inf."""
    r_bad = int(np.count_nonzero(~np.isfinite(R)))
    k_data = K.data if sp.issparse(K) else np.asarray(K).ravel()
    k_bad = int(np.count_nonzero(~np.isfinite(k_data)))
    if r_bad or k_bad:
        raise NonfiniteElementOutput(
            f"{context}: non-finite residual/tangent "
            f"(R bad entries={r_bad}, K bad entries={k_bad})")


def _trial_residual_norm(asm, Ut, U_prev, is_con):
    """Free-DOF residual 2-norm at a line-search trial; inf when nonfinite."""
    try:
        Rt, _ = asm(Ut, Ut - U_prev)
    except NonfiniteElementOutput:
        return np.inf
    Rt = np.asarray(Rt, dtype=float)
    Rt[is_con] = 0.0
    if not np.all(np.isfinite(Rt)):
        return np.inf
    return np.linalg.norm(Rt)


def _dirichlet_masks(dirichlet, ndof, context):
    con = np.array(sorted(dirichlet), dtype=int)
    if con.size and (con.min() < 0 or con.max() >= ndof):
        raise ValueError(f"{context}: Dirichlet DOF index out of range")
    is_con = np.zeros(ndof, dtype=bool)
    is_con[con] = True
    return con, is_con


def _monolithic_newton(step_asm, U, U_prev, con, is_con, *, max_newton, tol,
                       line_search, linear_solve, context, verbose=False,
                       progress=""):
    """Full Newton with backtracking, committed only below ``tol``.

    Returns ``(U, iterations, residual)``.  A step is accepted only when the
    free-DOF residual 2-norm is below ``tol`` within ``max_newton`` updates;
    the state reached by the final update is assessed too.  A line search
    that cannot reduce the residual in 30 bisections raises instead of
    applying a vanishing update.
    """
    iterations = 0
    for iteration in range(1, max_newton + 1):
        iterations = iteration
        R, K = step_asm(U, U - U_prev)
        R[is_con] = 0.0
        _assert_finite_system(context, R, K)
        rnorm = np.linalg.norm(R)
        if verbose:
            print(f"  {progress} iter {iteration:02d}: |R|={rnorm:.3e}",
                  flush=True)
        if rnorm < tol:
            return U, iterations, rnorm
        Kc, Rc = _apply_dirichlet_rows(K, R, con)
        du = linear_solve(Kc, -Rc)
        if not np.all(np.isfinite(du)):
            raise RuntimeError(f"{context}: non-finite linear solve update")
        alpha = 1.0
        if line_search:
            for _ls in range(30):
                Ut = U + alpha * du
                if _trial_residual_norm(step_asm, Ut, U_prev, is_con) < rnorm:
                    break
                alpha *= 0.5
            else:
                raise RuntimeError(
                    f"{context} iteration {iteration}: line search exhausted "
                    f"(no residual decrease in 30 bisections; residual "
                    f"{rnorm:.2e}, tol {tol:.2e})")
        U = U + alpha * du
    # Assess the state the last permitted update produced.  There is no
    # looser fallback threshold: a committed step always meets ``tol``.
    R, K = step_asm(U, U - U_prev)
    R[is_con] = 0.0
    _assert_finite_system(context, R, K)
    rfinal = np.linalg.norm(R)
    # ``not isfinite`` first: ``NaN >= tol`` is False and would hide divergence.
    if not np.isfinite(rfinal) or rfinal >= tol:
        raise RuntimeError(
            f"{context} did not converge in {max_newton} Newton iterations "
            f"(residual {rfinal:.2e}, tol {tol:.2e})")
    return U, iterations, rfinal


# --------------------------------------------------------------------------
# monolithic incremental solve
# --------------------------------------------------------------------------

def solve_steps(nodes, elems, dof_per_node, element_fn, dirichlet_fn, n_steps,
                max_newton=40, tol=1e-8, line_search=True, commit_fn=None,
                linear_solve=_scipy_solve, batch_fn=None, U0=None,
                verbose=False, external_fn=None, step_callback=None,
                state=None):
    """Incremental nonlinear solve with line search and step-wise history.

    ``dirichlet_fn(frac)`` returns the Dirichlet dict (global-dof -> value) at
    load fraction ``frac in (0, 1]``.  Each step holds ``U_prev`` (the start-of-
    step state) so the element sees ``DU = U - U_prev``.
    ``commit_fn(nodes, elems, dof_per_node, U, U_prev)`` (optional) is called
    after each step converges, to commit state variables (for example
    ``CompiledAbaqusElement.commit_step``).  ``batch_fn`` (optional) is a
    batched element evaluator (for example ``element_rk_batch``); when given,
    assembly uses one call for all elements.  ``U0`` (optional) is the initial
    state vector.  ``external_fn(step, frac, U, U_prev)`` (optional) returns
    ``(R_ext, K_ext)``; equilibrium is ``R_internal - R_ext = 0`` and ``K_ext``
    may be ``None`` for a dead load.  The callback supplies the actual load at
    each step; the driver does not multiply it by ``frac``.
    ``step_callback(step, frac, U, report)`` (optional) is called after state
    commit, with a copy of ``U`` and ``report = {'iterations': int,
    'residual_norm': float}``.

    ``state`` (optional, for example a
    :class:`~abaqus_ufl.fe.bounds.CompiledElementState`) receives
    ``state.begin_increment(kinc=step)`` before each step.  Compiled elements
    then see beginning-of-increment ``TIME = (step - 1) * dt`` with their own
    ``dt``, and ``KINC = step``.  ``TIME(2) = 0`` on the first step executes a
    generated UEL's initialization branch, which ignores SVARS seeded from
    Python.

    A step is committed only when the residual 2-norm over the free DOFs is
    below ``tol`` within ``max_newton`` Newton updates.  There is no looser
    fallback threshold, so choose ``tol`` for the problem's residual scale.
    With ``line_search``, a Newton direction that fails to reduce the residual
    in 30 bisections raises instead of being applied at a vanishing step.
    There is no cutback here; see :func:`solve_steps_staggered`.  Returns the
    converged ``U``; raises if a step fails.
    """
    nodes, elems, ndof = _stepping_mesh(nodes, elems, dof_per_node)
    n_steps = _check_positive_int(n_steps, "n_steps")
    asm = _make_assembler(nodes, elems, dof_per_node, element_fn, batch_fn)
    U = _initial_vector(U0, ndof)
    for step in range(1, n_steps + 1):
        frac = step / n_steps
        if state is not None:
            state.begin_increment(kinc=step)
        U_prev = U.copy()
        dirichlet = dirichlet_fn(frac)
        con, is_con = _dirichlet_masks(dirichlet, ndof, f"solve_steps: step {step}")
        for g, v in dirichlet.items():
            U[g] = v
        step_asm = _make_step_assembler(asm, external_fn, step, frac, U_prev)
        U, iterations, rfinal = _monolithic_newton(
            step_asm, U, U_prev, con, is_con, max_newton=max_newton, tol=tol,
            line_search=line_search, linear_solve=linear_solve,
            context=f"solve_steps: step {step}", verbose=verbose,
            progress=f"step {step:04d}/{n_steps}")
        if commit_fn is not None:
            commit_fn(nodes, elems, dof_per_node, U, U_prev)
        if step_callback is not None:
            step_callback(
                step, frac, U.copy(),
                {"iterations": iterations, "residual_norm": float(rfinal)})
        if verbose:
            print(f"  step {step:04d}/{n_steps}: committed", flush=True)
    return U


# --------------------------------------------------------------------------
# staggered (alternating) solve with adaptive, transactional stepping
# --------------------------------------------------------------------------

def _component_weights(scale, ndof, dof_per_node, name):
    """Expand a ``{component: weight}`` map to a per-DOF weight vector."""
    w = np.ones(ndof)
    if not scale:
        return w
    comp = np.arange(ndof) % dof_per_node
    for c, s in dict(scale).items():
        s = float(s)
        if not np.isfinite(s) or s <= 0.0:
            raise ValueError(f"{name}: component {c} weight must be > 0; got {s}")
        w[comp == c] = s
    return w


def _field_groups(field_components, dof_per_node):
    """Validate per-node component groups: nonempty, in range, disjoint."""
    try:
        groups = [list(g) for g in field_components]
    except TypeError as exc:
        raise ValueError("field_components must be a list of component lists") from exc
    if not groups:
        raise ValueError("field_components needs at least one field group")
    seen = set()
    for g in groups:
        if not g:
            raise ValueError("field_components: a field group is empty")
        for c in g:
            if (isinstance(c, bool) or not isinstance(c, (int, np.integer))
                    or not 0 <= c < dof_per_node):
                raise ValueError(
                    f"field_components: {c!r} is not a component index in "
                    f"[0, {dof_per_node})")
            if int(c) in seen:
                raise ValueError(
                    f"field_components: component {int(c)} belongs to more "
                    "than one field group")
            seen.add(int(c))
    return [[int(c) for c in g] for g in groups], sorted(seen)


def _field_dofs(field_components, n_node, dof_per_node):
    node = np.arange(n_node)
    return [np.sort(np.concatenate([node * dof_per_node + c for c in comps]))
            for comps in field_components]


def _component_residual_norms(R, dof_per_node):
    """L2 residual of each nodal component. Mixed-unit ‖R‖₂ hides this split."""
    return [float(np.linalg.norm(R[c::dof_per_node]))
            for c in range(dof_per_node)]


def _masked_scaled_norm(vec, weight, mask):
    """L2 norm of ``weight * vec`` on ``mask``. Empty mask is 0, not NaN."""
    if not np.any(mask):
        return 0.0
    v = np.asarray(vec, dtype=float)[mask] * np.asarray(weight, dtype=float)[mask]
    if not np.all(np.isfinite(v)):
        return np.inf
    return float(np.linalg.norm(v))


class StaggeredSolveFailed(RuntimeError):
    """A staggered inner block, outer sweep, or time step failed its contract."""

    def __init__(self, message, *, reason="", report=None):
        super().__init__(message)
        self.reason = reason
        self.report = report


class _StaggeredAttemptFailed(Exception):
    """Internal: one time-increment attempt failed and may be cut back."""

    def __init__(self, reason, R=None, extra="", rounds=0, iterations=0):
        super().__init__(reason)
        self.reason = reason
        self.R = R
        self.extra = extra
        self.rounds = rounds
        self.iterations = iterations


@dataclass
class StaggeredFieldReport:
    """One field's inner solve inside one outer sweep."""

    field: int
    iterations: int
    raw_residual: float
    scaled_residual: float
    scaled_update: float
    converged: bool
    reason: str = ""


@dataclass
class StaggeredSweepReport:
    """Certificate of one alternating-minimization sweep."""

    sweep: int
    fields: List[StaggeredFieldReport] = field(default_factory=list)
    raw_residual: float = np.inf
    scaled_residuals: Optional[List[float]] = None
    scaled_updates: Optional[List[float]] = None
    component_residuals: Optional[List[float]] = None
    accepted: bool = False
    reason: str = ""


@dataclass
class StaggeredAttempt:
    """One solve attempt of one (sub-)increment -- accepted or not."""

    step: int
    frac_start: float
    frac_end: float
    dt: float
    cutback: int
    sweeps: int
    raw_residual: float
    scaled_residual: float
    scaled_update: float
    converged: bool
    reason: str = ""
    projection_rounds: int = 0
    newton_method: str = "staggered"
    newton_iterations: int = 0


@dataclass
class StaggeredReport:
    """Per-solve record of staggered attempts and accepted nominal steps."""

    attempts: List[StaggeredAttempt] = field(default_factory=list)
    steps_accepted: int = 0
    cutbacks: int = 0
    stopped_early: bool = False

    def summary(self):
        return {
            "attempts": len(self.attempts),
            "accepted": self.steps_accepted,
            "cutbacks": self.cutbacks,
            "stopped_early": self.stopped_early,
        }


def solve_steps_staggered(nodes, elems, dof_per_node, element_fn, dirichlet_fn,
                          field_components, n_steps, max_stagger=60,
                          max_field_newton=15, tol=1e-7, line_search=True,
                          commit_fn=None, linear_solve=_scipy_solve,
                          batch_fn=None, field_bounds=None, *,
                          update_tol=None, box=None, step_callback=None,
                          verbose=False, U0=None, residual_scale=None,
                          field_scale=None, state=None, dt=1.0, max_cutbacks=0,
                          min_dt=None, cutback_factor=0.5, grow_after=2,
                          max_projection_rounds=0, attempt_callback=None,
                          accepted_callback=None, restart=None,
                          early_stop=None, return_report=False,
                          feasible_trials=False, bound_atol=1e-12,
                          newton_method="staggered",
                          cutback_on_element_rejection=False):
    """Staggered (alternating) incremental solve with adaptive stepping.

    Parameters after ``field_bounds`` are keyword-only, so later additions
    cannot rebind a positional call.

    Each field in ``field_components`` (list of per-node component-index lists,
    e.g. ``[[0,1],[2]]`` for u, phi) is solved in turn with the **other fields
    held fixed**, and the sweep is iterated to convergence.  This is robust
    through softening branches (for example phase-field fracture) where the
    monolithic Newton tangent is indefinite.  The groups must be disjoint,
    nonempty and in range, and every DOF that is not prescribed must belong to
    a group (otherwise ``ValueError``); a single group of all components is
    one monolithic Newton per increment.

    A field's sub-solve is the full system with every *other* DOF added to the
    Dirichlet set (frozen), so the Newton increment is nonzero only on that
    field's free DOFs.  ``batch_fn`` (optional) routes assembly through one
    batched call.  ``step_callback`` (optional) is the scientific output hook:
    it fires only when an accepted subincrement lands on a nominal output time,
    after ``commit_fn``.  ``commit_fn`` fires on every accepted subincrement
    (history).  ``attempt_callback`` (optional) fires for accepted and rejected
    attempts, after a rejected attempt has been rolled back.
    ``accepted_callback`` receives a complete restart dictionary after every
    accepted commit (and after nominal scientific output). ``restart`` restores
    one such dictionary, including primary/state/adaptive/report history.
    ``early_stop(U, frac, state) -> bool`` (optional) is evaluated once per
    accepted NOMINAL step, after ``step_callback`` and before
    ``accepted_callback``; a truthy return ends the stepping loop and sets
    ``report.stopped_early``.  ``verbose`` prints one line per accepted
    subincrement and the reason of every failure.

    Convergence is a **per-field** contract, not a mixed-unit ``||R||_2``:

    * each field's free DOFs must pass a scaled residual gate;
    * the outer sweep must also pass a scaled update gate (fixed-point);
    * an exhausted inner Newton or line search fails the time attempt;
    * the certificate is recomputed on a freshly assembled residual.

    An attempt is committed only when these gates hold below ``tol`` (and
    ``update_tol``); there is no looser fallback.  The residual after the LAST
    allowed inner-Newton update is assessed before the block is declared
    exhausted.  The tolerances are absolute on the weighted norms; they are
    not relative to the first residual.

    ``residual_scale`` / ``field_scale`` use the ``{component: weight}``
    convention of :func:`~abaqus_ufl.fe.bounds.solve_steps_bounded`:
    residuals are multiplied by ``residual_scale`` and updates divided by
    ``field_scale``.

    **Time and state.**  ``dt`` is the length of one nominal step.  When
    ``state`` is given (see :class:`~abaqus_ufl.fe.bounds.CompiledElementState`),
    each attempt first calls ``state.begin_increment(kinc, time_begin=t,
    dt=dt_sub)`` with the beginning time ``t`` of the attempted subincrement
    (the accumulated accepted time, starting at 0) and its own length
    ``dt_sub``.  ``KINC`` counts accepted subincrements from 1: a retried
    attempt keeps its ``KINC`` and an accepted one advances it.  Restart
    dictionaries carry ``time_end`` and ``next_kinc`` so a resumed solve
    passes the same ``TIME``/``DTIME``/``KINC``.

    **Feasible set.**  Two mutually exclusive declarations are accepted:

    ``field_bounds``   STATIC ``{component: (lo, hi)}`` commit-time clip;
                       components must exist and ranges be nonempty.
    ``box``            a :class:`~abaqus_ufl.fe.bounds.BoxSpec`, which adds the
                       dynamic ``irreversible`` / ``monotone_decreasing``
                       components.  The per-DOF bounds are derived ONCE from
                       the accepted state at the start of each attempted
                       subincrement, held fixed through every inner Newton,
                       sweep and active-mask repair round, recomputed from the
                       restored state after a rejected attempt, and
                       reconstructed from the restored ``U`` on restart.  They
                       are never advanced to the current iterate -- that would
                       ratchet the constraint and freeze the field.

    Either way the clip is followed by a fresh, scaled projected-residual
    certificate, so projection cannot silently invalidate the converged state
    or commit a non-KKT bound state.  On every bounded route a prescribed
    value outside the bounds (by more than 1e-12) raises ``ValueError``, and
    neither clipping nor snapping ever changes a prescribed value.

    Adaptive stepping (``max_cutbacks > 0``) snapshots ``U`` and ``state``
    (SVARS and the increment clock) before each trial, restores both on
    failure, multiplies the increment by ``cutback_factor``, and never steps
    over a nominal output time.  After ``grow_after`` consecutive easy
    increments the increment grows back toward the nominal one.  Default
    ``max_cutbacks=0`` takes exactly one increment per nominal step.  An
    element that raises :class:`~abaqus_ufl.fe.compiled_element.IncrementTooLarge`
    (the ``PNEWDT < 1`` signal) from assembly OR from ``commit_fn`` rejects
    the whole attempt transactionally, as does a nonfinite element output at
    an assembled iterate.  Any other exception that escapes an attempt
    propagates after ``U``, SVARS and the element clock are restored.

    ``max_projection_rounds`` (default 0) bounds the semi-smooth repair of a
    failed post-clamp certificate: when the clamped state fails the projected-
    residual gate *and* the clamp actually moved DOFs, the freshly clipped
    DOFs are frozen at their bounds and the staggered loop re-runs so the
    coupled fields re-equilibrate against the active set; the certificate is
    then recomputed as before.  A frozen DOF whose fresh residual has the
    wrong KKT sign is released for the next re-sweep.  Activations and releases
    both consume the same hard round budget; acceptance still requires the
    unchanged fresh projected-residual certificate.  With the default an
    active clamp conflict rejects the attempt.

    ``feasible_trials`` (default ``False``) keeps every trial of a bounded
    staggered solve inside its box.  By default the box is applied only after
    the sweeps converge, so between sweeps and inside each field's Newton
    iteration the fields are unconstrained, and an element may be evaluated
    far outside the box (for example at a damage value above 1, where an
    interpolation function defined on [0, 1] can make a coefficient change
    sign).  With ``feasible_trials`` each field block is solved by the
    projected semismooth Newton of
    :func:`~abaqus_ufl.fe.bounds._bounded_newton`: every assembly, including
    each line-search trial, is at a clipped state; bound-active DOFs drop out
    of the Newton system; convergence is on the projected residual.  The
    per-attempt bounds are the same accepted-start arrays; the sweep gate and
    the final certificate use the projected residual.  ``bound_atol`` scales
    ``field_scale`` into the activity tolerance of these blocks, as in
    :func:`~abaqus_ufl.fe.bounds.solve_steps_bounded`.  Requires ``box`` or
    ``field_bounds``.  The field reports of these blocks give the
    unprojected residual of the block's free DOFs as ``raw_residual``.

    ``newton_method="bounded_monolithic"`` instead solves all components
    together using the projected Newton/active-set method. It shares this
    driver's bounds, physical clock, rollback, callbacks and restart cursor.
    ``max_field_newton`` limits its coupled Newton iterations.  Acceptance
    requires a fresh finite weighted KKT residual below ``tol``; ``update_tol``,
    ``max_projection_rounds`` and ``feasible_trials`` apply only to staggered
    sweeps.  This method classifies activity exactly (zero activity tolerance)
    in its Newton iteration and in its certificate, so ``bound_atol`` does not
    apply to it.  The reported update is the actual last projected Newton
    update, recorded as an observation rather than an additional convergence
    condition.  Attempt records use ``newton_iterations`` and zero alternating
    ``sweeps`` in this mode.  Saved attempt records without these fields
    default to the staggered method.

    ``cutback_on_element_rejection=True`` makes the projected Newton solves
    (``bounded_monolithic``, and each field block under ``feasible_trials``)
    reject the attempt immediately when a line-search trial requests an
    element cutback; the transaction then restores state and subdivides.  With
    the default, those solves backtrack past a rejected trial.  Unprojected
    staggered blocks always treat an element cutback at any assembly,
    including line-search trials, as an attempt failure, so the flag does not
    change them.

    Returns ``U``, or ``(U, report)`` with a :class:`StaggeredReport` when
    ``return_report``.  Raises :class:`StaggeredSolveFailed` when an
    increment fails after ``max_cutbacks`` cutbacks or falls below ``min_dt``.
    """
    from .bounds import (BoxSpec, _bounded_newton, _is_component,
                         _restore_element_state, _snapshot_element_state,
                         projected_residual)

    if newton_method not in ("staggered", "bounded_monolithic"):
        raise ValueError("newton_method must be 'staggered' or 'bounded_monolithic'")

    nodes, elems, ndof = _stepping_mesh(nodes, elems, dof_per_node)
    n_steps = _check_positive_int(n_steps, "n_steps")
    n_node = len(nodes)
    field_components, covered = _field_groups(field_components, dof_per_node)
    fdofs = _field_dofs(field_components, n_node, dof_per_node)
    in_some_field = np.isin(np.arange(ndof) % dof_per_node, covered)
    asm = _make_assembler(nodes, elems, dof_per_node, element_fn, batch_fn)
    if restart is not None and U0 is not None:
        raise ValueError("U0 and restart are mutually exclusive")
    U = _initial_vector(restart["U"] if restart is not None else U0, ndof)
    w_res = _component_weights(residual_scale, ndof, dof_per_node,
                               "residual_scale")
    w_fld = _component_weights(field_scale, ndof, dof_per_node, "field_scale")
    w_upd = 1.0 / w_fld
    utol = float(tol if update_tol is None else update_tol)
    rtol = float(tol)
    if max_cutbacks < 0:
        raise ValueError(f"max_cutbacks must be >= 0; got {max_cutbacks}")
    if max_projection_rounds < 0:
        raise ValueError(f"max_projection_rounds must be >= 0; "
                         f"got {max_projection_rounds}")
    if not (0.0 < float(cutback_factor) < 1.0):
        raise ValueError(f"cutback_factor must be in (0, 1); got {cutback_factor}")
    nominal_span = 1.0 / n_steps
    dt = float(dt)
    if min_dt is None:
        min_dt = dt * (float(cutback_factor) ** max_cutbacks) if dt > 0.0 else 0.0
    min_dt = float(min_dt)
    min_span = (min_dt / (dt * n_steps)) if dt > 0.0 else 0.0

    # FEASIBLE SET.  ``field_bounds`` is the STATIC commit-time clip: each step
    # converges as the unconstrained staggered scheme, then bounded fields are
    # clamped once before the certificate.  ``box`` is the BoxSpec form, which
    # additionally supports the DYNAMIC step-start bounds (``irreversible`` /
    # ``monotone_decreasing``) a static map cannot express: ``d_{n+1} >= d_n``
    # is a per-DOF, per-subincrement quantity.  A static box alone lets a
    # damage field heal whenever its unconstrained equilibrium falls (for
    # example when the toughness rises at frozen history).
    if box is not None and field_bounds:
        raise ValueError(
            "solve_steps_staggered: field_bounds and box are mutually "
            "exclusive; the BoxSpec form supersedes the static map")
    if box is not None and not isinstance(box, BoxSpec):
        raise TypeError("solve_steps_staggered: box must be a BoxSpec")
    comp_of = np.arange(ndof) % dof_per_node
    static_lo = np.full(ndof, -np.inf)
    static_hi = np.full(ndof, np.inf)
    if field_bounds:
        for c, (c_lo, c_hi) in field_bounds.items():
            if not _is_component(c, dof_per_node):
                raise ValueError(f"field_bounds: {c!r} is not a component index "
                                 f"in [0, {dof_per_node})")
            if not c_lo <= c_hi:
                raise ValueError(f"field_bounds: component {c} has an empty "
                                 f"range ({c_lo}, {c_hi})")
            static_lo[comp_of == c] = c_lo
            static_hi[comp_of == c] = c_hi
    bounded = bool(field_bounds) or box is not None
    if feasible_trials and not bounded:
        raise ValueError("solve_steps_staggered: feasible_trials needs a box "
                         "or field_bounds")
    at_tol = float(bound_atol) * w_fld

    def _attempt_bounds(U_start):
        """Per-DOF (lo, hi) for ONE attempt, from the accepted start state."""
        if box is not None:
            return box.arrays(ndof, dof_per_node, U_start)
        return static_lo, static_hi

    def _clamp(Uv, lo, hi, prescribed):
        """Commit-time clip of the unprescribed DOFs; prescribed values stay."""
        if not bounded:
            return Uv
        out = np.clip(Uv, lo, hi)
        out[prescribed] = Uv[prescribed]
        return out

    def _assemble_free(Uv, U_prev, is_con, context=None, where=""):
        """Assemble with constrained rows zeroed.  An element cutback or a
        nonfinite residual/tangent fails the attempt (and may be cut back)."""
        try:
            R, K = asm(Uv, Uv - U_prev)
            R = np.asarray(R, dtype=float)
            R[is_con] = 0.0
            if context is not None:
                _assert_finite_system(context, R, K)
        except IncrementTooLarge as exc:
            raise _StaggeredAttemptFailed(
                f"element requested a smaller increment: {exc}") from exc
        except NonfiniteElementOutput as exc:
            raise _StaggeredAttemptFailed(
                "non-finite residual/tangent", extra=where) from exc
        return R, K

    def _trial_scaled_norm(Ut, U_prev, is_con, free):
        """Line-search metric: a nonfinite trial is an infinite residual."""
        try:
            Rt, _ = asm(Ut, Ut - U_prev)
        except IncrementTooLarge as exc:
            raise _StaggeredAttemptFailed(
                f"element requested a smaller increment: {exc}") from exc
        except NonfiniteElementOutput:
            return np.inf
        Rt = np.asarray(Rt, dtype=float)
        Rt[is_con] = 0.0
        return _masked_scaled_norm(Rt, w_res, free)

    def _fail_msg(step, R, reason, extra=""):
        R = np.zeros(ndof) if R is None else np.asarray(R, dtype=float)
        comp = _component_residual_norms(R, dof_per_node)
        i_dom = int(np.argmax(comp)) if comp else -1
        return (f"solve_steps_staggered: step {step} {reason}"
                f"{extra} (residual {float(np.linalg.norm(R)):.2e}; "
                f"per-component |R|={['%.3e' % v for v in comp]}; "
                f"dominant comp {i_dom})")

    def _abort(step, R, reason, extra="", report=None):
        msg = _fail_msg(step, R, reason, extra)
        if verbose:
            print(msg, flush=True)
        raise StaggeredSolveFailed(msg, reason=reason, report=report)

    def _try_increment(U_work, U_prev, frac_try, step, lo_arr, hi_arr):
        """Solve one trial increment. Raises ``_StaggeredAttemptFailed``."""
        bc = dirichlet_fn(frac_try)
        bc_con, is_bc = _dirichlet_masks(
            bc, ndof, f"solve_steps_staggered: step {step}")
        for g, v in bc.items():
            U_work[g] = v
        uncovered = np.flatnonzero(~is_bc & ~in_some_field)
        if uncovered.size:
            g = int(uncovered[0])
            raise ValueError(
                f"solve_steps_staggered: step {step}: free DOF {g} (node "
                f"{g // dof_per_node}, component {g % dof_per_node}) belongs "
                "to no field group; every DOF that is not prescribed must be "
                "solved by one group")
        # Validate prescribed values on every bounded route: a bounded solve
        # cannot satisfy an incompatible Dirichlet value, and accepting it
        # would solve a different boundary-value problem without saying so.
        if (bounded or newton_method == "bounded_monolithic") and bc_con.size:
            bad = bc_con[(U_work[bc_con] < lo_arr[bc_con] - 1e-12)
                         | (U_work[bc_con] > hi_arr[bc_con] + 1e-12)]
            if bad.size:
                g = int(bad[0])
                raise ValueError(
                    f"solve_steps_staggered: step {step}: prescribed value "
                    f"{U_work[g]:.6g} at dof {g} lies outside the box "
                    f"[{lo_arr[g]:.6g}, {hi_arr[g]:.6g}] -- an infeasible "
                    "Dirichlet value cannot be satisfied by a bounded solve")

        if newton_method == "bounded_monolithic":
            try:
                result = _bounded_newton(
                    asm, U_work, U_prev, bc_con, is_bc, lo_arr, hi_arr,
                    w_res=w_res, at_tol=0.0, max_newton=max_field_newton,
                    tol=rtol, linear_solve=linear_solve, line_search=line_search,
                    context=f"solve_steps_staggered bounded: step {step}",
                    verbose=verbose,
                    cutback_on_element_rejection=cutback_on_element_rejection)
            except NonfiniteElementOutput as exc:
                raise _StaggeredAttemptFailed("non-finite residual/tangent") from exc
            if not result["converged"]:
                raise _StaggeredAttemptFailed(
                    result["reason"], iterations=result["iterations"])
            # These are the residual and projected update actually assessed
            # by the converged Newton iteration, before any state commit.
            U_work = result["U"]
            R = result["residual_vector"]
            dU = result["last_update"]
            R_cert, _, _ = projected_residual(R, U_work, lo_arr, hi_arr, at_tol=0.0)
            kkt_norm = float(np.linalg.norm(w_res * R_cert))
            if (not np.all(np.isfinite(U_work)) or not np.all(np.isfinite(R))
                    or not np.all(np.isfinite(dU))
                    or not np.isfinite(kkt_norm) or kkt_norm >= rtol):
                raise _StaggeredAttemptFailed(
                    "bounded Newton certificate failed", R=R,
                    iterations=result["iterations"])
            field_reps, scl_R, scl_U = [], [], []
            for i_f, fd in enumerate(fdofs):
                free = np.zeros(ndof, dtype=bool)
                free[fd] = True
                free &= ~is_bc
                residual = _masked_scaled_norm(R_cert, w_res, free)
                update = _masked_scaled_norm(dU, w_upd, free)
                scl_R.append(residual)
                scl_U.append(update)
                field_reps.append(StaggeredFieldReport(
                    field=i_f, iterations=result["iterations"],
                    raw_residual=float(np.linalg.norm(R[free])),
                    scaled_residual=residual, scaled_update=update, converged=True))
            rfinal = float(np.linalg.norm(R))
            certificate = StaggeredSweepReport(
                sweep=0, fields=field_reps, raw_residual=rfinal,
                scaled_residuals=scl_R, scaled_updates=scl_U,
                component_residuals=_component_residual_norms(R, dof_per_node),
                accepted=True)
            return U_work, R, certificate, rfinal, 0

        R = None
        last_sweep = None
        pinned = np.zeros(ndof, dtype=bool)   # bound-activated this attempt
        n_rounds = 0
        for stag in range(max_stagger):
            U_sweep0 = U_work.copy()
            field_reps = []
            for i_f, fd in enumerate(fdofs):
                in_field = np.zeros(ndof, dtype=bool)
                in_field[fd] = True
                free = in_field & ~is_bc & ~pinned
                con = np.nonzero(is_bc | pinned | ~in_field)[0]
                is_con = np.zeros(ndof, dtype=bool)
                is_con[con] = True
                where = f" field {i_f}"
                if feasible_trials:
                    U_block0 = U_work.copy()
                    try:
                        res = _bounded_newton(
                            asm, U_work, U_prev, con, is_con, lo_arr, hi_arr,
                            w_res=w_res, at_tol=at_tol,
                            max_newton=max_field_newton, tol=rtol,
                            linear_solve=linear_solve,
                            line_search=line_search,
                            context=(f"solve_steps_staggered: step {step} "
                                     f"field {i_f}"),
                            cutback_on_element_rejection=(
                                cutback_on_element_rejection),
                            prescribed=is_bc)
                    except NonfiniteElementOutput as exc:
                        raise _StaggeredAttemptFailed(
                            "non-finite residual/tangent", extra=where) from exc
                    if not res["converged"]:
                        try:
                            R_fail, _ = _assemble_free(res["U"], U_prev, is_con)
                        except _StaggeredAttemptFailed:
                            R_fail = None
                        raise _StaggeredAttemptFailed(
                            res.get("reason", "bounded field solve failed"),
                            R=R_fail, extra=f"{where} after "
                            f"{res['iterations']} iterations")
                    U_work = res["U"]
                    R_block = res.get("residual_vector")
                    field_reps.append(StaggeredFieldReport(
                        field=i_f, iterations=int(res["iterations"]),
                        raw_residual=(float(np.linalg.norm(R_block[free]))
                                      if R_block is not None else float("nan")),
                        scaled_residual=float(res["residual"]),
                        scaled_update=_masked_scaled_norm(
                            U_work - U_block0, w_upd, free),
                        converged=True))
                    continue
                field_ok = False
                last_upd = 0.0
                last_raw = np.inf
                last_scl = np.inf
                n_it = 0
                context = f"solve_steps_staggered: step {step} field {i_f}"
                for iteration in range(1, max_field_newton + 1):
                    n_it = iteration
                    R, K = _assemble_free(U_work, U_prev, is_con, context, where)
                    last_raw = (float(np.linalg.norm(R[free]))
                                if np.any(free) else 0.0)
                    last_scl = _masked_scaled_norm(R, w_res, free)
                    if last_scl < rtol:
                        field_ok = True
                        break
                    Kc, Rc = _apply_dirichlet_rows(K, R, con)
                    du = np.asarray(linear_solve(Kc, -Rc), dtype=float).reshape(-1)
                    if du.size != ndof:
                        raise _StaggeredAttemptFailed(
                            "linear solve returned wrong size", R=R, extra=where)
                    du[is_con] = 0.0
                    if not np.all(np.isfinite(du)):
                        raise _StaggeredAttemptFailed(
                            "non-finite linear solve update", R=R, extra=where)
                    if not np.any(du[free]):
                        raise _StaggeredAttemptFailed(
                            "zero update with residual remaining", R=R,
                            extra=where)
                    alpha = 1.0
                    if line_search:
                        ls_ok = False
                        for _ls in range(30):
                            if _trial_scaled_norm(U_work + alpha * du, U_prev,
                                                  is_con, free) < last_scl:
                                ls_ok = True
                                break
                            alpha *= 0.5
                        if not ls_ok:
                            raise _StaggeredAttemptFailed(
                                "line-search exhausted", R=R, extra=where)
                    dU = alpha * du
                    last_upd = _masked_scaled_norm(dU, w_upd, free)
                    U_work = U_work + dU
                if not field_ok:
                    # ASSESS THE LAST ALLOWED UPDATE: a block that converged
                    # on its final permitted iteration is converged, not an
                    # exhausted block and a spurious cutback.
                    R, K = _assemble_free(U_work, U_prev, is_con, context, where)
                    last_raw = (float(np.linalg.norm(R[free]))
                                if np.any(free) else 0.0)
                    last_scl = _masked_scaled_norm(R, w_res, free)
                    field_ok = last_scl < rtol
                if not field_ok:
                    raise _StaggeredAttemptFailed(
                        "inner Newton exhausted", R=R,
                        extra=f"{where} after {n_it} iterations")
                field_reps.append(StaggeredFieldReport(
                    field=i_f, iterations=n_it, raw_residual=last_raw,
                    scaled_residual=last_scl, scaled_update=last_upd,
                    converged=True))

            R, K = _assemble_free(
                U_work, U_prev, is_bc | pinned,
                f"solve_steps_staggered: step {step} sweep {stag + 1}")
            rfinal = float(np.linalg.norm(R))
            dU_sw = U_work - U_sweep0
            R_gate = R
            if feasible_trials:
                R_gate, _, _ = projected_residual(
                    R, U_work, lo_arr, hi_arr, at_tol)
            scl_R, scl_U = [], []
            for fd in fdofs:
                free = np.zeros(ndof, dtype=bool)
                free[fd] = True
                free &= ~is_bc
                scl_R.append(_masked_scaled_norm(R_gate, w_res, free))
                scl_U.append(_masked_scaled_norm(dU_sw, w_upd, free))
            eta_R = max((r / rtol for r in scl_R), default=0.0)
            eta_U = max((u / utol for u in scl_U), default=0.0)
            last_sweep = StaggeredSweepReport(
                sweep=stag + 1, fields=field_reps, raw_residual=rfinal,
                scaled_residuals=scl_R, scaled_updates=scl_U,
                component_residuals=_component_residual_norms(R, dof_per_node),
                accepted=(eta_R < 1.0 and eta_U < 1.0))
            if last_sweep.accepted:
                U_clamped = _clamp(U_work, lo_arr, hi_arr, is_bc)
                activated = U_clamped != U_work
                U_work = U_clamped
                R, _ = _assemble_free(U_work, U_prev, is_bc)
                rfinal = float(np.linalg.norm(R))
                R_cert, at_lo, at_hi = projected_residual(
                    R, U_work, lo_arr, hi_arr, at_tol=0.0)
                post_scl_R = []
                for fd in fdofs:
                    free = np.zeros(ndof, dtype=bool)
                    free[fd] = True
                    free &= ~is_bc
                    post_scl_R.append(_masked_scaled_norm(
                        R_cert, w_res, free))
                last_sweep.raw_residual = rfinal
                last_sweep.scaled_residuals = post_scl_R
                last_sweep.component_residuals = _component_residual_norms(
                    R, dof_per_node)
                if all(r < rtol for r in post_scl_R):
                    return U_work, R, last_sweep, rfinal, n_rounds
                # Primal-dual active-mask update.  Newly clipped DOFs become
                # active.  Previously pinned DOFs whose *fresh* residual has
                # the wrong KKT sign are released: R < 0 at a lower bound or
                # R > 0 at an upper bound.  Fixed (lo == hi), prescribed, and
                # just-activated DOFs cannot be released in the same round.
                fixed = at_lo & at_hi
                released = (pinned & ~is_bc & ~fixed & ~activated
                            & ((at_lo & (R < 0.0))
                               | (at_hi & (R > 0.0))))
                if ((activated.any() or released.any())
                        and n_rounds < max_projection_rounds):
                    # Re-sweep against the updated active set.  The
                    # certificate stays the sole acceptance test, and this
                    # shared budget guarantees termination even if masks
                    # cycle.
                    pinned |= activated
                    pinned[released] = False
                    n_rounds += 1
                    continue
                last_sweep.accepted = False
                last_sweep.reason = "post-clamp projected residual failed"
                raise _StaggeredAttemptFailed(
                    "post-clamp projected residual failed", R=R,
                    extra=("; scaled |Rp|="
                           f"{['%.3e' % v for v in post_scl_R]}"),
                    rounds=n_rounds)
        extra = f" after {max_stagger} sweeps"
        if last_sweep is not None:
            extra += (
                f"; scaled |R|={['%.3e' % v for v in last_sweep.scaled_residuals]}"
                f" scaled |dU|={['%.3e' % v for v in last_sweep.scaled_updates]}")
        raise _StaggeredAttemptFailed("did not converge", R=R, extra=extra)

    if restart is None:
        report = StaggeredReport()
        time_end = 0.0
        span = nominal_span
        easy_streak = 0
        start_step = 1
        restart_frac = 0.0
        restart_cutbacks = 0
        kinc = 1
    else:
        saved = restart["report"]
        attempts = [StaggeredAttempt(**item) for item in saved["attempts"]]
        report = StaggeredReport(
            attempts=attempts,
            steps_accepted=int(saved["steps_accepted"]),
            cutbacks=int(saved["cutbacks"]),
            stopped_early=bool(saved.get("stopped_early", False)))
        time_end = float(restart["time_end"])
        span = float(restart["span"])
        easy_streak = int(restart["easy_streak"])
        start_step = int(restart["step"])
        restart_frac = float(restart["frac_done"])
        restart_cutbacks = int(restart["cutbacks_this_step"])
        # Dictionaries saved without ``next_kinc`` derive it from the record
        # of accepted subincrements, which is the same count.
        kinc = int(restart.get("next_kinc",
                               1 + sum(1 for a in attempts if a.converged)))
        if not (1 <= start_step <= n_steps and
                (start_step - 1) / n_steps - 1e-14 <= restart_frac
                <= start_step / n_steps + 1e-14):
            raise ValueError("restart nominal step/fraction cursor is inconsistent")
        _restore_element_state(state, restart.get("state_snapshot"))

    stop_now = False
    for step in range(start_step, n_steps + 1):
        if stop_now:
            break
        frac_target = step / n_steps
        frac_done = restart_frac if step == start_step else (step - 1) / n_steps
        cutbacks_this_step = restart_cutbacks if step == start_step else 0

        while frac_done < frac_target - 1e-15:
            span_try = min(span, frac_target - frac_done)
            if span_try + 1e-16 < min_span and span_try < frac_target - frac_done - 1e-16:
                _abort(step, None, "increment below minimum",
                       extra=f" span={span_try:.3e} min_span={min_span:.3e}",
                       report=report)
            # Never step over a nominal output time.
            frac_try = min(frac_done + span_try, frac_target)
            dt_sub = dt * (frac_try - frac_done) * n_steps
            U_attempt_start = U.copy()
            # The dynamic bound tracks the last ACCEPTED state, so it advances
            # only when an attempt is committed.  A rolled-back attempt
            # restores ``U`` and therefore recomputes the identical bound: a
            # cutback re-runs the SAME problem at a smaller increment, it does
            # not silently ratchet the constraint.
            lo_arr, hi_arr = _attempt_bounds(U_attempt_start)
            snap = _snapshot_element_state(state)
            try:
                if state is not None and hasattr(state, "begin_increment"):
                    state.begin_increment(kinc, time_begin=time_end, dt=dt_sub)
                U_new, R, last_sweep, rfinal, n_rounds = _try_increment(
                    U.copy(), U_attempt_start, frac_try, step, lo_arr, hi_arr)
                if commit_fn is not None:
                    # COMMIT INSIDE THE TRANSACTION.  A compiled UEL
                    # re-evaluates at the converged state during commit and can
                    # still answer PNEWDT < 1 there; that must reject the whole
                    # attempt (and roll SVARS back) rather than escape as a
                    # hard error after ``U`` was already advanced.
                    try:
                        commit_fn(nodes, elems, dof_per_node, U_new,
                                  U_attempt_start)
                    except IncrementTooLarge as exc:
                        raise _StaggeredAttemptFailed(
                            f"element requested a smaller increment: {exc}",
                            R=R, extra=" during state commit") from exc
            except _StaggeredAttemptFailed as exc:
                U = U_attempt_start
                _restore_element_state(state, snap)
                rec = StaggeredAttempt(
                    step=step, frac_start=frac_done, frac_end=frac_try,
                    dt=dt_sub, cutback=cutbacks_this_step,
                    sweeps=0, raw_residual=np.inf,
                    scaled_residual=np.inf, scaled_update=np.inf,
                    converged=False, reason=exc.reason,
                    projection_rounds=exc.rounds, newton_method=newton_method,
                    newton_iterations=exc.iterations)
                report.attempts.append(rec)
                if attempt_callback is not None:
                    attempt_callback(rec)
                cutbacks_this_step += 1
                report.cutbacks += 1
                easy_streak = 0
                if cutbacks_this_step > max_cutbacks:
                    extra = exc.extra
                    if max_cutbacks:
                        extra = (f" after {max_cutbacks} cutbacks"
                                 + (exc.extra or ""))
                    _abort(step, exc.R, exc.reason, extra=extra, report=report)
                span = max(span * float(cutback_factor), min_span)
                if verbose:
                    print(f"  step {step:04d}: cutback {cutbacks_this_step} "
                          f"({exc.reason}); new span={span:.3e}", flush=True)
                continue
            except BaseException:
                # Any other exception (an invalid boundary condition or field
                # group, an element or commit error) propagates, but only
                # after the attempt's state changes are undone.
                _restore_element_state(state, snap)
                raise

            rec = StaggeredAttempt(
                step=step, frac_start=frac_done, frac_end=frac_try,
                dt=dt_sub, cutback=cutbacks_this_step,
                sweeps=last_sweep.sweep if last_sweep else 0,
                raw_residual=rfinal,
                scaled_residual=(max(last_sweep.scaled_residuals)
                                 if last_sweep and last_sweep.scaled_residuals
                                 else rfinal),
                scaled_update=(max(last_sweep.scaled_updates)
                               if last_sweep and last_sweep.scaled_updates
                               else 0.0),
                converged=True, projection_rounds=n_rounds,
                newton_method=newton_method,
                newton_iterations=(last_sweep.fields[0].iterations
                                   if newton_method == "bounded_monolithic"
                                   and last_sweep and last_sweep.fields else 0))
            report.attempts.append(rec)
            if attempt_callback is not None:
                attempt_callback(rec)

            U = U_new
            frac_done = frac_try
            time_end += dt_sub
            kinc += 1
            if cutbacks_this_step == 0:
                easy_streak += 1
                if easy_streak >= int(grow_after):
                    span = min(span / float(cutback_factor), nominal_span)
                    easy_streak = 0
            else:
                easy_streak = 0

            on_nominal = abs(frac_done - frac_target) <= 1e-14
            if on_nominal:
                report.steps_accepted += 1
                if step_callback is not None:
                    step_callback(step, frac_target, U.copy(), {
                        "sweeps": last_sweep.sweep if last_sweep else 0,
                        "residual_norm": rfinal,
                        "component_residual_norms":
                            _component_residual_norms(R, dof_per_node),
                        "scaled_residual_norms":
                            list(last_sweep.scaled_residuals)
                            if last_sweep else [],
                        "scaled_update_norms":
                            list(last_sweep.scaled_updates)
                            if last_sweep else [],
                        "update_norm": rec.scaled_update,
                        "dt": dt_sub,
                        "cutbacks": cutbacks_this_step,
                        "newton_method": newton_method,
                        "newton_iterations": rec.newton_iterations,
                    })
            if early_stop is not None and on_nominal:
                stop_now = bool(early_stop(U, frac_target, state))
                report.stopped_early = stop_now
            if accepted_callback is not None:
                accepted_callback({
                    "U": U.copy(),
                    "state_snapshot": _snapshot_element_state(state),
                    "step": step,
                    "frac_done": frac_done,
                    "time_end": time_end,
                    "next_kinc": kinc,
                    "span": span,
                    "easy_streak": easy_streak,
                    "cutbacks_this_step": cutbacks_this_step,
                    "report": {
                        "attempts": [asdict(a) for a in report.attempts],
                        "steps_accepted": report.steps_accepted,
                        "cutbacks": report.cutbacks,
                        "stopped_early": report.stopped_early,
                    },
                })
            if verbose:
                comp = _component_residual_norms(R, dof_per_node)
                n_sw = last_sweep.sweep if last_sweep else 0
                work = (f"{rec.newton_iterations} coupled Newton iterations"
                        if newton_method == "bounded_monolithic" else f"{n_sw} sweeps")
                tag = "nominal" if on_nominal else "sub"
                print(f"  step {step:04d}/{n_steps} {tag} frac={frac_done:.6g}: "
                      f"{work}, |R|={rfinal:.3e}  "
                      f"|R_c|={['%.2e' % v for v in comp]}",
                      flush=True)
            if stop_now:
                break

    if return_report:
        return U, report
    return U


# --------------------------------------------------------------------------
# mixed-element (multi-material) assembly and solve
# --------------------------------------------------------------------------

def _prepare_groups(nodes, n_node, dof_per_node, groups):
    nodes = _stepping_nodes(nodes)
    n_node = _check_positive_int(n_node, "n_node")
    _check_positive_int(dof_per_node, "dof_per_node")
    if n_node > len(nodes):
        raise ValueError("n_node exceeds the number of node coordinates")
    prepared = []
    for i, g in enumerate(groups):
        elems = _stepping_elems(g["elems"], n_node, allow_empty=True)
        comps = np.asarray(g["comps"])
        if (comps.ndim != 1 or not comps.size
                or not np.issubdtype(comps.dtype, np.integer)
                or comps.min() < 0 or comps.max() >= dof_per_node
                or len(set(comps.tolist())) != comps.size):
            raise ValueError(
                f"group {i}: comps must be distinct component indices in "
                f"[0, {dof_per_node})")
        bfn, efn = g.get("batch_fn"), g.get("element_fn")
        if bfn is None and efn is None:
            raise ValueError(f"group {i} needs batch_fn or element_fn")
        prepared.append((elems, _dof_map(elems, comps, dof_per_node), bfn, efn))
    if not prepared:
        raise ValueError("groups must contain at least one element group")
    return nodes, n_node * dof_per_node, prepared


def _assemble_prepared_groups(nodes, ndof, prepared, U, DU):
    R = np.zeros(ndof)
    rows, cols, vals = [], [], []
    for i, (elems, gm, bfn, efn) in enumerate(prepared):
        nelem, ndofel = gm.shape
        if not nelem:
            continue
        shapes = ((nelem, ndofel), (nelem, ndofel, ndofel))
        if bfn is not None:
            R_all, K_all = _check_element_output(
                *bfn(nodes[elems], U[gm], DU[gm]), shapes, f"group {i} batch")
        else:
            R_all = np.empty(shapes[0])
            K_all = np.empty(shapes[1])
            for ei in range(nelem):
                R_all[ei], K_all[ei] = _check_element_output(
                    *efn(ei, nodes[elems[ei]], U[gm[ei]], DU[gm[ei]]),
                    ((ndofel,), (ndofel, ndofel)), f"group {i} element {ei}")
        np.add.at(R, gm.ravel(), R_all.ravel())
        rows.append(np.broadcast_to(gm[:, :, None], shapes[1]).ravel())
        cols.append(np.broadcast_to(gm[:, None, :], shapes[1]).ravel())
        vals.append(K_all.ravel())
    if not vals:
        return R, sp.csr_matrix((ndof, ndof))
    K = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(ndof, ndof)).tocsr()
    return R, K


def assemble_groups(nodes, n_node, dof_per_node, groups, U, DU=None, *,
                    du_is_total=False):
    """Assemble R, K for a MIXED-element / MIXED-DOF mesh (multi-material).

    The global layout keeps a UNIFORM ``dof_per_node`` per node (e.g. 3 = u_x,
    u_y, mu).  Each *group* drives one element type over its own element list
    and declares which per-node COMPONENTS that element touches via ``comps``:
    a coupled u-mu element uses ``comps=[0,1,2]``, a u-only element uses
    ``comps=[0,1]`` and never writes the mu rows/cols.  Two formulations with
    different DOF-per-node then coexist in one uniform global system *without*
    a per-node DOF map; the extra DOF nominally exists everywhere and is inert
    where no element provides its equation.  Such DOFs carry no equation and
    MUST be prescribed (Dirichlet) by the caller, else their rows are singular.

    ``groups``: list of dicts, each with ``elems`` (M, nne), ``comps``
    (distinct component indices), and either
    ``batch_fn(coords_all, U_all, DU_all) -> (R_all, K_all)`` or
    ``element_fn(ei, coords, U_e, DU_e) -> (R_e, K_e)``.  Element DOFs are
    node-major: ``[n0c0, n0c1, ..., n1c0, ...]`` over ``comps``.  Shared nodes
    and repeated connectivity accumulate by scatter-add.  ``DU`` must be given
    (``U - U_previous``), or ``du_is_total=True`` declares a from-zero
    evaluation.  A nonfinite element output raises
    :class:`~abaqus_ufl.fe.compiled_element.NonfiniteElementOutput`.
    """
    nodes, ndof, prepared = _prepare_groups(nodes, n_node, dof_per_node, groups)
    if DU is None and not du_is_total:
        raise ValueError("pass DU=U-U_previous, or du_is_total=True for a from-zero evaluation")
    if DU is not None and du_is_total:
        raise ValueError("specify DU or du_is_total, not both")
    U = np.asarray(U, dtype=float)
    DU = U if DU is None else np.asarray(DU, dtype=float)
    if U.shape != (ndof,) or DU.shape != (ndof,):
        raise ValueError("U and DU must be vectors of the global DOF size")
    return _assemble_prepared_groups(nodes, ndof, prepared, U, DU)


def solve_steps_multimaterial(nodes, n_node, dof_per_node, groups, dirichlet_fn,
                              n_steps, max_newton=40, tol=1e-8, line_search=True,
                              linear_solve=_scipy_solve, U0=None, state=None):
    """Incremental Newton solve for a MIXED-material mesh (see :func:`assemble_groups`).

    ``groups`` is the list passed to :func:`assemble_groups`; each group may also
    carry a ``commit(U, U_prev)`` callable (for example wrapping
    ``CompiledAbaqusElement.commit_group``) that is invoked after every
    converged step to commit that element type's state variables.  ``state``
    (optional) as in :func:`solve_steps`: ``state.begin_increment(kinc=step)``
    before each step.

    The acceptance rule is that of :func:`solve_steps`: a step is committed
    only when the free-DOF residual 2-norm is below ``tol`` within
    ``max_newton`` updates (the final update's state is assessed), and an
    exhausted 30-bisection line search raises instead of applying a vanishing
    update.
    """
    nodes, ndof, prepared = _prepare_groups(nodes, n_node, dof_per_node, groups)
    n_steps = _check_positive_int(n_steps, "n_steps")
    U = _initial_vector(U0, ndof)

    def asm(Uc, DUc):
        return _assemble_prepared_groups(nodes, ndof, prepared, Uc, DUc)

    for step in range(1, n_steps + 1):
        frac = step / n_steps
        if state is not None:
            state.begin_increment(kinc=step)
        U_prev = U.copy()
        dirichlet = dirichlet_fn(frac)
        con, is_con = _dirichlet_masks(
            dirichlet, ndof, f"solve_steps_multimaterial: step {step}")
        for g, v in dirichlet.items():
            U[g] = v
        U, _, _ = _monolithic_newton(
            asm, U, U_prev, con, is_con, max_newton=max_newton, tol=tol,
            line_search=line_search, linear_solve=linear_solve,
            context=f"solve_steps_multimaterial: step {step}")
        for g in groups:
            if g.get("commit") is not None:
                g["commit"](U, U_prev)
    return U

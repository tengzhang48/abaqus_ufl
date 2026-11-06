"""Bound-constrained (KKT) incremental solve for damage/phase-type fields.

Why this module exists
----------------------
A commit-time clip of a field (``solve_steps_staggered(field_bounds=...)``
without a certificate) is a useful safeguard against cumulative drift, but it
is not a constrained solver:

* the clipped state is generally **not** an equilibrium state -- clipping moves
  ``U`` without re-establishing ``R(U)=0``, so the committed state carries an
  unbalanced residual that the next step inherits;
* it cannot express ``d >= d_previous``, because the bound that matters for
  irreversibility is a *step-start* quantity, not a static constant.

This module solves the actual first-order (KKT) problem for a box

    lo <= U <= hi,   R(U) = 0 on the inactive set,
    R_i >= 0 where U_i = lo_i,   R_i <= 0 where U_i = hi_i

which is the stationarity condition of the incremental energy over the box when
the residual is the energy gradient ``R = dPsi/dU`` -- the convention every
abaqus_ufl element already uses (``R_e = -RHS``, ``K_e = dR_e/dU_e``).  The
signs above are what makes the difference physical rather than cosmetic: at a
lower bound the solver may only refuse to *decrease* the field.

Algorithm: primal-dual active set (PDAS) / projected semismooth Newton.  Each
iteration classifies every DOF as inactive or active-at-a-bound from the pair
``(U, R)``, solves the Newton system with the active DOFs frozen exactly like
Dirichlet rows, and projects the update back into the box.  Convergence is
measured on the **projected residual**, not ``|R|``: at a correctly signed
active bound the KKT residual is zero even though ``R`` is not.  A step is
committed only when the weighted projected residual is below ``tol``.

Anti-cycling.  Projecting inside a staggered sweep while masking residuals can
stall in active-set cycling, as can PDAS for contact.  When an active set
repeats without the KKT residual improving, this solver locks the union of the
disputed sets.  The lock never releases, so a detected cycle does not
"recover" -- the frozen DOFs hold ``du = 0``, the projected residual there
generally cannot reach ``tol``, and the attempt runs out ``max_newton`` and
becomes a *cutback*.  That is the intended outcome (a cutback, never a silent
acceptance), but it is a fail-safe, not a repair.  Cycling inside the FREE
subspace of a fixed active set is out of the lock's reach by construction;
``damping`` is the lever there.

Transactional state: an attempt that fails must leave nothing behind.  ``U``,
the element SVARS and the element increment clock are snapshotted before every
attempt and restored on failure, so a cutback re-runs from the accepted state
rather than from a half-converged one.  Only an accepted step commits.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Tuple

import numpy as np
import scipy.sparse as sp

from .compiled_element import IncrementTooLarge, NonfiniteElementOutput
from .stepping import (_apply_dirichlet_rows, _assert_finite_system,
                       _check_positive_int, _component_weights,
                       _dirichlet_masks, _initial_vector, _make_assembler,
                       _make_step_assembler, _scipy_solve, _stepping_mesh)

__all__ = [
    "BoxSpec", "BoundedAttempt", "BoundedReport", "CompiledElementState",
    "IncrementTooLarge", "active_set", "activity_equivalent_vi_certificate",
    "projected_residual", "solve_steps_bounded",
]


@dataclass(frozen=True)
class BoxSpec:
    """Feasible set for a bounded solve, declared per per-node component.

    ``bounds`` maps a component index to ``(lo, hi)`` -- e.g. ``{2: (0.0, 1.0)}``
    for a damage DOF, or ``{2: (0.0, 1.0), 3: (0.0, 1.0)}`` for a phase
    variable and a damage variable.  Use ``-np.inf``/``np.inf`` for a
    one-sided bound.

    ``irreversible`` lists components whose lower bound is raised to the
    **step-start value** of the field, i.e. ``lo_i = max(lo_static, U_prev_i)``.
    That is the discrete form of ``d_{n+1} >= d_n``; it is a per-DOF, per-step
    quantity and cannot be expressed as a static box.

    ``monotone_decreasing`` is the mirror image (``hi_i = min(hi_static,
    U_prev_i)``), for fields that may only decrease (e.g. a remaining-capacity
    variable).  It exists so an irreversibility direction is always an explicit
    declaration rather than an implicit convention.
    """

    bounds: Mapping[int, Tuple[float, float]]
    irreversible: Tuple[int, ...] = ()
    monotone_decreasing: Tuple[int, ...] = ()

    def __post_init__(self):
        for c, (lo, hi) in dict(self.bounds).items():
            if not (lo <= hi):
                raise ValueError(
                    f"BoxSpec: component {c} has empty box ({lo}, {hi})")
        overlap = set(self.irreversible) & set(self.monotone_decreasing)
        if overlap:
            raise ValueError(
                f"BoxSpec: components {sorted(overlap)} declared both "
                "irreversible and monotone_decreasing")

    def arrays(self, ndof, dof_per_node, U_prev=None):
        """Per-DOF ``(lo, hi)`` arrays for this step.

        ``U_prev`` (the accepted state at step start) is required whenever a
        monotone component is declared; passing the *current iterate* instead
        would ratchet the bound during the solve and freeze the field.
        """
        comp = np.arange(ndof) % dof_per_node
        lo = np.full(ndof, -np.inf)
        hi = np.full(ndof, np.inf)
        for c, (c_lo, c_hi) in dict(self.bounds).items():
            lo[comp == c] = c_lo
            hi[comp == c] = c_hi
        if self.irreversible or self.monotone_decreasing:
            if U_prev is None:
                raise ValueError(
                    "BoxSpec.arrays: U_prev is required for a monotone "
                    "component (the dynamic bound is the STEP-START value)")
            U_prev = np.asarray(U_prev, dtype=float)
            for c in self.irreversible:
                m = comp == c
                lo[m] = np.maximum(lo[m], U_prev[m])
            for c in self.monotone_decreasing:
                m = comp == c
                hi[m] = np.minimum(hi[m], U_prev[m])
        if np.any(lo > hi + 1e-14):
            bad = int(np.argmax(lo - hi))
            raise ValueError(
                f"BoxSpec.arrays: empty box at dof {bad} "
                f"(lo={lo[bad]:.6g} > hi={hi[bad]:.6g}); a step-start state "
                "outside the static box makes the step infeasible")
        return lo, hi


def projected_residual(R, U, lo, hi, at_tol=0.0):
    """KKT (projected) residual and the at-bound masks.

    Interior DOFs keep ``R``.  A DOF sitting at its lower bound is optimal when
    the residual pushes it *up* (``R >= 0``), so only the negative part is a
    violation; the upper bound is the mirror image.  ``at_tol`` is a per-DOF
    absolute tolerance for "sitting at" a bound -- pass a scaled tolerance for
    fields whose natural magnitude is not O(1).
    """
    R = np.asarray(R, dtype=float)
    U = np.asarray(U, dtype=float)
    at_lo = U <= lo + at_tol
    at_hi = U >= hi - at_tol
    Rp = R.copy()
    Rp[at_lo] = np.minimum(R[at_lo], 0.0)
    Rp[at_hi] = np.maximum(R[at_hi], 0.0)
    Rp[at_lo & at_hi] = 0.0            # a pinned DOF (lo == hi) has no equation
    return Rp, at_lo, at_hi


def active_set(R, U, lo, hi, at_tol=0.0):
    """DOFs whose bound is *pressing* -- at the bound with a residual that would
    drive the field out of the box.  These are frozen in the Newton system; the
    complement is the inactive set the linear solve runs on."""
    _, at_lo, at_hi = projected_residual(R, U, lo, hi, at_tol)
    R = np.asarray(R, dtype=float)
    return (at_lo & (R > 0.0)) | (at_hi & (R < 0.0)) | (at_lo & at_hi)


_CLOCK_FIELDS = ("time", "dt", "kinc")


def _copy_value(value):
    return np.array(value, copy=True) if isinstance(value, np.ndarray) else value


def _discard_trial(element):
    """Drop any cached trial evaluation of a rolled-back element."""
    for name in ("reject_step", "clear_batch_cache"):
        method = getattr(element, name, None)
        if callable(method):
            method()
            return


class CompiledElementState:
    """Transactional SVARS and increment clock of compiled elements.

    ``snapshot()``/``restore(snap)`` copy each element's committed ``svars``;
    ``restore`` also discards a cached trial (``reject_step``).
    ``clock()``/``restore_clock(clock)`` copy each element's increment clock:
    beginning-of-increment ``time``, ``dt`` and ``kinc`` of
    :class:`~abaqus_ufl.fe.compiled_element.CompiledAbaqusElement`.
    ``begin_increment(kinc, time_begin=None, dt=None)`` forwards to every
    element that has that method; the stepping drivers call it before each
    attempt with the attempt's beginning time and length.  Elements without
    SVARS, a clock or these methods (for example a Python element) are carried
    as ``None`` and skipped.
    """

    def __init__(self, elements):
        self.elements = list(elements)

    def snapshot(self):
        out = []
        for e in self.elements:
            svars = getattr(e, "svars", None)
            out.append(None if svars is None else np.array(svars, copy=True))
        return out

    def restore(self, snap):
        if len(snap) != len(self.elements):
            raise ValueError("SVARS snapshot does not match the element count")
        for e, s in zip(self.elements, snap):
            if s is not None:
                e.svars = np.array(s, dtype=float, copy=True)
            _discard_trial(e)

    def clock(self):
        return [{name: _copy_value(getattr(e, name, None)) for name in _CLOCK_FIELDS}
                for e in self.elements]

    def restore_clock(self, clock):
        if len(clock) != len(self.elements):
            raise ValueError("clock snapshot does not match the element count")
        for e, entry in zip(self.elements, clock):
            for name in _CLOCK_FIELDS:
                value = entry.get(name)
                if value is not None:
                    setattr(e, name, _copy_value(value))
            _discard_trial(e)

    def begin_increment(self, kinc, time_begin=None, dt=None):
        """Set ``KINC``, beginning-of-increment ``TIME`` and ``DTIME``.

        ``time_begin`` and ``dt`` default as in the element's own method
        (``(kinc - 1) * dt`` with the element's current ``dt``).
        """
        for e in self.elements:
            begin = getattr(e, "begin_increment", None)
            if callable(begin):
                begin(kinc, time_begin=time_begin, dt=dt)


def _snapshot_element_state(state):
    """SVARS plus the increment clock, so a rejected attempt (or a restart)
    is undone exactly.  A state without ``clock()`` snapshots SVARS only."""
    if state is None:
        return None
    clock = getattr(state, "clock", None)
    return state.snapshot(), (clock() if callable(clock) else None)


def _restore_element_state(state, snap):
    if state is None or snap is None:
        return
    svars, clock = snap
    elements = getattr(state, "elements", None)
    if elements is not None:
        if len(svars) != len(elements) or (clock is not None
                                           and len(clock) != len(elements)):
            raise ValueError("restart element-state slot count mismatch")
        for i, (e, saved) in enumerate(zip(elements, svars)):
            current = getattr(e, "svars", None)
            if (saved is None) != (current is None):
                raise ValueError(f"restart SVARS presence mismatch in element slot {i}")
            if saved is not None and np.shape(saved) != np.shape(current):
                raise ValueError(f"restart SVARS shape mismatch in element slot {i}")
    state.restore(svars)
    restore_clock = getattr(state, "restore_clock", None)
    if clock is not None and callable(restore_clock):
        restore_clock(clock)


@dataclass
class BoundedAttempt:
    """One solve attempt of one (sub-)increment -- accepted or not."""

    step: int
    frac_start: float
    frac_end: float
    dt: float
    cutback: int
    iterations: int
    kkt_residual: float
    converged: bool
    n_active_lower: int
    n_active_upper: int
    cycled: bool
    reason: str = ""


@dataclass
class BoundedReport:
    """Per-step record of everything the bounded solve did."""

    attempts: List[BoundedAttempt] = field(default_factory=list)
    steps_accepted: int = 0
    cutbacks: int = 0
    stopped_early: bool = False

    @property
    def converged_attempts(self):
        return [a for a in self.attempts if a.converged]

    def summary(self):
        return {
            "attempts": len(self.attempts),
            "accepted": self.steps_accepted,
            "cutbacks": self.cutbacks,
            "stopped_early": self.stopped_early,
            "max_kkt_residual": max(
                (a.kkt_residual for a in self.attempts if a.converged),
                default=float("nan")),
        }


def _kkt_norm(asm, U, U_prev, is_con, lo, hi, w, at_tol, *,
              cutback_on_element_rejection=False):
    """Weighted projected-residual norm at a trial point (line-search metric)."""
    try:
        R, _ = asm(U, U - U_prev)
    except IncrementTooLarge:
        if cutback_on_element_rejection:
            raise
        return np.inf              # a rejected trial is simply not a candidate
    except NonfiniteElementOutput:
        return np.inf
    R = np.asarray(R, dtype=float)
    if not np.all(np.isfinite(R)):
        return np.inf
    R[is_con] = 0.0
    Rp, _, _ = projected_residual(R, U, lo, hi, at_tol)
    return float(np.linalg.norm(w * Rp))


def activity_equivalent_vi_certificate(step_asm, U, U_prev, is_con, lo, hi,
                                       w_res, at_tol, tol, linear_solve,
                                       R=None, K=None):
    """Certify an exact-bound state via a nearby box-VI point.

    The ordinary projected residual and its signs remain authoritative.  This
    opt-in certificate solves one classical PDAS correction, requires that the
    *unprojected* correction is box-feasible and componentwise no larger than
    the already-declared activity tolerance, then reassembles at that nearby
    point and requires the unchanged weighted KKT gate.  It therefore cannot
    hide a resolvable growth direction or loosen the residual tolerance.

    Returns ``None`` when any condition fails; otherwise returns the certified
    KKT norm, correction, and projected residual at the nearby VI point.
    """
    U = np.asarray(U, dtype=float)
    U_prev = np.asarray(U_prev, dtype=float)
    if R is None or K is None:
        try:
            R, K = step_asm(U, U - U_prev)
        except (IncrementTooLarge, NonfiniteElementOutput):
            return None
    R = np.asarray(R, dtype=float).copy()
    R[is_con] = 0.0
    act = active_set(R, U, lo, hi, at_tol) & ~is_con
    con_all = np.nonzero(is_con | act)[0]
    Kc, Rc = _apply_dirichlet_rows(K, R, con_all)
    if not np.any(Rc):
        return None
    du = np.asarray(linear_solve(Kc, -Rc), dtype=float)
    if not np.all(np.isfinite(du)):
        return None
    du[con_all] = 0.0
    tol_arr = np.broadcast_to(np.asarray(at_tol, dtype=float), U.shape)
    if np.any(np.abs(du) > tol_arr):
        return None
    U_vi = U + du
    if np.any(U_vi < lo) or np.any(U_vi > hi):
        return None
    try:
        R_vi, _ = step_asm(U_vi, U_vi - U_prev)
    except (IncrementTooLarge, NonfiniteElementOutput):
        return None
    R_vi = np.asarray(R_vi, dtype=float)
    if not np.all(np.isfinite(R_vi)):
        return None
    R_vi[is_con] = 0.0
    Rp_vi, at_lo, at_hi = projected_residual(
        R_vi, U_vi, lo, hi, at_tol)
    Rp_vi[is_con] = 0.0
    norm = float(np.linalg.norm(w_res * Rp_vi))
    if norm >= tol:
        return None
    return dict(norm=norm, update=du, projected=Rp_vi,
                at_lo=at_lo, at_hi=at_hi)


def _bounded_newton(step_asm, U, U_prev, con, is_con, lo, hi, *, w_res, at_tol,
                    max_newton, tol, linear_solve, line_search, context,
                    damping=1.0, levenberg=0.0,
                    activity_vi_certificate=False, verbose=False,
                    cutback_on_element_rejection=False):
    """Projected semismooth-Newton / PDAS solve of one increment.

    Up to ``max_newton`` projected updates are applied; the state each update
    produces, including the last one, is assembled and accepted when its
    weighted KKT residual is below ``tol``.  ``iterations`` counts the
    assemblies that were assessed (updates + 1 on convergence).

    ``damping`` scales the Newton update before projection (``U += damping*du``).
    It is the lever for projected-map 2-cycles / frozen plateaus at a stable
    active set, where the anti-cycling lock is powerless because the cycling
    lives in the FREE subspace: with ``damping < 1`` the free fixed-point map
    contracts (spectral radius shrinks below one for the oscillatory mode)
    while the projected-KKT acceptance is unchanged.  Default ``1.0`` is the
    undamped method.

    ``levenberg`` augments the linear system with ``mu * diag(|K_ii|)`` after
    the Dirichlet rows are applied.  It is the lever for boundary rows whose
    diagonal is small relative to their off-diagonal coupling (for example a
    phase front sitting at an irreversibility floor): the coupled full-block
    direction then points OUT of the box for a row whose own equation points
    IN, the projection clips it back every iteration, and the outward residual
    persists.  The diagonal augmentation biases each row toward its own
    equation, letting the row lift off the floor.  ``levenberg=0`` (default)
    is the plain Newton system.

    ``activity_vi_certificate`` is an opt-in acceptance path for a sharp bound
    whose exact VI point lies within ``at_tol`` of the represented active state.
    It requires a feasible coupled correction and a fresh assembly below the
    unchanged ``tol``; the default ``False`` preserves classical PDAS behavior.

    ``cutback_on_element_rejection=True`` rejects the increment immediately
    when an element requests a smaller increment at a line-search trial.
    The default treats that trial as inadmissible and continues backtracking.
    A rejected return leaves rollback/subdivision to the caller.  A nonfinite
    line-search trial is always an infinite residual; a nonfinite residual or
    tangent at an iterate raises
    :class:`~abaqus_ufl.fe.compiled_element.NonfiniteElementOutput`.

    ``U`` is modified only through the returned value; the caller owns rollback.
    """
    U_initial = np.asarray(U, dtype=float)
    U = np.clip(U_initial, lo, hi)     # feasible start
    # A bound-adjacent initial/restart value can satisfy the projected KKT
    # test on the very first assembly.  Snap it before that early return, not
    # only after a Newton update, so downstream checks with ``at_tol=0`` see
    # the same active point rather than a numerically interior slack state.
    U = np.where((U <= lo + at_tol) & np.isfinite(lo), lo, U)
    U = np.where((U >= hi - at_tol) & np.isfinite(hi), hi, U)
    last_update = U - U_initial
    seen: Dict[bytes, float] = {}
    locked = np.zeros_like(is_con)
    cycled = False
    rnorm = np.inf
    iterations = 0
    at_lo = at_hi = np.zeros_like(is_con)

    # Pass ``max_newton + 1`` only assesses the state produced by the last
    # permitted update; it never computes another update.
    for iteration in range(1, max_newton + 2):
        final = iteration > max_newton
        if not final:
            iterations = iteration
        try:
            R, K = step_asm(U, U - U_prev)
        except IncrementTooLarge as exc:
            return dict(U=U, converged=False, iterations=iterations,
                        residual=np.inf, at_lo=at_lo, at_hi=at_hi,
                        cycled=cycled,
                        reason=f"element requested a smaller increment: {exc}")
        R = np.asarray(R, dtype=float)
        R[is_con] = 0.0
        _assert_finite_system(context, R, K)

        Rp, at_lo, at_hi = projected_residual(R, U, lo, hi, at_tol)
        Rp[is_con] = 0.0
        rnorm = float(np.linalg.norm(w_res * Rp))
        if verbose:
            print(f"  {context} iter {iteration:02d}: |Rkkt|={rnorm:.3e} "
                  f"active={int((at_lo | at_hi).sum())}", flush=True)
        if rnorm < tol:
            return dict(U=U, converged=True, iterations=iteration,
                        residual=rnorm, at_lo=at_lo, at_hi=at_hi, cycled=cycled,
                        residual_vector=R.copy(), last_update=last_update.copy())

        if activity_vi_certificate:
            certificate = activity_equivalent_vi_certificate(
                step_asm, U, U_prev, is_con, lo, hi, w_res, at_tol, tol,
                linear_solve, R=R, K=K)
            if certificate is not None:
                return dict(
                    U=U + certificate["update"], converged=True,
                    iterations=iteration,
                    residual=certificate["norm"], at_lo=at_lo, at_hi=at_hi,
                    cycled=cycled, activity_vi_certified=True,
                    certified_update_norm=float(np.max(
                        np.abs(certificate["update"]), initial=0.0)),
                    reason="activity-equivalent box-VI certificate")
        if final:
            break

        act = active_set(R, U, lo, hi, at_tol) & ~is_con
        # Anti-cycling: a repeated active set that did not improve the KKT
        # residual is the cycling signature.  Lock the union from then on --
        # the locked set only grows, so the iteration cannot revisit a cycle.
        key = np.packbits(act).tobytes()
        if key in seen and rnorm >= seen[key] * (1.0 - 1e-12):
            cycled = True
        seen[key] = min(seen.get(key, np.inf), rnorm)
        if cycled:
            locked |= act
            act = act | locked

        con_all = np.nonzero(is_con | act)[0]
        Kc, Rc = _apply_dirichlet_rows(K, R, con_all)
        if levenberg > 0.0:
            Kc = Kc + sp.diags(levenberg * np.abs(Kc.diagonal()))
        du = linear_solve(Kc, -Rc)
        if not np.all(np.isfinite(du)):
            return dict(U=U, converged=False, iterations=iterations,
                        residual=np.inf, at_lo=at_lo, at_hi=at_hi,
                        cycled=cycled, reason="non-finite linear solve update")
        du[con_all] = 0.0

        alpha = float(damping)
        ls_exhausted = False
        if line_search:
            for _ls in range(30):
                Ut = np.clip(U + alpha * du, lo, hi)
                try:
                    trial_norm = _kkt_norm(
                        step_asm, Ut, U_prev, is_con, lo, hi, w_res, at_tol,
                        cutback_on_element_rejection=cutback_on_element_rejection)
                except IncrementTooLarge as exc:
                    return dict(
                        U=U, converged=False, iterations=iterations,
                        residual=np.inf, at_lo=at_lo, at_hi=at_hi, cycled=cycled,
                        reason=f"element requested a smaller increment during line search: {exc}")
                if trial_norm < rnorm:
                    break
                alpha *= 0.5
            else:
                ls_exhausted = True
        if ls_exhausted:
            # All 30 bisections failed to reduce the KKT norm: ``du`` is a
            # certified non-descent direction.  Applying the 2^-30-scaled
            # update anyway would turn one bad direction into a zero-step
            # fixed point: every further iteration costs 31 assemblies while
            # the state crawls by ~2^-30|du|.  Reject the iteration instead:
            # ``U`` still holds the pre-iteration state, and the
            # non-converged return lets the caller's cutback machinery engage.
            reason = (f"line search exhausted: no KKT decrease in 30 "
                      f"bisections (|Rkkt|={rnorm:.3e}, "
                      f"|du|={float(np.linalg.norm(du)):.3e})")
            if verbose:
                print(f"  {context} iter {iteration:02d}: {reason}; "
                      f"rejecting iteration", flush=True)
            return dict(U=U, converged=False, iterations=iterations,
                        residual=rnorm, at_lo=at_lo, at_hi=at_hi,
                        cycled=cycled, reason=reason)
        # PROJECTION is part of the iteration, not a post-processing step: the
        # next assembly (and therefore the accepted state) is evaluated at the
        # projected point, so a converged solve is a KKT point of the assembled
        # system rather than a clipped copy of an unconstrained one.
        U_before_update = U
        U = np.clip(U + alpha * du, lo, hi)
        # SNAP to the bound inside the activity tolerance.  A DOF can
        # otherwise land a few ulps short of its bound and be judged "at the
        # bound" here (tolerance ``at_tol``) but "interior, with a large
        # residual" by any downstream check that uses a tighter tolerance --
        # i.e. convergence would rest on the slack rather than on the KKT
        # condition.  Snapping moves ``U`` by at most ``at_tol`` and makes the
        # accepted state exactly feasible and exactly active.
        U = np.where((U <= lo + at_tol) & np.isfinite(lo), lo, U)
        U = np.where((U >= hi - at_tol) & np.isfinite(hi), hi, U)
        last_update = U - U_before_update

    return dict(U=U, converged=False, iterations=iterations, residual=rnorm,
                at_lo=at_lo, at_hi=at_hi, cycled=cycled,
                reason=f"no KKT convergence in {max_newton} iterations")


def solve_steps_bounded(nodes, elems, dof_per_node, element_fn, dirichlet_fn,
                        n_steps, box, *, batch_fn=None, U0=None, commit_fn=None,
                        state=None, max_newton=40, tol=1e-8,
                        linear_solve=_scipy_solve, residual_scale=None,
                        field_scale=None, bound_atol=1e-12, line_search=True,
                        max_cutbacks=4, dt=1.0, external_fn=None,
                        step_callback=None, attempt_callback=None,
                        verbose=False, early_stop=None):
    """Incremental solve subject to a box (and optional irreversibility).

    Same contract as :func:`~abaqus_ufl.fe.stepping.solve_steps` -- element
    callable, ``dirichlet_fn(frac)``, optional ``batch_fn``/``U0``/
    ``external_fn``/``commit_fn`` -- with these additions:

    ``box``            :class:`BoxSpec` describing the feasible set.
    ``state``          transactional element state (see
                       :class:`CompiledElementState`): ``snapshot()``,
                       ``restore(snap)``, optional ``clock()``/
                       ``restore_clock()`` and optional
                       ``begin_increment(kinc, time_begin, dt)``.  Snapshot
                       and rollback of SVARS and clock happen per *attempt*,
                       so a failed attempt cannot leave committed state
                       behind.  Each attempt sets ``TIME`` to its beginning
                       time and ``DTIME`` to its own length; ``KINC`` counts
                       accepted subincrements from 1 and is kept on a retry.
    ``residual_scale`` ``{component: weight}`` applied before every residual
                       norm.  Coupled residuals can span many decades (for
                       example mechanics O(1e3) next to transport O(1e-11));
                       without explicit weights the small fields never enter
                       the norm and the gate is vacuous for them.
    ``field_scale``    ``{component: magnitude}`` -- the characteristic size of
                       each field, used to decide when a DOF is *at* a bound
                       (``bound_atol * scale``).  Activity detection must be
                       scale-aware for the same reason.
    ``max_cutbacks``   how many times a failed increment may be halved before
                       the step is declared failed.  Each cutback rolls ``U``
                       and the element state back and re-runs from the accepted
                       state at the reduced ``dt``.
    ``dt``             step time increment (the full step; sub-increments scale
                       it), reported in the attempt records.
    ``early_stop``     optional callable ``(U, frac, state) -> bool`` evaluated
                       once per accepted step, after ``step_callback``; a
                       truthy return ends the stepping loop and sets
                       ``report.stopped_early``.  The caller defines the
                       stopping criterion and the resulting claim boundary.
                       The same hook exists on
                       :func:`~abaqus_ufl.fe.stepping.solve_steps_staggered`.

    A step is committed only when its weighted KKT residual is below ``tol``.
    An element cutback (``IncrementTooLarge``) at an iterate, a line-search or
    commit evaluation, an exhausted line search, or Newton exhaustion is a
    cutback.  Returns ``(U, report)`` with ``report`` a :class:`BoundedReport`.
    Raises ``RuntimeError`` if a step cannot be completed within
    ``max_cutbacks``.
    """
    nodes, elems, ndof = _stepping_mesh(nodes, elems, dof_per_node)
    n_steps = _check_positive_int(n_steps, "n_steps")
    if not isinstance(box, BoxSpec):
        raise TypeError("box must be a BoxSpec")
    asm = _make_assembler(nodes, elems, dof_per_node, element_fn, batch_fn)
    w_res = _component_weights(residual_scale, ndof, dof_per_node, "residual_scale")
    w_fld = _component_weights(field_scale, ndof, dof_per_node, "field_scale")
    at_tol = float(bound_atol) * w_fld

    U = _initial_vector(U0, ndof)
    report = BoundedReport()
    time_end = 0.0
    kinc = 1

    for step in range(1, n_steps + 1):
        frac_done = (step - 1) / n_steps
        frac_target = step / n_steps
        span = frac_target - frac_done
        cutback = 0

        while frac_done < frac_target - 1e-15:
            frac_try = min(frac_done + span, frac_target)
            # dt of the SUB-INCREMENT: measured against a whole step, not
            # against what is left of it -- otherwise the last sub-increment
            # after a cutback reports (and hands the element) the full dt.
            dt_sub = float(dt) * (frac_try - frac_done) * n_steps
            U_attempt_start = U.copy()
            # The dynamic bound tracks the last ACCEPTED state, so it advances
            # only when an attempt is committed.  A rolled-back attempt restores
            # ``U`` and therefore recomputes the identical bound -- a cutback
            # re-runs the same problem at a smaller increment, it does not
            # silently ratchet the constraint.
            lo, hi = box.arrays(ndof, dof_per_node, U_attempt_start)
            snap = _snapshot_element_state(state)
            if state is not None and hasattr(state, "begin_increment"):
                state.begin_increment(kinc, time_begin=time_end, dt=dt_sub)

            dirichlet = dirichlet_fn(frac_try)
            con, is_con = _dirichlet_masks(
                dirichlet, ndof, f"solve_steps_bounded: step {step}")
            Ui = U.copy()
            for g, v in dirichlet.items():
                Ui[g] = v
            bad = np.nonzero(is_con & ((Ui < lo - 1e-12) | (Ui > hi + 1e-12)))[0]
            if bad.size:
                raise ValueError(
                    f"solve_steps_bounded: step {step}: prescribed value at "
                    f"dof {int(bad[0])} lies outside the box "
                    f"[{lo[bad[0]]:.6g}, {hi[bad[0]]:.6g}] -- an infeasible "
                    "Dirichlet value cannot be satisfied by a bounded solve")

            step_asm = _make_step_assembler(
                asm, external_fn, step, frac_try, U_attempt_start)
            res = _bounded_newton(
                step_asm, Ui, U_attempt_start, con, is_con, lo, hi,
                w_res=w_res, at_tol=at_tol, max_newton=max_newton, tol=tol,
                linear_solve=linear_solve, line_search=line_search,
                context=f"solve_steps_bounded: step {step}", verbose=verbose)

            if res["converged"] and commit_fn is not None:
                # A compiled element can request PNEWDT < 1 when evaluated
                # again during commit. Keep that evaluation in the attempt
                # transaction, before advancing U or recording acceptance.
                try:
                    commit_fn(nodes, elems, dof_per_node, res["U"], U_attempt_start)
                except IncrementTooLarge as exc:
                    res["converged"] = False
                    res["reason"] = f"element requested a smaller increment during commit: {exc}"

            rec = BoundedAttempt(
                step=step, frac_start=frac_done, frac_end=frac_try, dt=dt_sub,
                cutback=cutback, iterations=res["iterations"],
                kkt_residual=float(res["residual"]),
                converged=bool(res["converged"]),
                n_active_lower=int(res["at_lo"].sum()),
                n_active_upper=int(res["at_hi"].sum()),
                cycled=bool(res["cycled"]), reason=res.get("reason", ""))
            report.attempts.append(rec)
            if attempt_callback is not None:
                attempt_callback(rec)

            if not res["converged"]:
                # ROLLBACK: U, SVARS and the element clock return to the
                # accepted state.
                U = U_attempt_start
                _restore_element_state(state, snap)
                cutback += 1
                report.cutbacks += 1
                if cutback > max_cutbacks:
                    raise RuntimeError(
                        f"solve_steps_bounded: step {step} failed after "
                        f"{max_cutbacks} cutbacks (KKT residual "
                        f"{rec.kkt_residual:.2e}; {rec.reason})")
                span *= 0.5
                continue

            U = res["U"]
            frac_done = frac_try
            time_end += dt_sub
            kinc += 1

        report.steps_accepted += 1
        if step_callback is not None:
            step_callback(step, frac_target, U.copy(), report)
        if early_stop is not None and early_stop(U, frac_target, state):
            report.stopped_early = True
            break
        if verbose:
            print(f"  step {step:04d}/{n_steps}: committed "
                  f"(cutbacks={cutback})", flush=True)

    return U, report

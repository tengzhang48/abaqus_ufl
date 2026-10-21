# Lessons on constitutive integration

These principles apply when authoring a material law. They do not advertise
additional public solver or constitutive-model support.

## Choose a smooth local unknown

The obvious plastic increment can produce a singular or very stiff local
consistency equation. Reformulate with a standard, well-scaled unknown before
adding arbitrary damping. An implicit local update can avoid explicit
overshoot at resolved stress concentrations. Use bounded iteration, a
consistent Jacobian and finite real/imaginary residual checks. A branch on the
real state holds the branch fixed for an infinitesimal complex step; it does
not prove differentiability at the switching surface itself.

Check the returned state, local residual and tangent at the exact failing
nonlinear state. Changing a global linear solver cannot repair an invalid
constitutive update. Conversely, test matrix functions there before blaming
an eigenvalue degeneracy. See
[complex-step branching](complex-step-and-matrix-functions.md#cs-safe-branching-rules).

## Exploit exact regime solutions without changing the law

An exact integral for a hardening, softening or relaxation regime can provide
both a stable update and an independent material-point oracle. Preserve the
regime boundary, units and time level; an exact within-regime formula does not
automatically integrate a step crossing regimes exactly. Check accumulated
state against the analytic solution and refine steps at transitions.

When a criterion depends on stress state, evaluate it on the returned stress,
unless the defining theory explicitly uses the trial state. For multiplicative
plasticity, use one consistent updated plastic state in elastic reconstruction
and stress pull-back. Check F=Fe Fp and the proper volume factors on a
noncommuting example. A diagonal history cannot distinguish factor order.

## Objectivity and material isotropy are different tests

A superposed spatial rotation QF leaves C=F.T F unchanged. Passing
P(QF)=QP(F) therefore does not exercise a principal-basis reconstruction.
For a law claiming isotropy, also test QFQ.T with distinct principal values
and require the corresponding transformed stress. Transform tensor history
according to its configuration and leave scalar history unchanged. Do not
impose isotropy on a material that deliberately has directional structure.

Use eigenspace-invariant reconstruction, with `inv(V)` for unnormalized
eigenvectors, and test repeated as well as separated spectra. The
[matrix-function lessons](complex-step-and-matrix-functions.md#cs-safe-matrix-functions)
own the complex-step and degeneracy details.

## Lagging history changes temporal accuracy

Using start-of-increment history in a damage equation can improve robustness
but generally gives a first-order temporal approximation. An internal
benchmark found substantial finite-step damage bias despite converged Newton
steps. Check two or more time-step sizes and use Richardson extrapolation only
after verifying the asymptotic order. Compare different schemes at matched
small-step limits, not at one arbitrary steps-per-cycle setting.

## Protect accumulator arithmetic and extreme driving forces

When a theory guarantees a nonnegative work increment, group it as
`alpha_new = alpha_old + (W_new - W_old)`. Left-associated
`alpha_old + W_new - W_old` can round below the old accumulator even when
W_new equals W_old. Check monotonicity at representable neighboring values
as well as with a physically meaningful tolerance; do not clip away a real
negative increment caused by an incorrect law.

A tiny fatigue degradation floor can make the damage driving history enormous.
Then cancellation in `(1-d)*H` can dominate the result near d=1. A bounded or
rescaled source, or an intact-fraction variable with d+v=1, may improve
conditioning, but needs its own derivation and equivalent-residual check.
Lumped reaction terms can change bounds and accuracy. Sweep floor, mesh and
time step and inspect both local bounds and integral observables.

Clipping a converged unconstrained field is a safeguard, not a constrained
equilibrium certificate. A genuinely bounded/irreversible formulation needs
an appropriate projected/KKT residual and a transaction that restores old
state after rejection. This is a formulation lesson; the public FE hosts do
not supply bounded or staggered solves.

# Lessons on verification and evidence

The executable order lives in the skill's
[validation ladder](../../ai_skills/abaqus-ufl/references/validation-ladder.md).
The reasons tangent consistency needs independent physics evidence live in
[code generation and testing](codegen-and-testing.md#the-finite-difference-tangent-check-is-the-primary-automated-check).

## Separate model, setup and reported quantity

A discrepancy can mean that the equations lack a valid solution, that the
setup does not implement those equations, or that a reported observable does
not mean what its label says. Distinguish these questions before choosing a
test. A material point or analytic limit can isolate a constitutive question;
a deck/DOF/coverage contract can isolate setup; a matched-condition reference
comparison tests the result. A large FE run confounds all three.

Record solver status separately from validation disposition. A completed run
can have invalid setup; an intentionally stopped negative control can be
useful evidence without being an accepted simulation.

## State the reference assumed by a gate

A swelling law can be stress-free at c=c0 rather than c=0. Checking a hidden
zero-concentration reference can reject the correct law and encourage changing
the model to satisfy the wrong test. The stress-free gate must use the supplied
scalar/history values at F=I or an explicitly declared reference state, and
name that assumption when it fails. Keep the tangent probe separate when it
is deliberately hot or loaded.

The [verification API](../API_USAGE.md#verification-ladder) and
[swelling regression](../../tests/test_verify_contract.py) define the contract.
Changing the reference to hide a genuine stress error is not a model fix.

## Prove a gate ran and can fail

Silence can mean either success or a checker that never executed. Use explicit
results, finite checks and exit status; label unavailable gates as skips.
Reject a missing binary, missing output or iteration cap instead of presenting
zero work as equilibrium. A broken control should change the specific quantity
the gate claims to observe, and the same acceptance rule must reject it.
Perform mutations in temporary copies, preserving the healthy baseline.

Regressions prevent known defects from recurring; they do not discover every
new defect. An independent reference solution or analytic identity contributes
different evidence from another reading of the same implementation.

## A rate does not establish the limit

An error decreasing at the expected order can still converge to the wrong
limit if the reference shares a sign error, if a critical operator vanishes in
the manufactured state, or if the measured quantity ignores that operator.
Choose nonuniform manufactured fields that activate the disputed terms, and
pair rate tests with an independently derived value, conservation/direction
identity or exact limit. Refine space and time separately. A failure at an
extreme state can also be finite-difference noise: sweep perturbation size
against an analytic or independent derivative before editing the law.

## Compare the same observable at matched conditions

Match time, load, mesh, material state, units, sign and spatial support.
Classify reference-frame dead loads separately from follower loads on the
current surface; their agreement at small deformation does not imply agreement
at large rotation. A maximum of element averages is different from a
quadrature-point maximum.
Do not extrapolate beyond the reference's range and label it a discrepancy.
Recompute achieved stresses or balances from retained fields rather than
trusting the driver's own residual report alone. The
[output-bridge lesson](coupled-multiphysics.md#solver-completion-is-not-output-validation)
owns detailed coverage and identity checks.

For a prescribed fully damaged slit, the nodal maximum of damage tends to its
prescribed value regardless of loading, while tip history fields may be
singular. Integral or excess-damage measures, resolved with appropriate tip
and face grading, can be informative. Inspect the attainable values of a
node-index-derived metric before treating small changes as convergence.

## Validate long histories and nonunique patterns honestly

A small systematic spatial or temporal error can accumulate over many cycles
while Newton converges on every step. Use independent higher-order or enhanced
reference formulations and longer histories, including coarse meshes that make
the bias visible. Compare schemes in their small-step limits rather than at a
single shared step size. A short-cycle match does not establish fatigue life.

A bifurcation can select different shear-band multiplicities or orientations
under different perturbations. Check force response, peak/drop location and
the allowed band-angle family; report the perturbation source and observed
pattern. An image resemblance alone is not an acceptance criterion.

## A claim belongs to its retained artifact

Today's plotting script does not prove an older figure's colorbar, and today's
deck does not prove which conditions produced an archived run. Inspect the
actual retained image, source, log and data; record their generating revision
and configuration. Recompute important statistics from the evidence rather
than copying stale counts into prose. A truncated run is incomplete evidence,
not a different completed answer. Retire superseded numbers explicitly.

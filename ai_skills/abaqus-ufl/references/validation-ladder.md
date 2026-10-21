# Validation ladder

This reference owns the operational order. The
[task index](../../../docs/lessons/README.md) routes explanations, specific
[pitfalls](pitfalls.md) and [pivots](../../../docs/pivots.md).
Choose the cheapest gate that can expose the suspected failure.

| Gate | Execute and inspect | What it establishes |
|---|---|---|
| Evidence record | Read the actual checked deck, retained data/log and stated observable | Which configuration, quantity, runtime and oracle were already checked |
| Python oracle | Run the example's `check_reference.py`; assert the claimed regime was entered | An independent limit/invariant when the oracle is separately derived |
| Material consistency | Call `model.verify()` or `problem.verify()` at relevant states, with an explicit reference when needed | CS/FD consistency of implemented material methods; WeakForm.verify delegates to material verification |
| Element assembly | Run `check_assembled.py` where shipped, FD the element tangent at fixed old state, and test a patch/operator invariant | Residual/tangent, mapping, DOF/state and weak-form behavior for the checked element |
| Regeneration and build | Regenerate deterministically and compile every emitted source | Reproducibility and compiler acceptance, including fixed-form/interface declarations |
| Compiled point/element | Run `check_compiled.py` or a direct f2py wrapper with finite nonzero output and accepted PNEWDT | Executed generated-code interface, signs, state and tangents at that point/element |
| Small global solve | Run an appropriate `abaqus_ufl.fe` case, including history, failure and refinement checks | The particular host's assembly, boundary conditions, nonlinear solve and commit contract |
| Production/output bridge | Run the authorized Abaqus bundle and audit the output identities before comparing values | Abaqus call/procedure conventions and the checked solved-to-reported quantity bridge |
| Published response | Match the paper's geometry, material, loads, time/mesh limit and observable | Reproduction only within that matched, tested envelope |

## Commands and ownership

Use the nearest public example's commands and its [README](example-map.md).
The [example pipeline](../../../HOWTO_ADD_AN_EXAMPLE.md) owns required bundles;
the [API guide](../../../docs/API_USAGE.md#verification-ladder) owns verification
options. The [host guide](../../../abaqus_ufl/fe/README.md) owns serial/MPI usage.
A point/element f2py wrapper calls the generated boundary; mesh parsing,
assembly, Newton and time stepping belong in a host rather than that wrapper.

Keep integer flags as Fortran integers and derive hidden dimensions from the
actual arrays. Inspect actual/dummy argument order, tensor/Voigt/state layout,
normal and unsupported request flags, beginning TIME and later increments.
Stateful tests use fresh or explicitly reset old state; compiled test order
must not affect the result. Read the relevant
[interface pitfalls](pitfalls.md#abaqus-interface) for the symptom at hand.

## Acceptance and scope

Use an independent observable or identity and a broken control that the same
acceptance rule rejects. Compilation, Python/Fortran parity and Newton
convergence answer different questions. None alone proves the intended PDE,
branch coverage, physical accuracy or production output meaning. An absent
check is an explicit skip; a timeout, NaN or iteration cap is a failure, not a
completed result. A known open defect may use `xfail(strict=True)` against the
correct expected behavior; remove it when the fix passes.

Record quantity, oracle, runtime, state/parameter envelope, output bridge and
revision with the result. Use bounded claims such as implementation checked,
solver robustness open, or matched response reproduced; do not turn a small
smoke test into a full-model claim. The
[evidence lessons](../../../docs/lessons/verification.md) explain these boundaries.

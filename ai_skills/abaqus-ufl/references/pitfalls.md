# abaqus_ufl pitfalls

A pitfall is a specific mistake, its observable symptom and the check that
catches it. Explanations live in the linked lesson/API home. Start with the
[task index](../../../docs/lessons/README.md) when the failure is not localized.
Historical heading anchors are retained so existing links keep working.

## Abaqus interface

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Treat RHS as the host residual without converting signs | Newton goes uphill; FD gives the negative of K | Check K=dR/dU after R=-RHS; [AMATRX convention](../../../docs/lessons/abaqus-and-fortran.md#amatrx-sign-convention) |
| Ignore PNEWDT or nonfinite output | Zeroed work appears to be equilibrium | Assert PNEWDT>=1 and finite R/K/state before assembly or commit; [runtime contract](../../../abaqus_ufl/fe/README.md#exact-batches-and-material-state) |
| Assume RHS is 1D or pre-zeroed | Extra RHS columns retain garbage on early exits | Exercise NRHS>1 and unsupported requests; [interface lessons](../../../docs/lessons/abaqus-and-fortran.md#rhs-is-a-2d-array) |
| Recover old fields from U alone or silently use DU=U | Transient terms read zero history | Pass DU=U-U_previous, or explicitly declare a from-zero call; [old fields](../../../docs/lessons/abaqus-and-fortran.md#old-field-values-are-not-passed-to-the-uel) |
| Freeze TIME or pass end-of-increment time | First-increment initialization is never exercised, or seeded state is replaced | Check KINC and beginning TIME at first/later increments; [host API](../../../abaqus_ufl/fe/README.md#api) |
| Call the plane-strain UMAT with a 3D or plane-stress signature | PNEWDT requests a cutback | Verify NDI=3, NSHR=1, NTENS=4 and plane-strain F for ndim=2; [UMAT API](../../../docs/API_USAGE.md#generate_umat) |
| Keep a stale positional PROPS list | Finite Python response but NaN or wrong compiled response | Compare deck count/order with instance properties and probe identical inputs; [property contracts](../../../docs/lessons/abaqus-and-fortran.md#constantsn-must-match-the-umat-contract) |
| Use an unsupported scalar DOF/procedure or malformed deck block | Datacheck failure or an inactive scalar | Check active DOFs, procedure, trailing commas and keyword placement; [deck lessons](../../../docs/lessons/abaqus-and-fortran.md#abaqus-input-deck-compatibility) |

## DOF and state layout

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Zero F in plane strain | detF=0 and a cutback at zero displacement | Require F33=1; [kinematics](../../../docs/lessons/finite-element.md#initialize-f-to-identity-not-zeros-plane-strain) |
| Mix node-major, component-major or corner/midside field layouts | Wrong coupling or apparently singular rows | Assert the complete element DOF map, including mixed-degree midsides; [indexing](../../../docs/lessons/finite-element.md#element-local-dof-indexing-0-based-front-end-vs-1-based-arrays) |
| Build a separate Jacobian from lower-order scalar nodes | Distorted mixed-element patch fails | Use the geometry map consistently for every field; [mapping](../../../docs/lessons/finite-element.md#one-jacobian-for-all-fields-mixed-degree-elements) |
| Flatten tensor state in row-major order | A nonsymmetric old tensor returns transposed | Round-trip distinct scalar/vector/tensor values at multiple Gauss points; [state-order regression](../../../tests/test_reference_assembly_state_order.py) |
| Hardcode a topology's Gauss/state count in an adapter | Correct kernel, wrong host state or coordinates | Pin the active layout and check every adapter separately; [adapter lesson](../../../docs/lessons/finite-element.md#adapter-layout-and-accepted-state-are-part-of-the-element-contract) |
| Mutate committed state during a trial or failed commit | History advances after rejection | Seed nonzero old state, force failure and assert exact restoration; [state ownership](../../../docs/lessons/performance-and-parallel.md#exact-reuse-needs-complete-inputs-and-state-ownership) |
| Upgrade connectivity but retain old boundary/MPC/contact sets | Unconstrained midsides or overconstraints despite clean completion | Check geometric set membership, incompatible overlaps and the solved boundary invariant; [mesh-order lesson](../../../docs/lessons/finite-element.md#mesh-order-changes-require-new-set-membership-checks) |

## Material / WeakForm authoring

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Add a solved gradient field as local state | No transport coupling despite plausible constitutive values | Identify the field owner and global balance; [target/API guide](../../../docs/API_USAGE.md#choose-the-target) |
| Guess canonical storage/flux method names | Generation calls a missing method | Construct the WeakForm and read its detected equations; [guard regression](../../../tests/test_equation_method_guard.py) |
| Read class defaults instead of constructor overrides | Overrides have no effect in the compiled path | Compare PROPS with material-instance attributes; [adapter lesson](../../../docs/lessons/finite-element.md#adapter-layout-and-accepted-state-are-part-of-the-element-contract) |
| Use dynamic Python or an unsupported tensor helper | Verification works but generation rejects syntax | Check supported operations before rewriting the law; [DSL rules](../../../docs/API_USAGE.md#tensor-dsl-rules) |
| Hide a missing primitive with a second implementation | Reference and generated paths drift | Make the declaration express the model and compare the reached regime; [DSL lesson](../../../docs/lessons/codegen-and-testing.md#stop-and-surface-a-missing-primitive) |

## Fortran generation

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Ignore fixed-form column 72 | Compiled expression differs or continuation fails | Compile the wrapped source; [fixed form](../../../docs/lessons/abaqus-and-fortran.md#column-72-truncation) |
| Reuse case-insensitive names, or unpack into `_` | Duplicate declaration or invalid identifier | Inspect base names of indexed variables and compile; [namespace lesson](../../../docs/lessons/codegen-and-testing.md#single-letter-fortran-names-are-landmines) |
| Emit a complex `.real` temporary, loop counter or bound | IF/DO compile failure or dropped perturbation | Inspect generated declarations and run the intended branch; [generator constructs](../../../docs/lessons/complex-step-and-matrix-functions.md#generator-inserted-constructs-can-break-cs-even-when-the-source-is-safe) |
| Replay Python argument order after generation grouped state arguments | Compiles but tangent/Newton behavior is wrong | Compare call actuals to emitted dummies, especially history and dt; [history-order lesson](../../../docs/lessons/coupled-multiphysics.md#match-callercallee-ordering-for-multi-scalar-history-arguments) |
| Use inclusive slices with one extra endpoint | Adjacent state slots change | Check length end-start+1 using distinct neighboring state values; [slice lesson](../../../docs/lessons/codegen-and-testing.md#fortran-array-slices-are-inclusive--off-by-one-is-silent) |
| Trust compilation without a subroutine call | Runtime NaNs or an ABI mismatch remain invisible | Execute the generated boundary with nonzero outputs; [validation ladder](validation-ladder.md) |
| Repair generated output without its declaration/generator | Regeneration silently loses the fix | Require reproducible regeneration; [workflow](../../../docs/lessons/agent-workflow.md#make-progress-observable-and-claims-reproducible) |

## F-bar routing

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Force pure-mechanics F-bar onto a coupled formulation | Inconsistent swelling/storage coupling or rejected route | Select standard mixed/local-pressure theory explicitly; [coupled locking](../../../docs/lessons/coupled-multiphysics.md#f-bar-can-destabilize-coupled-problems-at-si-scale) |
| Use 3D volumetric scaling or contractions in 2D | F-bar tangent/value error | Check dimension-dependent exponent and in-plane contractions; [F-bar lesson](../../../docs/lessons/finite-element.md#f-bar-the-volumetric-exponent-is-dimension-dependent) |
| Leave auto-selected F-bar on in an intended locking control | The baseline is already stabilized | Inspect the generation summary and set formulation explicitly; [UEL API](../../../docs/API_USAGE.md#generate_uel) |

## Complex-step and matrix functions

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Strip the imaginary part from differentiated algebra | CS tangent misses a term while values look right | Compare derivatives at a reached nontrivial branch; [real-cast leak](../../../docs/lessons/complex-step-and-matrix-functions.md#defensive-real-casts-are-a-cs-leak) |
| Blame every real-part branch or repeated spectrum for a NaN | Solver tuning changes no root cause | Probe that exact state in Python and Fortran first; [branch rules](../../../docs/lessons/complex-step-and-matrix-functions.md#cs-safe-branching-rules) |
| Use gauge-dependent eigenvectors to reconstruct a matrix function | Rotated/repeated-spectrum derivatives fail | Prefer invariant/iterative matrix functions and test spectra explicitly; [matrix-function lessons](../../../docs/lessons/complex-step-and-matrix-functions.md#cs-safe-matrix-functions) |
| Reuse complex current fields between perturbations | Tangent depends on perturbation order | Reset current fields and keep old state real/frozen; [reset lesson](../../../docs/lessons/complex-step-and-matrix-functions.md#reset-all-complex-state-before-each-perturbation) |
| Interpret finite-difference noise as a model defect | An extreme-state mismatch changes strongly with FD step | Sweep FD step against an independent derivative; [rate/limit lesson](../../../docs/lessons/verification.md#a-rate-does-not-establish-the-limit) |
| Compare all nine columns with an intentionally in-plane tangent | Expected zero unused columns look like a bug | Compare active assembled columns and test the plane-strain contract; [active-work lesson](../../../docs/lessons/performance-and-parallel.md#perturb-only-active-columns) |

## UEL weak-form signs

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Copy a flux sign from a field name or sibling equation | Anti-diffusive localization despite tangent agreement | Derive storage*test-flux.grad(test), then assemble a nonuniform operator and a flipped control; [sign lesson](../../../docs/lessons/coupled-multiphysics.md#derive-the-weak-form-tangent-checks-cannot-prove-the-pde-sign) |
| Treat the scalar operator-sign screen as a proof for a mixed/anisotropic system | A legitimate mixed pair warns, or a cross-gradient error escapes | Inspect the full coupled block/dispersion spectrum; [sign conventions](../../../docs/lessons/coupled-multiphysics.md#the-flux-sign-depends-on-the-storage-convention-not-the-field-name) |
| Recover a stress-gradient flux from incomplete bilinear curvature | Error plateaus under refinement | Compare an auxiliary-field or higher-order reference; [gradient consistency](../../../docs/lessons/finite-element.md#gradient-of-stress-consistency) |

## Multiplicative-plasticity and finite-strain kinematics

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Reorder noncommuting plastic factors or pull back with stale state | Diagonal histories pass, rotated histories fail | Assert F=Fe Fp with the updated state and correct volume factors; [integration lesson](../../../docs/lessons/constitutive-integration.md#exploit-exact-regime-solutions-without-changing-the-law) |
| Substitute V.T for inv(V) with unnormalized eigenvectors | Off-axis spectral reconstruction is wrong | Test conjugated principal bases and repeated spectra; [isotropy lesson](../../../docs/lessons/constitutive-integration.md#objectivity-and-material-isotropy-are-different-tests) |
| Use spatial objectivity as the only isotropy check | Principal-basis errors stay hidden | Add the distinct QFQ.T test only when isotropy is claimed; [isotropy lesson](../../../docs/lessons/constitutive-integration.md#objectivity-and-material-isotropy-are-different-tests) |

## Finite-strain tensor-return shape trap

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Tensorize a scalar matrix-helper argument componentwise | Source compiles but gives NaNs, such as det33z(F(ii,jj)) | Inspect helper argument rank and execute the direct return; [guard regressions](../../../tests/test_translator_guards.py), [runtime regressions](../../../tests/test_runtime_guards.py) |

## Local return maps and nonlinear hardening

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Dampen a singular local unknown instead of reformulating it | Iteration stalls or drifts under refinement | Check the local consistency residual in a smooth, scaled unknown; [local integration](../../../docs/lessons/constitutive-integration.md#choose-a-smooth-local-unknown) |
| Evaluate a stress-state criterion on the wrong time level | Return-map values and the reported criterion disagree | Recompute from the returned stress and state; [regime solution lesson](../../../docs/lessons/constitutive-integration.md#exploit-exact-regime-solutions-without-changing-the-law) |
| Left-associate a monotone accumulator increment | An unchanged work value reduces the accumulator by an ulp | Check alpha_old+(W_new-W_old) near representable neighbors; [precision lesson](../../../docs/lessons/constitutive-integration.md#protect-accumulator-arithmetic-and-extreme-driving-forces) |

## Transcribing from the paper

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Use a constraint as a definition, a sibling's parameter, or one symbol for two roles | Internally consistent but wrong constitutive response | Cite definition equations and specific table units; [theory-first lessons](../../../docs/lessons/codegen-and-testing.md#validation-envelopes-must-cover-the-branch-being-claimed) |
| Call a tuned localization run a reproduction | Curves shift with loading model or missing setup | Match parameters, geometry, dead/follower loads, operating point and observable; [paper lessons](../../../docs/lessons/coupled-multiphysics.md#reproducing-multiphysics-papers) |

## Coupled dummy materials in Abaqus

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Omit density/specific heat or make an overlay mechanically significant | Coupled procedure rejects the deck or modifies the solution | Check dummy material completeness and overlay nonintrusiveness; [input scaffold API](../../../docs/API_USAGE.md#abaqus-input-scaffolding) |
| Assume the dummy/visualization field is the solved field | Plausible picture with the wrong scalar/history | Compare bridge values pointwise by stable identities; [output bridge](../../../docs/lessons/coupled-multiphysics.md#solver-completion-is-not-output-validation) |

## Non-Abaqus regression runtimes

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Promote a license-free host result to Abaqus evidence | Production procedure/output behavior is untested | State the executed interface and separate the production gate; [validation ladder](validation-ladder.md) |
| Pass a Fortran-order nonsymmetric tangent to a row-major consumer | Distributed Newton uses K.T | Check C-contiguity and solve a nonsymmetric control; [parallel lesson](../../../docs/lessons/performance-and-parallel.md#keep-replicated-batching-distinct-from-distributed-algebra) |
| Skip a collective after a rank-local cache hit or failure | MPI hangs, often with empty partitions | Require unanimous hits/error votes and run rank-local failure probes; [MPI tests](../../../tests/mpi_fe_contract.py) |
| Treat a KSP cap, exhausted line search or Newton cap as success | An unconverged state is committed | Force each failure path and assert no commit/unchanged state; [host API](../../../abaqus_ufl/fe/README.md#petsc-and-mpi) |
| Infer scaling or ABI safety from one tiny rank count | Larger runs slow down or fail | Match launcher/library ABI, measure rank/thread budgets and test intended ranks; [environment lesson](../../../docs/lessons/performance-and-parallel.md#diagnose-the-mpi-environment-before-the-algorithm) |

## Driving a generated element in Python (reference_assembly)

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Pass node-major coordinates where reference_assembly expects (ndim,nnode), or omit element config | Index error or wrong geometry | Match the API's coordinate orientation and active layout; [adapter lesson](../../../docs/lessons/finite-element.md#adapter-layout-and-accepted-state-are-part-of-the-element-contract) |
| Confuse raw RHS with the public host's R=-RHS | Hand-written Newton uses the wrong sign | FD one assembled column at fixed old state, perturbing U and DU together; [sign lesson](../../../docs/lessons/abaqus-and-fortran.md#amatrx-sign-convention) |

## Hard paper models

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Jump to the final figure or invent missing parameters | A large solve cannot localize the discrepancy | Build the smallest independent limit first and record unknown inputs; [framework rungs](../../../docs/lessons/coupled-multiphysics.md#hard-multiphysics-needs-framework-rungs-before-parameter-tuning) |
| Condense a variable with a gradient/global balance for convenience | Missing coupling or no well-defined transport | Derive the field's ownership and variational role; [locking/formulation lessons](../../../docs/lessons/coupled-multiphysics.md#locking--element-choice) |

## Debugging order

Route a failure through the [validation ladder](validation-ladder.md). Reproduce
the exact state and test a hypothesis before changing equations. Read existing
capability/evidence records before declaring a framework limitation. No symptom
alone uniquely identifies a physics, setup, generator or solver defect.

## Verification: the state and the regime, not just the method

| Mistake | Symptom | Check and explanation home |
|---|---|---|
| Assume c=0 is every model's stress-free reference | A valid swelling model fails only the reference gate | Supply c=c0 or a declared reference and read the failure assumption; [reference lesson](../../../docs/lessons/verification.md#state-the-reference-assumed-by-a-gate) |
| Check a benign state while claiming a nonlinear branch | Green suite, untested evolution | Assert branch activation and a branch-specific invariant; [validation envelopes](../../../docs/lessons/codegen-and-testing.md#validation-envelopes-must-cover-the-branch-being-claimed) |
| Trust a gate never observed to fail | Silence or a success string is accepted as evidence | Run a defect targeting that quantity and require rejection; [gate lesson](../../../docs/lessons/verification.md#prove-a-gate-ran-and-can-fail) |
| Use more quadrature to repair incompatible swelling strain | Damage bias survives refinement and quadrature changes | Compare compatible mixed/enhanced interpolation independently; [gradient lesson](../../../docs/lessons/finite-element.md#gradient-of-stress-consistency) |
| Compare lagged and implicit history at one coarse step size | Apparently different physics or fatigue life | Compare small-step limits and long histories; [time-level lesson](../../../docs/lessons/constitutive-integration.md#lagging-history-changes-temporal-accuracy) |
| Score a prescribed maximum or singular tip value | Metric saturates or changes without meaningful convergence | Use resolved integral/excess measures and inspect metric support; [observable lesson](../../../docs/lessons/verification.md#compare-the-same-observable-at-matched-conditions) |
| Clip the field and report unconstrained residual convergence | Bounds pass while constrained equilibrium fails | Recompute the appropriate projected/KKT certificate; [clipping pivot](../../../docs/pivots.md#p6-clipping-is-not-constrained-equilibrium) |

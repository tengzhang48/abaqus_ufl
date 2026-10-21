# Lessons on agent workflow

The [skill](../../ai_skills/abaqus-ufl/SKILL.md) routes tasks; this page owns
the workflow lessons and their rationale.

## Default to one agent and an independent oracle

An internal benchmark found that an implement/review/fix team consumed much
more effort without improving element accuracy over a single agent with its
own independent reference. Extra agents are not a correctness criterion.
Default to a single agent; seek independence in an analytic solution,
separately derived reference, higher-order discretization or experimental
observable. Another evaluation of the same formulas is consistency evidence.

Review should attack the validation envelope: which branches, rotated states,
gradients, near-singular limits and long histories have the checks not reached?
Turn a suspected gap into a counterexample and check, rather than adding an
unverified reviewer claim. If delegation is explicitly requested, assign
nonoverlapping ownership and reconcile assumptions and evidence at integration.
The [pivot](../pivots.md#p3-one-agent-with-an-independent-oracle) records the
change from the former team-as-default advice.

## Give the task a bounded contract

State the desired behavior, writable/read-only scope, theory/API constraints,
acceptance quantities and tolerances, and the required handoff. Read the
canonical repository skill, nearby API and example evidence before inventing
a second workflow. Preserve the user's existing authorization and scope;
repository access alone does not authorize disclosure, publication or a
different scientific model. Record missing theory or evidence as an open
question instead of silently choosing parameters.

## Keep one canonical implementation and index

Divergent package, lesson and installed-skill copies can hide an already fixed
bug. Record which checkout/revision is executing, keep the repository skill
with that source, and route through one index. Put an item's explanation in
one lesson, its symptom/check in one pitfall, and its change of approach in a
pivot. Link those homes instead of copying their prose. Preserve dated run
artifacts as evidence, with explicit status, rather than treating them as the
current implementation. The [pivot](../pivots.md#p4-one-canonical-copy) owns this
organization change.

## Make progress observable and claims reproducible

Use small independent gates before large runs. Report flushed progress for
long jobs and distinguish build time, a running solve, a failed attempt and a
completed result. A timeout does not identify the numerical cause. Record
revision, configuration, command, measured result, unrun checks and scope of
each claim. Re-run a disputed pass count instead of citing an earlier message.

Fix shared code at its authoritative generator or declaration, then enumerate
every consumer and early-exit obligation. The same defect can appear in
several wrappers with different text, so searching for one implementation's
pattern is insufficient. Check the retained artifact separately from its
current producer. The [verification lessons](verification.md#a-claim-belongs-to-its-retained-artifact)
explain why this matters for archived evidence.

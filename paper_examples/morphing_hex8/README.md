# Three-dimensional diffusion-driven morphing (Sec. 3.4, Fig. 6)

Grooved-sheet morphing with the Hex8 element-local-pressure gel UEL
(global u and mu; one condensed pressure per element).

- `pipeline_hex8.py` — the CURRENT declaration-to-Fortran pipeline
  (verify + generate with the neutral `pressuregel` prefix).
- `archived_submitted/chester_anand_local_pressure_hex8_pasta.for` —
  the EXACT generated UEL submitted with the completed Abaqus/Standard
  2022 production run, preserved unmodified (historical internal
  prefix; see the top-level naming note).
- `abaqus/` — the exact production deck, the C3D8 companion-mesh include,
  material properties, and portability-corrected run instructions.
- `reference_data/` — the reduced ODB record of the completed run
  (frame times, final fields, and maximum displacement).
- `figure/` — the Fig. 6 script and the six Abaqus/CAE exports.

At revision `e0985d9`, `tests/test_uel_local_pressure_contract.py` contains 15
current-generator tests covering the Quad4 and Hex8 condensed Jacobians, the
pressure-response discriminator, a frozen-Schur failure control, local-solve
idempotence, and UEL call contracts. Separately, a July 28 audit compiled the
exact archived Figure 6 source and obtained a complete black-box
condensed-Jacobian relative error of `6.494e-10`. The later dispatch, cutback,
and state-commit guards are properties of the current generators and are not
attributed to the archived production source.

The C3D8 mesh is not an output-only, negligible-stiffness overlay. It shares
nodes with the UEL mesh and contributes to mechanical equilibrium with initial
shear modulus 200, compared with 800 for the UEL. The reported Figure 6
deformation is therefore from the combined UEL/native-element model. A
companion-mesh stiffness-sensitivity audit remains pending.

# Lessons: Abaqus Interface & Fortran Codegen

Hard-won traps and best practices for writing (or generating) Abaqus user
subroutines (UEL/UMAT), compiling fixed-form Fortran, and authoring input
decks that the Abaqus parser will actually accept. Each item states the trap,
why it bites, and the fix.

---

## Abaqus UEL/UMAT interface

### AMATRX sign convention

Abaqus defines `AMATRX` as the Jacobian of the system, i.e.
`AMATRX = -dRHS/dU`. If your residual mixes signs (for example, a stress term
enters with `-` but a flux term enters with `+`), each tangent block must take
the *opposite* sign from its residual contribution.

```
RHS -= c_dot * N * wdetJ       ->  AMATRX += dc_dot/du * N * dN/dX * wdetJ
RHS += j_R . dN/dX * wdetJ     ->  AMATRX -= dj_R/du * dN/dX * dN/dX * wdetJ
                                              ^^^^^^ negative!
```

**Fix:** Build a sign table for every residual term *before* writing assembly
code. Track `rhs_sign` and `amatrx_sign = -rhs_sign` per term.

### RHS is a 2D array

Abaqus declares `RHS(MLVARX, *)`, not `RHS(NDOFEL)`. Writing `RHS(i)` compiles
but accesses the wrong memory.

**Fix:** Always index as `RHS(row, 1)`.

### Old field values are not passed to the UEL

Abaqus does not provide previous-step field values. The UEL receives
`U(NDOFEL)` (total DOFs at the end of the increment) and `DU(MLVARX, *)`
(incremental DOFs). Reconstruct old values yourself:

```fortran
      u_old(i, a) = U(idx) - DU(idx, 1)
```

**Fix:** Parse current and old DOFs in the same loop.

### DTIME = 0 on the initial stiffness evaluation

Abaqus calls the UEL with `DTIME = 0` to evaluate the initial stiffness before
taking any time increment. Any `1/DTIME` in your code (e.g. a backward-Euler
rate `c_dot = (c - c_old) / dt`) then raises a floating-point exception and
crashes Abaqus.

**Fix:** Floor the time step before dividing:

```fortran
      dt_safe = DTIME
      IF (dt_safe .LT. 1.0d-14) dt_safe = 1.0d-14
```

### Do not assume RHS/AMATRX are pre-zeroed

Abaqus warns that `RHS` and `AMATRX` may contain dirty memory depending on the
nonlinear solver settings (Full Newton vs. Quasi-Newton). Relying on zero
initialization causes sporadic wrong results that are nearly impossible to
debug.

**Fix:** Explicitly zero both arrays at the top of the UEL:

```fortran
      DO i = 1, NDOFEL
        RHS(i, 1) = 0.0d0
        DO j = 1, NDOFEL
          AMATRX(i, j) = 0.0d0
        END DO
      END DO
```

### Guard against degenerate (inverted) elements

During Newton iterations in large-deformation problems, elements can
temporarily invert (`det(J_xi) <= 0`). Proceeding with the Gauss loop makes
`1/det(J_xi)` produce NaN, which propagates silently through the element
stiffness and eventually crashes the solver with an unhelpful message.

**Fix:** Check the Jacobian mapping immediately and ask Abaqus to cut the step
instead of continuing:

```fortran
      CALL map_grad_2d(dshxi, coords, nNode, dsh, detJ, Jinv, stat)
      IF (stat .EQ. 0) THEN
        PNEWDT = 0.25d0    ! tell Abaqus to cut the time step
        RETURN
      END IF
```

### Use PNEWDT for time-step control

`PNEWDT` is the UEL's lever on the increment size. Use it aggressively:

- `PNEWDT = 0.25` for degenerate elements
- `PNEWDT = 0.5` for large strain increments
- `PNEWDT = 1.5` when converging easily (lets Abaqus grow the step)

### DOF numbering for mixed elements

On a `*User Element` card, DOF labels must be listed explicitly. For mixed
elements where not all nodes carry all DOFs (e.g. pressure only on corner
nodes), the labels must match exactly what the UEL expects.

```
*User Element, Type=U2, Nodes=8, Coordinates=2,
 Properties=9, Variables=1
1, 2, 11, 12
5, 1, 2, 12
```

Here DOFs 1,2 = displacement, 11 = pressure, 12 = a second scalar field (e.g.
chemical potential). The first data line is the DOF list applying from nodal
position 1 (corners: all four DOFs); each subsequent line starts with a nodal
position and gives the DOF list that applies from that position on
(`5, 1, 2, 12`: midside nodes 5-8 carry 1, 2, 12 only). Abaqus does NOT infer
varying per-node DOF sets from the element connectivity; without the second
data line every node gets the first list.

**Fix:** Document the DOF layout with a table mapping DOF index -> node -> field.

---

## Fortran fixed-form compilation

### Column-72 truncation

Abaqus compiles `.for` files as fixed-form Fortran: anything past **column 72
is silently ignored**. A long expression like

```fortran
      cR0 = ((DCMPLX(1.0d+00, 0.0d0) - DCMPLX(props(9), 0.0d0)) / DCMPLX(props(6), 0.0d0))
```

is truncated to

```fortran
      cR0 = ((DCMPLX(1.0d+00, 0.0d0) - DCMPLX(props(9), 0.0d0)) / DCM
```

producing phantom syntax errors or, worse, silently wrong compiled math.

**Fix:** Wrap every line at column 72 with a `&` continuation marker in column 6.
For generated code, run the output through a line-wrapping utility.

**Codegen gotcha:** A line wrapper that detects comment lines by checking
whether the *stripped* line starts with `c` will also match Fortran variable
names beginning with `c` (e.g. `cR0`, `cdot_out`) and wrongly skip them. Detect
comments by checking **column 1 of the original line**, before stripping
whitespace.

### No Fortran 2003+ features (BLOCK, etc.)

`BLOCK ... END BLOCK` for local declarations is a Fortran 2008 feature. Abaqus
traditionally compiles with F77/F90 flags unless the user edits
`abaqus_v6.env`, so `BLOCK` fails to compile with no obvious error.

**Fix:** Declare all variables at the top of the subroutine. No `BLOCK`, no
variable-length arrays, no F2003+ features.

### Use standard intrinsics

`CDLOG`, `CDEXP`, `CDSQRT` are GNU extensions, not standard Fortran. The Intel
compiler used by Abaqus may not support them.

**Fix:** Use the standard generic intrinsics `LOG`, `EXP`, `SQRT`, which accept
both REAL and COMPLEX arguments in standard Fortran.

### Always IMPLICIT NONE

Without `IMPLICIT NONE`, undeclared variables starting with `i`-`n` default to
INTEGER and `a`-`h`, `o`-`z` default to `REAL*4`. A typo like `detj` instead of
`detJ` then silently creates a single-precision variable.

**Fix:** Put `IMPLICIT NONE` in every subroutine.

### Fortran is case-insensitive

Uppercase and lowercase identifiers are the same symbol. A generated
declaration like `DOUBLE PRECISION :: J` collides with an integer loop
variable `j`, and compilers such as `gfortran` reject the code.

**Fix:** Never use single uppercase letters (`I`, `J`, `K`, `L`, `M`, `N`) as
variable names — they collide with lowercase loop indices.

---

## Abaqus input-deck compatibility

Input-deck parsing is stricter and more version-sensitive than the docs
suggest; the behaviors below were observed with the Abaqus 2022 parser.

### constants=N must match the UMAT contract

`*User Material, constants=N` easily drifts from the generated UMAT's `NPROPS`
and property order. Supplying **too many** values shifts later properties
relative to what the UMAT expects; supplying **too few** triggers
`INVALID DATA ASSOCIATED WITH THIS USER DEFINED MATERIAL DEFINITION`, or the
UMAT reads undefined entries.

**Fix:** Treat the generated `.for` header and the property order in your
generator (e.g. the `props` dict) as the single source of truth. Use 8 values
per line (standard Abaqus layout); the last line gets whatever remains; no
trailing comma on the final line.

```
*User Material, constants=38
** props 1- 8
 10.0,    20.0,    1.2,     0.0,     0.9,     0.05,    3.5,     20.0,
** props 9-16
 1.0e-8,  0.0,     0.0,     1.0,     0.25,    0.0,     0.0,     1.0,
** props 17-24
 1.0,     1.0,     1.0,     0.0,     0.0,     0.0,     10.0,    0.5,
** props 25-32
 0.3,     1.0,     0.0,     0.0,     0.9,     1.0,     0.0,     0.8,
** props 33-38
 1.0,     0.0,     0.0,     0.0,     0.0,     0.0
```

### *Depvar needs a trailing comma

A bare integer on the line after `*Depvar` parses, but the expanded deck shows
no value and downstream `*Initial Conditions` fails with `INSUFFICIENT DATA
CARDS`.

**Fix:** Add a trailing comma (matching working Abaqus examples, where
single-value data lines consistently use one):

```
*Depvar
 9,
```

### *Initial Conditions respects the 8-value line limit

Putting all SDV values plus the elset name on one line overflows the limit:

```
SOLID, 0.8, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0   <- 10 values
```

Abaqus reads the first 8 data values, discards the rest, and reports
`INSUFFICIENT DATA CARDS` because only 8 of 9 SDVs were provided.

**Fix:** Split across lines (elset + up to 8 values on the first line, the
remainder on the next). No trailing comma on the last line:

```
*Initial Conditions, type=Solution
 SOLID, 0.8, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0
 0.0
```

### Comments use ** not bare *

`* This is a comment` is parsed as an incomplete keyword, producing cascading
fatal errors (`AMBIGUOUS KEYWORD`, then `NODE LABEL IS NOT AN INTEGER`, ...).

**Fix:** Every comment line starts with `**`. No exceptions.

### *User Element Variables must be >= 1

Emitting `Variables=0` for a UEL with no state variables is rejected:
`INTEGER VALUE FOR PARAMETER VARIABLE MUST BE BETWEEN 1 AND 2147483647`.

**Fix:** Enforce a minimum of 1 in the generator:

```python
if config.variables > 0:
    fields.append("Variables={}".format(config.variables))
else:
    fields.append("Variables=1")  # Abaqus minimum is 1
```

### No blank lines between keyword blocks

Inserting a blank line between consecutive `*Element`, `*Elset`, or `*Node`
blocks throws `BLANK LINE IN ELEMENT DATA`.

**Fix:** Write keyword blocks continuously, using `**` separator lines rather
than bare blank lines.

### *Conductivity belongs inside *Material

`*Conductivity` emitted as a standalone keyword is rejected; Abaqus requires it
inside a `*Material` definition.

### *Elset must be explicit for *Initial Conditions

Defining an element set only via `elset=SOLID` on `*Element` and then
referencing `SOLID` in `*Initial Conditions` can fail to resolve the implicit
set for initial-condition data cards.

**Fix:** Define the set explicitly:

```
*Element, type=C3D8, elset=SOLID
 1, 1, 2, 3, 4, 5, 6, 7, 8
**
*Elset, elset=SOLID
 1
```

### Quick checklist

| Item | Correct | Wrong |
|------|---------|-------|
| Comments | `** comment` | `* comment` |
| `constants` | Matches actual value count | Off-by-one |
| `*Depvar` | `9,` (with comma) | `9` (bare) |
| Material data | 8 per line, comma on continuations | Trailing comma on last line |
| `*Initial Conditions` | Split to <= 8 values/line | All values on one line |
| `*User Element` vars | `Variables=1` minimum | `Variables=0` |
| Blank lines | Never between keyword data | Blank line between blocks |
| `*Conductivity` | Inside `*Material` | Standalone |
| Element sets | Explicit `*Elset` | Implicit via `elset=` only |

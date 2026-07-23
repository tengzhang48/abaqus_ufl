# How to add a validated example

An example is only released once it is **verified** and **license-clear**. The
template makes the verified part mechanical.

## Steps

1. **Copy the template.**
   ```bash
   cp -r examples/_template examples/<your_name>
   ```

2. **Author the model** in `examples/<your_name>/build.py`: rename the class,
   fill in `props` and `stress_PK1` (UMAT) or the weak form (UEL). Then:
   ```bash
   cd examples/<your_name> && python build.py
   ```
   `build.py` runs `model.verify()` first — a complex-step-vs-finite-difference
   tangent check. **If it fails, stop; the generated code would be wrong.** On
   success it writes `<your_name>_umat.for`.

3. **Write the deck** `abaqus/job.inp`: keep it to one element for a clean
   material-point comparison. The `*User Material, constants=` values must match
   `props` in the same order.

4. **Generate the reference** (independent Python integration of the same path):
   ```bash
   python abaqus/generate_reference.py    # writes abaqus/reference.json
   ```

5. **Run in Abaqus and compare** (needs Abaqus/Standard on PATH):
   ```bash
   cd abaqus && bash run.sh
   python ../../../tools/compare_results.py job_extracted.json reference.json
   ```

6. **Write `README.md`** from the template header. Set the **Verification** line
   to the level you actually reached — do not label a smoke run "validated".

7. **Add one row** to `examples/README.md` with the honest verification level.

## Before it can ship (release gates)

- [ ] Verification level in the README matches reality.
- [ ] No third-party source without its license (put any vendored code under a
      clearly-labeled `third_party/` with the exact notice, or leave it out).
- [ ] No PDFs, ODB/large outputs, `__pycache__`, dev notes (`*REVIEW*`,
      `*_PLAN*`, `*_LOG*`, `*NOTE*`), or absolute HPC paths in scripts.
- [ ] Added to the `examples/README.md` manifest; private models are **not**.

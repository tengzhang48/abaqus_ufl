# abaqus_ufl Public Example Map

This file is a positive release map, not the development lab's capability
inventory.

Before reusing an example:

1. read `examples/README.md`;
2. confirm that the folder exists in this checkout;
3. read its theory, non-scope, and evidence record; and
4. reuse only the rung that the example actually passed.

Do not infer that a missing research model is a package defect, and do not add
an internal, unpublished, license-unclear, planned, or not-yet-ported path to
this map.

## Pipeline reference

| Path | Reuse for | Scope |
|---|---|---|
| `examples/_template` | Complete minimal Python-to-compiled-UMAT pipeline | Working elastic demonstration; replace its physics and model-specific checks. |
| `HOWTO_ADD_AN_EXAMPLE.md` | Pipeline responsibilities and release gates | Canonical public contract; filenames remain example-owned. |
| `tools/` | Simple Abaqus run/extract/compare support | One-element/simple output only; difficult output bridges remain example-owned. |

## Released scientific examples

The authoritative allowlist is the table in `examples/README.md`. Keep any
scientific example entries here synchronized with that table after their
folders land in the public repository. Do not pre-advertise v1 targets.

## Authoritative docs

- `README.md`
- `HOWTO_ADD_AN_EXAMPLE.md`
- `docs/API_USAGE.md`
- `examples/README.md`

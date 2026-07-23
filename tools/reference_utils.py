"""Utilities for Abaqus validation reference generation.

The references here are Python material-point values for deliberately
homogeneous one-element Abaqus decks. Small-strain UMAT examples use the
package's compression-positive material API, then convert back to Abaqus'
tension-positive Voigt convention for comparison against ODB field output.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]


def load_symbol(relative_path: str, symbol: str) -> Any:
    """Load ``symbol`` from a Python file relative to the repository root."""
    path = ROOT / relative_path
    module_name = "abaqus_ufl_ref_" + "_".join(path.with_suffix("").parts[-3:])
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, symbol)


def _plain(value: Any) -> Any:
    """Convert NumPy scalars/arrays to JSON-serializable Python values."""
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        if abs(value.imag) > 1.0e-12:
            raise ValueError(f"Reference value has non-negligible imaginary part: {value}")
        return float(value.real)
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    return value


def tensor_to_abq_voigt_from_cp(tensor: np.ndarray) -> list[float]:
    """Compression-positive tensor -> Abaqus tension-positive Voigt."""
    return [
        -float(tensor[0, 0]),
        -float(tensor[1, 1]),
        -float(tensor[2, 2]),
        -float(tensor[0, 1]),
        -float(tensor[0, 2]),
        -float(tensor[1, 2]),
    ]


def tensor_to_abq_voigt(tensor: np.ndarray) -> list[float]:
    """Standard tensor -> Abaqus Voigt for stress-like output."""
    return [
        float(tensor[0, 0]),
        float(tensor[1, 1]),
        float(tensor[2, 2]),
        float(tensor[0, 1]),
        float(tensor[0, 2]),
        float(tensor[1, 2]),
    ]


def engineering_strain_voigt_x(final_stretch_x: float) -> list[float]:
    """Abaqus logarithmic strain Voigt for the one-axis NLGEOM decks."""
    return [math.log(final_stretch_x), 0.0, 0.0, 0.0, 0.0, 0.0]


def mises_from_tensor(stress: np.ndarray) -> float:
    """Von Mises equivalent stress from a 3x3 stress tensor."""
    mean = np.trace(stress) / 3.0
    dev = stress - mean * np.eye(3)
    return float(np.sqrt(1.5 * np.trace(dev.T @ dev)))


def statev_from_state(state: dict[str, Any]) -> list[float]:
    """Flatten state variables in the UMAT generator's STATEV order."""
    values: list[float] = []
    for value in state.values():
        arr = np.asarray(value)
        if arr.ndim == 0:
            values.append(float(arr))
        elif arr.shape == (3,):
            values.extend(float(arr[i]) for i in range(3))
        elif arr.shape == (3, 3):
            for j in range(3):
                for i in range(3):
                    values.append(float(arr[i, j]))
        else:
            raise ValueError(f"Unsupported state variable shape: {arr.shape}")
    return values


def small_strain_path(
    model_cls: type,
    *,
    final_stretch_x: float,
    ninc: int,
) -> dict[str, Any]:
    """Run a homogeneous small-strain UMAT material-point path.

    The Abaqus decks use ``nlgeom=YES`` and a single prescribed x stretch, so
    the total x strain passed to UMAT is the logarithmic strain. The material
    API is compression-positive; Abaqus field output is tension-positive.
    """
    model = model_cls()
    sigma = np.zeros((3, 3), dtype=float)
    strain = np.zeros((3, 3), dtype=float)
    state = {
        name: np.array(value, dtype=float, copy=True)
        for name, value in getattr(model, "state_vars", {}).items()
    }
    dt = 1.0 / ninc
    previous_abq_strain = 0.0

    params = list(model.stress_update.__signature__.parameters) if hasattr(
        model.stress_update, "__signature__") else None
    if params is None:
        import inspect
        params = list(inspect.signature(model.stress_update).parameters)

    for inc in range(1, ninc + 1):
        stretch_x = 1.0 + (final_stretch_x - 1.0) * (float(inc) / float(ninc))
        current_abq_strain = math.log(stretch_x)
        dstrain = np.diag([-(current_abq_strain - previous_abq_strain), 0.0, 0.0])
        previous_abq_strain = current_abq_strain
        args: list[Any] = [sigma, strain, dstrain]
        for name in params[3:]:
            if name == "dt":
                args.append(dt)
            elif name.endswith("_old"):
                state_name = name[:-4]
                args.append(state[state_name])
            else:
                raise ValueError(f"Unsupported stress_update argument: {name}")

        result = model.stress_update(*args)
        if isinstance(result, tuple):
            sigma, state_update = result
            state.update(state_update)
        else:
            sigma = result
        sigma = np.asarray(sigma, dtype=float)
        strain = strain + dstrain

    stress_abq = tensor_to_abq_voigt_from_cp(sigma)
    return {
        "stress_abq": stress_abq,
        "mises": mises_from_tensor(sigma),
        "strain_abq": engineering_strain_voigt_x(final_stretch_x),
        "statev": statev_from_state(state),
    }


def finite_strain_uniaxial_reference(
    model_cls: type,
    *,
    final_stretch_x: float,
) -> dict[str, Any]:
    """Compute Cauchy stress for the finite-strain UMAT one-element decks."""
    model = model_cls()
    F = np.diag([final_stretch_x, 1.0, 1.0])
    P = np.asarray(model.stress_PK1(F), dtype=float)
    J = float(np.linalg.det(F))
    sigma = (P @ F.T) / J
    return {
        "stress_abq": tensor_to_abq_voigt(sigma),
        "mises": mises_from_tensor(sigma),
    }


def stress_checks(ref: dict[str, Any], *, rtol: float = 1.0e-5) -> dict[str, Any]:
    return {
        "stress_voigt": {
            "path": "fields.S.data",
            "expected": ref["stress_abq"],
            "rtol": rtol,
            "atol": 1.0e-10,
        },
        "mises": {
            "path": "fields.S.mises",
            "expected": ref["mises"],
            "rtol": rtol,
            "atol": 1.0e-10,
        },
    }


def strain_checks(ref: dict[str, Any], *, rtol: float = 1.0e-5) -> dict[str, Any]:
    return {
        "strain_voigt": {
            "path": "fields.LE.data",
            "expected": ref["strain_abq"],
            "rtol": rtol,
            "atol": 1.0e-10,
        },
    }


def statev_checks(ref: dict[str, Any], *, rtol: float = 1.0e-5) -> dict[str, Any]:
    if not ref["statev"]:
        return {}
    return {
        "statev": {
            "path": "fields.SDV.data",
            "expected": ref["statev"],
            "rtol": rtol,
            "atol": 1.0e-10,
        },
    }


def write_reference(case_dir: Path, description: str, checks: dict[str, Any]) -> None:
    payload = {
        "description": description,
        "checks": checks,
    }
    path = case_dir / "reference.json"
    path.write_text(json.dumps(_plain(payload), indent=2) + "\n")
    print(f"Wrote {path}")

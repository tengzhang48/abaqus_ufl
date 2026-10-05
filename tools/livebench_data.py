"""Record numerical results from an already compiled verification run."""

import numpy as np


def comparison_curve(identifier, title, x_label, y_label, x, observed, reference,
                     description, *, reference_label="Closed form", rtol=1e-9,
                     atol=1e-11, extra_series=()):
    """Gate every plotted sample against its independent reference."""
    x, observed, reference = (np.asarray(v, dtype=float) for v in (x, observed, reference))
    if not (x.ndim == 1 and x.shape == observed.shape == reference.shape
            and len(x) >= 2 and np.isfinite([x, observed, reference]).all()):
        raise AssertionError("invalid response samples: " + identifier)
    if not np.allclose(observed, reference, rtol=rtol, atol=atol):
        raise AssertionError("plotted response disagrees with reference: " + identifier)
    return {
        "kind": "curve", "id": identifier, "title": title,
        "x_label": x_label, "y_label": y_label, "description": description,
        "x": x.tolist(), "x_scale": "linear", "y_scale": "linear",
        "series": [
            {"label": "Generated Fortran", "role": "computed", "values": observed.tolist()},
            {"label": reference_label, "role": "reference", "values": reference.tolist()},
            *extra_series,
        ],
        "comparison": {"max_abs_error": float(np.max(abs(observed - reference))),
                       "rtol": rtol, "atol": atol},
    }


def mesh_field(identifier, title, nodes, elems, values, description, *, displacements=None):
    result = {
        "kind": "mesh", "id": identifier, "title": title, "description": description,
        "nodes": np.asarray(nodes).tolist(), "elements": np.asarray(elems).tolist(),
        "values": np.asarray(values).tolist(),
    }
    if displacements is not None:
        result["displacements"] = np.asarray(displacements).tolist()
    return result

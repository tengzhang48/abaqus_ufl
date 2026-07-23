"""Generic ODB extractor for abaqus_ufl validation cases.

Usage:
    abaqus python extract_odb.py --config extract_config.json

The config JSON specifies:
    - odb_path: path to .odb file (default: job.odb)
    - frame: "last" or integer index (default: "last")
    - fields: list of field output names, e.g. ["S", "SDV", "E", "RF", "U"]
    - element_index: which element value to extract (default: 0, for single-element tests)
    - ip_index: which integration point (default: 0)
    - node_index: which node (for nodal output, default: 0)
    - node_label: specific node label to extract (for nodal output, overrides node_index)

Output: writes extracted values as JSON to stdout or --output file.
"""
from __future__ import print_function
import json
import sys
import os

try:
    from odbAccess import openOdb
except ImportError:
    # Allow dry-run testing outside Abaqus Python
    print("ERROR: odbAccess not available. Run this script with 'abaqus python'.", file=sys.stderr)
    sys.exit(1)


def parse_args():
    args = {"config": "extract_config.json", "output": None}
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--config" and i + 1 < len(sys.argv):
            args["config"] = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == "--output" and i + 1 < len(sys.argv):
            args["output"] = sys.argv[i + 1]
            i += 2
        else:
            i += 1
    return args


def load_config(path):
    with open(path, "r") as f:
        return json.load(f)


def abaqus_path(path):
    """Return a native string path for Abaqus Python 2 ODB APIs."""
    try:
        unicode_type = unicode
    except NameError:
        unicode_type = str
    if isinstance(path, unicode_type) and unicode_type is not str:
        return path.encode("utf-8")
    return path


def find_nodal_value(values, node_label=None, node_index=0):
    """Find a nodal value by label or index.

    If node_label is provided, searches for the value with matching nodeLabel.
    Otherwise falls back to node_index.
    """
    if node_label is not None:
        for v in values:
            if v.nodeLabel == node_label:
                return v
        return None
    if node_index >= len(values):
        return None
    return values[node_index]


def optional_float(value):
    try:
        return float(value)
    except Exception:
        return None


def numeric_data(data):
    if data is None:
        return None
    try:
        return [float(x) for x in data]
    except TypeError:
        return [float(data)]


def extract_field_value(field_output, element_index=0, ip_index=0, node_index=0,
                        node_label=None):
    """Extract scalar/tensor/vector data from a field output.

    For element output (integration point), uses element_index and ip_index.
    For nodal output, uses node_label (preferred) or node_index.
    """
    values = field_output.values
    if not values:
        return None

    # Determine if this is nodal or elemental
    if hasattr(values[0], "nodeLabel"):
        # Nodal output
        v = find_nodal_value(values, node_label=node_label, node_index=node_index)
        if v is None:
            return None
        data = v.data
        result = {"data": numeric_data(data)}
        if hasattr(v, "mises"):
            value = optional_float(v.mises)
            if value is not None:
                result["mises"] = value
        return result
    else:
        # Element / integration point output
        # Abaqus fieldOutput.values is a flat list of all integration points
        # for all elements. For a single-element test, we index directly.
        # For multi-element tests, the caller must know the flat index.
        # If ip_index is explicitly provided, use it as a flat index into
        # the values list. Otherwise fall back to element_index (works for
        # single-element single-IP tests where values[0] is the only entry).
        idx = ip_index if ip_index is not None else element_index
        if idx >= len(values):
            return None
        v = values[idx]
        data = v.data
        result = {"data": numeric_data(data)}
        if hasattr(v, "mises"):
            value = optional_float(v.mises)
            if value is not None:
                result["mises"] = value
        if hasattr(v, "maxPrincipal"):
            result["max_principal"] = optional_float(v.maxPrincipal)
        if hasattr(v, "minPrincipal"):
            result["min_principal"] = optional_float(v.minPrincipal)
        if hasattr(v, "maxInPlanePrincipal"):
            result["max_inplane_principal"] = optional_float(v.maxInPlanePrincipal)
        if hasattr(v, "minInPlanePrincipal"):
            result["min_inplane_principal"] = optional_float(v.minInPlanePrincipal)
        if hasattr(v, "outOfPlanePrincipal"):
            result["outofplane_principal"] = optional_float(v.outOfPlanePrincipal)
        return result


def sdv_sort_key(name):
    try:
        return int(name[3:])
    except Exception:
        return name


def extract_sdv_fields(frame, element_index=0, ip_index=0, node_index=0, node_label=None):
    names = [name for name in frame.fieldOutputs.keys() if name.startswith("SDV")]
    names.sort(key=sdv_sort_key)
    if not names:
        return None
    values = []
    for name in names:
        extracted = extract_field_value(
            frame.fieldOutputs[name],
            element_index=element_index,
            ip_index=ip_index,
            node_index=node_index,
            node_label=node_label,
        )
        if extracted and extracted.get("data") is not None:
            values.extend(extracted["data"])
    return {"data": values}


def extract_frame(odb, step_name=None, frame_spec="last"):
    """Get a specific frame from the ODB.

    frame_spec: "last" or integer index.
    If step_name is None, uses the last step.
    """
    steps_list = list(odb.steps.values())
    if step_name is None:
        step = steps_list[-1]
    else:
        step = odb.steps[step_name]

    if frame_spec == "last":
        frame_index = len(step.frames) - 1
    else:
        frame_index = int(frame_spec)
    frame = step.frames[frame_index]
    return frame, frame_index


def main():
    args = parse_args()
    config = load_config(args["config"])

    odb_path = abaqus_path(config.get("odb_path", "job.odb"))
    frame_spec = config.get("frame", "last")
    fields = [abaqus_path(field) for field in config.get("fields", ["S"])]
    element_index = config.get("element_index", 0)
    ip_index = config.get("ip_index", None)
    node_index = config.get("node_index", 0)
    node_label = config.get("node_label", None)
    step_name = config.get("step_name", None)
    if step_name is not None:
        step_name = abaqus_path(step_name)

    odb = openOdb(odb_path, readOnly=True)
    try:
        frame, frame_index = extract_frame(odb, step_name=step_name, frame_spec=frame_spec)
        steps_list = list(odb.steps.values())
        step_keys = list(odb.steps.keys())
        results = {
            "odb": odb_path,
            "step": step_name if step_name else step_keys[-1],
            "frame_index": frame_index,
            "frame_value": float(frame.frameValue),
            "fields": {}
        }

        for field_name in fields:
            if field_name not in frame.fieldOutputs:
                if field_name == "SDV":
                    results["fields"][field_name] = extract_sdv_fields(
                        frame,
                        element_index=element_index,
                        ip_index=ip_index,
                        node_index=node_index,
                        node_label=node_label,
                    )
                else:
                    results["fields"][field_name] = None
                continue
            field_output = frame.fieldOutputs[field_name]
            results["fields"][field_name] = extract_field_value(
                field_output,
                element_index=element_index,
                ip_index=ip_index,
                node_index=node_index,
                node_label=node_label
            )
    finally:
        odb.close()

    json_str = json.dumps(results, indent=2)
    if args["output"]:
        with open(args["output"], "w") as f:
            f.write(json_str)
        print("Wrote extracted results to", args["output"])
    else:
        print(json_str)


if __name__ == "__main__":
    main()

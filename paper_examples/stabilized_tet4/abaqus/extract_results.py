"""Extract the Figure 14 response from an Abaqus ODB.

Run with the Abaqus Python interpreter, for example:

    abq2022 python extract_results.py scovazzi_fig14_n4.odb
"""

from __future__ import print_function

import argparse
import json
import math
import os

from odbAccess import openOdb


def scalar(data):
    if isinstance(data, (tuple, list)):
        return float(data[0])
    try:
        return float(data)
    except TypeError:
        return float(data[0])


def values_by_node(frame, field_name):
    if field_name not in frame.fieldOutputs:
        available = ", ".join(sorted(frame.fieldOutputs.keys()))
        raise RuntimeError("{} is absent; available fields: {}".format(
            field_name, available))
    result = {}
    for value in frame.fieldOutputs[field_name].values:
        result[int(value.nodeLabel)] = value.data
    return result


def center_label(odb):
    containers = [odb.rootAssembly]
    containers.extend(odb.rootAssembly.instances.values())
    node_set = None
    available = []
    for container in containers:
        available.extend(container.nodeSets.keys())
        if "PUNCH_CENTER" in container.nodeSets:
            node_set = container.nodeSets["PUNCH_CENTER"]
            break
    if node_set is None:
        raise RuntimeError("PUNCH_CENTER is absent; available sets: {}".format(
            ", ".join(sorted(set(available)))))

    labels = []
    for item in node_set.nodes:
        if hasattr(item, "label"):
            labels.append(int(item.label))
        else:
            for node in item:
                labels.append(int(node.label))
    if len(labels) != 1:
        raise RuntimeError("PUNCH_CENTER must contain exactly one node")
    return labels[0]


def extract(odb_path, output_json, output_csv):
    odb = openOdb(path=odb_path, readOnly=True)
    try:
        if "Compression" not in odb.steps:
            raise RuntimeError("ODB does not contain the Compression step")
        step = odb.steps["Compression"]
        node_label = center_label(odb)
        rows = []

        for frame in step.frames:
            displacement = values_by_node(frame, "U")
            thetat = values_by_node(frame, "NT11")
            u_center = displacement[node_label]
            theta_center = scalar(thetat[node_label])

            u_magnitudes = []
            for data in displacement.values():
                u_magnitudes.append(math.sqrt(sum(float(item) ** 2 for item in data)))
            theta_values = []
            for data in thetat.values():
                theta_values.append(scalar(data))

            rows.append({
                "step_time": float(frame.frameValue),
                "pressure_N_per_mm2": 320.0 * float(frame.frameValue),
                "center_u3_mm": float(u_center[2]),
                "center_thetat": theta_center,
                "min_thetat": min(theta_values),
                "max_thetat": max(theta_values),
                "max_displacement_magnitude_mm": max(u_magnitudes),
            })

        final = rows[-1]
        summary = {
            "odb": os.path.abspath(odb_path),
            "step": "Compression",
            "frames": len(rows),
            "punch_center_node": node_label,
            "paper_center_abs_u3_mm_approx": 0.70,
            "final": final,
        }
        with open(output_json, "w") as stream:
            json.dump(summary, stream, indent=2, sort_keys=True)
            stream.write("\n")

        columns = [
            "step_time", "pressure_N_per_mm2", "center_u3_mm",
            "center_thetat", "min_thetat", "max_thetat",
            "max_displacement_magnitude_mm",
        ]
        with open(output_csv, "w") as stream:
            stream.write(",".join(columns) + "\n")
            for row in rows:
                stream.write(",".join("{:.16g}".format(row[key]) for key in columns) + "\n")
        return summary
    finally:
        odb.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("odb")
    parser.add_argument("--json", default="scovazzi_fig14_results.json")
    parser.add_argument("--csv", default="scovazzi_fig14_history.csv")
    args = parser.parse_args()
    summary = extract(args.odb, args.json, args.csv)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

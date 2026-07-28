#!/usr/bin/env python3
"""Generate the native Abaqus deck for Scovazzi et al. Figure 14.

Gmsh creates an optimized unstructured Tet4 mesh of the classical quarter
block. The top face is partitioned at x=y=0.5 so the punch pressure covers
exactly one quarter of that face. The resulting nodes and elements are written
directly to a standalone Abaqus input file; Gmsh is not needed to run it.
"""

from __future__ import print_function

import argparse
import json
from pathlib import Path

import gmsh
import numpy as np


MU = 80.194
LAM = 400889.806
C_TAU_U = 2.0
C_TAU_THETA = 0.1
PRESSURE = 320.0
DEFAULT_SEED = 14032023


def _element_connectivity(dim, entity_tags, expected_name, expected_nodes):
    blocks = []
    for entity_tag in entity_tags:
        element_types, _, node_blocks = gmsh.model.mesh.getElements(dim, entity_tag)
        for element_type, node_block in zip(element_types, node_blocks):
            properties = gmsh.model.mesh.getElementProperties(element_type)
            name, number_of_nodes = properties[0], properties[3]
            if name == expected_name and number_of_nodes == expected_nodes:
                blocks.append(np.asarray(node_block, dtype=np.int64).reshape(-1, expected_nodes))
            elif len(node_block):
                raise RuntimeError(
                    "unexpected {}D Gmsh element type: {} ({} nodes)".format(
                        dim, name, number_of_nodes))
    if not blocks:
        raise RuntimeError("Gmsh generated no {} elements".format(expected_name))
    return np.vstack(blocks)


def _classify_top_patch_surface():
    # OpenCASCADE bounding boxes are padded by about 1e-7 for this geometry.
    tolerance = 1.0e-6
    top_surfaces = []
    patch_surfaces = []
    for _, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        if abs(zmin - 1.0) < tolerance and abs(zmax - 1.0) < tolerance:
            top_surfaces.append(tag)
            if (xmin >= -tolerance and ymin >= -tolerance
                    and xmax <= 0.5 + tolerance and ymax <= 0.5 + tolerance):
                patch_surfaces.append(tag)
    if not top_surfaces:
        raise RuntimeError("could not identify the z=1 top surface")
    if not patch_surfaces:
        raise RuntimeError("could not identify the 0.5 x 0.5 punch patch")
    return top_surfaces, patch_surfaces


def make_gmsh_mesh(n, seed):
    """Return optimized unstructured Tet4 mesh data for target size 1/n."""
    if n < 2 or n % 2:
        raise ValueError("n must be an even integer >= 2")

    h = 1.0 / n
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads1D", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", 1)
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", 1)
        gmsh.option.setNumber("Mesh.RandomSeed", seed)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)

        gmsh.model.add("scovazzi_fig14")
        boxes = []
        for x0 in (0.0, 0.5):
            for y0 in (0.0, 0.5):
                boxes.append((3, gmsh.model.occ.addBox(x0, y0, 0.0, 0.5, 0.5, 1.0)))
        gmsh.model.occ.fragment([boxes[0]], boxes[1:], removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()

        top_surfaces, patch_surfaces = _classify_top_patch_surface()
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")

        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        node_tags = np.asarray(node_tags, dtype=np.int64)
        coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
        order = np.argsort(node_tags)
        node_tags = node_tags[order]
        points = coordinates[order]
        tag_to_label = {int(tag): index + 1 for index, tag in enumerate(node_tags)}

        volume_entities = [tag for _, tag in gmsh.model.getEntities(3)]
        tet_tags = _element_connectivity(3, volume_entities, "Tetrahedron 4", 4)
        triangle_tags = _element_connectivity(2, patch_surfaces, "Triangle 3", 3)

        tets = np.asarray([
            [tag_to_label[int(tag)] for tag in element] for element in tet_tags
        ], dtype=np.int64)
        punch_triangles = np.asarray([
            [tag_to_label[int(tag)] for tag in element] for element in triangle_tags
        ], dtype=np.int64)

        p0 = points[tets[:, 0] - 1]
        matrices = np.stack(
            (points[tets[:, 1] - 1] - p0,
             points[tets[:, 2] - 1] - p0,
             points[tets[:, 3] - 1] - p0),
            axis=2,
        )
        negative = np.linalg.det(matrices) < 0.0
        if np.any(negative):
            saved = tets[negative, 1].copy()
            tets[negative, 1] = tets[negative, 2]
            tets[negative, 2] = saved

        # Keep a valid outward (+z) orientation. SFM3D3 pressure P follows the
        # element normal, so the deck applies a negative magnitude below.
        punch_areas = []
        for triangle in punch_triangles:
            xyz = points[triangle - 1]
            normal = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
            if normal[2] < 0.0:
                triangle[1], triangle[2] = triangle[2], triangle[1]
                normal *= -1.0
            if normal[2] <= 0.0:
                raise RuntimeError("degenerate punch surface triangle")
            punch_areas.append(0.5 * normal[2])

        tolerance = max(1.0e-9, h * 1.0e-7)
        x_symm = list(np.flatnonzero(np.abs(points[:, 0]) < tolerance) + 1)
        y_symm = list(np.flatnonzero(np.abs(points[:, 1]) < tolerance) + 1)
        bottom = list(np.flatnonzero(np.abs(points[:, 2]) < tolerance) + 1)
        top = list(np.flatnonzero(np.abs(points[:, 2] - 1.0) < tolerance) + 1)
        center_candidates = np.flatnonzero(
            np.linalg.norm(points - np.array([0.0, 0.0, 1.0]), axis=1) < tolerance)
        if len(center_candidates) != 1:
            raise RuntimeError("expected exactly one punch-center node")

        return {
            "points": points,
            "tets": tets,
            "punch_triangles": punch_triangles,
            "punch_areas": np.asarray(punch_areas),
            "x_symm": x_symm,
            "y_symm": y_symm,
            "bottom": bottom,
            "top": top,
            "punch_center": int(center_candidates[0] + 1),
            "top_surface_entities": top_surfaces,
            "punch_surface_entities": patch_surfaces,
            "gmsh_version": gmsh.__version__,
        }
    finally:
        gmsh.finalize()


def tet_metrics(points, tets):
    p0 = points[tets[:, 0] - 1]
    matrices = np.stack(
        (points[tets[:, 1] - 1] - p0,
         points[tets[:, 2] - 1] - p0,
         points[tets[:, 3] - 1] - p0),
        axis=2,
    )
    volumes = np.linalg.det(matrices) / 6.0
    if np.min(volumes) <= 0.0:
        raise RuntimeError("generated mesh contains an inverted or zero-volume Tet4")
    if not np.isclose(np.sum(volumes), 1.0, rtol=1.0e-9, atol=1.0e-9):
        raise RuntimeError("tetrahedra do not fill the unit cube")
    if len(np.unique(tets)) != len(points):
        raise RuntimeError("at least one generated node is unused")

    face_nodes = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
    minimum_angles = np.empty(len(tets), dtype=float)
    maximum_angles = np.empty(len(tets), dtype=float)
    for element_index, labels in enumerate(tets):
        xyz = points[labels - 1]
        normals = []
        for opposite, face in enumerate(face_nodes):
            a, b, c = (xyz[index] for index in face)
            normal = np.cross(b - a, c - a)
            if np.dot(normal, xyz[opposite] - a) > 0.0:
                normal *= -1.0
            normal /= np.linalg.norm(normal)
            normals.append(normal)
        angles = []
        for first in range(4):
            for second in range(first + 1, 4):
                cosine = np.clip(np.dot(normals[first], normals[second]), -1.0, 1.0)
                angles.append(180.0 - np.degrees(np.arccos(cosine)))
        minimum_angles[element_index] = min(angles)
        maximum_angles[element_index] = max(angles)
    return volumes, minimum_angles, maximum_angles


def append_nset(lines, name, labels, width=16):
    lines.append("*Nset, nset={}".format(name))
    for start in range(0, len(labels), width):
        lines.append(", ".join(str(value) for value in labels[start:start + width]))


def make_h_property_bins(volumes, maximum_bins=16):
    """Group elements by the paper's h=volume^(1/3) definition."""
    h_values = volumes ** (1.0 / 3.0)
    order = np.argsort(h_values)
    groups = np.array_split(order, min(maximum_bins, len(order)))
    bins = []
    for index, group in enumerate(groups, start=1):
        bins.append({
            "name": "UEL_H_{:02d}".format(index),
            "element_labels": list(np.asarray(group, dtype=int) + 1),
            "h_elem_mm": float(np.mean(h_values[group])),
        })
    return bins, h_values


def build_deck(mesh, n, seed, h_property_bins, output):
    points = mesh["points"]
    tets = mesh["tets"]
    punch_triangles = mesh["punch_triangles"]
    h = 1.0 / n

    lines = [
        "*Heading",
        "Scovazzi et al. (2023) Figure 14 - unstructured mixed Tet4 block",
        "** Generated by abaqus_block_fig14/generate_inp.py; do not edit manually.",
        "** Unit system: N, mm. Quarter block: 1 x 1 x 1 mm.",
        "** Gmsh target h={:.12g}, seed={}.".format(h, seed),
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Physical Constants, absolute zero=-273.15",
        "*Node",
    ]
    for label, xyz in enumerate(points, start=1):
        lines.append("{}, {:.12g}, {:.12g}, {:.12g}".format(label, *xyz))

    lines.extend([
        "** Mixed volume elements: (U1,U2,U3,thetat) at every node.",
        "*User Element, type=U1, nodes=4, coordinates=3, properties=5, variables=1, unsymm",
        "1, 2, 3, 11",
        "*Element, type=U1, elset=MIXED_TET4",
    ])
    for label, connectivity in enumerate(tets, start=1):
        lines.append("{}, {}, {}, {}, {}".format(label, *connectivity))
    for property_bin in h_property_bins:
        lines.append("*Elset, elset={}".format(property_bin["name"]))
        labels = property_bin["element_labels"]
        for start in range(0, len(labels), 16):
            lines.append(", ".join(str(value) for value in labels[start:start + 16]))
        lines.extend([
            "*UEL Property, elset={}".format(property_bin["name"]),
            "{:.12g}, {:.12g}, {:.12g}, {:.12g}, {:.12g}".format(
                MU, LAM, C_TAU_U, C_TAU_THETA, property_bin["h_elem_mm"]),
        ])
    lines.extend([
        "** Near-zero C3D4T overlay makes the UEL mesh and NT11 visible in the ODB.",
        "** Its stress and strain are dummy visualization values and are not requested.",
        "*Element, type=C3D4T, elset=VISUALIZATION",
    ])
    overlay_offset = 1000000
    for local_label, connectivity in enumerate(tets, start=1):
        lines.append("{}, {}, {}, {}, {}".format(
            overlay_offset + local_label, *connectivity))
    lines.extend([
        "*Material, name=VISUALIZATION_DUMMY",
        "*Elastic",
        "1.0e-9, 0.3",
        "*Conductivity",
        "1.0e-12",
        "*Density",
        "1.0e-12",
        "*Specific Heat",
        "1.0",
        "*Solid Section, elset=VISUALIZATION, material=VISUALIZATION_DUMMY",
        ",",
        "** Surface elements carry the follower punch pressure.",
        "*Element, type=SFM3D3, elset=PUNCH_SURFACE_ELEMENTS",
    ])
    surface_offset = 2000000
    for local_label, connectivity in enumerate(punch_triangles, start=1):
        lines.append("{}, {}, {}, {}".format(
            surface_offset + local_label, *connectivity))
    lines.append("*Surface Section, elset=PUNCH_SURFACE_ELEMENTS")

    append_nset(lines, "X_SYMM", mesh["x_symm"])
    append_nset(lines, "Y_SYMM", mesh["y_symm"])
    append_nset(lines, "BOTTOM", mesh["bottom"])
    append_nset(lines, "TOP", mesh["top"])
    append_nset(lines, "PUNCH_CENTER", [mesh["punch_center"]])
    append_nset(lines, "ALL_NODES", list(range(1, len(points) + 1)))

    lines.extend([
        "** Classical quarter-block constraints.",
        "*Boundary",
        "X_SYMM, 1, 1, 0.0",
        "Y_SYMM, 2, 2, 0.0",
        "BOTTOM, 3, 3, 0.0",
        "TOP, 1, 2, 0.0",
        "*Amplitude, name=LOAD_RAMP, definition=TABULAR",
        "0.0, 0.0, 1.0, 1.0",
        "*Step, name=Compression, nlgeom=YES, inc=10000, unsymm=YES",
        "*Coupled Temperature-displacement, deltmx=0.05",
        "0.01, 1.0, 1.0e-8, 0.05",
        "*Dload, amplitude=LOAD_RAMP",
        "PUNCH_SURFACE_ELEMENTS, P, {:.12g}".format(-PRESSURE),
        "*Restart, write, frequency=0",
        "*Output, field, frequency=5",
        "*Node Output",
        "U, NT, RF",
        "*Output, history, frequency=1",
        "*Node Output, nset=PUNCH_CENTER",
        "U3, NT11, RF3",
        "*End Step",
        "",
    ])
    output.write_text("\n".join(lines))


def generate(n, seed, output):
    mesh = make_gmsh_mesh(n, seed)
    volumes, minimum_angles, maximum_angles = tet_metrics(mesh["points"], mesh["tets"])
    punch_area = float(np.sum(mesh["punch_areas"]))
    if not np.isclose(punch_area, 0.25, rtol=1.0e-9, atol=1.0e-9):
        raise RuntimeError("punch area is not 0.25 mm2: {}".format(punch_area))

    h_property_bins, h_values = make_h_property_bins(volumes)
    build_deck(mesh, n, seed, h_property_bins, output)
    metadata = {
        "benchmark": "Scovazzi et al. (2023), Section 6.3, Figure 14",
        "input_file": output.name,
        "mesher": "Gmsh {} optimized Delaunay Tet4".format(mesh["gmsh_version"]),
        "n_divisions_equivalent": n,
        "target_h_mm": 1.0 / n,
        "seed": seed,
        "nodes": len(mesh["points"]),
        "mixed_tets": len(mesh["tets"]),
        "visualization_tets": len(mesh["tets"]),
        "punch_triangles": len(mesh["punch_triangles"]),
        "punch_center_node": mesh["punch_center"],
        "boundary_node_counts": {
            "x_symm": len(mesh["x_symm"]),
            "y_symm": len(mesh["y_symm"]),
            "bottom": len(mesh["bottom"]),
            "top": len(mesh["top"]),
        },
        "material_N_per_mm2": {"mu": MU, "lambda": LAM},
        "stabilization": {"c_tau_u": C_TAU_U, "c_tau_theta": C_TAU_THETA},
        "stabilization_h_mm": {
            "definition": "reference_tet_volume^(1/3), grouped into property bins",
            "bins": len(h_property_bins),
            "minimum": float(np.min(h_values)),
            "maximum": float(np.max(h_values)),
            "mean": float(np.mean(h_values)),
        },
        "pressure_N_per_mm2": PRESSURE,
        "paper_center_abs_u3_mm_approx": 0.70,
        "mesh_quality": {
            "total_volume_mm3": float(np.sum(volumes)),
            "min_tet_volume_mm3": float(np.min(volumes)),
            "max_tet_volume_mm3": float(np.max(volumes)),
            "mean_tet_volume_mm3": float(np.mean(volumes)),
            "min_dihedral_angle_deg": float(np.min(minimum_angles)),
            "max_dihedral_angle_deg": float(np.max(maximum_angles)),
            "tets_outside_abaqus_angle_guidance": int(np.count_nonzero(
                (minimum_angles < 10.0) | (maximum_angles > 160.0))),
            "punch_area_mm2": punch_area,
        },
    }
    metadata_path = output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata, metadata_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    output = args.output
    if output is None:
        output = here / "scovazzi_block_fig14_n{}.inp".format(args.n)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata, metadata_path = generate(args.n, args.seed, output)
    quality = metadata["mesh_quality"]
    print("Wrote {}".format(output))
    print("Wrote {}".format(metadata_path))
    print("nodes={} tets={} pressure_triangles={}".format(
        metadata["nodes"], metadata["mixed_tets"], metadata["punch_triangles"]))
    print("volume={:.12g} pressure_area={:.12g}".format(
        quality["total_volume_mm3"], quality["punch_area_mm2"]))
    print("dihedral range={:.3f}..{:.3f} deg; outside guidance={}".format(
        quality["min_dihedral_angle_deg"], quality["max_dihedral_angle_deg"],
        quality["tets_outside_abaqus_angle_guidance"]))


if __name__ == "__main__":
    main()

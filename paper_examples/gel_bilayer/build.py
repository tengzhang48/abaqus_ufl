"""Build a three-field Quad8 Chester-Anand gel comparison case.

This example targets the paper's 2D swell-induced bending deck, but uses the
standard Abaqus-UFL mixed formulation:

    u  - quadratic displacement on all 8 nodes
    p  - bilinear pressure on the 4 corner nodes
    mu - quadratic chemical potential on all 8 nodes

The original supplemental deck is a 4-node UEL mesh.  This script converts the
rubber and gel elements to serendipity Quad8 connectivity so the interface
remains conforming.
"""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import abaqus_ufl as au
from abaqus_ufl.core.tensor import det, exp, inv, log


HERE = Path(__file__).resolve().parent
CHESTER_DIR = HERE.parent
SOURCE_INP = CHESTER_DIR / "gelUEL" / "gelUEL" / "code" / "SwellInducedBending.inp"
VISUALIZATION_ELEMENT_OFFSET = 10000


class ChesterAnandGelMaterial(au.Material):
    """Neo-Hookean + Flory-Huggins gel material for the mixed u-p-mu UEL."""

    props = dict(
        G=1.0,
        K=100.0,
        chi=0.1,
        D=5.0e-9,
        mu0=0.0,
        Omega=1.0e-4,
        Rgas=8.314,
        theta=298.0,
        phi0=0.5,
    )

    def stress_PK1(self, F, p, mu):
        return self.G * (F - inv(F).T) + p * inv(F).T

    def pressure_resid(self, F, p, mu):
        J = det(F)
        Je = exp(p / self.K)
        phi = self.phi0 * Je / J
        RT = self.Rgas * self.theta
        return (
            mu
            - self.mu0
            - RT * (log(1.0 - phi) + phi + self.chi * phi**2)
            + (phi / self.phi0) * self.Omega * p
        )

    def solvent_flux(self, F, p, mu, grad_mu):
        J = det(F)
        Je = exp(p / self.K)
        phi = self.phi0 * Je / J
        C = F.T @ F
        Cinv = inv(C)
        cR0 = (1.0 - self.phi0) / self.Omega
        cR = cR0 + (self.phi0 - phi) / (self.Omega * phi)
        M = self.D * cR / (self.Rgas * self.theta)
        return -M * (Cinv @ grad_mu)

    def solvent_storage(self, F, F_old, p, p_old, dt):
        J = det(F)
        J_old = det(F_old)
        Je = exp(p / self.K)
        Je_old = exp(p_old / self.K)
        return (J / Je - J_old / Je_old) / (self.Omega * dt)


class ChesterAnandUPMuQuad8(au.WeakForm):
    """Plane-strain mixed Quad8 gel: quadratic u/mu and corner pressure."""

    material = ChesterAnandGelMaterial
    ndim = 2

    def define_fields(self):
        self.u = au.VectorField("u", degree=2)
        self.p = au.ScalarField("p", degree=1)
        self.mu = au.ScalarField("mu", degree=2)

    def momentum_equation(self, v, F, p, mu):
        return self.material.stress_PK1(F, p, mu)

    def pressure_equation(self, q, F, p, mu):
        return self.material.pressure_resid(F, p, mu)

    def transport_equation(self, w, F, p, mu, grad_mu, F_old, p_old, dt):
        c_dot = self.material.solvent_storage(F, F_old, p, p_old, dt)
        j_R = self.material.solvent_flux(F, p, mu, grad_mu)
        return c_dot, j_R


def _parse_nodes(lines: list[str]) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    in_nodes = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower().startswith("*node") and not nodes:
            in_nodes = True
            continue
        if in_nodes and stripped.startswith("*"):
            break
        if in_nodes and stripped and not stripped.startswith("**"):
            parts = [p.strip() for p in stripped.split(",")]
            if len(parts) >= 3:
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]))
    return nodes


def _parse_elements(lines: list[str], keyword: str) -> list[tuple[int, tuple[int, int, int, int]]]:
    elements: list[tuple[int, tuple[int, int, int, int]]] = []
    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower() == keyword.lower():
            in_block = True
            continue
        if in_block and stripped.startswith("*"):
            break
        if in_block and stripped and not stripped.startswith("**"):
            parts = [int(p.strip()) for p in stripped.split(",") if p.strip()]
            elements.append((parts[0], tuple(parts[1:5])))
    return elements


def _to_quad8(
    nodes: dict[int, tuple[float, float]],
    elements: list[tuple[int, tuple[int, int, int, int]]],
    next_node: int,
    edge_to_mid: dict[tuple[int, int], int],
) -> tuple[list[tuple[int, tuple[int, ...]]], int]:
    converted: list[tuple[int, tuple[int, ...]]] = []

    def midpoint(a: int, b: int) -> int:
        nonlocal next_node
        key = tuple(sorted((a, b)))
        if key in edge_to_mid:
            return edge_to_mid[key]
        xa, ya = nodes[a]
        xb, yb = nodes[b]
        edge_to_mid[key] = next_node
        nodes[next_node] = ((xa + xb) * 0.5, (ya + yb) * 0.5)
        next_node += 1
        return edge_to_mid[key]

    for elem_id, conn4 in elements:
        n1, n2, n3, n4 = conn4
        converted.append(
            (
                elem_id,
                (
                    n1,
                    n2,
                    n3,
                    n4,
                    midpoint(n1, n2),
                    midpoint(n2, n3),
                    midpoint(n3, n4),
                    midpoint(n4, n1),
                ),
            )
        )
    return converted, next_node


def _ids_to_lines(ids: list[int], per_line: int = 16) -> list[str]:
    lines = []
    for i in range(0, len(ids), per_line):
        lines.append(", ".join(f"{v:8d}" for v in ids[i : i + per_line]) + "\n")
    return lines


def _element_lines(keyword: str, elements: list[tuple[int, tuple[int, ...]]]) -> list[str]:
    out = [keyword + "\n"]
    for elem_id, conn in elements:
        out.append(", ".join(str(v) for v in (elem_id, *conn)) + "\n")
    return out


def _nodes_from_elements(elements: list[tuple[int, tuple[int, ...]]]) -> list[int]:
    return sorted({node for _, conn in elements for node in conn})


def _nodes_on_boundary(
    nodes: dict[int, tuple[float, float]], *, x: float | None = None, y: float | None = None
) -> list[int]:
    tol = 1.0e-10
    selected = []
    for node_id, (xn, yn) in nodes.items():
        if x is not None and abs(xn - x) > tol:
            continue
        if y is not None and abs(yn - y) > tol:
            continue
        selected.append(node_id)
    return sorted(selected)


def _validate_deck_labels(deck: str) -> dict[str, int]:
    """Reject duplicate labels and undefined element-connectivity nodes."""
    node_labels: list[int] = []
    element_labels: list[int] = []
    connectivity_nodes: set[int] = set()
    block = None

    for raw_line in deck.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("**"):
            continue
        if line.startswith("*"):
            keyword = line.split(",", 1)[0].lower()
            block = "node" if keyword == "*node" else "element" if keyword == "*element" else None
            continue
        if block not in {"node", "element"}:
            continue

        values = [value.strip() for value in line.split(",") if value.strip()]
        if block == "node":
            node_labels.append(int(values[0]))
        else:
            element_labels.append(int(values[0]))
            connectivity_nodes.update(int(value) for value in values[1:])

    duplicate_nodes = sorted(
        label for label in set(node_labels) if node_labels.count(label) > 1
    )
    if duplicate_nodes:
        raise ValueError(f"duplicate Abaqus node labels: {duplicate_nodes}")

    duplicate_elements = sorted(
        label for label in set(element_labels) if element_labels.count(label) > 1
    )
    if duplicate_elements:
        raise ValueError(f"duplicate Abaqus element labels: {duplicate_elements}")

    missing_nodes = sorted(connectivity_nodes.difference(node_labels))
    if missing_nodes:
        raise ValueError(f"element connectivity references undefined nodes: {missing_nodes}")

    return {
        "nodes": len(node_labels),
        "elements": len(element_labels),
        "connectivity_nodes": len(connectivity_nodes),
    }


def _abaqus_data_lines(values: tuple[str, ...], per_line: int = 8) -> str:
    """Format keyword data without exceeding Abaqus's eight-value line limit."""
    return "".join(
        ", ".join(values[start : start + per_line]) + "\n"
        for start in range(0, len(values), per_line)
    )


def _write_material_include(path: Path) -> None:
    header = """*******************************************************************
** Material properties for the generated u-p-mu Quad8 UEL.
** PROPS: G, K, chi, D, mu0, Omega, Rgas, theta, phi0
*******************************************************************
"""
    values = (
        "1.e6",
        "1.e8",
        "0.1",
        "5.e-9",
        "<mu0>",
        "1.e-4",
        "8.3145",
        "<theta>",
        "<phi0>",
    )
    path.write_text(header + _abaqus_data_lines(values), encoding="utf-8")


def _write_converted_input(path: Path) -> None:
    lines = SOURCE_INP.read_text(encoding="utf-8").splitlines()
    nodes = _parse_nodes(lines)
    source_node_ids = set(nodes)
    rubber4 = _parse_elements(lines, "*Element, type=CPE4H")
    gel4 = _parse_elements(lines, "*Element, type=U1")

    edge_to_mid: dict[tuple[int, int], int] = {}
    next_node = max(nodes) + 1
    rubber8, next_node = _to_quad8(nodes, rubber4, next_node, edge_to_mid)
    gel8, next_node = _to_quad8(nodes, gel4, next_node, edge_to_mid)
    rigid_ref_node = next_node
    temp_nodes = tuple(range(rigid_ref_node + 1, rigid_ref_node + 5))
    temp_element = max(eid for eid, _ in (*rubber8, *gel8)) + 1
    gel_visualization8 = [
        (VISUALIZATION_ELEMENT_OFFSET + eid, conn) for eid, conn in gel8
    ]
    max_x = max(x for x, _ in nodes.values())
    min_x = min(x for x, _ in nodes.values())
    max_y = max(y for _, y in nodes.values())

    real_nodes = sorted(n for n in nodes if n < 900000)
    rubber_nodes = _nodes_from_elements(rubber8)
    gel_nodes = _nodes_from_elements(gel8)
    top_nodes = _nodes_on_boundary(nodes, y=max_y)
    left_nodes = _nodes_on_boundary(nodes, x=min_x)
    right_nodes = _nodes_on_boundary(nodes, x=max_x)
    slider_bottom_node = min(right_nodes, key=lambda node: nodes[node][1])
    slider_top_node = max(right_nodes, key=lambda node: nodes[node][1])
    slider_nodes = [
        node
        for node in right_nodes
        if node not in (slider_bottom_node, slider_top_node)
    ]
    right_contact_nodes = [
        node for node in right_nodes if node in source_node_ids
    ]

    out: list[str] = [
        "*Heading\n",
        "Plane strain swell-induced bending, generated u-p-mu Quad8 variant\n",
        "** Converted from gelUEL/gelUEL/code/SwellInducedBending.inp\n",
        "** The generated UEL uses p on DOF 11 and mu on DOF 12.\n",
        "** Stored mu = physical mu - initMU, so initial DOF 12 is zero.\n",
        "** GelVisualization duplicates the gel UEL with negligible-stiffness\n",
        "** CPE8R elements so Abaqus/CAE displays both physical mesh layers.\n",
        "** Its stress and strain are visualization artifacts, not UEL results.\n",
        "*Parameter\n",
        "phi0 = 0.999\n",
        "theta = 25.0 + 273.0\n",
        "initMU = -14382.975\n",
        "** Effective reference for the shifted computational mu field.\n",
        "mu0 = -initMU\n",
        "tf = 3600.0*6.0\n",
        "G = 50.e6\n",
        "** A small rubber compressibility avoids singular CPE8H pressure modes.\n",
        "nuRubber = 0.499\n",
        "** Abaqus neo-Hooke D1=2/K0 with K0=2G(1+nu)/(3(1-2nu)).\n",
        "D = 3.0*(1.0-2.0*nuRubber)/(G*(1.0+nuRubber))\n",
        "C10 = 0.5*G\n",
        "*Node\n",
    ]
    for node_id in real_nodes:
        x, y = nodes[node_id]
        out.append(f"{node_id}, {x:.12g}, {y:.12g}\n")

    out.extend(_element_lines("*Element, type=CPE8H, elset=RubberSet", rubber8))
    out.append("*User Element, Type=U1, Nodes=8, Coordinates=2, Properties=9, Variables=9, Unsymm\n")
    out.append("1, 2, 11, 12\n")
    out.append("5, 1, 2, 12\n")
    out.extend(_element_lines("*Element, type=U1, elset=GelSet", gel8))
    out.extend(
        _element_lines(
            "*Element, type=CPE8R, elset=GelVisualization",
            gel_visualization8,
        )
    )

    out.append("*Nset, nset=RubberSet\n")
    out.extend(_ids_to_lines(rubber_nodes))
    out.append("*Elset, elset=RubberSet\n")
    out.extend(_ids_to_lines([eid for eid, _ in rubber8]))
    out.append("*Nset, nset=GelSet\n")
    out.extend(_ids_to_lines(gel_nodes))
    out.append("*Elset, elset=GelSet\n")
    out.extend(_ids_to_lines([eid for eid, _ in gel8]))
    out.append("*Nset, nset=Nall\n")
    out.extend(_ids_to_lines(real_nodes))
    out.append("*Nset, nset=Top\n")
    out.extend(_ids_to_lines(top_nodes))
    out.append("*Nset, nset=Left\n")
    out.extend(_ids_to_lines(left_nodes))
    out.append("*Nset, nset=Right\n")
    out.extend(_ids_to_lines(right_nodes))
    out.append("*Nset, nset=Right1\n")
    out.extend(_ids_to_lines(slider_nodes))
    out.append(
        "** Match the original Quad4 node-to-surface contact sampling. The six Quad8\n"
        "** midside nodes remain in Right1 so the Slider MPC keeps the edge straight,\n"
        "** but putting them in both Slider and contact overconstrains the edge.\n"
    )
    out.append("*Nset, nset=RightContact\n")
    out.extend(_ids_to_lines(right_contact_nodes))
    out.append(f"*Nset, nset=n1\n{slider_bottom_node}\n")
    out.append(f"*Nset, nset=n2\n{slider_top_node}\n")
    out.append("*Nset, nset=PinnedCorner\n272\n")
    out.append("*Surface, Type=Node, Name=RightSurf\nRightContact\n")
    out.append("*Surface, Type=Node, Name=TopSurf\nTop\n")

    out.extend(
        [
            "*Node\n",
            f"{temp_nodes[0]}, 0.0, 0.0\n",
            f"{temp_nodes[1]}, 1.0e-5, 0.0\n",
            f"{temp_nodes[2]}, 1.0e-5, 1.0e-5\n",
            f"{temp_nodes[3]}, 0.0, 1.0e-5\n",
            "*Nset, nset=tempElement\n",
            ",".join(str(node) for node in temp_nodes) + "\n",
            "*Element, Type=CPE4T, elset=tempElement\n",
            f"{temp_element}," + ",".join(str(node) for node in temp_nodes) + "\n",
            "*Node\n",
            f"{rigid_ref_node},0.0,0.e-3\n",
            "*Nset, Nset=RigidRef\n",
            f"{rigid_ref_node}\n",
            "*Surface, type=SEGMENTS, name=RigidSurf\n",
            "START,0.0,   0.e-3\n",
            "LINE, 0.0, -50.e-3\n",
            "*Rigid Body, ref node=RigidRef, analytical surface=RigidSurf\n",
            "*uel property, elset=GelSet\n",
            "*Include, Input=ElasticGelProps_upmu_quad8.inp\n",
            "*Solid section, elset=RubberSet, material=RubberMaterial\n",
            "*Material, name=RubberMaterial\n",
            "*Hyperelastic, neo hooke\n",
            "<C10>,<D>\n",
            "** Display-only overlay. Do not interpret its S or E fields.\n",
            "*Solid section, elset=GelVisualization, material=VisualizationMaterial\n",
            "*Material, name=VisualizationMaterial\n",
            "*Elastic\n",
            "1.e-20, 0.3\n",
            "*Solid section, elset=tempElement, material=Material-2\n",
            "*Material, name=Material-2\n",
            "*Elastic\n",
            "1.e-20\n",
            "*Conductivity\n",
            "1.0\n",
            "*Density\n",
            "1.0\n",
            "*Specific heat\n",
            "1.0\n",
            "*MPC\n",
            "Slider, right1, n1, n2\n",
            "*Contact Pair, Interaction=int, Type=node to surface\n",
            "RightSurf,RigidSurf\n",
            "TopSurf,RigidSurf\n",
            "*Surface Interaction, name=int\n",
            "*Friction\n",
            "0.0\n",
            "*Surface Behavior, pressure-overclosure=hard\n",
            "** Computational mu=0 is the initial dry physical chemical potential.\n",
            "** Abaqus temperature initial conditions address DOF 11, which is p here.\n",
            "*Initial conditions, type=temperature\n",
            "tempElement, 0.0\n",
            "** This decay is physical mu(t)-initMU: zero initially, -initMU finally.\n",
            "*Amplitude, name=chemProfile, definition=decay\n",
            "<mu0>,<initMU>,0.0,300.0\n",
            "*Step, Name=Swell, nlgeom=yes, inc=50000\n",
            "** DELTMX is omitted because Abaqus also applies it to pressure DOF 11.\n",
            "*Coupled temperature-displacement, creep=none\n",
            "10.0,<tf>,1.e-6,100.0\n",
            "*Controls, Parameters=Line Search\n",
            "10,1.0,0.0001,0.25,0.10\n",
            "*Controls, Parameters=Time Incrementation\n",
            ",,,,,,,10,,,,,,\n",
            "*Boundary\n",
            "Left,1,1\n",
            "PinnedCorner,1,2\n",
            "RigidRef,encastre\n",
            "*Boundary, amplitude=chemProfile\n",
            "Top,12,12,1.0\n",
            "*Boundary\n",
            "tempElement,encastre\n",
            "tempElement,11,11,0.0\n",
            "*Output, field, number interval=100, time marks=no\n",
            "*node output, nset=Nall\n",
            "u\n",
            "*node output, nset=tempElement\n",
            "u\n",
            "*node output, nset=RigidRef\n",
            "u\n",
            "*element output, elset=RubberSet\n",
            "le\n",
            "*End Step\n",
        ]
    )
    deck = "".join(out)
    counts = _validate_deck_labels(deck)
    path.write_text(deck, encoding="utf-8")
    print(
        "Validated Abaqus labels: "
        f"{counts['nodes']} nodes, {counts['elements']} elements"
    )


def main() -> None:
    problem = ChesterAnandUPMuQuad8()
    print("Verifying Chester-Anand u-p-mu Quad8 material tangents...")
    if not problem.verify(verbose=False):
        raise RuntimeError("material tangent verification failed")

    au.generate_uel(
        problem,
        str(HERE / "chester_anand_upmu_quad8.for"),
        element="Quad8",
        formulation="standard",
        mat_prefix="chesteranandupmu",
    )
    _write_material_include(HERE / "ElasticGelProps_upmu_quad8.inp")
    _write_converted_input(HERE / "SwellInducedBending_upmu_quad8.inp")
    print("Generated chester_anand_upmu_quad8.for")
    print("Generated ElasticGelProps_upmu_quad8.inp")
    print("Generated SwellInducedBending_upmu_quad8.inp")


if __name__ == "__main__":
    main()

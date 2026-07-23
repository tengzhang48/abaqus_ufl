"""Abaqus input-file scaffold generation for UEL jobs.

This module intentionally implements a small, flat-mesh Abaqus keyword
scanner, not a full Abaqus input parser.  It is meant for the common workflow
where Abaqus/CAE exports nodes, elements, elsets, nsets, and sections in one
flat ``mesh.inp`` file, then abaqus_ufl derives the UEL include files from
that mesh without editing it in place.
"""

from dataclasses import dataclass, field
from pathlib import Path
import hashlib
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..core.fields import ScalarField, VectorField
from .element_config import ElementConfig, ELEMENT_CONFIGS
from .umat_gen import _nstate_per_gp, _state_var_info


MECH_DOFS = frozenset({1, 2, 3})
VALID_DOFS = frozenset({1, 2, 3, 8, 9, 11, 12})
EXTRA_ELEMENT_SUFFIX = {
    frozenset({11}): "T",
    frozenset({11, 12}): "T",
    frozenset({8}): "P",
    frozenset({9}): "E",
}

_SECTION_KEYWORDS = {
    "solid section",
    "shell section",
    "surface section",
    "beam section",
}

_DUMMY_VARIANTS = {
    ("CPE", 4): {
        "fully": "CPE4",
        "reduced": "CPE4R",
        "incompatible": "CPE4I",
        "hybrid": "CPE4H",
    },
    ("CPE", 8): {
        "fully": "CPE8",
        "reduced": "CPE8R",
        "hybrid": "CPE8H",
    },
    ("CPS", 4): {
        "fully": "CPS4",
        "reduced": "CPS4R",
        "incompatible": "CPS4I",
    },
    ("CPS", 8): {
        "fully": "CPS8",
        "reduced": "CPS8R",
    },
    ("CAX", 4): {
        "fully": "CAX4",
        "reduced": "CAX4R",
        "incompatible": "CAX4I",
        "hybrid": "CAX4H",
    },
    ("CAX", 8): {
        "fully": "CAX8",
        "reduced": "CAX8R",
        "hybrid": "CAX8H",
    },
    ("C3D", 8): {
        "fully": "C3D8",
        "reduced": "C3D8R",
        "incompatible": "C3D8I",
        "hybrid": "C3D8H",
    },
    ("C3D", 20): {
        "fully": "C3D20",
        "reduced": "C3D20R",
        "hybrid": "C3D20H",
    },
}


@dataclass(frozen=True)
class KeywordBlock:
    """One Abaqus keyword block and its raw line range."""

    keyword: str
    params: Dict[str, object]
    data_lines: List[str]
    start_line: int
    end_line: int
    raw_lines: List[str] = field(default_factory=list)


@dataclass
class UELModelConfig:
    """Configuration shared by generated UEL input scaffold files.

    Geometry is derived from ``ElementConfig`` to avoid drift between the
    Fortran generator and the input deck generator.
    """

    jtype: str
    element_config: ElementConfig
    dofs: Tuple[int, ...]
    n_properties: int
    variables: int
    target_elset: str
    n_iproperties: int = 0
    node_dofs: Optional[Tuple[Tuple[int, ...], ...]] = None
    element_family: Optional[str] = None
    dummy_elem_start: Optional[int] = None
    extra_node_start: Optional[int] = None
    extra_elem_id: Optional[int] = None
    dummy_element_type: Optional[str] = None
    dummy_element_variant: str = "fully"
    dummy_elastic_modulus: float = 1.0e-20
    dummy_poisson: float = 0.3
    dummy_density: Optional[float] = None
    dummy_expose_surfaces: bool = False
    uel_elset_suffix: str = "_uel"

    @classmethod
    def from_element(
        cls,
        jtype: str,
        element: str,
        dofs: Sequence[int],
        n_properties: int,
        variables: int,
        target_elset: str,
        **kwargs,
    ) -> "UELModelConfig":
        key = element.lower()
        if key not in ELEMENT_CONFIGS:
            raise ValueError(
                "Unknown element '{}'. Use one of: {}".format(
                    element, sorted(ELEMENT_CONFIGS)))
        return cls(
            jtype=jtype,
            element_config=ELEMENT_CONFIGS[key],
            dofs=tuple(dofs),
            n_properties=n_properties,
            variables=variables,
            target_elset=target_elset,
            **kwargs,
        )

    @classmethod
    def from_weakform(
        cls,
        weakform,
        jtype: str,
        element: str,
        target_elset: str,
        variables: Optional[int] = None,
        n_properties: Optional[int] = None,
        scalar_dofs: Optional[Dict[str, int]] = None,
        **kwargs,
    ) -> "UELModelConfig":
        """Build scaffold config from the same WeakForm used by generate_uel.

        Vector fields map to mechanical Abaqus DOFs ``1..ndim``. Scalar fields
        map to caller-provided labels when ``scalar_dofs`` is given; otherwise
        they receive consecutive labels starting at 11 in declaration order.
        """
        key = element.lower()
        if key not in ELEMENT_CONFIGS:
            raise ValueError(
                "Unknown element '{}'. Use one of: {}".format(
                    element, sorted(ELEMENT_CONFIGS)))
        element_config = ELEMENT_CONFIGS[key]
        if element_config.ndim != weakform.ndim:
            raise ValueError(
                "WeakForm ndim={} does not match element '{}' ndim={}"
                .format(weakform.ndim, element, element_config.ndim))
        dofs = _dofs_from_weakform(weakform, scalar_dofs=scalar_dofs)
        node_dofs = _node_dofs_from_weakform(
            weakform,
            element_config,
            scalar_dofs=scalar_dofs,
        )
        if n_properties is None:
            n_properties = len(weakform._mat.props_names)
        if variables is None:
            variables = _nstate_per_gp(_state_var_info(weakform._mat))
            variables *= element_config.n_gauss_points
        return cls(
            jtype=jtype,
            element_config=element_config,
            dofs=dofs,
            node_dofs=node_dofs,
            n_properties=n_properties,
            variables=variables,
            target_elset=target_elset,
            **kwargs,
        )

    @property
    def n_nodes(self) -> int:
        return self.element_config.n_nodes

    @property
    def n_coords(self) -> int:
        return self.element_config.ndim

    @property
    def n_int(self) -> int:
        return self.element_config.n_gauss_points

    @property
    def family(self) -> str:
        if self.element_family is not None:
            return self.element_family.upper()
        if self.n_coords == 3:
            return "C3D"
        return "CPE"

    @property
    def uel_elset_name(self) -> str:
        return self.target_elset + self.uel_elset_suffix

    def __post_init__(self) -> None:
        if not self.jtype.upper().startswith("U"):
            raise ValueError("jtype must be an Abaqus UEL type such as 'U1'")
        if not self.dofs:
            raise ValueError("dofs must not be empty")
        if len(set(self.dofs)) != len(self.dofs):
            raise ValueError(
                "dofs must be unique; got {}".format(list(self.dofs)))
        if self.node_dofs is not None:
            if len(self.node_dofs) != self.n_nodes:
                raise ValueError(
                    "node_dofs length must match n_nodes; got {} for {} nodes"
                    .format(len(self.node_dofs), self.n_nodes))
            for node_index, dofs in enumerate(self.node_dofs, start=1):
                if not dofs:
                    raise ValueError(
                        "node_dofs for node {} is empty".format(node_index))
                if tuple(sorted(dofs)) != tuple(dofs):
                    raise ValueError(
                        "node_dofs for node {} must be sorted in Abaqus node "
                        "DOF order".format(node_index))
                if len(set(dofs)) != len(dofs):
                    raise ValueError(
                        "node_dofs for node {} must be unique; got {}"
                        .format(node_index, list(dofs)))
        active_dofs = self.active_dofs()
        unknown = sorted(set(active_dofs) - VALID_DOFS)
        if unknown:
            raise ValueError("Unsupported Abaqus DOF labels: {}".format(unknown))
        if not self.target_elset:
            raise ValueError("target_elset is required")
        if self.n_properties < 0 or self.n_iproperties < 0 or self.variables < 0:
            raise ValueError("property and variable counts must be non-negative")
        nonmech = frozenset(set(active_dofs) - MECH_DOFS)
        if nonmech and nonmech not in EXTRA_ELEMENT_SUFFIX:
            raise ValueError(
                "Unsupported non-mechanical DOF combination: {}".format(
                    sorted(nonmech)))
        self.resolve_dummy_element_type()

    def active_dofs(self) -> Tuple[int, ...]:
        if self.node_dofs is None:
            return self.dofs
        return tuple(sorted(set(dof for node in self.node_dofs for dof in node)))

    def needs_extra_element(self) -> bool:
        return bool(set(self.active_dofs()) - MECH_DOFS)

    def resolve_dummy_element_type(self) -> str:
        if self.dummy_element_type is not None:
            return self.dummy_element_type.upper()
        key = (self.family, self.n_nodes)
        variants = _DUMMY_VARIANTS.get(key)
        if variants is None:
            raise ValueError(
                "No dummy element table entry for family={} n_nodes={}".format(
                    self.family, self.n_nodes))
        variant = self.dummy_element_variant.lower()
        if variant not in variants:
            raise ValueError(
                "Dummy variant '{}' is invalid for {}{}. Valid variants: "
                "{}".format(
                    self.dummy_element_variant, self.family, self.n_nodes,
                    sorted(variants)))
        return variants[variant]

    def resolve_extra_element_type(self) -> Optional[str]:
        nonmech = frozenset(set(self.active_dofs()) - MECH_DOFS)
        if not nonmech:
            return None
        suffix = EXTRA_ELEMENT_SUFFIX[nonmech]
        return "{}{}{}".format(self.family, self.n_nodes, suffix)


def _dofs_from_weakform(
    weakform,
    scalar_dofs: Optional[Dict[str, int]] = None,
) -> Tuple[int, ...]:
    scalar_dofs = scalar_dofs or {}
    dofs: List[int] = []
    next_scalar_dof = 11
    for name in weakform.field_names:
        field = weakform.fields[name]
        if isinstance(field, VectorField):
            if name != "u":
                raise ValueError(
                    "Only vector field 'u' can be mapped to Abaqus "
                    "mechanical DOFs automatically; pass an explicit "
                    "UELModelConfig for custom vector fields.")
            dofs.extend(range(1, weakform.ndim + 1))
        elif isinstance(field, ScalarField):
            if name in scalar_dofs:
                dof = scalar_dofs[name]
            else:
                dof = next_scalar_dof
                next_scalar_dof += 1
            dofs.append(dof)
        else:
            raise ValueError("Unsupported field type for '{}'".format(name))
    return tuple(dofs)


def _node_dofs_from_weakform(
    weakform,
    element_config: ElementConfig,
    scalar_dofs: Optional[Dict[str, int]] = None,
) -> Tuple[Tuple[int, ...], ...]:
    """Return Abaqus active DOFs for each element node.

    Mixed-order UELs need this because Abaqus's ``*User Element`` card defines
    active DOFs by node ranges. For example, a Quad8 with quadratic
    displacement but bilinear temperature has DOF 11 only on the four corner
    nodes, while all eight nodes keep mechanical DOFs 1 and 2.
    """
    scalar_dofs = scalar_dofs or {}
    scalar_labels: Dict[str, int] = {}
    next_scalar_dof = 11
    for name in weakform.field_names:
        field = weakform.fields[name]
        if isinstance(field, ScalarField):
            if name in scalar_dofs:
                scalar_labels[name] = scalar_dofs[name]
            else:
                scalar_labels[name] = next_scalar_dof
                next_scalar_dof += 1

    node_dofs: List[Tuple[int, ...]] = []
    for node_index in range(1, element_config.n_nodes + 1):
        is_corner = node_index <= element_config.n_corner_nodes
        dofs: List[int] = []
        for name in weakform.field_names:
            field = weakform.fields[name]
            field_has_node = is_corner or field.degree >= 2
            if not field_has_node:
                continue
            if isinstance(field, VectorField):
                if name != "u":
                    raise ValueError(
                        "Only vector field 'u' can be mapped to Abaqus "
                        "mechanical DOFs automatically; pass an explicit "
                        "UELModelConfig for custom vector fields.")
                dofs.extend(range(1, weakform.ndim + 1))
            elif isinstance(field, ScalarField):
                dofs.append(scalar_labels[name])
            else:
                raise ValueError("Unsupported field type for '{}'".format(name))
        node_dofs.append(tuple(dofs))
    return tuple(node_dofs)


@dataclass
class ScaffoldReport:
    num_elements: int
    num_nodes: int
    uel_elset_name: str
    target_elset_name: str
    dummy_elset_name: str
    dummy_element_type: str
    dummy_elem_start: int
    extra_element_emitted: bool
    extra_element_type: Optional[str]
    extra_node_start: Optional[int]
    extra_elem_id: Optional[int]
    mesh_hash: str
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    dry_run: bool = False


def _parse_keyword(line: str) -> Tuple[str, Dict[str, object]]:
    body = line.strip()[1:].strip()
    parts = [p.strip() for p in body.split(",")]
    keyword = parts[0].lower()
    params = {}
    for part in parts[1:]:
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            params[key.strip().lower()] = value.strip()
        else:
            params[part.lower()] = True
    return keyword, params


def scan_keywords(mesh_path: Path) -> List[KeywordBlock]:
    """Scan a flat Abaqus input file into keyword blocks."""
    lines = Path(mesh_path).read_text().splitlines(keepends=True)
    blocks = []
    current_start = None
    current_keyword = None
    current_params = None
    current_data = []
    current_raw = []

    def flush(end_line: int) -> None:
        if current_keyword is None:
            return
        blocks.append(
            KeywordBlock(
                keyword=current_keyword,
                params=current_params or {},
                data_lines=list(current_data),
                start_line=current_start or 1,
                end_line=end_line,
                raw_lines=list(current_raw),
            )
        )

    for idx, line in enumerate(lines, start=1):
        stripped = line.lstrip()
        if stripped.startswith("**") or not stripped.strip():
            if current_keyword is not None:
                current_raw.append(line)
            continue
        if stripped.startswith("*"):
            flush(idx - 1)
            current_keyword, current_params = _parse_keyword(stripped)
            current_start = idx
            current_data = []
            current_raw = [line]
        elif current_keyword is not None:
            current_data.append(line)
            current_raw.append(line)

    flush(len(lines))
    return blocks


def _block_line_ranges(blocks: Iterable[KeywordBlock]) -> set:
    line_numbers = set()
    for block in blocks:
        line_numbers.update(range(block.start_line, block.end_line + 1))
    return line_numbers


def _parse_ints(lines: Iterable[str]) -> List[int]:
    values = []
    for line in lines:
        clean = line.split("**", 1)[0].strip()
        if not clean:
            continue
        for token in clean.replace(",", " ").split():
            values.append(int(token))
    return values


def _expand_elset(block: KeywordBlock) -> List[int]:
    values = _parse_ints(block.data_lines)
    if block.params.get("generate"):
        if len(values) % 3 != 0:
            raise ValueError(
                "Generated elset at line {} must use start,end,step "
                "triples".format(block.start_line))
        out = []
        for i in range(0, len(values), 3):
            start, end, step = values[i], values[i + 1], values[i + 2]
            if step == 0:
                raise ValueError(
                    "Generated elset at line {} has zero step".format(
                        block.start_line))
            out.extend(range(start, end + (1 if step > 0 else -1), step))
        return out
    return values


def _parse_element_connectivity(
    block: KeywordBlock,
    expected_n_nodes: Optional[int] = None,
) -> Dict[int, List[int]]:
    if block.params.get("generate"):
        raise ValueError(
            "*Element, generate is not supported by the UEL scaffold scanner "
            "(line {}). Expand the element block first.".format(
                block.start_line))
    if expected_n_nodes is not None:
        tokens = _parse_ints(block.data_lines)
        record_len = expected_n_nodes + 1
        if len(tokens) % record_len != 0:
            raise ValueError(
                "*Element block at line {} has {} integer tokens, which is "
                "not divisible by element record length {}. This usually "
                "means the element type does not match the requested UEL "
                "node count, or the deck uses unsupported connectivity "
                "syntax.".format(block.start_line, len(tokens), record_len))
        conn = {}
        for i in range(0, len(tokens), record_len):
            elem_id = tokens[i]
            conn[elem_id] = tokens[i + 1:i + record_len]
        return conn

    conn = {}
    for line in block.data_lines:
        values = _parse_ints([line])
        if values:
            conn[values[0]] = values[1:]
    return conn


def _parse_node_ids(block: KeywordBlock) -> List[int]:
    ids = []
    for line in block.data_lines:
        clean = line.split("**", 1)[0].strip()
        if clean:
            ids.append(int(clean.replace(",", " ").split()[0]))
    return ids


def _find_blocks(blocks: Iterable[KeywordBlock], keyword: str) -> List[KeywordBlock]:
    keyword = keyword.lower()
    return [b for b in blocks if b.keyword == keyword]


def _elset_name(block: KeywordBlock) -> Optional[str]:
    value = block.params.get("elset")
    return str(value) if value is not None else None


def _names_match(a: Optional[str], b: str) -> bool:
    return a is not None and a.lower() == b.lower()


def _reject_assembly(blocks: Sequence[KeywordBlock]) -> None:
    assembly_keywords = {"part", "instance", "end instance", "assembly", "end assembly"}
    for block in blocks:
        if block.keyword in assembly_keywords:
            raise ValueError(
                "Assembly/instance mesh syntax is not supported in v1. "
                "Found *{} at line {}. Export a flat mesh from Abaqus/CAE."
                .format(block.keyword.title(), block.start_line))


def _reject_includes(blocks: Sequence[KeywordBlock]) -> None:
    for block in blocks:
        if block.keyword == "include":
            input_name = block.params.get("input", "(missing input=...)")
            raise ValueError(
                "Nested *Include files are not supported by the UEL scaffold "
                "scanner in v1. Found include '{}' at line {}. Pass the flat "
                "mesh/job deck that contains the target *Node, *Element, and "
                "*Elset blocks directly.".format(input_name, block.start_line))


def _all_elset_names(blocks: Sequence[KeywordBlock]) -> List[str]:
    names = []
    for block in blocks:
        if block.keyword in {"elset", "element"}:
            name = _elset_name(block)
            if name is not None:
                names.append(name)
    return sorted(set(names), key=str.lower)


def resolve_target_elset(
    blocks: Sequence[KeywordBlock],
    name: str,
    expected_n_nodes: Optional[int] = None,
) -> Tuple[List[KeywordBlock], List[int], Dict[int, List[int]]]:
    """Resolve a target elset to element blocks, IDs, and connectivity."""
    element_blocks = _find_blocks(blocks, "element")
    inline_blocks = [
        block for block in element_blocks
        if _names_match(_elset_name(block), name)
    ]
    if inline_blocks:
        conn = {}
        ids = []
        for block in inline_blocks:
            block_conn = _parse_element_connectivity(block, expected_n_nodes)
            conn.update(block_conn)
            ids.extend(block_conn)
        return inline_blocks, ids, conn

    elset_blocks = [
        block for block in _find_blocks(blocks, "elset")
        if _names_match(_elset_name(block), name)
    ]
    if not elset_blocks:
        raise ValueError(
            "Target elset '{}' was not found. Available elsets: {}"
            .format(name, ", ".join(_all_elset_names(blocks)) or "(none)"))

    target_ids = []
    for block in elset_blocks:
        target_ids.extend(_expand_elset(block))
    target_set = set(target_ids)

    matched_blocks = []
    conn = {}
    missing = set(target_ids)
    for block in element_blocks:
        try:
            block_conn = _parse_element_connectivity(block, expected_n_nodes)
        except ValueError:
            # A real Abaqus deck may contain retained element blocks of a
            # different type. Fall back to one-record-per-line parsing to
            # discover element IDs; if the target elset actually hits such a
            # block, generate() will reject it later on node-count mismatch.
            block_conn = _parse_element_connectivity(block)
        hit = target_set.intersection(block_conn)
        if hit:
            matched_blocks.append(block)
            for elem_id in sorted(hit):
                conn[elem_id] = block_conn[elem_id]
            missing.difference_update(hit)

    if missing:
        raise ValueError(
            "Target elset '{}' references element IDs not present in any "
            "*Element block: {}".format(name, sorted(missing)[:10]))
    if not matched_blocks:
        raise ValueError("Target elset '{}' resolved to no elements".format(name))
    return matched_blocks, target_ids, conn


def _collect_mesh_ids(
    blocks: Sequence[KeywordBlock],
    expected_n_nodes: Optional[int] = None,
) -> Tuple[List[int], List[int]]:
    elem_ids = []
    node_ids = []
    for block in blocks:
        if block.keyword == "element":
            try:
                conn = _parse_element_connectivity(block, expected_n_nodes)
            except ValueError:
                conn = _parse_element_connectivity(block)
            elem_ids.extend(conn.keys())
        elif block.keyword == "node":
            node_ids.extend(_parse_node_ids(block))
    return elem_ids, node_ids


def _target_section_blocks(
    blocks: Sequence[KeywordBlock],
    target_elset: str,
) -> List[KeywordBlock]:
    out = []
    for block in blocks:
        if block.keyword in _SECTION_KEYWORDS and _names_match(
            _elset_name(block), target_elset
        ):
            out.append(block)
    return out


def _target_related_elsets(
    blocks: Sequence[KeywordBlock],
    target_elset: str,
    target_ids: Sequence[int],
) -> Tuple[List[KeywordBlock], List[Tuple[str, List[int]]]]:
    """Find elsets that must move after the generated UEL element block.

    Abaqus input decks are order-sensitive. If an elset references elements
    that are stripped from mesh_nodes_only.inp and recreated later as UELs, the
    elset definition must move into uel_scaffold.inp after those UEL elements
    exist. V1 supports elsets wholly contained in the target UEL set; mixed
    retained/UEL elsets need splitting and are rejected clearly.
    """
    target = set(target_ids)
    strip_blocks: List[KeywordBlock] = []
    reemit: List[Tuple[str, List[int]]] = []

    for block in blocks:
        if block.keyword != "elset":
            continue
        ids = _expand_elset(block)
        overlap = target.intersection(ids)
        if not overlap:
            continue
        name = _elset_name(block)
        if not name:
            continue
        if not set(ids).issubset(target):
            raise ValueError(
                "Elset '{}' mixes target UEL elements and retained elements; "
                "v1 scaffold generation cannot split mixed elsets".format(name))

        strip_blocks.append(block)
        if not _names_match(name, target_elset):
            reemit.append((name, ids))

    return strip_blocks, reemit


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _format_keyword_user_element(config: UELModelConfig) -> List[str]:
    fields = [
        "*User Element",
        "Type={}".format(config.jtype.upper()),
        "Nodes={}".format(config.n_nodes),
        "Coordinates={}".format(config.n_coords),
        "Properties={}".format(config.n_properties),
    ]
    variables = config.variables if config.variables > 0 else 1
    fields.append("Variables={}".format(variables))
    if config.n_iproperties:
        fields.append("Iproperties={}".format(config.n_iproperties))
    header = ", ".join(fields)
    lines = [header + "\n"]
    if config.node_dofs is None:
        if config.dofs and config.dofs != tuple(sorted(config.dofs)):
            raise ValueError("dofs must be sorted in Abaqus node DOF order")
        lines.append(", ".join(str(d) for d in config.dofs) + "\n")
    else:
        previous: Optional[Tuple[int, ...]] = None
        for node_index, dofs in enumerate(config.node_dofs, start=1):
            if dofs == previous:
                continue
            if node_index == 1:
                row = list(dofs)
            else:
                row = [node_index] + list(dofs)
            lines.append(", ".join(str(d) for d in row) + "\n")
            previous = dofs
    return lines


def _format_element_block(
    elem_type: str,
    elset: str,
    connectivity: Sequence[Tuple[int, List[int]]],
) -> List[str]:
    lines = ["*Element, type={}, elset={}\n".format(elem_type, elset)]
    for elem_id, nodes in connectivity:
        row = [elem_id] + list(nodes)
        lines.append(", ".join(str(v) for v in row) + "\n")
    return lines


def _format_filtered_element_block(
    block: KeywordBlock,
    connectivity: Sequence[Tuple[int, List[int]]],
) -> List[str]:
    if not block.raw_lines:
        raise ValueError(
            "Element block at line {} has no raw keyword line".format(
                block.start_line))
    lines = [block.raw_lines[0]]
    for elem_id, nodes in connectivity:
        row = [elem_id] + list(nodes)
        lines.append(", ".join(str(v) for v in row) + "\n")
    return lines


def _format_elset(name: str, ids: Sequence[int]) -> List[str]:
    lines = ["*Elset, elset={}\n".format(name)]
    for i in range(0, len(ids), 16):
        lines.append(", ".join(str(v) for v in ids[i:i + 16]) + "\n")
    return lines


def _format_nset(name: str, ids: Sequence[int]) -> List[str]:
    lines = ["*Nset, nset={}\n".format(name)]
    for i in range(0, len(ids), 16):
        lines.append(", ".join(str(v) for v in ids[i:i + 16]) + "\n")
    return lines


def _extra_nodes(config: UELModelConfig, start: int) -> List[Tuple[int, List[float]]]:
    if config.n_coords == 2:
        coords = [
            [0.0, 0.0],
            [1.0e-5, 0.0],
            [1.0e-5, 1.0e-5],
            [0.0, 1.0e-5],
            [0.5e-5, 0.0],
            [1.0e-5, 0.5e-5],
            [0.5e-5, 1.0e-5],
            [0.0, 0.5e-5],
        ]
    else:
        coords = [
            [0.0, 0.0, 0.0],
            [1.0e-5, 0.0, 0.0],
            [1.0e-5, 1.0e-5, 0.0],
            [0.0, 1.0e-5, 0.0],
            [0.0, 0.0, 1.0e-5],
            [1.0e-5, 0.0, 1.0e-5],
            [1.0e-5, 1.0e-5, 1.0e-5],
            [0.0, 1.0e-5, 1.0e-5],
        ]
    if len(coords) < config.n_nodes:
        raise ValueError(
            "Extra element node template supports {} nodes, but config needs {}"
            .format(len(coords), config.n_nodes))
    return [(start + i, coords[i]) for i in range(config.n_nodes)]


def _write_mesh_nodes_only(
    mesh_path: Path,
    out_path: Path,
    remove_blocks: Sequence[KeywordBlock],
    target_elset: str,
    mesh_hash: str,
    replacements: Optional[Dict[int, List[str]]] = None,
) -> None:
    source_lines = Path(mesh_path).read_text().splitlines(keepends=True)
    removed = _block_line_ranges(remove_blocks)
    replacements = replacements or {}
    header = [
        "** =====================================================================\n",
        "** AUTO-GENERATED by abaqus_ufl. DO NOT EDIT.\n",
        "** Source: {} (sha256: {})\n".format(Path(mesh_path).name, mesh_hash[:12]),
        "** Stripped: target UEL element blocks for elset '{}'\n".format(target_elset),
        "** =====================================================================\n",
    ]
    kept = []
    for idx, line in enumerate(source_lines, start=1):
        if idx in replacements:
            kept.extend(replacements[idx])
        if idx in removed:
            continue
        kept.append(line)
    out_path.write_text("".join(header + kept))


def _target_element_block_rewrites(
    target_blocks: Sequence[KeywordBlock],
    target_ids: Sequence[int],
    target_elset: str,
    expected_n_nodes: int,
) -> Tuple[List[KeywordBlock], Dict[int, List[str]]]:
    target = set(target_ids)
    remove_blocks: List[KeywordBlock] = []
    replacements: Dict[int, List[str]] = {}

    for block in target_blocks:
        block_conn = _parse_element_connectivity(block, expected_n_nodes)
        non_target = [
            (elem_id, nodes)
            for elem_id, nodes in sorted(block_conn.items())
            if elem_id not in target
        ]
        remove_blocks.append(block)
        if not non_target:
            continue

        inline_elset = _elset_name(block)
        if inline_elset is not None and not _names_match(inline_elset, target_elset):
            raise ValueError(
                "Target elset '{}' selects only part of *Element block at "
                "line {}, which also defines inline elset '{}'. This is "
                "ambiguous because preserving the inline elset would mix "
                "retained elements and generated UEL elements. Split the "
                "element block in the source deck or remove the inline elset."
                .format(target_elset, block.start_line, inline_elset))

        replacements[block.start_line] = _format_filtered_element_block(
            block, non_target)

    return remove_blocks, replacements


def _write_uel_scaffold(
    out_path: Path,
    config: UELModelConfig,
    connectivity: Sequence[Tuple[int, List[int]]],
    dummy_connectivity: Sequence[Tuple[int, List[int]]],
    reemit_elsets: Sequence[Tuple[str, List[int]]],
    extra_node_start: Optional[int],
    extra_elem_id: Optional[int],
) -> None:
    lines = ["** Generated by abaqus_ufl. Do not edit manually.\n"]
    lines.extend(_format_keyword_user_element(config))

    elem_ids = [elem_id for elem_id, _ in connectivity]
    lines.extend(_format_element_block(config.jtype.upper(), config.target_elset,
                                       connectivity))
    if config.uel_elset_name.lower() != config.target_elset.lower():
        lines.extend(_format_elset(config.uel_elset_name, elem_ids))
    for name, ids in reemit_elsets:
        lines.extend(_format_elset(name, ids))

    dummy_type = config.resolve_dummy_element_type()
    lines.extend(_format_element_block(dummy_type, "elDummy", dummy_connectivity))
    if config.dummy_expose_surfaces:
        lines.extend([
            "*Surface, name=elDummy_surface, type=ELEMENT\n",
            "elDummy\n",
        ])

    if config.needs_extra_element():
        extra_type = config.resolve_extra_element_type()
        if extra_node_start is None or extra_elem_id is None:
            raise ValueError("extra node/element IDs are required")
        extra_nodes = _extra_nodes(config, extra_node_start)
        lines.append("*Node\n")
        for node_id, coords in extra_nodes:
            lines.append(
                ", ".join([str(node_id)] + ["{:.12g}".format(c) for c in coords])
                + "\n")
        lines.extend(_format_nset("extraElement", [n for n, _ in extra_nodes]))
        lines.extend(_format_element_block(
            extra_type or "",
            "extraElement",
            [(extra_elem_id, [n for n, _ in extra_nodes])],
        ))

    out_path.write_text("".join(lines))


def write_job_inp(
    out_dir: Path,
    heading: str = "abaqus_ufl generated UEL job",
    mesh_include: str = "mesh_nodes_only.inp",
    scaffold_include: str = "uel_scaffold.inp",
    materials_include: str = "materials.inp",
    steps_include: str = "steps.inp",
    filename: str = "job.inp",
) -> Path:
    """Write the top-level Abaqus include driver for a scaffolded UEL job."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / filename
    includes = [
        mesh_include,
        scaffold_include,
        materials_include,
        steps_include,
    ]
    lines = ["*Heading\n", heading + "\n"]
    for include in includes:
        lines.append("*Include, input={}\n".format(include))
    path.write_text("".join(lines))
    return path


def ensure_coupled_dummy_material(
    materials_path: Path,
    density: float = 1.0,
    specific_heat: float = 1.0,
) -> Path:
    """Ensure an auxiliary coupled thermal/scalar material is Abaqus-valid.

    Scaffolded UEL decks use an ``extraElement`` dummy Abaqus element to expose
    non-mechanical DOFs such as temperature or scalar fields. In a coupled
    temperature-displacement step, Abaqus requires the dummy material with
    ``*Conductivity`` to also define ``*Density`` and ``*Specific Heat`` before
    it accepts the input file. This helper is idempotent and only modifies
    files containing the scaffold's ``extra_material`` with conductivity.
    """
    path = Path(materials_path)
    text = path.read_text()
    if "*Material, name=extra_material" not in text or "*Conductivity" not in text:
        return path

    lines = text.splitlines()
    has_density = any(line.strip().lower() == "*density" for line in lines)
    has_specific_heat = any(
        line.strip().lower() == "*specific heat" for line in lines)
    if has_density and has_specific_heat:
        return path

    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        lower = line.strip().lower()
        if lower == "*conductivity" and not has_density:
            out.extend(["*Density", "{:.12g},".format(density)])
            has_density = True
        out.append(line)
        if lower == "*conductivity" and not has_specific_heat:
            if i + 1 < len(lines):
                i += 1
                out.append(lines[i])
            out.extend(["*Specific Heat", "{:.12g},".format(specific_heat)])
            has_specific_heat = True
        i += 1

    path.write_text("\n".join(out) + "\n")
    return path


def generate(
    config: UELModelConfig,
    mesh_path: Path,
    out_dir: Path,
    dry_run: bool = False,
) -> ScaffoldReport:
    """Generate UEL include files from a flat Abaqus mesh input file."""
    mesh_path = Path(mesh_path)
    out_dir = Path(out_dir)
    blocks = scan_keywords(mesh_path)
    _reject_assembly(blocks)
    _reject_includes(blocks)
    mesh_hash = _sha256(mesh_path)

    target_blocks, elem_ids, conn_by_id = resolve_target_elset(
        blocks, config.target_elset, expected_n_nodes=config.n_nodes)
    ordered_conn = [(elem_id, conn_by_id[elem_id]) for elem_id in sorted(conn_by_id)]
    for elem_id, nodes in ordered_conn:
        if len(nodes) != config.n_nodes:
            raise ValueError(
                "Element {} in target elset '{}' has {} nodes; config expects {}"
                .format(elem_id, config.target_elset, len(nodes), config.n_nodes))

    elem_ids_all, node_ids_all = _collect_mesh_ids(
        blocks, expected_n_nodes=config.n_nodes)
    max_elem = max(elem_ids_all) if elem_ids_all else 0
    max_node = max(node_ids_all) if node_ids_all else 0
    dummy_start = config.dummy_elem_start or (max_elem + 1)
    dummy_end = dummy_start + len(ordered_conn) - 1
    existing_elems = set(elem_ids_all)
    if any(elem_id in existing_elems for elem_id in range(dummy_start, dummy_end + 1)):
        raise ValueError("dummy_elem_start collides with existing element IDs")

    extra_node_start = None
    extra_elem_id = None
    if config.needs_extra_element():
        extra_node_start = config.extra_node_start or (max_node + 1)
        extra_nodes = range(extra_node_start, extra_node_start + config.n_nodes)
        if any(node_id in set(node_ids_all) for node_id in extra_nodes):
            raise ValueError("extra_node_start collides with existing node IDs")
        extra_elem_id = config.extra_elem_id or (dummy_end + 1)
        if extra_elem_id in existing_elems or dummy_start <= extra_elem_id <= dummy_end:
            raise ValueError("extra_elem_id collides with existing element IDs")

    dummy_conn = [
        (dummy_start + i, list(nodes))
        for i, (_, nodes) in enumerate(ordered_conn)
    ]
    section_blocks = _target_section_blocks(blocks, config.target_elset)
    target_elset_blocks, reemit_elsets = _target_related_elsets(
        blocks, config.target_elset, elem_ids)
    target_remove_blocks, replacements = _target_element_block_rewrites(
        target_blocks, elem_ids, config.target_elset, config.n_nodes)
    remove_blocks = target_remove_blocks + section_blocks + target_elset_blocks

    warnings = []
    if (config.dummy_elastic_modulus > 1.0e-15
            and not config.dummy_expose_surfaces
            and config.dummy_density is None
            and config.dummy_element_variant == "fully"):
        warnings.append(
            "dummy_elastic_modulus is non-negligible but no dummy capability "
            "is enabled")

    report = ScaffoldReport(
        num_elements=len(ordered_conn),
        num_nodes=len(set(node_ids_all)),
        uel_elset_name=config.uel_elset_name,
        target_elset_name=config.target_elset,
        dummy_elset_name="elDummy",
        dummy_element_type=config.resolve_dummy_element_type(),
        dummy_elem_start=dummy_start,
        extra_element_emitted=config.needs_extra_element(),
        extra_element_type=config.resolve_extra_element_type(),
        extra_node_start=extra_node_start,
        extra_elem_id=extra_elem_id,
        mesh_hash=mesh_hash,
        warnings=warnings,
        dry_run=dry_run,
    )

    if dry_run:
        return report

    out_dir.mkdir(parents=True, exist_ok=True)
    _write_mesh_nodes_only(
        mesh_path,
        out_dir / "mesh_nodes_only.inp",
        remove_blocks,
        config.target_elset,
        mesh_hash,
        replacements=replacements,
    )
    _write_uel_scaffold(
        out_dir / "uel_scaffold.inp",
        config,
        ordered_conn,
        dummy_conn,
        reemit_elsets,
        extra_node_start,
        extra_elem_id,
    )
    return report


__all__ = [
    "KeywordBlock",
    "UELModelConfig",
    "ScaffoldReport",
    "scan_keywords",
    "resolve_target_elset",
    "write_job_inp",
    "generate",
]

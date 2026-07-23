"""
Element configuration for UEL code generation.

Each ElementConfig describes the geometry, shape functions, quadrature,
and Jacobian mapping for one element type.  The generator reads these
to emit the correct Fortran declarations, subroutine calls, and loop
bounds — without hard-coding any element-specific details.

Built-in configs are in ELEMENT_CONFIGS.  Users may also construct
their own and pass them directly to ``generate_uel``.
"""

from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class ElementLayout:
    """Nodal DOF layout induced by an ElementConfig and field set."""

    corner_fields: list
    midside_fields: list
    ndof_corner: int
    ndof_midside: int
    n_corner: int
    n_midside: int
    n_nodes: int
    ndofel: int
    dof_map: dict
    node_dof_start: list
    field_nodes: dict
    field_shape: dict


@dataclass(frozen=True)
class ElementConfig:
    """Complete description of an element for UEL generation.

    Attributes
    ----------
    name : str
        Short identifier, e.g. ``'quad8'``, ``'hex20'``.
    ndim : int
        Spatial dimension (2 or 3).
    n_nodes : int
        Total nodes on the element (high-order count).
    n_corner_nodes : int
        Corner (vertex) nodes — used for low-order (degree-1) fields.

    gauss_subroutine : str
        Name of the Fortran quadrature subroutine, e.g.
        ``'gauss_2d_3x3'``.
    n_gauss_points : int
        Number of integration points returned by *gauss_subroutine*.
    gauss_xi_args : str
        Fortran snippet for the xi arguments to the shape subroutine,
        e.g. ``'xi_gp(kk,1), xi_gp(kk,2)'``.

    shape_subroutine : str
        High-order shape-function subroutine, e.g. ``'shape_quad8'``.
    shape_linear_subroutine : str
        Low-order (linear/corner) shape-function subroutine, e.g.
        ``'shape_quad4'``.

    jac_subroutine : str
        Jacobian-mapping subroutine, e.g. ``'map_grad_2d'``.
    jac_inv_size : int
        Dimension of Jinv (2 for 2-D, 3 for 3-D).

    supports_fbar : bool
        Whether an F-bar formulation exists for this element.
    plane_strain : bool
        If True the generated code sets F(3,3)=1 (2-D plane strain).

    template_files : list of str
        Filenames of shipped Fortran templates to include.
    """

    name: str
    ndim: int
    n_nodes: int
    n_corner_nodes: int

    # Quadrature
    gauss_subroutine: str
    n_gauss_points: int
    gauss_xi_args: str

    # Shape functions (high-order = all nodes, linear = corners only)
    shape_subroutine: str
    shape_linear_subroutine: str

    # Jacobian / mapping
    jac_subroutine: str
    jac_inv_size: int

    # Element capabilities
    supports_fbar: bool
    plane_strain: bool

    # Templates to include
    template_files: List[str] = field(default_factory=list)

    # -- Derived names used in generated Fortran ---------------------

    @property
    def sh_name(self):
        """Value-shape array name (high-order), e.g. ``'sh8'``."""
        return f'sh{self.n_nodes}'

    @property
    def dshxi_name(self):
        """Parametric-gradient shape array, e.g. ``'dshxi8'``."""
        return f'dshxi{self.n_nodes}'

    @property
    def dsh_name(self):
        """Physical-gradient shape array, e.g. ``'dsh8'``."""
        return f'dsh{self.n_nodes}'

    @property
    def sh_linear_name(self):
        """Value-shape array name (linear), e.g. ``'sh4'``."""
        return f'sh{self.n_corner_nodes}'

    @property
    def dshxi_linear_name(self):
        """Parametric-gradient shape array (linear), e.g. ``'dshxi4'``."""
        return f'dshxi{self.n_corner_nodes}'

    @property
    def dsh_linear_name(self):
        """Physical-gradient shape array (linear), e.g. ``'dsh4'``."""
        return f'dsh{self.n_corner_nodes}'

    @property
    def coords_name(self):
        """Coordinate array name, e.g. ``'coords_2d'`` or ``'coords_3d'``."""
        return f'coords_{self.ndim}d'

    # -- Helpers for emitting Fortran --------------------------------

    def shape_call(self, prefix=''):
        """Return the CALL line for the high-order shape subroutine."""
        return (f'{prefix}CALL {self.shape_subroutine}({self.gauss_xi_args},\n'
                f'     &{" " * max(1, 19 - len(self.shape_subroutine))}'
                f'{self.sh_name}, {self.dshxi_name})')

    def shape_linear_call(self, prefix=''):
        """Return the CALL line for the linear shape subroutine."""
        return (f'{prefix}CALL {self.shape_linear_subroutine}({self.gauss_xi_args},\n'
                f'     &{" " * max(1, 19 - len(self.shape_linear_subroutine))}'
                f'{self.sh_linear_name}, {self.dshxi_linear_name})')

    def sh_for_degree(self, degree):
        """Return (sh_name, dsh_name, n_nodes) for a given field degree."""
        if degree >= 2:
            return self.sh_name, self.dsh_name, self.n_nodes
        return self.sh_linear_name, self.dsh_linear_name, self.n_corner_nodes

    def field_node_indices(self, field_obj):
        """Return zero-based node indices that carry a field."""
        if field_obj.degree >= 2:
            return list(range(self.n_nodes))
        return list(range(self.n_corner_nodes))

    def field_shape_name(self, field_obj):
        """Return the shape-function family name for a field."""
        if field_obj.degree >= 2:
            return self.shape_subroutine.replace('shape_', '')
        return self.shape_linear_subroutine.replace('shape_', '')

    def build_layout(self, fields, field_order):
        """Build a complete nodal DOF layout for the declared fields."""
        corner_fields = []
        midside_fields = []
        for name in field_order:
            field_obj = fields[name]
            corner_fields.append(field_obj)
            if field_obj.degree >= 2:
                midside_fields.append(field_obj)

        ndof_corner = sum(f.ndof_per_node for f in corner_fields)
        ndof_midside = sum(f.ndof_per_node for f in midside_fields)
        n_midside = self.n_nodes - self.n_corner_nodes
        ndofel = (self.n_corner_nodes * ndof_corner
                  + n_midside * ndof_midside)

        dof_map = {}
        node_dof_start = []
        idx = 0
        for node in range(self.n_nodes):
            node_dof_start.append(idx)
            is_corner = node < self.n_corner_nodes
            node_fields = corner_fields if is_corner else midside_fields
            for field_obj in node_fields:
                n = field_obj.ndof_per_node
                dof_map[(node, field_obj.name)] = list(range(idx, idx + n))
                idx += n

        if idx != ndofel:
            raise RuntimeError(
                f"Element layout error for {self.name}: expected {ndofel} "
                f"DOFs, assigned {idx}")

        field_nodes = {}
        field_shape = {}
        for name in field_order:
            field_obj = fields[name]
            field_nodes[name] = self.field_node_indices(field_obj)
            field_shape[name] = self.field_shape_name(field_obj)

        return ElementLayout(
            corner_fields=corner_fields,
            midside_fields=midside_fields,
            ndof_corner=ndof_corner,
            ndof_midside=ndof_midside,
            n_corner=self.n_corner_nodes,
            n_midside=n_midside,
            n_nodes=self.n_nodes,
            ndofel=ndofel,
            dof_map=dof_map,
            node_dof_start=node_dof_start,
            field_nodes=field_nodes,
            field_shape=field_shape,
        )


# =====================================================================
# Built-in registry
# =====================================================================

ELEMENT_CONFIGS = {
    'quad8': ElementConfig(
        name='quad8',
        ndim=2,
        n_nodes=8,
        n_corner_nodes=4,
        gauss_subroutine='gauss_2d_3x3',
        n_gauss_points=9,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2)',
        shape_subroutine='shape_quad8',
        shape_linear_subroutine='shape_quad4',
        jac_subroutine='map_grad_2d',
        jac_inv_size=2,
        supports_fbar=False,
        plane_strain=True,
        template_files=['shape_quad8.for'],
    ),
    'quad8r': ElementConfig(
        name='quad8r',
        ndim=2,
        n_nodes=8,
        n_corner_nodes=4,
        gauss_subroutine='gauss_2d_2x2',
        n_gauss_points=4,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2)',
        shape_subroutine='shape_quad8',
        shape_linear_subroutine='shape_quad4',
        jac_subroutine='map_grad_2d',
        jac_inv_size=2,
        supports_fbar=False,
        plane_strain=True,
        template_files=['shape_quad8.for'],
    ),
    'quad4': ElementConfig(
        name='quad4',
        ndim=2,
        n_nodes=4,
        n_corner_nodes=4,
        gauss_subroutine='gauss_2d_2x2',
        n_gauss_points=4,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2)',
        shape_subroutine='shape_quad4',
        shape_linear_subroutine='shape_quad4',
        jac_subroutine='map_grad_2d',
        jac_inv_size=2,
        supports_fbar=True,
        plane_strain=True,
        template_files=['shape_quad4.for'],
    ),
    'hex8': ElementConfig(
        name='hex8',
        ndim=3,
        n_nodes=8,
        n_corner_nodes=8,
        gauss_subroutine='gauss_hex8',
        n_gauss_points=8,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2), xi_gp(kk,3)',
        shape_subroutine='shape_hex8',
        shape_linear_subroutine='shape_hex8',
        jac_subroutine='map_grad_3d',
        jac_inv_size=3,
        supports_fbar=True,
        plane_strain=False,
        template_files=['shape_hex8.for', 'gauss_hex.for', 'face_hex.for'],
    ),
    'tet4': ElementConfig(
        name='tet4',
        ndim=3,
        n_nodes=4,
        n_corner_nodes=4,
        gauss_subroutine='gauss_tet4',
        n_gauss_points=4,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2), xi_gp(kk,3)',
        shape_subroutine='shape_tet4',
        shape_linear_subroutine='shape_tet4',
        jac_subroutine='map_grad_3d',
        jac_inv_size=3,
        supports_fbar=False,
        plane_strain=False,
        template_files=['shape_tet4.for'],
    ),
    'tet4r': ElementConfig(
        name='tet4r',
        ndim=3,
        n_nodes=4,
        n_corner_nodes=4,
        gauss_subroutine='gauss_tet1',
        n_gauss_points=1,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2), xi_gp(kk,3)',
        shape_subroutine='shape_tet4',
        shape_linear_subroutine='shape_tet4',
        jac_subroutine='map_grad_3d',
        jac_inv_size=3,
        supports_fbar=False,
        plane_strain=False,
        template_files=['shape_tet4.for'],
    ),
    'hex20': ElementConfig(
        name='hex20',
        ndim=3,
        n_nodes=20,
        n_corner_nodes=8,
        gauss_subroutine='gauss_hex20',
        n_gauss_points=27,
        gauss_xi_args='xi_gp(kk,1), xi_gp(kk,2), xi_gp(kk,3)',
        shape_subroutine='shape_hex20',
        shape_linear_subroutine='shape_hex8',
        jac_subroutine='map_grad_3d',
        jac_inv_size=3,
        supports_fbar=False,
        plane_strain=False,
        template_files=['shape_hex20.for', 'shape_hex8.for', 'gauss_hex.for',
                        'face_hex.for'],
    ),
}

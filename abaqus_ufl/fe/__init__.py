"""Optional serial runtime for small standard Quad4 UEL verification meshes.

Install ``abaqus_ufl[fe]`` for SciPy, Meson, and Ninja. A Fortran compiler is
also required. This API supports uniform degree-one nodal fields, with history
carried by the previous nodal solution. Material/local SVARS are excluded.
"""

from .compiled_element import CompiledAbaqusElement, build_compiled_uel
from .driver import assemble, newton_solve
from .mesh import boundary_nodes, structured_quad_mesh

__all__ = ["CompiledAbaqusElement", "build_compiled_uel", "assemble",
           "newton_solve", "boundary_nodes", "structured_quad_mesh"]

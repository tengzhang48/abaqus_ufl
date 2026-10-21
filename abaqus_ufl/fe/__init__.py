"""Optional FE runtimes for standard Quad4 UEL verification meshes.

Install ``abaqus_ufl[fe]`` for SciPy, Meson, and Ninja. A Fortran compiler is
also required. This API supports uniform degree-one nodal fields, with history
carried by the previous nodal solution; material SVARS require explicit commit.
PETSc/MPI imports are lazy, so the serial host needs neither MPI library.
"""

from .compiled_element import CompiledAbaqusElement, IncrementTooLarge, build_compiled_uel
from .driver import assemble, assemble_batched, newton_solve
from .mesh import boundary_nodes, structured_quad_mesh
from .petsc_backend import make_petsc_solver, make_fieldsplit_solver
from .petsc_mpi import (partition_structured_quad, partition_rcb, partition_unstructured,
                        solve_steps_mpi, solve_steps_mpi_ghost, solve_steps_mpi_local)

__all__ = ["CompiledAbaqusElement", "IncrementTooLarge", "build_compiled_uel", "assemble", "assemble_batched",
           "newton_solve", "boundary_nodes", "structured_quad_mesh",
           "make_petsc_solver", "make_fieldsplit_solver", "partition_structured_quad",
           "partition_rcb", "partition_unstructured", "solve_steps_mpi",
           "solve_steps_mpi_ghost", "solve_steps_mpi_local"]

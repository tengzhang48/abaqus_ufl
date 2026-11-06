"""Optional FE runtimes for standard Quad4 UEL verification meshes.

Install ``abaqus_ufl[fe]`` for SciPy, Meson, and Ninja. A Fortran compiler is
also required. The compiled adapter supports uniform degree-one nodal fields,
with history carried by the previous nodal solution; material SVARS require
explicit commit. PETSc/MPI imports are lazy, so the serial host needs neither
MPI library.

``newton_solve`` solves one increment. The serial stepping drivers
(``solve_steps``, ``solve_steps_staggered``, ``solve_steps_multimaterial``)
and the bounded solver (``solve_steps_bounded``, with ``BoxSpec``) advance
many increments, commit a step only below its tolerance, and, where stated,
cut back rejected increments transactionally with ``CompiledElementState``.
"""

from .compiled_element import (CompiledAbaqusElement, IncrementTooLarge,
                               NonfiniteElementOutput, build_compiled_uel)
from .driver import assemble, assemble_batched, newton_solve
from .mesh import boundary_nodes, structured_quad_mesh
from .stepping import (StaggeredAttempt, StaggeredFieldReport, StaggeredReport,
                       StaggeredSolveFailed, StaggeredSweepReport,
                       assemble_groups, global_dof, solve_steps,
                       solve_steps_multimaterial, solve_steps_staggered)
from .bounds import (BoundedAttempt, BoundedReport, BoxSpec,
                     CompiledElementState, active_set,
                     activity_equivalent_vi_certificate, projected_residual,
                     solve_steps_bounded)
from .petsc_backend import make_petsc_solver, make_fieldsplit_solver
from .petsc_mpi import (partition_structured_quad, partition_rcb, partition_unstructured,
                        solve_steps_mpi, solve_steps_mpi_ghost, solve_steps_mpi_local)

__all__ = ["CompiledAbaqusElement", "IncrementTooLarge", "NonfiniteElementOutput",
           "build_compiled_uel", "assemble", "assemble_batched",
           "newton_solve", "boundary_nodes", "structured_quad_mesh",
           "solve_steps", "solve_steps_staggered", "solve_steps_multimaterial",
           "assemble_groups", "global_dof", "StaggeredSolveFailed",
           "StaggeredReport", "StaggeredAttempt", "StaggeredSweepReport",
           "StaggeredFieldReport", "solve_steps_bounded", "BoxSpec",
           "BoundedAttempt", "BoundedReport", "CompiledElementState",
           "active_set", "activity_equivalent_vi_certificate",
           "projected_residual",
           "make_petsc_solver", "make_fieldsplit_solver", "partition_structured_quad",
           "partition_rcb", "partition_unstructured", "solve_steps_mpi",
           "solve_steps_mpi_ghost", "solve_steps_mpi_local"]

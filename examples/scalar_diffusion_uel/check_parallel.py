"""Compiled Quad4 mechanics and heat diffusion on serial and MPI hosts.

Run from a shared working directory:
    mpiexec -n 2 python examples/scalar_diffusion_uel/check_parallel.py
The default .parallel-build folder retains generated source and build logs.
"""

import argparse
import json
import math
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from build import HeatDiffusionProblem
from abaqus_ufl.fe import (CompiledAbaqusElement, boundary_nodes, newton_solve,
                           solve_steps_mpi, structured_quad_mesh)
from abaqus_ufl.fe.mpi_batch import MPIBatchedElement, build_module


def check(build_dir='.parallel-build'):
    from petsc4py import PETSc
    comm = PETSc.COMM_WORLD
    problem = HeatDiffusionProblem()
    module = build_module(problem, build_dir)
    element = CompiledAbaqusElement(module, problem, cache_batches=True, dt=0.1)
    batched = MPIBatchedElement(module, problem, cache_batches=True, dt=0.1)
    root_solve = batched.root_linear_solver()
    nodes, elems = structured_quad_mesh(4, 2)
    steps = 2
    ndof = 3 * len(nodes)

    def patch_bc(frac):
        bc = {int(3*i+2): 5.0*frac for i in range(len(nodes))}
        for i in boundary_nodes(nodes):
            bc[int(3*i)] = 0.02 * frac * nodes[i, 0]
            bc[int(3*i+1)] = -0.01 * frac * nodes[i, 1]
        return bc

    def clock(step, frac):
        element.begin_increment(step)

    serial, replicated = np.zeros(ndof), np.zeros(ndof)
    timings = {}
    start = time.monotonic()
    for step in range(1, steps+1):
        element.begin_increment(step)
        serial = newton_solve(nodes, elems, 3, element.element_rk, patch_bc(step/steps),
                              U0=serial, U_previous=serial, batch_fn=element.element_rk_batch, tol=1e-12)
    timings['serial_patch_s'] = time.monotonic() - start
    start = time.monotonic()
    for step in range(1, steps+1):
        batched.begin_increment(step)
        replicated = newton_solve(nodes, elems, 3, batched.element_rk, patch_bc(step/steps),
                                  U0=replicated, U_previous=replicated,
                                  batch_fn=batched.element_rk_batch,
                                  linear_solve=root_solve, tol=1e-12)
    timings['mpi_batch_patch_s'] = time.monotonic() - start
    start = time.monotonic()
    distributed, info = solve_steps_mpi(
        nodes, elems, 3, element.element_rk, patch_bc, steps,
        pc='lu', solver='mumps', tol=1e-12, begin_step=clock)
    timings['distributed_patch_s'] = time.monotonic() - start
    exact = np.column_stack((0.02*nodes[:, 0], -0.01*nodes[:, 1], np.full(len(nodes), 5.0))).ravel()

    # A transient mode exercises DU and the previous accepted nodal solution.
    initial = np.zeros(ndof)
    initial[2::3] = np.sin(math.pi * nodes[:, 0])
    bc = {int(3*i+c): 0.0 for i in range(len(nodes)) for c in (0, 1)}
    for i in np.flatnonzero((nodes[:, 0] == 0) | (nodes[:, 0] == 1)):
        bc[int(3*i+2)] = 0.0
    serial_heat = initial.copy()
    for step in range(1, steps+1):
        element.begin_increment(step)
        serial_heat = newton_solve(nodes, elems, 3, element.element_rk, bc,
                                   U0=serial_heat, U_previous=serial_heat, tol=1e-12)
    distributed_heat, _ = solve_steps_mpi(
        nodes, elems, 3, element.element_rk, lambda frac: bc, steps,
        pc='lu', solver='mumps', tol=1e-12, begin_step=clock, u0=initial)
    h = 1/4
    cosine = math.cos(math.pi*h)
    eigenvalue = 6*(1-cosine)/(h*h*(2+cosine))
    heat_exact = initial[2::3] / (1 + 0.1*problem._mat.k/problem._mat.rho_cp*eigenvalue)**steps
    errors = {
        'mpi_batch_vs_serial': float(np.max(abs(replicated-serial))),
        'distributed_vs_serial': float(np.max(abs(distributed-serial))),
        'affine_patch': float(np.max(abs(distributed-exact))),
        'heat_vs_serial': float(np.max(abs(distributed_heat-serial_heat))),
        'heat_vs_discrete_oracle': float(np.max(abs(distributed_heat[2::3]-heat_exact))),
    }
    if max(errors.values()) > 1e-10:
        raise AssertionError(errors)
    result = dict(ranks=comm.getSize(), errors=errors, timings=timings,
                  cache_hits=batched.batch_cache_hits, batch_calls=batched.batch_calls,
                  ghost_dofs=info['n_ghost_total'])
    if comm.getRank() == 0:
        print(json.dumps(result, sort_keys=True))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', default='.parallel-build')
    check(parser.parse_args().build_dir)

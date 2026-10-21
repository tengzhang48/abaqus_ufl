"""Time the public Ogden Quad4 with four versus nine F perturbations.

Both use freshly generated compiled UELs, identical nontrivial plane-strain
inputs and uncached batches. Compilation is excluded from kernel timings.
"""

import argparse
import importlib.util
import json
from pathlib import Path
import statistics
import sys
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from abaqus_ufl.fe import CompiledAbaqusElement, build_compiled_uel, structured_quad_mesh
from abaqus_ufl.generators import uel_gen


def benchmark(build_dir=None, batches=200):
    if batches < 1:
        raise ValueError('batches must be positive')
    spec = importlib.util.spec_from_file_location('_benchmark_ogden', ROOT / 'examples/ogden_umat/build.py')
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    problem = example.OgdenPlaneStrainProblem()
    nodes, elems = structured_quad_mesh(8, 8)
    gm = (elems[:, :, None] * 2 + np.arange(2)).reshape(len(elems), 8)
    U = (nodes @ np.array([[0.13, 0.03], [0.09, -0.06]])).ravel()
    args = (nodes[elems], U[gm], 0.4*U[gm])
    original = uel_gen._generate_cs_engine
    timings, results = {}, []
    for ndim in (3, 2):
        folder = None if build_dir is None else Path(build_dir) / str(ndim)
        with patch.object(uel_gen, '_generate_cs_engine',
                          lambda wf, prefix: original(wf, prefix, ndim=ndim)):
            module = build_compiled_uel(problem, workdir=folder)
        element = CompiledAbaqusElement(module, problem, dt=0.1)
        results.append(element.element_rk_batch(*args))
        samples = []
        for _ in range(5):
            started = time.perf_counter()
            for repeat in range(batches):
                element.element_rk_batch(*args)
            samples.append((time.perf_counter()-started) / (batches*len(elems)))
        timings[str(ndim)] = statistics.median(samples) * 1e6
    residual_error, tangent_error = [float(np.max(abs(a-b))) for a, b in zip(*results)]
    if residual_error > 1e-12 or tangent_error > 1e-11:
        raise AssertionError('four/nine-perturbation results disagree')
    return dict(model='public Ogden Quad4', elements=len(elems), batches=batches,
                samples=5, microseconds_per_element=timings,
                speedup=timings['3']/timings['2'], max_residual_error=residual_error,
                max_tangent_error=tangent_error)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir')
    parser.add_argument('--batches', type=int, default=200)
    args = parser.parse_args()
    print(json.dumps(benchmark(args.build_dir, args.batches), sort_keys=True))

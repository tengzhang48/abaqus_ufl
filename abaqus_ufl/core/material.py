"""
Material base class for constitutive models.

Key design: Material methods operate at a single Gauss point. No FEM
concepts (shape functions, assembly) appear here. The Material class
is purely constitutive.

Method signatures define the tangent block structure:
  - Arguments that are field variables (F, p, mu, grad_mu) -> tangent blocks
  - Arguments that are history (F_old, p_old) -> no tangent (frozen)
  - Arguments that are parameters (dt) -> no tangent
  - Arguments matching state_vars names with '_old' suffix -> internal state

Internal variables (state_vars):
  Subclasses may define a state_vars dict mapping variable names to
  initial values. When present:
  - stress_PK1(self, F, ep_old, Fp_old, dt) receives state as args
  - Return tuple: (P, {'ep': ep_new, 'Fp': Fp_new})
  - state_old is always real; state_new is complex during CS tangent
  - verify() tests at both default and user-specified state

  state_vars iteration order determines STATEV layout in generated
  Fortran. Python 3.7+ preserves dict insertion order; for explicit
  control use collections.OrderedDict.

The framework detects which methods are defined and what arguments they
take to determine tangent blocks automatically.
"""

import inspect
import warnings
import numpy as np
from .defs import FIELD_ARGS, HISTORY_ARGS, PARAM_ARGS, METHOD_INFO


# ===================================================================
# Material base class
# ===================================================================

class Material:
    """
    Base class for constitutive models.

    Subclass and define:
      props = dict(G=1.0, K=100.0, ...)
      def stress_PK1(self, F): ...

    For materials with internal variables:
      state_vars = dict(ep=0.0, Fp=np.eye(3))
      def stress_PK1(self, F, ep_old, Fp_old, dt):
          ...
          return P, {'ep': ep_new, 'Fp': Fp_new}

    Then:
      model = MyMaterial(G=0.5)
      model.verify()
    """

    # Subclasses override these
    props = {}
    state_vars = {}

    def __init__(self, **kwargs):
        # Set default props, then apply overrides
        for name, default in self.__class__.props.items():
            setattr(self, name, kwargs.pop(name, default))

        if kwargs:
            raise TypeError(
                f"Unknown parameters: {list(kwargs.keys())}. "
                f"Valid props: {list(self.__class__.props.keys())}"
            )

        # Build set of state variable argument names (xxx_old)
        self._state_arg_names = set()
        if self.state_vars:
            for sv_name in self.state_vars:
                self._state_arg_names.add(f'{sv_name}_old')

        # Detect which material methods are defined. WeakForm refreshes this
        # metadata after user fields are declared so names such as T, c, eta,
        # and grad_T are recognized without hard-coded branches.
        self._field_args = dict(FIELD_ARGS)
        self._history_args = set(HISTORY_ARGS)
        self._param_args = set(PARAM_ARGS)
        self._detect_methods(warn_unknown=False)

    def refresh_metadata(self, field_args=None, history_args=None,
                         param_args=None):
        """Refresh signature metadata with problem-specific argument names."""
        self._field_args = dict(field_args or FIELD_ARGS)
        self._history_args = set(history_args or HISTORY_ARGS)
        self._param_args = set(param_args or PARAM_ARGS)
        self._detect_methods(warn_unknown=True)

    def _detect_methods(self, warn_unknown=True):
        """Detect material methods and classify their arguments."""
        self._methods = {}
        for method_name, info in METHOD_INFO.items():
            method = getattr(self, method_name, None)
            if method is not None and callable(method):
                sig = inspect.signature(method)
                params = [p for p in sig.parameters.keys()]

                # Classify each parameter
                field_params = [p for p in params
                                if p in self._field_args]
                history_params = [p for p in params
                                  if p in self._history_args]
                param_params = [p for p in params
                                if p in self._param_args]
                state_params = [p for p in params
                                if p in self._state_arg_names]

                # Warn about unknown argument names
                known = (set(self._field_args) | self._history_args
                         | self._param_args | self._state_arg_names)
                unknown = [p for p in params if p not in known]
                if unknown and warn_unknown:
                    warnings.warn(
                        f"{self.__class__.__name__}.{method_name}() "
                        f"has unknown argument(s): {unknown}. "
                        f"Recognized: {sorted(known)}. "
                        f"Unknown args will be ignored for tangent "
                        f"detection.",
                        stacklevel=2,
                    )

                self._methods[method_name] = {
                    'callable': method,
                    'all_params': params,
                    'field_params': field_params,
                    'history_params': history_params,
                    'param_params': param_params,
                    'state_params': state_params,
                    'output_shape': info['output_shape'],
                    'description': info['description'],
                }

    @property
    def has_state(self):
        """True if material defines state variables."""
        return bool(self.state_vars)

    @property
    def props_array(self):
        """Property values as ordered array (for Fortran props(*))."""
        return np.array([getattr(self, name)
                         for name in self.__class__.props])

    @property
    def props_names(self):
        """Property names in order."""
        return list(self.__class__.props.keys())

    def defined_methods(self):
        """Return list of material method names that are defined."""
        return list(self._methods.keys())

    def tangent_blocks(self):
        """
        Return dict describing all nonzero tangent blocks.

        Each entry: (method_name, field_arg) -> {
            'output_shape': ..., 'input_shape': ..., 'kind': ...
        }
        """
        blocks = {}
        for method_name, info in self._methods.items():
            for field in info['field_params']:
                finfo = self._field_args[field]
                blocks[(method_name, field)] = {
                    'output_shape': info['output_shape'],
                    'input_shape': finfo['shape'],
                    'kind': finfo['kind'],
                }
        return blocks

    def _make_state(self, **kwargs):
        """
        Build a state dict with default test values, overridden by kwargs.

        Includes field/history variables (F, p, mu, etc.) and internal
        state variables (as xxx_old, e.g. ep_old, Fp_old). Scalar state
        variables are stored as plain Python floats; arrays are copied.
        """
        state = {
            'F': np.array([[1.3, 0.1, 0.0],
                           [0.05, 1.2, 0.0],
                           [0.0, 0.0, 1.1]]),
            'p': 5.0,
            'mu': -100.0,
            'grad_mu': np.array([100.0, -50.0, 30.0]),
            'F_old': np.array([[1.25, 0.09, 0.0],
                               [0.04, 1.15, 0.0],
                               [0.0, 0.0, 1.05]]),
            'p_old': 4.5,
            'dt': 0.01,
            # Non-differentiated parameters (PARAM_ARGS). Representative values
            # so methods that take total time / Gauss-point reference position
            # can be verified; they are not perturbed in the tangent.
            'time': 1.0,
            'X': np.array([0.5, 0.5, 0.0]),
        }

        for name, info in self._field_args.items():
            if name in state:
                continue
            kind = info.get('kind')
            if kind == 'matrix':
                state[name] = np.eye(3)
            elif kind == 'vector':
                state[name] = np.zeros(3)
            else:
                state[name] = 0.0

        for name in self._history_args:
            if name in state:
                continue
            base = name[:-4] if name.endswith('_old') else None
            if base and base in self._field_args:
                kind = self._field_args[base].get('kind')
                if kind == 'matrix':
                    state[name] = np.eye(3)
                elif kind == 'vector':
                    state[name] = np.zeros(3)
                else:
                    state[name] = 0.0
            else:
                state[name] = 0.0

        # Add state variable defaults (as xxx_old)
        for sv_name, sv_init in self.state_vars.items():
            if np.isscalar(sv_init):
                state[f'{sv_name}_old'] = float(sv_init)
            else:
                state[f'{sv_name}_old'] = np.asarray(sv_init).copy()

        state.update(kwargs)
        return state

    def _call_method(self, method_name, state):
        """
        Call a material method with the right arguments from state dict.

        For methods returning (output, state_dict), returns only output.
        """
        info = self._methods[method_name]
        args = {p: state[p] for p in info['all_params']}
        result = info['callable'](**args)
        # Unpack tuple returns
        if isinstance(result, tuple):
            return result[0]
        return result

    def _call_method_full(self, method_name, state):
        """
        Call a material method and return the full result.

        For methods returning (output, state_dict), returns the tuple.
        For methods returning just output, returns (output, {}).
        """
        info = self._methods[method_name]
        args = {p: state[p] for p in info['all_params']}
        result = info['callable'](**args)
        if isinstance(result, tuple):
            return result
        return result, {}

    def _call_method_as_func(self, method_name):
        """
        Return a callable(**state) for use with the tangent engine.

        The returned function accepts keyword arguments matching the
        method's parameter names. Unpacks tuple returns.
        """
        info = self._methods[method_name]
        method = info['callable']
        param_names = info['all_params']

        def func(**kw):
            args = {p: kw[p] for p in param_names}
            result = method(**args)
            if isinstance(result, tuple):
                return result[0]
            return result

        return func

    def verify(self, state=None, tol=1e-6, verbose=True):
        """
        Run all verification checks.

        Checks:
          1. Stress-free reference (P=0 at F=I, p=0) — if stress_PK1
          2. CS vs FD for all tangent blocks
          3. Major symmetry of dP/dF (if stress_PK1 defined)

        For materials with state_vars, also tests at nonzero state
        to verify the algorithmic tangent from return mapping.

        Args:
            state: optional dict overriding default test state.
                   For materials with state_vars, can include
                   {'ep': 0.05, 'Fp': some_Fp} to test at nonzero
                   internal state.
            tol: tolerance for CS vs FD comparison
            verbose: print detailed results

        Returns:
            True if all checks pass
        """
        try:
            from .verify import run_verification
        except ImportError:
            from verify import run_verification

        # Build state dict
        if state is not None:
            # Convert user-friendly state names to _old suffix
            full_state = self._make_state()
            for k, v in state.items():
                if k in self.state_vars:
                    # User passed 'ep' -> set 'ep_old'
                    full_state[f'{k}_old'] = np.asarray(v).copy()
                else:
                    full_state[k] = v
            return run_verification(
                self, state=full_state, tol=tol, verbose=verbose)
        else:
            return run_verification(
                self, state=None, tol=tol, verbose=verbose)


class SmallStrainMaterial(Material):
    """Base class for small-strain incremental UMAT materials.

    User-facing stresses and strains are compression-positive. The Abaqus UMAT
    wrapper converts to and from Abaqus' tensile-positive `STRESS`, `STRAN`,
    and `DSTRAN` convention.

    Subclasses define:

        def stress_update(self, sigma_old, strain_old, dstrain,
                          state_old..., dt):
            return sigma_new, state_new

    All tensor arguments are 3x3 symmetric tensors in the compression-positive
    convention.
    """

    stress_convention = "compression_positive"

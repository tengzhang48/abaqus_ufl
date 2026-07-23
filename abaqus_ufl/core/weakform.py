"""
WeakForm base class.

Users subclass WeakForm and define:
  - material: a Material class (not instance)
  - define_fields(): declare VectorField/ScalarField on self
  - equation methods: momentum_equation, pressure_equation,
    transport_equation, phase_equation, pressure_transport_equation

The framework infers from these:
  - DOF layout per node (corner vs midside, which fields)
  - Shape function sets per field
  - Total NDOFEL
  - Residual structure (how each equation method maps to the weak form)
  - Tangent block structure (which derivatives are needed)

Equation return conventions:
  momentum_equation(v, F, ...)  -> P (tensor)    => integral P : Grad(N_a) dV
  pressure_equation(q, F, ...)  -> r_p (scalar)  => integral r_p * N_a dV
  transport-like equations return (storage, flux) and assemble
    integral storage*N_a dV - integral flux . Grad(N_a) dV

The first argument of each equation method is the test function name
(v, q, w) — it identifies which field's equation this is, but is not
used in computation (it's there for documentation/readability).
"""

import inspect
import numpy as np
from .fields import ScalarField, VectorField, LocalScalar
from .defs import FIELD_ARGS, HISTORY_ARGS, PARAM_ARGS
from ..generators.element_config import ELEMENT_CONFIGS


# Known equation methods and their structure
#
# Sign convention:
#   rhs_sign: the sign with which the material return enters RHS
#   amatrx_sign = -rhs_sign (because AMATRX = -dRHS/dU)
#
# For transport_equation (mixed): two sub-terms with different signs
#   storage (c_dot): enters as -c_dot*N  -> rhs_sign = -1, amatrx_sign = +1
#   flux (j_R):      enters as +j_R.dN   -> rhs_sign = +1, amatrx_sign = -1

EQUATION_INFO = {
    'momentum_equation': {
        'test_field_rank': 1,       # test function is vector
        'return_type': 'tensor',    # returns P (3x3)
        'assembly': 'grad',         # integral P : Grad(N) dV
        'rhs_sign': -1,             # RHS -= P : Grad(N) * wdetJ
        'description': 'Momentum balance (PK1 stress)',
    },
    'pressure_equation': {
        'test_field_rank': 0,       # test function is scalar
        'return_type': 'scalar',    # returns r_p
        'assembly': 'value',        # integral r_p * N dV
        'rhs_sign': -1,             # RHS -= r_p * N * wdetJ
        'description': 'Pressure constraint',
    },
    'transport_equation': {
        'test_field_rank': 0,       # test function is scalar
        'return_type': 'tuple',     # returns (c_dot, j_R)
        'assembly': 'mixed',        # storage + flux with different signs
        'description': 'Transport (storage + flux)',
        'sub_terms': {
            'storage': {
                'assembly': 'value',
                'rhs_sign': -1,     # RHS -= c_dot * N * wdetJ
                'material_method': 'solvent_storage',
            },
            'flux': {
                'assembly': 'grad',
                'rhs_sign': +1,     # RHS += j_R . Grad(N) * wdetJ
                'material_method': 'solvent_flux',
            },
        },
    },
    'phase_equation': {
        'test_field_rank': 0,
        'return_type': 'tuple',
        'assembly': 'mixed',
        'description': 'Phase-field balance (storage/value + flux)',
        'sub_terms': {
            'storage': {
                'assembly': 'value',
                'rhs_sign': -1,
                'material_method': 'phase_storage',
            },
            'flux': {
                'assembly': 'grad',
                'rhs_sign': +1,
                'material_method': 'phase_flux',
            },
        },
    },
    'pressure_transport_equation': {
        'test_field_rank': 0,
        'return_type': 'tuple',
        'assembly': 'mixed',
        'description': 'Pressure/diffusion balance (storage + flux)',
        'sub_terms': {
            'storage': {
                'assembly': 'value',
                'rhs_sign': -1,
                'material_method': 'pressure_storage',
            },
            'flux': {
                'assembly': 'grad',
                'rhs_sign': +1,
                'material_method': 'pressure_flux',
            },
        },
    },
    'species_transport_equation': {
        'test_field_rank': 0,
        'return_type': 'tuple',
        'assembly': 'mixed',
        'description': 'Species/concentration balance (storage + flux)',
        'sub_terms': {
            'storage': {
                'assembly': 'value',
                'rhs_sign': -1,
                'material_method': 'species_storage',
            },
            'flux': {
                'assembly': 'grad',
                'rhs_sign': +1,
                'material_method': 'species_flux',
            },
        },
    },
}

# Map test function name to field name
# The first arg of each equation method is the test function
# Convention: v->u, q->p, w->mu
TEST_TO_FIELD = {
    'v': 'u',
    'q': 'p',
    'w': 'mu',
    'r': 'p',
}


class WeakForm:
    """
    Base class for weak form definitions.

    Subclass and define:
      material = MyMaterialClass
      def define_fields(self): ...
      def momentum_equation(self, v, F, ...): ...
      [optional] def pressure_equation(self, q, F, ...): ...
      [optional] def transport_equation(self, w, F, ...): ...

    Then:
      problem = MyProblem(G=0.5, K=50.0)  # kwargs passed to Material
      problem.verify()
    """

    material = None  # subclass overrides with a Material class
    ndim = 2         # spatial dimension (2D default, 3D later)

    def __init__(self, ndim=None, **mat_kwargs):
        if ndim is not None:
            self.ndim = ndim

        # Instantiate material
        if self.material is None:
            raise TypeError("WeakForm subclass must define 'material' class attribute")
        self._mat = self.material(**mat_kwargs)

        # Call define_fields to populate field declarations
        self._fields = {}        # name -> Field object
        self._field_order = []   # field names in declaration order
        self._local_fields = {}  # name -> LocalScalar object
        self._local_field_order = []
        self._layout_config = ELEMENT_CONFIGS['quad8']
        self.define_fields()

        # Finalize field dimensions
        for name in self._field_order:
            field = self._fields[name]
            field.ndim = self.ndim
            if isinstance(field, VectorField):
                field.ndof_per_node = self.ndim

        self._build_argument_metadata()
        self._mat.refresh_metadata(self._field_arg_defs,
                                   self._history_args,
                                   self._param_args)

        # Compute DOF layout
        self._compute_dof_layout()

        # Detect equation methods
        self._equations = {}
        self._detect_equations()

    def define_fields(self):
        """Override in subclass to declare fields."""
        raise NotImplementedError("Subclass must implement define_fields()")

    def __setattr__(self, name, value):
        """Intercept field assignments in define_fields()."""
        if isinstance(value, (ScalarField, VectorField)):
            if not hasattr(self, '_fields'):
                # During __init__ before _fields exists
                super().__setattr__('_fields', {})
                super().__setattr__('_field_order', [])
            self._fields[value.name] = value
            self._field_order.append(value.name)
        elif isinstance(value, LocalScalar):
            if not hasattr(self, '_local_fields'):
                super().__setattr__('_local_fields', {})
                super().__setattr__('_local_field_order', [])
            self._local_fields[value.name] = value
            self._local_field_order.append(value.name)
        super().__setattr__(name, value)

    # =================================================================
    # DOF layout computation
    # =================================================================

    def _set_element_config(self, element_config):
        """Set the element config used for layout-sensitive metadata."""
        self._layout_config = element_config
        self._compute_dof_layout(element_config)

    def _compute_dof_layout(self, element_config=None):
        """Compute the DOF layout using ElementConfig topology rules."""
        cfg = element_config or self._layout_config
        self._layout_config = cfg
        layout = cfg.build_layout(self._fields, self._field_order)
        self._corner_fields = layout.corner_fields
        self._midside_fields = layout.midside_fields
        self._ndof_corner = layout.ndof_corner
        self._ndof_midside = layout.ndof_midside
        self._n_corner = layout.n_corner
        self._n_midside = layout.n_midside
        self._n_nodes = layout.n_nodes
        self._ndofel = layout.ndofel
        self._dof_map = layout.dof_map
        self._node_dof_start = layout.node_dof_start
        self._field_nodes = layout.field_nodes
        self._field_shape = layout.field_shape

    def _build_argument_metadata(self):
        """Build recognized field/history names from declared fields."""
        field_args = dict(FIELD_ARGS)
        history_args = set(HISTORY_ARGS)

        for name in self._field_order:
            field = self._fields[name]
            if isinstance(field, ScalarField):
                field_args[name] = {
                    'shape': (),
                    'kind': 'scalar',
                    'field': name,
                }
                field_args[f'grad_{name}'] = {
                    'shape': (3,),
                    'kind': 'vector',
                    'field': name,
                    'gradient': True,
                }
                history_args.add(f'{name}_old')
            elif isinstance(field, VectorField) and name != 'u':
                field_args[name] = {
                    'shape': (self.ndim,),
                    'kind': 'vector',
                    'field': name,
                }
                field_args[f'grad_{name}'] = {
                    'shape': (self.ndim, self.ndim),
                    'kind': 'matrix',
                    'field': name,
                    'gradient': True,
                }
                history_args.add(f'{name}_old')

        for name in self._local_field_order:
            field_args[name] = {
                'shape': (),
                'kind': 'scalar',
                'field': name,
                'local': True,
            }
            history_args.add(f'{name}_old')

        self._field_arg_defs = field_args
        self._history_args = history_args
        self._param_args = set(PARAM_ARGS)

    def _test_field_for(self, test_func_name):
        """Resolve a test-function argument name to a declared field name."""
        for fname in self._field_order:
            if getattr(self._fields[fname], 'test', None) == test_func_name:
                return fname
        mapped = TEST_TO_FIELD.get(test_func_name)
        if mapped in self._fields:
            return mapped
        return test_func_name

    # =================================================================
    # Equation detection
    # =================================================================

    def _detect_equations(self):
        """Detect which equation methods are defined and their signatures."""
        for eq_name, eq_info in EQUATION_INFO.items():
            method = getattr(self, eq_name, None)
            if method is not None and callable(method):
                sig = inspect.signature(method)
                params = list(sig.parameters.keys())

                # First param is the test function name
                test_func_name = params[0]
                test_field_name = self._test_field_for(test_func_name)

                # Remaining params are the field variables
                field_params = params[1:]

                # Classify params
                field_vars = [p for p in field_params
                              if p in self._field_arg_defs]
                history_vars = [p for p in field_params
                                if p in self._history_args]
                param_vars = [p for p in field_params
                              if p in self._param_args]

                self._equations[eq_name] = {
                    'callable': method,
                    'test_field': test_field_name,
                    'test_func': test_func_name,
                    'all_params': field_params,
                    'field_vars': field_vars,
                    'history_vars': history_vars,
                    'param_vars': param_vars,
                    'return_type': eq_info['return_type'],
                    'assembly': eq_info['assembly'],
                    'description': eq_info['description'],
                }

    # =================================================================
    # Public API
    # =================================================================

    @property
    def fields(self):
        """Return dict of field name -> Field object."""
        return dict(self._fields)

    @property
    def local_fields(self):
        """Return dict of element-local field name -> LocalScalar."""
        return dict(self._local_fields)

    @property
    def local_field_names(self):
        """Return local field names in declaration order."""
        return list(self._local_field_order)

    @property
    def field_names(self):
        """Return field names in declaration order."""
        return list(self._field_order)

    @property
    def ndofel(self):
        """Total number of DOFs per element."""
        return self._ndofel

    @property
    def equations(self):
        """Return dict of equation info."""
        return dict(self._equations)

    def dof_indices(self, field_name, node):
        """Return global DOF indices for a field at a given node."""
        return self._dof_map.get((node, field_name), [])

    def tangent_blocks(self):
        """
        Determine all tangent blocks needed for the element Jacobian.

        Each block is identified by:
          (equation_name, field_variable) -> {shape info, signs}

        For transport equations (mixed assembly), blocks are split into
        storage sub-blocks and flux sub-blocks with different signs.

        The amatrx_sign tells the assembly: AMATRX += amatrx_sign * T * N * N * wdetJ
        """
        blocks = {}
        for eq_name, eq_info in self._equations.items():
            eq_meta = EQUATION_INFO[eq_name]

            for field_var in eq_info['field_vars']:
                finfo = self._field_arg_defs.get(field_var)
                if finfo is None:
                    continue

                if eq_meta.get('assembly') == 'mixed':
                    # Transport equation: split into storage and flux
                    sub_terms = eq_meta['sub_terms']

                    # Determine which sub-term this field_var belongs to
                    # by checking which material method uses it
                    storage_method = sub_terms['storage']['material_method']
                    flux_method = sub_terms['flux']['material_method']

                    mat = self._mat
                    storage_info = mat._methods.get(storage_method, {})
                    flux_info = mat._methods.get(flux_method, {})

                    storage_fields = storage_info.get('field_params', [])
                    flux_fields = flux_info.get('field_params', [])

                    # A field_var can appear in storage, flux, or both
                    if field_var in storage_fields:
                        key = (eq_name + '_storage', field_var)
                        blocks[key] = {
                            'equation': eq_name,
                            'sub_term': 'storage',
                            'test_field': eq_info['test_field'],
                            'wrt_field': field_var,
                            'wrt_kind': finfo['kind'],
                            'row_assembly': 'value',
                            'col_assembly': 'grad' if finfo['kind'] == 'matrix' else ('grad' if finfo['kind'] == 'vector' else 'value'),
                            'amatrx_sign': +1,
                            'material_method': storage_method,
                        }

                    if field_var in flux_fields:
                        key = (eq_name + '_flux', field_var)
                        blocks[key] = {
                            'equation': eq_name,
                            'sub_term': 'flux',
                            'test_field': eq_info['test_field'],
                            'wrt_field': field_var,
                            'wrt_kind': finfo['kind'],
                            'row_assembly': 'grad',
                            'col_assembly': 'grad' if finfo['kind'] in ('matrix', 'vector') else 'value',
                            'amatrx_sign': -1,
                            'material_method': flux_method,
                        }
                else:
                    # Simple equation (momentum or constraint)
                    rhs_sign = eq_meta.get('rhs_sign', -1)
                    row_asm = eq_meta['assembly']

                    if finfo['kind'] == 'matrix':
                        col_asm = 'grad'
                    elif finfo['kind'] == 'vector':
                        col_asm = 'grad'
                    else:
                        col_asm = 'value'

                    blocks[(eq_name, field_var)] = {
                        'equation': eq_name,
                        'sub_term': None,
                        'test_field': eq_info['test_field'],
                        'wrt_field': field_var,
                        'wrt_kind': finfo['kind'],
                        'row_assembly': row_asm,
                        'col_assembly': col_asm,
                        'amatrx_sign': -rhs_sign,
                        'material_method': None,
                    }
        # Optionally drop selected (equation, wrt_field) tangent blocks from the
        # element Jacobian while keeping the RESIDUAL fully coupled. This realises
        # a targeted "staggered stiffness": e.g. dropping the mechanics-damage
        # block K_uphi removes the large g'(phi)*sigma0 coupling that inflates the
        # displacement correction in a phase-field front, without changing the
        # converged solution (the residual is unchanged). Declare on the WeakForm
        # subclass, e.g.:  drop_tangent_coupling = [('momentum_equation', 'phi')]
        drop = set(getattr(self, 'drop_tangent_coupling', None) or [])
        if drop:
            blocks = {
                k: v for k, v in blocks.items()
                if (v['equation'], v['wrt_field']) not in drop
            }

        # Sort by key for deterministic argument order in generated code
        from collections import OrderedDict
        return OrderedDict(sorted(blocks.items()))

    def summary(self):
        """Print a summary of the weak form definition."""
        print(f"\n{'='*60}")
        print(f"  WeakForm: {self.__class__.__name__}")
        print(f"  Material: {self._mat.__class__.__name__}")
        print(f"  Dimension: {self.ndim}D")
        print(f"{'='*60}")

        print(f"\n  Fields ({len(self._fields)}):")
        for name in self._field_order:
            f = self._fields[name]
            ftype = 'Vector' if isinstance(f, VectorField) else 'Scalar'
            shape = self._field_shape[name]
            nodes = len(self._field_nodes[name])
            print(f"    {name:8s}: {ftype:6s} degree={f.degree} "
                  f"({nodes} nodes, {shape}), "
                  f"{f.ndof_per_node} DOF/node")

        cfg = self._layout_config
        print(f"\n  DOF layout ({cfg.name}):")
        print(f"    Corner nodes (1-{self._n_corner}):  "
              f"{self._ndof_corner} DOFs each")
        if self._n_midside:
            first_mid = self._n_corner + 1
            print(f"    Higher-order nodes ({first_mid}-{self._n_nodes}): "
                  f"{self._ndof_midside} DOFs each")
        print(f"    NDOFEL = {self._ndofel}")

        if self._local_fields:
            print(f"\n  Local scalars ({len(self._local_fields)}):")
            for name in self._local_field_order:
                f = self._local_fields[name]
                print(f"    {name:8s}: storage={f.storage}, "
                      f"interpolation={f.interpolation}, "
                      f"condensed={f.condensed}")

        print(f"\n  Equations ({len(self._equations)}):")
        for eq_name, eq_info in self._equations.items():
            print(f"    {eq_name}:")
            print(f"      test field: {eq_info['test_field']} "
                  f"(test func: {eq_info['test_func']})")
            print(f"      depends on: {eq_info['field_vars']}")
            print(f"      returns: {eq_info['return_type']} "
                  f"({eq_info['assembly']})")

        blocks = self.tangent_blocks()
        print(f"\n  Tangent blocks ({len(blocks)}):")
        for (eq, var), info in blocks.items():
            print(f"    d({eq})/d({var})")

        n_perturbations = 0
        for (eq, var), info in blocks.items():
            kind = info['wrt_kind']
            if kind == 'matrix':
                n_perturbations += 9
            elif kind == 'vector':
                n_perturbations += 3
            else:
                n_perturbations += 1
        print(f"\n  Total complex perturbations per Gauss point: "
              f"{n_perturbations}")
        print(f"{'='*60}")

    def verify(self, **mat_kwargs):
        """
        Verify the material tangents used by this weak form.

        Delegates to Material.verify() with an appropriate test state.
        """
        return self._mat.verify(**mat_kwargs)

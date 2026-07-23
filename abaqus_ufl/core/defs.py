"""
Shared constants for Material method introspection.

Separated from material.py to avoid circular imports with verify.py.
"""

# Arguments that are field variables (differentiated in tangent)
FIELD_ARGS = {
    'F':       {'shape': (3, 3), 'kind': 'matrix'},
    'sigma_old': {'shape': (3, 3), 'kind': 'matrix'},
    'strain_old': {'shape': (3, 3), 'kind': 'matrix'},
    'dstrain': {'shape': (3, 3), 'kind': 'matrix'},
    'p':       {'shape': (),     'kind': 'scalar'},
    'mu':      {'shape': (),     'kind': 'scalar'},
    'grad_mu': {'shape': (3,),   'kind': 'vector'},
    'B':       {'shape': (3,),   'kind': 'vector'},
}

# Arguments that are history variables (not differentiated)
HISTORY_ARGS = {'F_old', 'p_old'}

# Arguments that are parameters (not differentiated)
# 'X' is the Gauss-point reference position (vector); 'time' the total time.
PARAM_ARGS = {'dt', 'time', 'X'}

# Known material methods and their return structure
METHOD_INFO = {
    'stress_PK1': {
        'output_shape': (3, 3),
        'description': 'First Piola-Kirchhoff stress',
    },
    'stress_update': {
        'output_shape': (3, 3),
        'description': 'Small-strain Cauchy stress update',
    },
    'pressure_resid': {
        'output_shape': (),
        'description': 'Pressure constraint residual',
    },
    'solvent_flux': {
        'output_shape': (3,),
        'description': 'Referential solvent flux',
    },
    'solvent_storage': {
        'output_shape': (),
        'description': 'Concentration rate',
    },
    'phase_flux': {
        'output_shape': (3,),
        'description': 'Phase-field gradient flux',
    },
    'phase_storage': {
        'output_shape': (),
        'description': 'Phase-field value/storage term',
    },
    'pressure_flux': {
        'output_shape': (3,),
        'description': 'Pressure/diffusion gradient flux',
    },
    'pressure_storage': {
        'output_shape': (),
        'description': 'Pressure/diffusion storage term',
    },
    'species_flux': {
        'output_shape': (3,),
        'description': 'Species/concentration gradient flux',
    },
    'species_storage': {
        'output_shape': (),
        'description': 'Species/concentration storage term',
    },
    'h_field': {
        'output_shape': (3,),
        'description': 'Lagrangian h-field H = dW/dB (magneto-mechanics)',
    },
    'volumetric_PK1': {
        'output_shape': (3, 3),
        'description': 'Volumetric PK1 part (selective reduced integration)',
    },
}

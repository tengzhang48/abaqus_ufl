"""
Field declarations for the WeakForm.

Fields describe the function spaces (polynomial degree, rank) for each
unknown in the weak form. The framework uses field declarations to:
  1. Determine DOF layout per node (which fields live on which nodes)
  2. Select shape function sets (degree=2 -> Quad8, degree=1 -> Quad4)
  3. Compute NDOFEL (total DOFs per element)
  4. Map between the flat Abaqus Uall array and per-field nodal arrays

For Quad8 (2D, 8-node serendipity):
  degree=2 fields: all 8 nodes (quadratic interpolation)
  degree=1 fields: corner nodes 1-4 only (bilinear interpolation)
  
Node numbering (Quad8):
    4-----7-----3        corners: 1,2,3,4
    |           |        midsides: 5,6,7,8
    8           6
    |           |
    1-----5-----2

Corner nodes carry DOFs for ALL fields (degree=1 and degree=2).
Midside nodes carry DOFs for degree=2 fields only.
"""


class ScalarField:
    """
    Scalar field declaration (1 DOF per node).

    Args:
        name: field name (used in equation method signatures)
        degree: polynomial degree (1=bilinear, 2=quadratic)
        test: optional weak-form test function name for this field
    """
    def __init__(self, name, degree=2, test=None):
        if degree not in (1, 2):
            raise ValueError(
                f"ScalarField '{name}' has unsupported degree {degree}. "
                "Only degree 1 and 2 are currently supported.")
        self.name = name
        self.degree = degree
        self.test = test
        self.rank = 0          # scalar
        self.ndof_per_node = 1
        self.ndim = None       # set by WeakForm after init

    def __repr__(self):
        return f"ScalarField('{self.name}', degree={self.degree})"


class VectorField:
    """
    Vector field declaration (ndim DOFs per node).

    Args:
        name: field name (used in equation method signatures)
        degree: polynomial degree (1=bilinear, 2=quadratic)
        test: optional weak-form test function name for this field
    """
    def __init__(self, name, degree=2, test=None):
        if degree not in (1, 2):
            raise ValueError(
                f"VectorField '{name}' has unsupported degree {degree}. "
                "Only degree 1 and 2 are currently supported.")
        self.name = name
        self.degree = degree
        self.test = test
        self.rank = 1          # vector
        self.ndof_per_node = None  # set to ndim by WeakForm
        self.ndim = None           # set by WeakForm after init

    def __repr__(self):
        return f"VectorField('{self.name}', degree={self.degree})"


class LocalScalar:
    """
    Element-local scalar declaration.

    Local scalars are metadata for variables such as a constant pressure
    condensed inside a UEL. They are not nodal unknowns, do not contribute
    to NDOFEL, and are typically stored in SVARS by a formulation-specific
    generator.

    Args:
        name: argument name used in material/equation methods
        storage: backing storage convention, currently informational
        condensed: whether the variable is eliminated locally
        initial: default initial value
    """
    def __init__(self, name, storage='SVARS', interpolation='constant',
                 condensed=True, initial=0.0):
        self.name = name
        self.storage = storage
        self.interpolation = interpolation
        self.condensed = condensed
        self.initial = initial
        self.rank = 0
        self.ndof_per_node = 0
        self.ndim = None
        self.is_local = True

    def __repr__(self):
        return (f"LocalScalar('{self.name}', storage='{self.storage}', "
                f"interpolation='{self.interpolation}', "
                f"condensed={self.condensed})")

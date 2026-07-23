"""
UINTER code generator.

Given a SurfaceInteraction instance with traction() [and optionally flux()],
generates a complete self-contained .for file that works as an Abaqus UINTER
subroutine.

The generated file contains:
  1. tensor_ops.for   - shared math utilities (real + complex)
  2. traction sub     - COMPLEX*16, translated from Python
  3. flux sub         - COMPLEX*16, translated from Python (if defined)
  4. CS tangent engine - computes DDSDDR [+ DDFDDT, DDSDDT, DDFDDR]
  5. UINTER wrapper   - standard Abaqus UINTER signature

State variables:
  Same pattern as umat_gen.py. Materials may define state_vars dict.
  STATEV read/write uses the standard STATEV TIME(2)==0 initialization pattern.
  State_old is real and never perturbed during CS calls.

Usage:
    model = LinearCohesive(K_n=1e6, K_t=1e5)
    generate_uinter(model, 'linear_cohesive.for')

    model = HealingCohesive(K_n=1e6, K_t=1e5, sigma_n=1e3, ...)
    generate_uinter(model, 'healing_cohesive.for')
"""

import inspect
import textwrap
import ast
import os
import numpy as np
from collections import OrderedDict

# Reuse AST translator and state-variable helpers from umat_gen
# In the package, these would be imported from a shared module.
# For standalone usage, we import directly from umat_gen.
try:
    from .umat_gen import (
        _state_var_info,
        _nstate_per_gp,
        _fortran_dims,
        _fortran_call,
        _wrap_fortran_source,
        FortranTranslator,
    )
except ImportError:
    # Standalone: try sibling import
    try:
        from umat_gen import (
            _state_var_info,
            _nstate_per_gp,
            _fortran_dims,
            _fortran_call,
            _wrap_fortran_source,
            FortranTranslator,
        )
    except ImportError:
        raise ImportError(
            "uinter_gen.py requires umat_gen.py for AST translation. "
            "Ensure umat_gen.py is in the same directory or package."
        )


# ===================================================================
# SurfaceInteraction base class
# ===================================================================

class SurfaceInteraction:
    """
    Base class for user-defined surface interaction models.

    Subclasses must define:
        props: dict mapping property names to default values
        traction(self, rdisp, ...): returns stress vector [+ state dict]

    Optionally define:
        state_vars: dict mapping state variable names to initial values
        flux(self, rdisp, temp, ...): returns [q_slave, q_master]
    """

    props = {}
    state_vars = {}

    def __init__(self, **kwargs):
        # Store property names in insertion order
        self.props_names = list(self.__class__.props.keys())

        # Set property values (default or user-supplied)
        for name, default in self.__class__.props.items():
            setattr(self, name, kwargs.get(name, default))

        # Detect methods
        self._methods = {}
        if hasattr(self, 'traction') and callable(self.traction):
            self._methods['traction'] = {
                'callable': self.traction,
                'field_params': self._extract_field_params('traction'),
            }
        if hasattr(self, 'flux') and callable(self.flux):
            self._methods['flux'] = {
                'callable': self.flux,
                'field_params': self._extract_field_params('flux'),
            }

        # Detect helper methods (private, start with _)
        self._helpers = {}
        for name in dir(self.__class__):
            if name.startswith('_') and not name.startswith('__'):
                method = getattr(self.__class__, name)
                if callable(method) and name != '_extract_field_params':
                    self._helpers[name] = method

    def _extract_field_params(self, method_name):
        """Extract non-self, non-state, non-dt parameter names."""
        method = getattr(self, method_name)
        sig = inspect.signature(method)
        params = []
        sv_names = set(self.__class__.state_vars.keys())
        for pname in sig.parameters:
            if pname == 'self':
                continue
            if pname == 'dt':
                continue
            # Check if it's a state_old parameter
            base = pname.replace('_old', '')
            if base in sv_names:
                continue
            params.append(pname)
        return params

    @property
    def has_state(self):
        return bool(self.__class__.state_vars)

    @property
    def has_flux(self):
        return 'flux' in self._methods


# ===================================================================
# Traction subroutine generation
# ===================================================================

def _generate_traction_subroutine(interaction, prefix):
    """
    Generate the COMPLEX*16 traction subroutine.

    Signature (stateless):
        SUBROUTINE prefix_traction(rdisp_z, ndir, props, stress_z)

    Signature (with state):
        SUBROUTINE prefix_traction(rdisp_z, ndir,
            d_old, H_old, dt, props,
            stress_z, d_new, H_new)

    Args follow the convention:
        1. Field inputs (COMPLEX):  rdisp_z
        2. Dimension:               ndir
        3. State old (COMPLEX):     d_old, H_old, ...
        4. Real parameters:         dt, props
        5. Output (COMPLEX):        stress_z
        6. State new (COMPLEX):     d_new, H_new, ...
    """
    sv_info = _state_var_info(interaction)
    has_state = bool(sv_info)
    nprops = len(interaction.props_names)

    # --- Build argument list ---
    args = ['rdisp_z', 'ndir']
    if has_state:
        for name in sv_info:
            args.append(f'{name}_old')
        args.append('dt')
    args.append('props')
    args.append('stress_z')
    if has_state:
        for name in sv_info:
            args.append(f'{name}_new')

    lines = []
    sig = f'      SUBROUTINE {prefix}_traction({", ".join(args)})'
    lines.append(sig)
    lines.append('      IMPLICIT NONE')

    # --- Argument declarations ---
    lines.append('      INTEGER, INTENT(IN) :: ndir')
    lines.append(
        '      DOUBLE COMPLEX, INTENT(IN) :: rdisp_z(ndir)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN) :: '
                f'{name}_old{dims}')
        lines.append(
            '      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.append(
        f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')
    lines.append(
        '      DOUBLE COMPLEX, INTENT(OUT) :: stress_z(ndir)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(OUT) :: '
                f'{name}_new{dims}')

    lines.append('')

    # --- AST translation of the traction method ---
    # This uses the FortranTranslator from umat_gen.py.
    # The translator needs to be extended to handle:
    #   - rdisp[i] → rdisp_z(i+1)  (0-indexed Python to 1-indexed Fortran)
    #   - Vector outputs instead of tensor outputs
    #
    # For Phase 1, we generate a placeholder that calls the
    # translator. The full integration with the AST translator
    # requires extending FortranTranslator to recognize 'rdisp'
    # as a vector field argument.

    lines.append('C     --- Translated from Python traction() ---')
    lines.append('C     TODO: AST translation integration')
    lines.append('C     The body below is a placeholder.')
    lines.append('C     Full implementation will use')
    lines.append('C     FortranTranslator from umat_gen.py')
    lines.append('C     extended for vector field arguments.')
    lines.append('')

    # Placeholder: zero traction
    lines.append('      stress_z(1) = DCMPLX(0.0d0, 0.0d0)')
    lines.append('      IF (ndir .GE. 2) THEN')
    lines.append(
        '        stress_z(2) = DCMPLX(0.0d0, 0.0d0)')
    lines.append('      END IF')
    lines.append('      IF (ndir .GE. 3) THEN')
    lines.append(
        '        stress_z(3) = DCMPLX(0.0d0, 0.0d0)')
    lines.append('      END IF')

    if has_state:
        for name in sv_info:
            lines.append(f'      {name}_new = {name}_old')

    lines.append('')
    lines.append('      RETURN')
    lines.append(f'      END SUBROUTINE {prefix}_traction')

    return '\n'.join(lines)


# ===================================================================
# Flux subroutine generation (optional)
# ===================================================================

def _generate_flux_subroutine(interaction, prefix):
    """
    Generate the COMPLEX*16 flux subroutine (if flux() is defined).

    Signature:
        SUBROUTINE prefix_flux(rdisp_z, temp_z, ndir,
            [state_old], dt, props, flux_z)
    """
    if not interaction.has_flux:
        return ''

    sv_info = _state_var_info(interaction)
    has_state = bool(sv_info)
    nprops = len(interaction.props_names)

    args = ['rdisp_z', 'temp_z', 'ndir']
    if has_state:
        for name in sv_info:
            args.append(f'{name}_old')
        args.append('dt')
    args.append('props')
    args.append('flux_z')

    lines = []
    sig = f'      SUBROUTINE {prefix}_flux({", ".join(args)})'
    lines.append(sig)
    lines.append('      IMPLICIT NONE')

    lines.append('      INTEGER, INTENT(IN) :: ndir')
    lines.append(
        '      DOUBLE COMPLEX, INTENT(IN) :: rdisp_z(ndir)')
    lines.append(
        '      DOUBLE COMPLEX, INTENT(IN) :: temp_z(2)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX, INTENT(IN) :: '
                f'{name}_old{dims}')
        lines.append(
            '      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.append(
        f'      DOUBLE PRECISION, INTENT(IN) :: props({nprops})')
    lines.append(
        '      DOUBLE COMPLEX, INTENT(OUT) :: flux_z(2)')

    lines.append('')
    lines.append('C     --- Translated from Python flux() ---')
    lines.append('C     TODO: AST translation integration')
    lines.append('')
    lines.append('      flux_z(1) = DCMPLX(0.0d0, 0.0d0)')
    lines.append('      flux_z(2) = DCMPLX(0.0d0, 0.0d0)')
    lines.append('')
    lines.append('      RETURN')
    lines.append(f'      END SUBROUTINE {prefix}_flux')

    return '\n'.join(lines)


# ===================================================================
# Complex-step tangent engine
# ===================================================================

def _generate_cs_tangent(prefix, nprops, sv_info=None,
                         has_flux=False):
    """
    Generate the complex-step tangent engine for UINTER.

    Computes up to 4 Jacobian blocks:
      Phase 1 (perturb RDISP): DDSDDR, DDFDDR (if flux)
      Phase 2 (perturb TEMP):  DDSDDT, DDFDDT (if flux)

    State_old is DOUBLE PRECISION, reset before each perturbation.
    """
    if sv_info is None:
        sv_info = OrderedDict()
    has_state = bool(sv_info)

    # --- Build argument list ---
    args = ['RDISP', 'NDIR']
    if has_flux:
        args.append('TEMP')
    if has_state:
        for name in sv_info:
            args.append(f'{name}_old')
        args.append('dt')
    args.extend(['props', 'DDSDDR'])
    if has_flux:
        args.extend(['DDFDDT', 'DDSDDT', 'DDFDDR'])

    lines = []
    sig = f'      SUBROUTINE {prefix}_cs_tangent({", ".join(args)})'
    lines.append(sig)
    lines.append('      IMPLICIT NONE')

    # --- Declarations ---
    lines.append(
        '      INTEGER, INTENT(IN) :: NDIR')
    lines.append(
        '      DOUBLE PRECISION, INTENT(IN) :: RDISP(NDIR)')
    if has_flux:
        lines.append(
            '      DOUBLE PRECISION, INTENT(IN) :: TEMP(2)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE PRECISION, INTENT(IN) :: '
                f'{name}_old{dims}')
        lines.append(
            '      DOUBLE PRECISION, INTENT(IN) :: dt')
    lines.append(
        f'      DOUBLE PRECISION, INTENT(IN) :: '
        f'props({nprops})')
    lines.append(
        '      DOUBLE PRECISION, INTENT(OUT) :: '
        'DDSDDR(NDIR,NDIR)')
    if has_flux:
        lines.append(
            '      DOUBLE PRECISION, INTENT(OUT) :: '
            'DDFDDT(2,2)')
        lines.append(
            '      DOUBLE PRECISION, INTENT(OUT) :: '
            'DDSDDT(NDIR,2)')
        lines.append(
            '      DOUBLE PRECISION, INTENT(OUT) :: '
            'DDFDDR(2,NDIR)')

    lines.append('')
    lines.append(
        '      DOUBLE PRECISION, PARAMETER :: CS_H = 1.0d-10')
    lines.append('      DOUBLE COMPLEX :: rdisp_z(3)')
    lines.append('      DOUBLE COMPLEX :: stress_z(3)')
    if has_flux:
        lines.append('      DOUBLE COMPLEX :: temp_z(2)')
        lines.append('      DOUBLE COMPLEX :: flux_z(2)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
        lines.append('      DOUBLE PRECISION :: dt_safe')
    lines.append('      INTEGER :: j, i')
    lines.append('')

    # --- dt safety ---
    if has_state:
        lines.append('      dt_safe = dt')
        lines.append(
            '      IF (dt_safe .LT. 1.0d-14) '
            'dt_safe = 1.0d-14')
        lines.append('')

    # --- Initialize outputs ---
    lines.append('      DDSDDR = 0.0d0')
    if has_flux:
        lines.append('      DDFDDT = 0.0d0')
        lines.append('      DDSDDT = 0.0d0')
        lines.append('      DDFDDR = 0.0d0')
    lines.append('')

    # ============================================================
    # Phase 1: Perturb RDISP components → DDSDDR [+ DDFDDR]
    # ============================================================
    lines.append(
        'C     --- Phase 1: perturb RDISP ---')
    lines.append('      DO j = 1, NDIR')
    lines.append('')

    # Reset rdisp_z to real
    lines.append(
        'C       Reset rdisp to real')
    lines.append('        DO i = 1, NDIR')
    lines.append(
        '          rdisp_z(i) = DCMPLX(RDISP(i), 0.0d0)')
    lines.append('        END DO')

    # Reset state_old to real
    if has_state:
        lines.append(
            'C       Reset state_old to real')
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'        {name}_old_z = '
                    f'DCMPLX({name}_old, 0.0d0)')
            elif shape == 'vector':
                lines.append(
                    f'        CALL real2complex3('
                    f'{name}_old, {name}_old_z)')
            elif shape == 'tensor':
                lines.append(
                    f'        CALL real2complex33('
                    f'{name}_old, {name}_old_z)')

    # Perturb component j
    lines.append('')
    lines.append(
        '        rdisp_z(j) = rdisp_z(j) '
        '+ DCMPLX(0.0d0, CS_H)')
    lines.append('')

    # Call traction
    call_args = ['rdisp_z', 'NDIR']
    if has_state:
        for name in sv_info:
            call_args.append(f'{name}_old_z')
        call_args.append('dt_safe')
    else:
        pass
    call_args.extend(['props', 'stress_z'])
    if has_state:
        for name in sv_info:
            call_args.append(f'{name}_new_z')

    lines.append(
        _fortran_call('        ', f'{prefix}_traction',
                      call_args))

    # Extract DDSDDR column
    lines.append(
        '        DO i = 1, NDIR')
    lines.append(
        '          DDSDDR(i,j) = AIMAG(stress_z(i)) / CS_H')
    lines.append('        END DO')

    # Extract DDFDDR column (if flux defined)
    if has_flux:
        lines.append('')
        lines.append(
            'C       Also evaluate flux for DDFDDR')

        # Reset temp_z to real (no perturbation in Phase 1)
        lines.append(
            '        temp_z(1) = DCMPLX(TEMP(1), 0.0d0)')
        lines.append(
            '        temp_z(2) = DCMPLX(TEMP(2), 0.0d0)')

        flux_args = ['rdisp_z', 'temp_z', 'NDIR']
        if has_state:
            for name in sv_info:
                flux_args.append(f'{name}_old_z')
            flux_args.append('dt_safe')
        flux_args.extend(['props', 'flux_z'])

        lines.append(
            _fortran_call('        ', f'{prefix}_flux',
                          flux_args))
        lines.append(
            '        DO i = 1, 2')
        lines.append(
            '          DDFDDR(i,j) = '
            'AIMAG(flux_z(i)) / CS_H')
        lines.append('        END DO')

    lines.append('')
    lines.append('      END DO')
    lines.append('')

    # ============================================================
    # Phase 2: Perturb TEMP → DDSDDT, DDFDDT (if flux defined)
    # ============================================================
    if has_flux:
        lines.append(
            'C     --- Phase 2: perturb TEMP ---')
        lines.append('      DO j = 1, 2')
        lines.append('')

        # Reset all to real
        lines.append(
            '        DO i = 1, NDIR')
        lines.append(
            '          rdisp_z(i) = '
            'DCMPLX(RDISP(i), 0.0d0)')
        lines.append('        END DO')
        lines.append(
            '        temp_z(1) = DCMPLX(TEMP(1), 0.0d0)')
        lines.append(
            '        temp_z(2) = DCMPLX(TEMP(2), 0.0d0)')

        if has_state:
            for name, entry in sv_info.items():
                shape = entry['shape']
                if shape == 'scalar':
                    lines.append(
                        f'        {name}_old_z = '
                        f'DCMPLX({name}_old, 0.0d0)')
                elif shape == 'vector':
                    lines.append(
                        f'        CALL real2complex3('
                        f'{name}_old, {name}_old_z)')
                elif shape == 'tensor':
                    lines.append(
                        f'        CALL real2complex33('
                        f'{name}_old, {name}_old_z)')

        # Perturb TEMP(j)
        lines.append('')
        lines.append(
            '        temp_z(j) = temp_z(j) '
            '+ DCMPLX(0.0d0, CS_H)')
        lines.append('')

        # Call traction (for DDSDDT)
        lines.append(
            _fortran_call('        ', f'{prefix}_traction',
                          call_args))
        lines.append(
            '        DO i = 1, NDIR')
        lines.append(
            '          DDSDDT(i,j) = '
            'AIMAG(stress_z(i)) / CS_H')
        lines.append('        END DO')

        # Call flux (for DDFDDT)
        lines.append(
            _fortran_call('        ', f'{prefix}_flux',
                          flux_args))
        lines.append(
            '        DO i = 1, 2')
        lines.append(
            'C         Note: DDFDDT = -d(FLUX)/d(TEMP)')
        lines.append(
            '          DDFDDT(i,j) = '
            '-AIMAG(flux_z(i)) / CS_H')
        lines.append('        END DO')

        lines.append('')
        lines.append('      END DO')
        lines.append('')

    lines.append('      RETURN')
    lines.append(
        f'      END SUBROUTINE {prefix}_cs_tangent')

    return '\n'.join(lines)


# ===================================================================
# UINTER wrapper generation
# ===================================================================

def _generate_uinter_wrapper(prefix, nprops, sv_info=None,
                              has_flux=False):
    """
    Generate the Abaqus UINTER wrapper subroutine.

    Handles:
      - STATEV read/write with TIME(2)==0 initialization
      - Call to CS tangent engine
      - Call to real traction (and flux) for actual STRESS update
      - LOPENCLOSE / LSDI flag management
      - PNEWDT safety
    """
    if sv_info is None:
        sv_info = OrderedDict()
    has_state = bool(sv_info)
    nstatv = _nstate_per_gp(sv_info) if has_state else 0

    lines = []

    # --- Subroutine signature ---
    lines.append('      SUBROUTINE UINTER(STRESS,DDSDDR,')
    lines.append(
        '     &   DVISCOUS,DSTRUCTURAL,FLUX,DDFDDT,')
    lines.append('     &   DDSDDT,DDFDDR,')
    lines.append(
        '     &   STATEV,SED,SFD,SPD,SVD,SCD,PNEWDT,')
    lines.append('     &   RDISP,DRDISP,')
    lines.append(
        '     &   TEMP,DTEMP,PREDEF,DPRED,')
    lines.append(
        '     &   TIME,DTIME,FREQR,CINAME,SLNAME,MSNAME,')
    lines.append(
        '     &   PROPS,COORDS,ALOCALDIR,DROT,AREA,')
    lines.append('     &   CHRLNGTH,')
    lines.append(
        '     &   NODE,NDIR,NSTATV,NPRED,NPROPS,MCRD,')
    lines.append(
        '     &   KSTEP,KINC,KIT,LINPER,LOPENCLOSE,')
    lines.append('     &   LSTATE,LSDI,LPRINT)')
    lines.append('')
    lines.append("      INCLUDE 'ABA_PARAM.INC'")
    lines.append('')
    lines.append('      CHARACTER*80 CINAME,SLNAME,MSNAME')
    lines.append('')
    lines.append(
        '      DIMENSION STRESS(NDIR),DDSDDR(NDIR,NDIR),')
    lines.append(
        '     &   FLUX(2),DDFDDT(2,2),')
    lines.append(
        '     &   DDSDDT(NDIR,2),DDFDDR(2,NDIR),')
    lines.append(
        '     &   STATEV(NSTATV),')
    lines.append(
        '     &   RDISP(NDIR),DRDISP(NDIR),')
    lines.append(
        '     &   TEMP(2),DTEMP(2),')
    lines.append(
        '     &   PREDEF(2,NPRED),DPRED(2,NPRED),')
    lines.append(
        '     &   TIME(2),PROPS(NPROPS),COORDS(MCRD),')
    lines.append(
        '     &   ALOCALDIR(3,3),DROT(2,2),')
    lines.append(
        '     &   DVISCOUS(NDIR,NDIR),')
    lines.append(
        '     &   DSTRUCTURAL(NDIR,NDIR)')
    lines.append('')

    # --- Local variables ---
    lines.append('C     --- Local variables ---')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE PRECISION :: {name}_old{dims}')
            lines.append(
                f'      DOUBLE PRECISION :: {name}_new{dims}')
    lines.append('      DOUBLE COMPLEX :: stress_z(3)')
    if has_flux:
        lines.append('      DOUBLE COMPLEX :: flux_z(2)')
    lines.append('      DOUBLE COMPLEX :: rdisp_z(3)')
    if has_flux:
        lines.append('      DOUBLE COMPLEX :: temp_z(2)')
    if has_state:
        for name, entry in sv_info.items():
            dims = _fortran_dims(entry)
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_old_z{dims}')
            lines.append(
                f'      DOUBLE COMPLEX :: {name}_new_z{dims}')
        lines.append('      DOUBLE PRECISION :: dt_safe')
    lines.append('      INTEGER :: i')
    lines.append('')

    # ============================================================
    # State variable read
    # ============================================================
    if has_state:
        lines.append(
            'C     ========================================')
        lines.append(
            'C     Read state variables')
        lines.append(
            'C     ========================================')
        lines.append(
            '      IF (TIME(2) .EQ. 0.0d0) THEN')

        # Initialize from defaults
        offset = 0
        for name, entry in sv_info.items():
            size = entry['size']
            init = entry['init']
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'        {name}_old = {init}d0')
            elif shape == 'vector':
                for k in range(3):
                    lines.append(
                        f'        {name}_old({k+1}) = '
                        f'{float(init[k])}d0')
            elif shape == 'tensor':
                for col in range(3):
                    for row in range(3):
                        lines.append(
                            f'        {name}_old({row+1},'
                            f'{col+1}) = '
                            f'{float(init[row, col])}d0')
            offset += size

        lines.append('      ELSE')

        # Read from STATEV
        offset = 0
        for name, entry in sv_info.items():
            size = entry['size']
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'        {name}_old = '
                    f'STATEV({offset + 1})')
            elif shape == 'vector':
                for k in range(3):
                    lines.append(
                        f'        {name}_old({k+1}) = '
                        f'STATEV({offset + k + 1})')
            elif shape == 'tensor':
                idx = 0
                for col in range(3):
                    for row in range(3):
                        lines.append(
                            f'        {name}_old({row+1},'
                            f'{col+1}) = '
                            f'STATEV({offset + idx + 1})')
                        idx += 1
            offset += size

        lines.append('      END IF')
        lines.append('')

    # ============================================================
    # dt safety
    # ============================================================
    if has_state:
        lines.append('      dt_safe = DTIME')
        lines.append(
            '      IF (dt_safe .LT. 1.0d-14) '
            'dt_safe = 1.0d-14')
        lines.append('')

    # ============================================================
    # Call CS tangent engine
    # ============================================================
    lines.append(
        'C     ========================================')
    lines.append(
        'C     Compute Jacobian via complex-step')
    lines.append(
        'C     ========================================')

    cs_args = ['RDISP', 'NDIR']
    if has_flux:
        cs_args.append('TEMP')
    if has_state:
        for name in sv_info:
            cs_args.append(f'{name}_old')
        cs_args.append('dt_safe')
    cs_args.extend(['PROPS', 'DDSDDR'])
    if has_flux:
        cs_args.extend(['DDFDDT', 'DDSDDT', 'DDFDDR'])

    lines.append(
        _fortran_call('      ', f'{prefix}_cs_tangent',
                      cs_args))
    lines.append('')

    # ============================================================
    # Call real traction for STRESS and state update
    # ============================================================
    lines.append(
        'C     ========================================')
    lines.append(
        'C     Compute real traction + state update')
    lines.append(
        'C     ========================================')

    # Convert rdisp to complex with zero imaginary
    lines.append('      DO i = 1, NDIR')
    lines.append(
        '        rdisp_z(i) = DCMPLX(RDISP(i), 0.0d0)')
    lines.append('      END DO')

    if has_state:
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'      {name}_old_z = '
                    f'DCMPLX({name}_old, 0.0d0)')
            elif shape == 'vector':
                lines.append(
                    f'      CALL real2complex3('
                    f'{name}_old, {name}_old_z)')
            elif shape == 'tensor':
                lines.append(
                    f'      CALL real2complex33('
                    f'{name}_old, {name}_old_z)')

    # Call traction with real inputs (complex with Im=0)
    real_args = ['rdisp_z', 'NDIR']
    if has_state:
        for name in sv_info:
            real_args.append(f'{name}_old_z')
        real_args.append('dt_safe')
    real_args.extend(['PROPS', 'stress_z'])
    if has_state:
        for name in sv_info:
            real_args.append(f'{name}_new_z')

    lines.append(
        _fortran_call('      ', f'{prefix}_traction',
                      real_args))
    lines.append('')

    # Extract real part of stress
    lines.append('      DO i = 1, NDIR')
    lines.append(
        '        STRESS(i) = DBLE(stress_z(i))')
    lines.append('      END DO')
    lines.append('')

    # ============================================================
    # Call real flux (if defined)
    # ============================================================
    if has_flux:
        lines.append(
            'C     --- Compute real flux ---')
        lines.append(
            '      temp_z(1) = DCMPLX(TEMP(1), 0.0d0)')
        lines.append(
            '      temp_z(2) = DCMPLX(TEMP(2), 0.0d0)')

        flux_args = ['rdisp_z', 'temp_z', 'NDIR']
        if has_state:
            for name in sv_info:
                flux_args.append(f'{name}_old_z')
            flux_args.append('dt_safe')
        flux_args.extend(['PROPS', 'flux_z'])

        lines.append(
            _fortran_call('      ', f'{prefix}_flux',
                          flux_args))
        lines.append('      FLUX(1) = DBLE(flux_z(1))')
        lines.append('      FLUX(2) = DBLE(flux_z(2))')
        lines.append('')

    # ============================================================
    # State variable write
    # ============================================================
    if has_state:
        lines.append(
            'C     ========================================')
        lines.append(
            'C     Write state variables')
        lines.append(
            'C     ========================================')

        # Extract real parts of state_new
        for name, entry in sv_info.items():
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'      {name}_new = DBLE({name}_new_z)')
            elif shape == 'vector':
                for k in range(3):
                    lines.append(
                        f'      {name}_new({k+1}) = '
                        f'DBLE({name}_new_z({k+1}))')
            elif shape == 'tensor':
                for col in range(3):
                    for row in range(3):
                        lines.append(
                            f'      {name}_new({row+1},'
                            f'{col+1}) = '
                            f'DBLE({name}_new_z({row+1},'
                            f'{col+1}))')

        lines.append('')

        # Write to STATEV
        offset = 0
        for name, entry in sv_info.items():
            size = entry['size']
            shape = entry['shape']
            if shape == 'scalar':
                lines.append(
                    f'      STATEV({offset + 1}) = '
                    f'{name}_new')
            elif shape == 'vector':
                for k in range(3):
                    lines.append(
                        f'      STATEV({offset + k + 1}) = '
                        f'{name}_new({k+1})')
            elif shape == 'tensor':
                idx = 0
                for col in range(3):
                    for row in range(3):
                        lines.append(
                            f'      STATEV({offset + idx + 1}'
                            f') = {name}_new({row+1},'
                            f'{col+1})')
                        idx += 1
            offset += size

        lines.append('')

    # ============================================================
    # Contact status flags
    # ============================================================
    lines.append(
        'C     ========================================')
    lines.append(
        'C     Contact status management')
    lines.append(
        'C     ========================================')
    lines.append(
        'C     Default: set LOPENCLOSE based on normal gap')
    lines.append(
        'C     RDISP(1) > 0 = penetration (closed)')
    lines.append(
        'C     RDISP(1) < 0 = open gap')
    lines.append(
        '      IF (RDISP(1) .LT. -CHRLNGTH*1.0d-3) THEN')
    lines.append(
        'C       Open')
    lines.append(
        '        IF (LOPENCLOSE .EQ. 1) LSDI = 1')
    lines.append(
        '        LOPENCLOSE = 0')
    lines.append('      ELSE')
    lines.append(
        'C       Closed')
    lines.append(
        '        IF (LOPENCLOSE .EQ. 0) LSDI = 1')
    lines.append(
        '        LOPENCLOSE = 1')
    lines.append('      END IF')
    lines.append('')

    lines.append('      RETURN')
    lines.append('      END SUBROUTINE UINTER')

    return '\n'.join(lines)


# ===================================================================
# Main entry point
# ===================================================================

def generate_uinter(interaction, output_path, prefix=None):
    """
    Generate a complete .for file for Abaqus UINTER.

    Args:
        interaction: a SurfaceInteraction instance with traction()
        output_path: path to write the .for file
        prefix: prefix for generated subroutine names
                (default: lowercase class name)
    """
    if 'traction' not in interaction._methods:
        raise ValueError(
            "SurfaceInteraction must define traction()")

    if prefix is None:
        prefix = interaction.__class__.__name__.lower()

    nprops = len(interaction.props_names)
    sv_info = _state_var_info(interaction)
    has_state = bool(sv_info)
    has_flux = interaction.has_flux

    # Read tensor_ops.for template
    template_dir = os.path.join(
        os.path.dirname(__file__), 'templates')
    tensor_ops_path = os.path.join(
        template_dir, 'tensor_ops.for')
    if not os.path.exists(tensor_ops_path):
        tensor_ops_path = os.path.join(
            os.path.dirname(__file__), '..', '..', 'tensor_ops.for')
    if not os.path.exists(tensor_ops_path):
        # Fallback for standalone usage
        tensor_ops_path = 'tensor_ops.for'

    tensor_ops_src = ''
    if os.path.exists(tensor_ops_path):
        with open(tensor_ops_path, 'r') as f:
            tensor_ops_src = f.read()

    # Generate sections
    traction_src = _generate_traction_subroutine(
        interaction, prefix)
    flux_src = _generate_flux_subroutine(
        interaction, prefix) if has_flux else ''
    cs_src = _generate_cs_tangent(
        prefix, nprops, sv_info, has_flux)
    uinter_src = _generate_uinter_wrapper(
        prefix, nprops, sv_info, has_flux)

    # Header
    nstatv = _nstate_per_gp(sv_info) if has_state else 0
    hdr = []
    hdr.append(f'C{"=" * 70}')
    hdr.append(f'C     {os.path.basename(output_path)}')
    hdr.append('C')
    hdr.append('C     Generated by abaqus_ufl uinter_gen v0.1.0')
    hdr.append(
        f'C     Interaction: {interaction.__class__.__name__}')
    hdr.append(
        f'C     Props: '
        f'{", ".join(f"{k}={getattr(interaction, k)}" for k in interaction.props_names)}')
    hdr.append(f'C     NPROPS = {nprops}')
    if has_state:
        hdr.append(f'C     NSTATV = {nstatv}')
        hdr.append(
            f'C     State vars: '
            f'{", ".join(sv_info.keys())}')
    if has_flux:
        hdr.append('C     Thermal coupling: YES')
    else:
        hdr.append('C     Thermal coupling: NO')
    hdr.append('C')
    hdr.append(
        'C     *SURFACE INTERACTION, USER, '
        f'PROPERTIES={nprops}')
    if has_state:
        hdr.append(f'C     DEPVAR={nstatv}')
    hdr.append(
        'C     This file is self-contained.')
    hdr.append(f'C{"=" * 70}')
    header = '\n'.join(hdr)

    # Assemble
    sections = [
        header,
        '',
        'C' + '=' * 70,
        'C     Section 1: Abaqus UINTER interface',
        'C' + '=' * 70,
        uinter_src,
        '',
        'C' + '=' * 70,
        f'C     Section 2: Traction: {prefix}_traction',
        'C' + '=' * 70,
        traction_src,
    ]

    if has_flux:
        sections.extend([
            '',
            'C' + '=' * 70,
            f'C     Section 3: Flux: {prefix}_flux',
            'C' + '=' * 70,
            flux_src,
        ])

    sections.extend([
        '',
        'C' + '=' * 70,
        'C     Section 4: Complex-step tangent engine',
        'C' + '=' * 70,
        cs_src,
    ])

    if tensor_ops_src:
        sections.extend([
            '',
            'C' + '=' * 70,
            'C     Section 5: Tensor operations',
            'C' + '=' * 70,
            tensor_ops_src,
        ])

    full_src = '\n'.join(sections)

    # Apply column-72 wrapping
    full_src = _wrap_fortran_source(full_src)

    with open(output_path, 'w') as f:
        f.write(full_src)

    subs = [f'{prefix}_traction',
            f'{prefix}_cs_tangent', 'UINTER']
    if has_flux:
        subs.insert(1, f'{prefix}_flux')

    print(f"  Generated: {output_path}")
    print(
        f"  Interaction: {interaction.__class__.__name__}")
    props_str = ', '.join(interaction.props_names)
    print(f"  Props:     {nprops} ({props_str})")
    if has_state:
        sv_str = ', '.join(
            f'{n}({v["size"]})' for n, v in sv_info.items())
        print(f"  State:     {nstatv} ({sv_str})")
    print(f"  Flux:      {'YES' if has_flux else 'NO'}")
    print(f"  Subroutines: {', '.join(subs)}")

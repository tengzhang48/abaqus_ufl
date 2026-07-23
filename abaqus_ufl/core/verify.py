"""
Verification suite for Material instances.

Called by model.verify(). Runs:
  1. Stress-free reference state: P(F=I, p=0) ≈ 0
  2. CS vs FD for every tangent block
  3. Major symmetry of dP/dF: A_{iJkL} = A_{kLiJ}

Handles materials with or without state_vars. For materials with
state_vars, the CS perturbation converts state_old to complex with
zero imaginary part (state is real, never perturbed). The FD
perturbation uses real state directly.

CS_H = 1e-10: validated across all tangent blocks in prior work.
"""

import numpy as np

from ._cs_state import complexify_state as _complexify_state

CS_H = 1e-10


class VerificationError(AssertionError):
    """Raised when Material.verify() detects a failed check."""
    pass


class OperatorSignWarning(UserWarning):
    """A storage/flux pair implies an INDEFINITE single-field operator.

    The generator assembles residual r = storage*eta - flux.grad(eta), so
    the single-field AMATRX block is (ds/df)*M - grad(N).Dflux.grad(N).
    If ds/df and the eigenvalues of -d(flux)/d(grad f) have mixed signs,
    the gradient term is anti-diffusive (the flipped phase_flux trap,
    2026-06-12): solves diverge exactly when the field localizes, while
    verify()'s CS-vs-FD checks stay green.
    """
    pass


# Storage/flux method pairs as assembled by EQUATION_INFO. Used by the
# operator-sign screening below.
_STORAGE_FLUX_PAIRS = [
    ('phase_storage', 'phase_flux'),
    ('species_storage', 'species_flux'),
    ('solvent_storage', 'solvent_flux'),
    ('pressure_storage', 'pressure_flux'),
]


def _quad4_operator_matrices(h):
    """Q4 mass M(4,4) and Laplacian L(4,4) on a square of side h."""
    import itertools
    g = 1.0 / np.sqrt(3.0)
    M = np.zeros((4, 4))
    L = np.zeros((4, 4))
    for xi, eta in itertools.product((-g, g), (-g, g)):
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dN = 0.25 * np.array([[-(1 - eta), -(1 - xi)],
                              [(1 - eta), -(1 + xi)],
                              [(1 + eta), (1 + xi)],
                              [-(1 + eta), (1 - xi)]]) * (2.0 / h)
        detJ = (h / 2.0) ** 2
        M += np.outer(N, N) * detJ
        L += dN @ dN.T * detJ
    return M, L


def _check_operator_signs(material, state, verbose=True):
    """Screen storage/flux pairs for the anti-diffusive sign pairing.

    Assembles the single-field AMATRX block K = (ds/df)*M - Df*L (the
    package convention r = storage*eta - flux.grad(eta)) at the resolving
    scale h = l_eff/10, l_eff = sqrt(|Df|/|ds/df|), and warns if K is
    INDEFINITE — the signature of an anti-diffusive gradient term (the
    flipped-flux trap, 2026-06-12/13). This is the operator-level test
    (same logic as tests/test_phase_flux_operator.py), assembled rather
    than inferred from a sign heuristic: the earlier heuristic version
    used a shared rate-inflated tolerance and a dt->0 evaluation that
    HID pure-diffusion (rate + gradient, no comparable reaction
    stiffness) sign errors such as those in a coupled two-field
    diffusion equation (2026-06-13). Block assembly at l_eff/10
    self-calibrates the two
    terms, so a single relative tolerance suffices.

    Warns, never fails. A consistently-negated equation (storage
    AND flux both negated) gives an all-negative block (NOT indefinite)
    and is correctly silent; a non-convex reaction at a spinodal with
    the correct diffusion sign likewise gives a definite block.

    Pairs whose flux has zero or multiple grad arguments are skipped, as
    are pairs whose flux gradient-field is absent from the storage
    signature (the phase half of a mixed CH split). KNOWN LIMITATION:
    the mu-definition half of a mixed split (storage depends on its own
    primal field, flux = +kappa*grad(primal), which is CORRECT there)
    can false-warn if the verify state sits at a spinodal (psi'' < 0);
    signatures cannot distinguish that from a degradation-coupled
    storage. Assess mixed pairs with the coupled dispersion pattern
    (test_phase_flux_operator.py, mixed-split bulk-phase-stability test)
    and dismiss if the dispersion is stable.
    """
    import warnings as _warnings

    for sm, fm in _STORAGE_FLUX_PAIRS:
        if sm not in material._methods or fm not in material._methods:
            continue
        fparams = material._methods[fm]['all_params']
        grads = [p for p in fparams if p.startswith('grad_')]
        if len(grads) != 1:
            continue
        field = grads[0][5:]
        sparams = material._methods[sm]['all_params']
        if field not in sparams or field not in state \
                or grads[0] not in state:
            continue

        try:
            s_val = float(np.real(
                _cs_tangent(material, sm, field, state)))
            Df = np.real(_cs_tangent(material, fm, grads[0], state))
        except Exception:
            continue   # screening only — never block verify() itself

        # dominant gradient-diffusivity eigenvalue (2D in-plane block)
        lam_Df = np.linalg.eigvalsh(0.5 * (Df + Df.T)[:2, :2])
        d_dom = lam_Df[np.argmax(np.abs(lam_Df))]
        if abs(d_dom) < 1e-300 or abs(s_val) < 1e-300:
            if verbose:
                print(f"  Operator sign ({sm}/{fm}, '{field}'): "
                      f"degenerate (s={s_val:.2e}, Df={d_dom:.2e}), skip")
            continue

        l_eff = np.sqrt(abs(d_dom) / abs(s_val))
        M, L = _quad4_operator_matrices(l_eff / 10.0)
        K = s_val * M - d_dom * L
        w = np.linalg.eigvalsh(0.5 * (K + K.T))
        scale = float(np.max(np.abs(w)))
        indefinite = (w.max() > 1e-9 * scale) and (w.min() < -1e-9 * scale)
        if indefinite:
            _warnings.warn(
                f"{material.__class__.__name__}: storage/flux pair "
                f"({sm}, {fm}) assembles an INDEFINITE single-field "
                f"operator for field '{field}' at h=l_eff/10 "
                f"(d{sm}/d{field} = {s_val:.3e}, dominant "
                f"d{fm}/d{grads[0]} = {d_dom:.3e}; same sign => "
                f"anti-diffusive). Likely a flipped {fm} sign. The "
                f"assembled convention is r = storage*eta - "
                f"flux.grad(eta).", OperatorSignWarning, stacklevel=3)
        elif verbose:
            print(f"  Operator sign ({sm}/{fm}, field '{field}'): OK")

# FIELD_ARGS is imported from material.py (or defs.py).
# To avoid circular imports, we accept it from the material itself.


# ===================================================================
# Core tangent computation
# ===================================================================


def _cs_tangent(material, method_name, field_var, state, h=CS_H):
    """
    Compute tangent d(output)/d(field_var) via complex step.

    Args:
        material:   Material instance
        method_name: e.g. 'stress_PK1'
        field_var:   e.g. 'F', 'p', 'mu', 'grad_mu'
        state:       dict with all method arguments (real)
        h:           CS perturbation size

    Returns:
        tangent array with shape = output_shape + input_shape
    """
    func = material._call_method_as_func(method_name)
    info = material._methods[method_name]
    out_shape = info['output_shape']

    val = state[field_var]

    if np.isscalar(val) or (isinstance(val, np.ndarray) and val.ndim == 0):
        # Scalar perturbation: one evaluation
        cs_state = _complexify_state(state)
        cs_state[field_var] = complex(float(val), h)
        result = func(**cs_state)
        return np.asarray(result).imag / h

    else:
        # Array perturbation: loop over components
        val = np.asarray(val)
        in_shape = val.shape
        tangent = np.zeros(out_shape + in_shape)

        for idx in np.ndindex(in_shape):
            cs_state = _complexify_state(state)
            perturbed = np.asarray(state[field_var], dtype=complex).copy()
            perturbed[idx] += 1j * h
            cs_state[field_var] = perturbed
            result = func(**cs_state)
            result = np.asarray(result)

            if out_shape == ():
                tangent[idx] = result.imag / h
            else:
                tangent_slice = tuple(
                    [slice(None)] * len(out_shape) + list(idx))
                tangent[tangent_slice] = result.imag / h

        return tangent


def _fd_tangent(material, method_name, field_var, state, eps=1e-7):
    """
    Compute tangent d(output)/d(field_var) via central finite differences.

    Same interface as _cs_tangent but uses real perturbations.
    """
    func = material._call_method_as_func(method_name)
    info = material._methods[method_name]
    out_shape = info['output_shape']

    val = state[field_var]

    if np.isscalar(val) or (isinstance(val, np.ndarray) and val.ndim == 0):
        s_p = dict(state); s_p[field_var] = float(val) + eps
        s_m = dict(state); s_m[field_var] = float(val) - eps
        r_p = np.asarray(func(**s_p))
        r_m = np.asarray(func(**s_m))
        return (r_p - r_m) / (2.0 * eps)

    else:
        val = np.asarray(val)
        in_shape = val.shape
        tangent = np.zeros(out_shape + in_shape)

        for idx in np.ndindex(in_shape):
            s_p = dict(state)
            s_m = dict(state)
            v_p = val.copy().astype(float)
            v_m = val.copy().astype(float)
            v_p[idx] += eps
            v_m[idx] -= eps
            s_p[field_var] = v_p
            s_m[field_var] = v_m
            r_p = np.asarray(func(**s_p))
            r_m = np.asarray(func(**s_m))

            tangent_slice = tuple(
                [slice(None)] * len(out_shape) + list(idx))
            tangent[tangent_slice] = (r_p - r_m) / (2.0 * eps)

        return tangent


# ===================================================================
# Verification checks
# ===================================================================

def _check_stress_free(material, verbose=True):
    """Check P = 0 at F = I (stress-free reference)."""
    if 'stress_PK1' not in material._methods:
        return True

    state = material._make_state(F=np.eye(3))
    # Set pressure and chemical potential to zero if present
    info = material._methods['stress_PK1']
    if 'p' in info['all_params']:
        state['p'] = 0.0
    if 'mu' in info['all_params']:
        state['mu'] = 0.0

    P = material._call_method('stress_PK1', state)
    norm = np.max(np.abs(P))

    ok = norm < 1e-10
    if verbose:
        status = "PASS" if ok else "FAIL"
        print(f"  Stress-free reference: ||P(F=I)|| = {norm:.2e} "
              f"[{status}]")
    return ok


def _check_tangent_block(material, method_name, field_var, state,
                         tol=1e-6, atol=1e-12, verbose=True):
    """Check CS vs FD for one tangent block."""
    cs = _cs_tangent(material, method_name, field_var, state)
    fd = _fd_tangent(material, method_name, field_var, state)

    err = np.max(np.abs(cs - fd))
    norm = np.max(np.abs(cs))

    if norm > atol:
        rel_err = err / norm
        ok = rel_err < tol
    else:
        # Near-zero tangent: use absolute tolerance
        rel_err = err
        ok = err < atol

    if verbose:
        status = "PASS" if ok else "FAIL"
        print(f"  d({method_name})/d({field_var}): "
              f"err = {err:.2e}, norm = {norm:.2e} [{status}]")
    return ok


def _check_major_symmetry(material, state, tol=1e-6, verbose=True):
    """Check A_{iJkL} = A_{kLiJ} for dP/dF tangent."""
    if 'stress_PK1' not in material._methods:
        return True

    info = material._methods['stress_PK1']
    if 'F' not in info['field_params']:
        return True

    dPdF = _cs_tangent(material, 'stress_PK1', 'F', state)

    max_asym = 0.0
    for i in range(3):
        for J in range(3):
            for k in range(3):
                for L in range(3):
                    d = abs(dPdF[i, J, k, L] - dPdF[k, L, i, J])
                    if d > max_asym:
                        max_asym = d

    norm = np.max(np.abs(dPdF)) + 1e-30
    rel_asym = max_asym / norm
    ok = rel_asym < tol

    if verbose:
        status = "PASS" if ok else "FAIL"
        print(f"  Major symmetry dP/dF: "
              f"max asym/norm = {rel_asym:.2e} [{status}]")
    return ok


# ===================================================================
# Main entry point
# ===================================================================

def run_verification(material, state=None, tol=1e-6, verbose=True):
    """
    Run all verification checks for a Material instance.

    Checks:
      1. Stress-free reference (P=0 at F=I)
      2. CS vs FD for every tangent block
      3. Major symmetry of dP/dF

    Args:
        material: Material instance
        state: optional state dict (default: material._make_state())
        tol: relative tolerance for CS vs FD comparison
        verbose: print results

    Returns:
        True if all checks pass
    """
    if state is None:
        state = material._make_state()

    all_pass = True

    if verbose:
        name = material.__class__.__name__
        print(f"Verifying {name}:")

    # 1. Stress-free reference
    ok = _check_stress_free(material, verbose)
    all_pass = all_pass and ok

    # 2. CS vs FD for each tangent block
    blocks = material.tangent_blocks()
    for (method_name, field_var), binfo in blocks.items():
        ok = _check_tangent_block(
            material, method_name, field_var, state, tol,
            verbose=verbose)
        all_pass = all_pass and ok

    # 2b. Operator-sign screening for storage/flux pairs (warns only)
    _check_operator_signs(material, state, verbose)

    # 3. Major symmetry — only for stateless materials that expect it.
    #    Plasticity algorithmic tangents are generally non-symmetric, and a
    #    material may declare ``symmetric_tangent = False`` when its stress is
    #    deliberately not an energy derivative (e.g. the stabilized u-theta
    #    formulation of Scovazzi et al. 2023, Eq. 23).
    if not material.has_state and getattr(material, "symmetric_tangent", True):
        ok = _check_major_symmetry(material, state, tol, verbose)
        all_pass = all_pass and ok
    elif verbose:
        reason = ("non-symmetric for materials with state variables"
                  if material.has_state
                  else "material declares symmetric_tangent = False")
        print(f"  Major symmetry dP/dF: skipped ({reason})")

    if verbose:
        if all_pass:
            print(f"  ALL CHECKS PASSED")
        else:
            print(f"  SOME CHECKS FAILED")

    if not all_pass:
        failed = []
        # Re-run checks silently to capture which ones failed
        ok = _check_stress_free(material, verbose=False)
        if not ok:
            failed.append("stress-free reference")
        blocks = material.tangent_blocks()
        for (method_name, field_var), binfo in blocks.items():
            ok = _check_tangent_block(
                material, method_name, field_var, state, tol,
                verbose=False)
            if not ok:
                failed.append(f"d({method_name})/d({field_var})")
        if not material.has_state:
            ok = _check_major_symmetry(material, state, tol, verbose=False)
            if not ok:
                failed.append("major symmetry dP/dF")
        report = (
            f"Verification failed for {material.__class__.__name__}. "
            f"Failed checks ({len(failed)}): " + ", ".join(failed)
        )
        raise VerificationError(report)

    return True

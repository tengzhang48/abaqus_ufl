"""Decorator for hand-written Fortran helper sidecars.

Materials may mark a method with ``@au.fortran_helper(...)``.  The Python
body still runs during ``verify()`` and Python reference paths, but the
Fortran generator emits a CALL to the named subroutine instead of
translating the method body.
"""


def fortran_helper(*, inputs, outputs, subroutine=None, pass_props=False):
    """Mark a Python method as a stub for a hand-written Fortran helper.

    Parameters
    ----------
    inputs : list[tuple[str, str]]
        [(name, kind), ...] where kind describes the Fortran shape,
        e.g. ``'scalar'``, ``'vector(3)'``, ``'tensor(3,3)'``.
    outputs : list[tuple[str, str]]
        Same shape vocabulary.  **Multiple outputs are supported** and are
        emitted as additional trailing ``INTENT(OUT)`` arguments.  The
        Python body must ``return`` a tuple with matching arity/shape.
    subroutine : str, optional
        Fortran subroutine name.  Defaults to ``<material_class_lower>_<method>``.
    pass_props : bool, optional
        If ``True``, the generator appends ``props`` to the emitted
        ``CALL`` argument list (before the output arguments).  This lets
        the sidecar subroutine access material properties via the standard
        ``DOUBLE PRECISION, INTENT(IN) :: props(NPROPS)`` array.

    Examples
    --------
    Single-output helper::

        class MyModel(au.SmallStrainMaterial):
            @au.fortran_helper(
                inputs=[('s', 'tensor(3,3)')],
                outputs=[('q', 'scalar')],
                subroutine='my_q_mises')
            def _q_mises(self, s):
                # Python body runs in verify()
                return np.sqrt(1.5) * np.linalg.norm(s)

    Multi-output helper::

        class MyModel(au.SmallStrainMaterial):
            @au.fortran_helper(
                inputs=[('sigma', 'tensor(3,3)')],
                outputs=[('s', 'tensor(3,3)'), ('p', 'scalar')])
            def _split_dev_vol(self, sigma):
                p = np.trace(sigma) / 3.0
                return sigma - p * np.eye(3), p
    """
    def deco(func):
        func.__fortran_helper__ = {
            'inputs': inputs,
            'outputs': outputs,
            'subroutine_override': subroutine,
            'pass_props': pass_props,
        }
        return func
    return deco

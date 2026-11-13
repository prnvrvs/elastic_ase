"""Calls into spglib that behave the same on spglib 2.4 to 3.x.

Silences the deprecation warning ("Set OLD_ERROR_HANDLING to false ...")
on spglib >= 2.7 and handles failures uniformly across versions.
"""
import warnings


def call(name, *args, **kwargs):
    """Return ``spglib.<name>(*args, **kwargs)``, or ``None`` if spglib fails."""
    import spglib

    error = getattr(getattr(spglib, "error", None), "SpglibError", None) or getattr(spglib, "SpglibError", None)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*OLD_ERROR_HANDLING.*", category=DeprecationWarning)
        try:
            return getattr(spglib, name)(*args, **kwargs)
        except Exception as exc:
            if error is not None and isinstance(exc, error):
                return None
            raise

"""Compatibility shim for legacy imports.

Older code/tests import `xespresso.gui.pages.codes_config`.
This shim re-exports the real implementation from `qtgui.pages.codes_config`.
"""
try:
    # Re-export everything from the qtgui implementation
    from qtgui.pages.codes_config import *  # noqa: F401,F403
except Exception:  # pragma: no cover - best-effort shim for headless/test envs
    # Minimal fallback so imports don't fail in very constrained environments
    def _codes_config_stub():
        """Fallback stub for codes_config module."""
        return None

    __all__ = ["_codes_config_stub"]

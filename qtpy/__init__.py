"""Lightweight compatibility shim for projects that import `qtpy`.

If the real `qtpy` package is available it will be used. Otherwise this
module exposes the `QtWidgets`, `QtCore` and `QtGui` submodules from
`PySide6` and registers them in `sys.modules` so imports like
`from qtpy.QtWidgets import QWidget` continue to work.
"""
from __future__ import annotations

import sys

try:
    # Prefer the real qtpy if available
    import qtpy as _real_qtpy  # type: ignore
    # Re-export standard submodules from the real package
    from qtpy import QtWidgets, QtCore, QtGui  # type: ignore
except Exception:
    # Fallback to PySide6 when qtpy is not installed
    try:
        from PySide6 import QtWidgets, QtCore, QtGui  # type: ignore
    except Exception:
        # Let import error propagate to caller with original traceback
        raise

    # Make submodules importable as `qtpy.QtWidgets`, etc.
    sys.modules.setdefault(__name__, sys.modules.get(__name__))
    sys.modules.setdefault(f"{__name__}.QtWidgets", QtWidgets)
    sys.modules.setdefault(f"{__name__}.QtCore", QtCore)
    sys.modules.setdefault(f"{__name__}.QtGui", QtGui)

    # Expose names on the package namespace for `from qtpy import QtWidgets` style
    globals()["QtWidgets"] = QtWidgets
    globals()["QtCore"] = QtCore
    globals()["QtGui"] = QtGui

__all__ = ["QtWidgets", "QtCore", "QtGui"]

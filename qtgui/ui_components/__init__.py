"""UI components package.

Re-exports helper functions for sidebar, menu, toolbar and statusbar so callers
can continue to import from `qtgui.ui_components`.
"""
from .sidebar import create_sidebar
from .menu import setup_menu
from .toolbar import setup_toolbar
from .statusbar import setup_statusbar

__all__ = [
    "create_sidebar",
    "setup_menu",
    "setup_toolbar",
    "setup_statusbar",
]

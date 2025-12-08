"""
Utility functions for the xespresso PyQt GUI.

This module provides common utilities used across the GUI pages.
"""

from .validation import validate_path, validate_path_under_base, safe_path_exists, safe_makedirs

# ASE_ESPRESSO_COMMAND template for Quantum ESPRESSO execution
# LAUNCHER, PACKAGE, PARALLEL, and PREFIX are placeholders replaced by xespresso
ASE_ESPRESSO_COMMAND_TEMPLATE = "LAUNCHER PACKAGE.x PARALLEL -in PREFIX.PACKAGEi > PREFIX.PACKAGEo"

# Optional imports - visualization and database modules may have optional dependencies
try:
    from .visualization import create_structure_figure, get_structure_info_text, MATPLOTLIB_AVAILABLE
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    create_structure_figure = None
    get_structure_info_text = None

try:
    from .database import DatabaseSaveWidget, DatabaseLoadWidget, ASE_DB_AVAILABLE
except ImportError:
    ASE_DB_AVAILABLE = False
    DatabaseSaveWidget = None
    DatabaseLoadWidget = None

__all__ = [
    'validate_path',
    'validate_path_under_base',
    'safe_path_exists',
    'safe_makedirs',
    'ASE_ESPRESSO_COMMAND_TEMPLATE',
    'create_structure_figure',
    'get_structure_info_text',
    'MATPLOTLIB_AVAILABLE',
    'DatabaseSaveWidget',
    'DatabaseLoadWidget',
    'ASE_DB_AVAILABLE',
]

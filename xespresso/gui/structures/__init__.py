"""
Structures module for xespresso GUI.

This module provides structure handling functionality for the GUI,
including loading, exporting, and managing atomic structures.

The key principle: structure modules handle loading/export operations,
while GUI pages coordinate the user interaction.
"""

from xespresso.gui.structures.base import BaseStructureHandler
from xespresso.gui.structures.loader import (
    StructureLoader,
    load_structure_from_file,
    load_structure_from_upload
)
from xespresso.gui.structures.exporter import (
    StructureExporter,
    export_structure
)

__all__ = [
    'BaseStructureHandler',
    'StructureLoader',
    'StructureExporter',
    'load_structure_from_file',
    'load_structure_from_upload',
    'export_structure',
]

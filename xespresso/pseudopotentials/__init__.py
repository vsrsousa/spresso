"""
pseudopotentials module

This module provides utilities for managing pseudopotential configurations
on different machines. It allows users to:
- Define available pseudopotentials and their paths
- Store pseudopotential information for different elements
- Configure pseudopotential settings per machine
- Detect available pseudopotentials on a machine
- Support multiple pseudopotential libraries on the same machine
- List available pseudopotentials for calculations
"""

from .config import Pseudopotential, PseudopotentialsConfig
from .manager import (
    PseudopotentialsManager,
    create_pseudopotentials_config,
    load_pseudopotentials_config,
    DEFAULT_PSEUDOPOTENTIALS_DIR
)
from .detector import (
    detect_pseudopotentials,
    detect_pseudopotentials_remote,
    detect_upf_files,
    extract_element_from_filename
)

__all__ = [
    'Pseudopotential',
    'PseudopotentialsConfig',
    'PseudopotentialsManager',
    'create_pseudopotentials_config',
    'load_pseudopotentials_config',
    'detect_pseudopotentials',
    'detect_pseudopotentials_remote',
    'detect_upf_files',
    'extract_element_from_filename',
    'DEFAULT_PSEUDOPOTENTIALS_DIR',
]

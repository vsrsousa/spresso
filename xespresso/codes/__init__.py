"""
codes module

This module provides utilities for managing Quantum ESPRESSO code configurations
on different machines. It allows users to:
- Define available QE executables and their paths
- Store QE version information
- Configure code settings per machine
- Detect available codes on a machine
- Support multiple QE versions on the same machine
- List available modules on remote systems
"""

from .manager import (
    CodesManager,
    create_codes_config,
    load_codes_config,
    detect_qe_codes,
    add_version_to_config,
)
from .config import Code, CodesConfig

__all__ = [
    'CodesManager',
    'Code',
    'CodesConfig',
    'create_codes_config',
    'load_codes_config',
    'detect_qe_codes',
    'add_version_to_config',
]

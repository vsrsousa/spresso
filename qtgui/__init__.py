"""
PySide6 GUI for xespresso - Quantum ESPRESSO Configuration Interface

This application provides a user-friendly interface for:
- Configuring machines (local/remote execution environments)
- Setting up Quantum ESPRESSO codes
- Viewing and selecting molecular structures
- Configuring calculations and workflows
- Submitting computational jobs

Two versions are available:
- Full version (default): Full-featured but may have stability issues
- Simple version (--simple flag): Simplified, more stable interface

Run with:
    python -m qtgui          # Full version
    python -m qtgui --simple # Simple stable version
"""

__version__ = "2.0.0"
__author__ = "xespresso team"

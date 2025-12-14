"""
PySide6 GUI for xespresso - Quantum ESPRESSO Configuration Interface

This application provides a user-friendly interface for:
- Configuring machines (local/remote execution environments)
- Setting up Quantum ESPRESSO codes
- Viewing and selecting molecular structures
- Configuring calculations and workflows
- Submitting computational jobs

Two versions are available:
- Full version (default): Opens the Session Manager to launch isolated workflow windows
- Simple version (--simple flag): Simplified, more stable interface

Run with:
    python -m qtgui             # Session Manager (full features)
    python -m qtgui --workspace # Open a single workflow window directly
    python -m qtgui --simple    # Simple stable version
"""

__version__ = "2.0.0"
__author__ = "xespresso team"

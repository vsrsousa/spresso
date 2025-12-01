"""
PySide6 entry point for xespresso GUI.

Run with: python -m qtgui
For the simplified stable version: python -m qtgui --simple
"""

import sys


def main():
    """Main entry point that chooses between simple and full GUI."""
    if '--simple' in sys.argv or '--stable' in sys.argv:
        # Use simplified, more stable version
        from qtgui.main_app_simple import main as simple_main
        simple_main()
    else:
        # Use full version (may have stability issues)
        from qtgui.main_app import main as full_main
        full_main()


if __name__ == "__main__":
    main()

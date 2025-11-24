#!/usr/bin/env python3
"""
Launcher script for xespresso Streamlit GUI.

Usage:
    xespresso-gui
    
Or:
    python -m xespresso.gui
"""

import sys
import os
from pathlib import Path

def main():
    """Launch the xespresso Streamlit GUI."""
    try:
        import streamlit.web.cli as stcli
    except ImportError:
        print("Error: Streamlit is not installed.")
        print("Please install it with: pip install streamlit")
        sys.exit(1)
    
    # Get the path to the streamlit app
    gui_dir = Path(__file__).parent
    app_path = gui_dir / "streamlit_app.py"
    
    if not app_path.exists():
        print(f"Error: Could not find streamlit_app.py at {app_path}")
        sys.exit(1)
    
    # Launch streamlit
    sys.argv = ["streamlit", "run", str(app_path)]
    sys.exit(stcli.main())

if __name__ == "__main__":
    main()

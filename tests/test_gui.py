"""
Tests for the GUI module.

These tests verify that the GUI components can be imported and
basic functionality works without actually running the Streamlit app.
"""

import pytest
import sys
from pathlib import Path


def test_gui_module_import():
    """Test that GUI module can be imported."""
    try:
        from gui import streamlit_app
        assert True
    except ImportError as e:
        pytest.fail(f"Failed to import GUI module: {e}")


def test_gui_dependencies():
    """Test that required GUI dependencies are available."""
    try:
        import streamlit
        import plotly.graph_objects as go
        assert True
    except ImportError as e:
        pytest.skip(f"GUI dependencies not installed: {e}")


def test_structure_visualization_function():
    """Test the 3D structure plot creation function."""
    try:
        import plotly.graph_objects as go
        from ase.build import bulk
        from gui.streamlit_app import create_3d_structure_plot
        
        # Create a simple structure
        atoms = bulk('Fe', 'bcc', a=2.87)
        
        # Create plot
        fig = create_3d_structure_plot(atoms)
        
        # Verify plot was created
        assert fig is not None
        assert len(fig.data) > 0  # Should have at least one trace
        
    except ImportError as e:
        pytest.skip(f"Required dependencies not available: {e}")


def test_display_structure_info_function():
    """Test that display_structure_info function exists."""
    try:
        from gui.streamlit_app import display_structure_info
        assert callable(display_structure_info)
    except ImportError as e:
        pytest.skip(f"GUI module not available: {e}")


def test_gui_launcher_exists():
    """Test that GUI launcher script exists."""
    launcher_path = Path(__file__).parent.parent / "scripts" / "xespresso-gui"
    assert launcher_path.exists(), "GUI launcher script not found"


def test_gui_main_module_exists():
    """Test that __main__.py exists for python -m execution."""
    main_path = Path(__file__).parent.parent / "xespresso" / "gui" / "__main__.py"
    assert main_path.exists(), "GUI __main__.py not found"


def test_gui_app_file_exists():
    """Test that streamlit_app.py exists."""
    app_path = Path(__file__).parent.parent / "xespresso" / "gui" / "streamlit_app.py"
    assert app_path.exists(), "streamlit_app.py not found"


def test_gui_readme_exists():
    """Test that GUI README exists."""
    readme_path = Path(__file__).parent.parent / "xespresso" / "gui" / "README.md"
    assert readme_path.exists(), "GUI README.md not found"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

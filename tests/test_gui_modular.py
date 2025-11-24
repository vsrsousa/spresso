"""
Tests for the modularized GUI structure.

These tests verify that the modular GUI components are properly structured
and can be imported (when dependencies are available).
"""

import pytest
import sys
from pathlib import Path


def test_gui_module_structure():
    """Test that GUI module structure exists."""
    gui_path = Path(__file__).parent.parent / "xespresso" / "gui"
    assert gui_path.exists(), "GUI directory should exist"
    
    # Check main files
    assert (gui_path / "__init__.py").exists()
    assert (gui_path / "streamlit_app.py").exists()
    assert (gui_path / "__main__.py").exists()
    
    # Check pages directory
    pages_path = gui_path / "pages"
    assert pages_path.exists(), "Pages directory should exist"
    assert (pages_path / "__init__.py").exists()
    
    # Check utils directory  
    utils_path = gui_path / "utils"
    assert utils_path.exists(), "Utils directory should exist"
    assert (utils_path / "__init__.py").exists()


def test_page_modules_exist():
    """Test that all page modules exist."""
    pages_path = Path(__file__).parent.parent / "xespresso" / "gui" / "pages"
    
    expected_modules = [
        "machine_config.py",
        "codes_config.py",
        "structure_viewer.py",
        "calculation_setup.py",
        "workflow_builder.py",
        "job_submission.py",
        "results_postprocessing.py"
    ]
    
    for module in expected_modules:
        assert (pages_path / module).exists(), f"{module} should exist"


def test_utility_modules_exist():
    """Test that all utility modules exist."""
    utils_path = Path(__file__).parent.parent / "xespresso" / "gui" / "utils"
    
    expected_modules = [
        "validation.py",
        "visualization.py",
        "connection.py",
        "dry_run.py"
    ]
    
    for module in expected_modules:
        assert (utils_path / module).exists(), f"{module} should exist"


def test_documentation_exists():
    """Test that documentation files exist."""
    gui_path = Path(__file__).parent.parent / "xespresso" / "gui"
    
    assert (gui_path / "README.md").exists(), "README should exist"
    assert (gui_path / "ARCHITECTURE.md").exists(), "ARCHITECTURE doc should exist"


def test_original_backup_exists():
    """Test that original streamlit_app was backed up."""
    gui_path = Path(__file__).parent.parent / "xespresso" / "gui"
    assert (gui_path / "streamlit_app_original.py").exists(), "Original backup should exist"


def test_validation_module_syntax():
    """Test that validation module has correct syntax."""
    import ast
    
    validation_path = Path(__file__).parent.parent / "xespresso" / "gui" / "utils" / "validation.py"
    
    with open(validation_path, 'r') as f:
        code = f.read()
    
    # This will raise SyntaxError if invalid
    ast.parse(code)


def test_machine_config_module_syntax():
    """Test that machine_config module has correct syntax."""
    import ast
    
    machine_config_path = Path(__file__).parent.parent / "xespresso" / "gui" / "pages" / "machine_config.py"
    
    with open(machine_config_path, 'r') as f:
        code = f.read()
    
    # This will raise SyntaxError if invalid
    ast.parse(code)


def test_dry_run_module_syntax():
    """Test that dry_run module has correct syntax."""
    import ast
    
    dry_run_path = Path(__file__).parent.parent / "xespresso" / "gui" / "utils" / "dry_run.py"
    
    with open(dry_run_path, 'r') as f:
        code = f.read()
    
    # This will raise SyntaxError if invalid
    ast.parse(code)


def test_main_app_reduced_size():
    """Test that main app was actually reduced in size."""
    gui_path = Path(__file__).parent.parent / "xespresso" / "gui"
    
    # Read both files
    with open(gui_path / "streamlit_app.py", 'r') as f:
        new_lines = len(f.readlines())
    
    with open(gui_path / "streamlit_app_original.py", 'r') as f:
        original_lines = len(f.readlines())
    
    # New version should be smaller
    assert new_lines < original_lines, f"New app should be smaller (was {original_lines}, now {new_lines})"
    
    # Should be at least 30% reduction
    reduction = (original_lines - new_lines) / original_lines
    assert reduction >= 0.30, f"Should have at least 30% reduction (got {reduction*100:.1f}%)"


def test_xespresso_integration_documented():
    """Test that xespresso integrations are documented."""
    arch_doc = Path(__file__).parent.parent / "xespresso" / "gui" / "ARCHITECTURE.md"
    
    with open(arch_doc, 'r') as f:
        content = f.read()
    
    # Check that key integrations are mentioned
    assert "test_ssh_connection" in content, "SSH connection testing integration should be documented"
    assert "write_input" in content, "Input file generation integration should be documented"
    assert "Machine" in content, "Machine class integration should be documented"
    assert "detect_qe_codes" in content, "Code detection integration should be documented"


def test_modular_structure_imports():
    """Test that modular structure can be imported (module-level, not functions)."""
    # This just tests that the module files can be found
    # Actual function imports may fail without dependencies
    
    try:
        import xespresso.gui.pages
        import xespresso.gui.utils
        # If we get here, module structure is correct
        assert True
    except ImportError as e:
        # ImportError is OK if it's about missing dependencies (ase, streamlit, etc)
        # but not if it's about module structure
        if "xespresso.gui" in str(e) and "ase" not in str(e).lower():
            pytest.fail(f"Module structure issue: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

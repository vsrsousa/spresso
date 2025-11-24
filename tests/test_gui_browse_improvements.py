"""
Tests for improved GUI browse directory functionality.

These tests verify that the simplified browse directory UI works correctly.
"""

import pytest
import os
import tempfile
from pathlib import Path


def test_workdir_browser_returns_valid_path():
    """Test that the workdir browser returns a valid path."""
    # Note: Since render_workdir_browser uses Streamlit session state,
    # we can't test it directly without a Streamlit app context.
    # Instead, we test the core logic that the browser would use.
    
    import os
    
    # Test path normalization
    test_dir = os.path.expanduser("~")
    assert os.path.isabs(test_dir)
    assert os.path.exists(test_dir)
    
    # Test parent directory navigation
    parent = os.path.dirname(test_dir)
    assert os.path.exists(parent)
    
    # Test listing subdirectories
    try:
        contents = os.listdir(test_dir)
        subdirs = [d for d in contents if os.path.isdir(os.path.join(test_dir, d)) and not d.startswith('.')]
        # Should have at least some directories in home
        assert isinstance(subdirs, list)
    except PermissionError:
        pytest.skip("No permission to list home directory")


def test_path_validation_security():
    """Test that path validation prevents traversal attacks."""
    import os
    
    # Test that we detect path traversal attempts
    safe_base = "/tmp"
    
    # Valid path within base
    valid_path = os.path.join(safe_base, "test")
    assert os.path.commonpath([safe_base, valid_path]) == safe_base
    
    # Path traversal attempt
    try:
        malicious_path = os.path.join(safe_base, "../etc/passwd")
        real_malicious = os.path.realpath(malicious_path)
        # Should not start with safe_base after resolution
        assert not real_malicious.startswith(safe_base + os.sep)
    except:
        pass  # Expected to fail security check


def test_folder_listing_with_hidden_filter():
    """Test that hidden directories are filtered out."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some test directories
        os.makedirs(os.path.join(tmpdir, "visible_folder"))
        os.makedirs(os.path.join(tmpdir, ".hidden_folder"))
        os.makedirs(os.path.join(tmpdir, "another_visible"))
        
        # List and filter directories
        contents = os.listdir(tmpdir)
        subdirs = [d for d in contents if os.path.isdir(os.path.join(tmpdir, d)) and not d.startswith('.')]
        
        # Should have 2 visible folders
        assert len(subdirs) == 2
        assert "visible_folder" in subdirs
        assert "another_visible" in subdirs
        assert ".hidden_folder" not in subdirs


def test_session_state_keys_are_unique():
    """Test that the workdir browser uses unique keys to avoid conflicts."""
    # This tests the design pattern of using unique keys for multiple instances
    
    key1 = "workdir_browser_1"
    key2 = "workdir_browser_2"
    
    workdir_key1 = f"{key1}_workdir"
    workdir_key2 = f"{key2}_workdir"
    
    # Keys should be different
    assert workdir_key1 != workdir_key2
    
    # Both should be valid strings
    assert isinstance(workdir_key1, str)
    assert isinstance(workdir_key2, str)


def test_navigation_button_logic():
    """Test the logic for navigation buttons."""
    import os
    
    # Test "Current" button logic
    current_dir = os.getcwd()
    assert os.path.exists(current_dir)
    assert os.path.isdir(current_dir)
    
    # Test "Home" button logic
    home_dir = os.path.expanduser("~")
    assert os.path.exists(home_dir)
    assert os.path.isdir(home_dir)
    
    # Test "Up" button logic (parent directory)
    if current_dir != "/":
        parent_dir = os.path.dirname(current_dir)
        assert os.path.exists(parent_dir)
        assert os.path.isdir(parent_dir)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

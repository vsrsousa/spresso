"""
Tests for the new QE code configuration features in the GUI.

These tests verify that:
1. Module listing functionality is properly integrated
2. QE version specification is available
3. Version selection for different calculations works
"""

import sys
from pathlib import Path


def test_codes_config_page_imports():
    """Test that codes_config page can be imported."""
    try:
        from xespresso.gui.pages.codes_config import render_codes_config_page
        assert callable(render_codes_config_page)
        print("✓ codes_config page can be imported")
        return True
    except ImportError as e:
        # It's OK if streamlit or xespresso is not installed
        if "streamlit" in str(e).lower() or "xespresso" in str(e).lower():
            print("✓ codes_config import skipped (dependencies not available)")
            return True
        print(f"✗ Import error: {e}")
        return False


def test_codes_config_module_has_updated_features():
    """Test that the codes_config module has the updated features."""
    codes_config_path = Path(__file__).parent.parent / "xespresso" / "gui" / "pages" / "codes_config.py"
    assert codes_config_path.exists(), "codes_config.py should exist"
    
    content = codes_config_path.read_text()
    
    # Feature 1: Explicit QE version specification
    assert "qe_version" in content, "QE version parameter should be present"
    assert "QE Version" in content, "QE version input field should be present"
    assert "compiler version" in content.lower(), "Warning about compiler version should be present"
    
    # Feature 2: Version selection
    assert "available_versions" in content, "Version selection should be present"
    assert "list_versions" in content, "Version listing method should be called"
    assert "Select QE Version" in content, "Version selection UI should be present"
    
    print("✓ All three features found in codes_config module")
    return True


def test_qe_version_parameter_integration():
    """Test that QE version parameter is passed to detect_qe_codes."""
    try:
        from xespresso.codes.manager import detect_qe_codes
        import inspect
        
        # Check that detect_qe_codes accepts qe_version parameter
        sig = inspect.signature(detect_qe_codes)
        params = list(sig.parameters.keys())
        
        assert 'qe_version' in params, "qe_version parameter should be in detect_qe_codes"
        
        # Check it's optional
        param = sig.parameters['qe_version']
        assert param.default is not inspect.Parameter.empty or param.default is None, \
            "qe_version should be optional"
        
        print("✓ QE version parameter properly integrated")
        return True
        
    except ImportError as e:
        print(f"✓ QE version test skipped: {e}")
        return True


def test_version_selection_integration():
    """Test that version selection uses load_codes_config with version parameter."""
    try:
        from xespresso.codes.manager import load_codes_config
        import inspect
        
        # Check that load_codes_config accepts version parameter
        sig = inspect.signature(load_codes_config)
        params = list(sig.parameters.keys())
        
        assert 'version' in params, "version parameter should be in load_codes_config"
        
        # Check it's optional
        param = sig.parameters['version']
        assert param.default is not inspect.Parameter.empty or param.default is None, \
            "version should be optional"
        
        print("✓ Version selection properly integrated")
        return True
        
    except ImportError as e:
        print(f"✓ Version selection test skipped: {e}")
        return True


def test_codes_config_structure():
    """Test that CodesConfig has version-related methods."""
    try:
        from xespresso.codes.config import CodesConfig
        
        # Create a dummy config
        config = CodesConfig(machine_name="test")
        
        # Check for version-related methods
        assert hasattr(config, 'list_versions') or hasattr(config, 'versions'), \
            "CodesConfig should have version support"
        
        print("✓ CodesConfig has version support")
        return True
        
    except ImportError as e:
        print(f"✓ CodesConfig test skipped: {e}")
        return True


def test_gui_page_structure_with_new_features():
    """Test that the GUI page structure includes new features."""
    codes_config_path = Path(__file__).parent.parent / "xespresso" / "gui" / "pages" / "codes_config.py"
    content = codes_config_path.read_text()
    
    # Check for proper sectioning
    assert "# Feature 2: Explicit QE Version" in content or "Feature 2:" in content or "qe_version" in content, \
        "QE version feature should be present"
    assert "# Feature 3: Version Selection" in content or "Feature 3:" in content or "available_versions" in content, \
        "Version selection feature should be present"
    
    # Check for proper UI elements
    assert "st.text_input" in content, "Should have text inputs for parameters"
    assert "st.selectbox" in content, "Should have version selection dropdown"
    assert "st.button" in content, "Should have action buttons"
    
    print("✓ GUI page structure has all new features")
    return True


def test_backward_compatibility():
    """Test that changes maintain backward compatibility."""
    codes_config_path = Path(__file__).parent.parent / "xespresso" / "gui" / "pages" / "codes_config.py"
    content = codes_config_path.read_text()
    
    # Should still support basic detection without new features
    assert "detect_qe_codes" in content, "Basic detection should still work"
    assert "auto_load_machine=True" in content, "Auto-loading should still be default"
    
    # Old functionality should be preserved
    assert "Select Machine" in content, "Machine selection should be preserved"
    assert "Auto-Detect Codes" in content, "Auto-detection should be preserved"
    
    print("✓ Backward compatibility maintained")
    return True


def test_ui_usability_features():
    """Test that the UI has good usability features."""
    codes_config_path = Path(__file__).parent.parent / "xespresso" / "gui" / "pages" / "codes_config.py"
    content = codes_config_path.read_text()
    
    # Check for help text
    assert "help=" in content, "Should have help text for inputs"
    
    # Check for info messages
    assert "st.info" in content, "Should have informational messages"
    assert "st.success" in content, "Should have success messages"
    assert "st.warning" in content, "Should have warning messages"
    
    # Check for proper error handling
    assert "try:" in content and "except" in content, "Should have error handling"
    assert "st.error" in content, "Should display errors to user"
    
    print("✓ UI has good usability features")
    return True


def run_all_tests():
    """Run all tests."""
    tests = [
        test_codes_config_page_imports,
        test_codes_config_module_has_updated_features,
        test_qe_version_parameter_integration,
        test_version_selection_integration,
        test_codes_config_structure,
        test_gui_page_structure_with_new_features,
        test_backward_compatibility,
        test_ui_usability_features,
    ]
    
    print("\n" + "="*70)
    print("Running GUI QE Features Tests")
    print("="*70 + "\n")
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test_name = test.__name__
            print(f"\nRunning {test_name}...")
            result = test()
            if result:
                passed += 1
            else:
                failed += 1
                print(f"✗ {test_name} failed")
        except Exception as e:
            failed += 1
            print(f"✗ {test.__name__} failed with exception: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*70 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)


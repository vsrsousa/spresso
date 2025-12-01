"""
Tests for the simplified Qt GUI.

These tests verify the simplified GUI components work correctly.
"""

import pytest
import os


class TestSimpleSessionState:
    """Tests for SimpleSessionState."""
    
    def test_simple_state_creation(self):
        """Test creating a simple session state."""
        from qtgui.main_app_simple import SimpleSessionState
        
        state = SimpleSessionState()
        assert state is not None
    
    def test_simple_state_defaults(self):
        """Test default values."""
        from qtgui.main_app_simple import SimpleSessionState
        
        state = SimpleSessionState()
        
        assert state['current_structure'] is None
        assert state['working_directory'] == os.path.expanduser("~")
        assert state['workflow_config'] == {}
    
    def test_simple_state_get_set(self):
        """Test getting and setting values."""
        from qtgui.main_app_simple import SimpleSessionState
        
        state = SimpleSessionState()
        
        state['test_key'] = 'test_value'
        assert state['test_key'] == 'test_value'
        assert state.get('test_key') == 'test_value'
        assert state.get('nonexistent', 'default') == 'default'
    
    def test_simple_state_set_method(self):
        """Test the set method."""
        from qtgui.main_app_simple import SimpleSessionState
        
        state = SimpleSessionState()
        
        state.set('key', 'value')
        assert state.get('key') == 'value'


class TestSimpleMainWindow:
    """Tests for SimpleMainWindow."""
    
    def test_main_window_import(self):
        """Test importing SimpleMainWindow."""
        from qtgui.main_app_simple import SimpleMainWindow
        assert SimpleMainWindow is not None
    
    def test_main_function_exists(self):
        """Test main function exists."""
        from qtgui.main_app_simple import main
        assert callable(main)


class TestSimpleTabs:
    """Tests for the tab widgets."""
    
    def test_structure_tab_import(self):
        """Test importing StructureTab."""
        from qtgui.main_app_simple import StructureTab
        assert StructureTab is not None
    
    def test_calculation_tab_import(self):
        """Test importing CalculationTab."""
        from qtgui.main_app_simple import CalculationTab
        assert CalculationTab is not None
    
    def test_file_browser_tab_import(self):
        """Test importing FileBrowserTab."""
        from qtgui.main_app_simple import FileBrowserTab
        assert FileBrowserTab is not None


class TestEntryPoint:
    """Tests for the entry point."""
    
    def test_main_module_import(self):
        """Test __main__ module imports correctly."""
        from qtgui import __main__
        assert hasattr(__main__, 'main')
        assert callable(__main__.main)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

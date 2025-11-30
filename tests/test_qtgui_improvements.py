"""
Tests for the improved Qt GUI functionality.

These tests verify:
- Enhanced SessionState with multiple sessions, save/load
- ConfigurationDialog non-blocking behavior
- Session management features
"""

import pytest
import os
import json
import tempfile
from pathlib import Path


class TestSessionState:
    """Tests for the SessionState class."""
    
    def test_session_state_singleton(self):
        """Test that SessionState is a singleton."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state1 = SessionState()
        state2 = SessionState()
        assert state1 is state2
    
    def test_session_state_defaults(self):
        """Test default session state values."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        assert state['current_structure'] is None
        assert state['current_machine'] is None
        assert state['current_machine_name'] is None
        assert state['current_codes'] is None
        assert state['selected_code_version'] is None
        assert state['workflow_config'] == {}
        assert state['working_directory'] == os.path.expanduser("~")
        assert state['session_name'] == "Default Session"
    
    def test_session_state_get_set(self):
        """Test getting and setting state values."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        state['test_key'] = 'test_value'
        assert state['test_key'] == 'test_value'
        assert state.get('test_key') == 'test_value'
        assert state.get('nonexistent', 'default') == 'default'
    
    def test_session_state_contains(self):
        """Test the contains check."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        assert 'working_directory' in state
        assert 'nonexistent_key' not in state
    
    def test_session_name(self):
        """Test getting session name."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        assert state.get_session_name() == "Default Session"
        state['session_name'] = "My Custom Session"
        assert state.get_session_name() == "My Custom Session"
    
    def test_session_id(self):
        """Test getting current session ID."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        assert state.get_current_session_id() == "default"
    
    def test_session_reset(self):
        """Test resetting session to defaults."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        state['working_directory'] = '/some/custom/path'
        state['custom_key'] = 'custom_value'
        
        state.reset()
        
        assert state['working_directory'] == os.path.expanduser("~")
        assert state['session_name'] == "Default Session"
    
    def test_session_listeners(self):
        """Test state change listeners."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        
        callback_called = {'count': 0}
        
        def listener():
            callback_called['count'] += 1
        
        state.add_listener(listener)
        state['test_key'] = 'value1'
        
        assert callback_called['count'] == 1
        
        state['test_key'] = 'value2'
        assert callback_called['count'] == 2
        
        state.remove_listener(listener)
        state['test_key'] = 'value3'
        assert callback_called['count'] == 2  # Should not increase


class TestSessionPersistence:
    """Tests for session save/load functionality."""
    
    def test_session_save_load(self):
        """Test saving and loading a session."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Set some values
            state['working_directory'] = '/test/path'
            state['session_name'] = 'Test Session'
            
            # Save session
            state.save_session()
            
            # Verify file was created
            session_file = os.path.join(tmpdir, "default.json")
            assert os.path.exists(session_file)
            
            # Verify content
            with open(session_file, 'r') as f:
                saved_data = json.load(f)
            
            assert saved_data['working_directory'] == '/test/path'
            assert saved_data['session_name'] == 'Test Session'
    
    def test_create_session(self):
        """Test creating a new session."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Create new session
            new_id = state.create_session("New Session")
            
            assert new_id is not None
            assert new_id.startswith("session_")
            assert state.get_current_session_id() == new_id
            assert state.get_session_name() == "New Session"
    
    def test_list_sessions(self):
        """Test listing sessions."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Initially empty
            sessions = state.list_sessions()
            assert len(sessions) == 0
            
            # Create sessions
            state.create_session("Session 1")
            state.create_session("Session 2")
            
            sessions = state.list_sessions()
            assert len(sessions) >= 2


class TestConfigurationDialog:
    """Tests for the ConfigurationDialog class."""
    
    def test_dialog_import(self):
        """Test that ConfigurationDialog can be imported."""
        try:
            from qtgui.dialogs import ConfigurationDialog
            assert ConfigurationDialog is not None
        except ImportError as e:
            pytest.fail(f"Failed to import ConfigurationDialog: {e}")
    
    def test_dialog_attributes(self):
        """Test that ConfigurationDialog has expected attributes."""
        from qtgui.dialogs import ConfigurationDialog
        
        assert hasattr(ConfigurationDialog, 'configuration_changed')
        assert hasattr(ConfigurationDialog, 'show_machine_tab')
        assert hasattr(ConfigurationDialog, 'show_codes_tab')
        assert hasattr(ConfigurationDialog, 'show_pseudopotentials_tab')


class TestMainApp:
    """Tests for the main application."""
    
    def test_main_function_exists(self):
        """Test that main function exists."""
        from qtgui.main_app import main
        assert callable(main)
    
    def test_main_window_import(self):
        """Test that MainWindow can be imported."""
        from qtgui.main_app import MainWindow
        assert MainWindow is not None
    
    def test_session_state_import(self):
        """Test that session_state can be imported."""
        from qtgui.main_app import session_state
        assert session_state is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

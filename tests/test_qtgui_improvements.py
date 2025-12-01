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
            
            # Verify file was created with session name (not session ID)
            session_file = os.path.join(tmpdir, "Test Session.json")
            assert os.path.exists(session_file), f"Expected 'Test Session.json' but found: {os.listdir(tmpdir)}"
            
            # Verify content
            with open(session_file, 'r') as f:
                saved_data = json.load(f)
            
            assert saved_data['working_directory'] == '/test/path'
            assert saved_data['session_name'] == 'Test Session'
            # Verify session ID is stored inside the file
            assert '_session_id' in saved_data
            assert saved_data['_session_id'] == 'default'
    
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
            # Set the sessions dir and clear existing sessions loaded from default
            state._sessions_dir = tmpdir
            state._sessions = {}
            
            # Initially empty
            sessions = state.list_sessions()
            assert len(sessions) == 0
            
            # Create sessions
            state.create_session("Session 1")
            state.create_session("Session 2")
            
            sessions = state.list_sessions()
            assert len(sessions) >= 2
    
    def test_get_sessions_dir(self):
        """Test getting the sessions directory."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            assert state.get_sessions_dir() == tmpdir
    
    def test_list_session_files(self):
        """Test listing session files from the filesystem."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Create some test session files with portable paths
            # Include _session_id inside the file (as per the new behavior)
            session1_data = {'session_name': 'Test Session 1', '_session_id': 'session_001', 'working_directory': tempfile.gettempdir()}
            session2_data = {'session_name': 'Test Session 2', '_session_id': 'session_002', 'working_directory': os.path.expanduser('~')}
            
            with open(os.path.join(tmpdir, 'Test Session 1.json'), 'w') as f:
                json.dump(session1_data, f)
            
            with open(os.path.join(tmpdir, 'Test Session 2.json'), 'w') as f:
                json.dump(session2_data, f)
            
            # Create sessions_index.json which should be ignored
            with open(os.path.join(tmpdir, 'sessions_index.json'), 'w') as f:
                json.dump({}, f)
            
            # List session files
            session_files = state.list_session_files()
            
            # Should find 2 session files (excluding sessions_index.json)
            assert len(session_files) == 2
            
            # Check that session names are extracted correctly
            session_names = [s[1] for s in session_files]
            assert 'Test Session 1' in session_names
            assert 'Test Session 2' in session_names
            
            # Check that session IDs are extracted correctly from inside the files
            session_ids = [s[0] for s in session_files]
            assert 'session_001' in session_ids
            assert 'session_002' in session_ids
    
    def test_list_session_files_empty_dir(self):
        """Test listing session files in empty directory."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            session_files = state.list_session_files()
            assert session_files == []
    
    def test_load_session_from_file(self):
        """Test loading a session directly from a file."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Create a test session file with _session_id inside
            session_data = {
                'session_name': 'Test Loaded Session',
                '_session_id': 'my_custom_session_id',
                'working_directory': '/test/dir',
                'session_created': '2024-01-01T00:00:00',
                'session_modified': '2024-01-01T00:00:00'
            }
            
            session_file = os.path.join(tmpdir, 'Test Loaded Session.json')
            with open(session_file, 'w') as f:
                json.dump(session_data, f)
            
            # Load the session
            result = state.load_session_from_file(session_file)
            
            assert result is True
            assert state.get_session_name() == 'Test Loaded Session'
            assert state['working_directory'] == '/test/dir'
            # Session ID should be read from inside the file
            assert state.get_current_session_id() == 'my_custom_session_id'
    
    def test_load_session_from_file_nonexistent(self):
        """Test loading a session from a non-existent file."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Use a path within the temp directory that doesn't exist
            nonexistent_path = os.path.join(tmpdir, 'nonexistent', 'file.json')
            result = state.load_session_from_file(nonexistent_path)
            
            assert result is False
    
    def test_load_session_from_file_invalid_format(self):
        """Test loading a session from a file with invalid format."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SessionState()
            state._sessions_dir = tmpdir
            
            # Create an invalid session file (not JSON)
            session_file = os.path.join(tmpdir, 'invalid.json')
            with open(session_file, 'w') as f:
                f.write('not valid json')
            
            result = state.load_session_from_file(session_file)
            
            assert result is False


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

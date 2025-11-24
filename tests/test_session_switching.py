"""
Tests for session switching and session-specific working directory functionality.
"""

import pytest
import os

# Mock streamlit session_state for testing
class MockSessionState(dict):
    """Mock Streamlit session state."""
    def __getattr__(self, key):
        return self.get(key)
    
    def __setattr__(self, key, value):
        self[key] = value


@pytest.fixture
def mock_st(monkeypatch):
    """Mock streamlit for testing."""
    import sys
    from unittest.mock import MagicMock
    
    # Create mock streamlit module
    mock_streamlit = MagicMock()
    # Create a fresh session state for each test
    mock_streamlit.session_state = MockSessionState()
    
    # Mock the streamlit import
    sys.modules['streamlit'] = mock_streamlit
    
    yield mock_streamlit
    
    # Cleanup
    if 'streamlit' in sys.modules:
        del sys.modules['streamlit']
    # Clear session state after each test
    if hasattr(mock_streamlit, 'session_state'):
        mock_streamlit.session_state.clear()


def test_create_new_session(mock_st):
    """Test creating a new session."""
    from gui.utils.session_manager import create_new_session, get_current_session_id
    
    # Create first session
    session_id = create_new_session()
    
    # Verify session was created
    assert session_id == 'session_1'
    assert get_current_session_id() == session_id
    assert '_active_sessions' in mock_st.session_state
    assert session_id in mock_st.session_state._active_sessions
    
    # Verify session has default working directory
    assert 'working_directory' in mock_st.session_state
    assert mock_st.session_state.working_directory == os.path.expanduser("~")
    
    # Create second session
    session_id2 = create_new_session()
    
    # Verify second session
    assert session_id2 == 'session_2'
    assert get_current_session_id() == session_id2
    # Should have reset working directory to default
    assert mock_st.session_state.working_directory == os.path.expanduser("~")


def test_switch_session_preserves_working_directory(mock_st):
    """Test that switching sessions preserves each session's working directory."""
    from gui.utils.session_manager import (
        create_new_session, switch_session, get_current_session_id
    )
    
    # Create first session
    session1_id = create_new_session()
    
    # Set a specific working directory for session 1
    test_dir1 = "/tmp/session1_workdir"
    mock_st.session_state.working_directory = test_dir1
    
    # Create second session (this should save session 1's state)
    session2_id = create_new_session()
    
    # Set a different working directory for session 2
    test_dir2 = "/tmp/session2_workdir"
    mock_st.session_state.working_directory = test_dir2
    
    # Switch back to session 1
    switch_session(session1_id)
    
    # Verify we're in session 1
    assert get_current_session_id() == session1_id
    
    # Verify session 1's working directory was restored
    assert mock_st.session_state.working_directory == test_dir1
    
    # Switch to session 2
    switch_session(session2_id)
    
    # Verify we're in session 2
    assert get_current_session_id() == session2_id
    
    # Verify session 2's working directory was restored
    assert mock_st.session_state.working_directory == test_dir2


def test_rename_session(mock_st):
    """Test renaming a session."""
    from gui.utils.session_manager import (
        create_new_session, rename_session, get_active_sessions
    )
    
    # Create a session
    session_id = create_new_session()
    
    # Verify default name
    sessions = get_active_sessions()
    assert sessions[session_id]['name'] == 'Session 1'
    
    # Rename session
    rename_session(session_id, 'My Custom Session')
    
    # Verify name was changed
    sessions = get_active_sessions()
    assert sessions[session_id]['name'] == 'My Custom Session'


def test_close_session(mock_st):
    """Test closing a session."""
    from gui.utils.session_manager import (
        create_new_session, close_session, get_active_sessions, get_current_session_id
    )
    
    # Create two sessions
    session1_id = create_new_session()
    session2_id = create_new_session()
    
    # Verify both exist
    sessions = get_active_sessions()
    assert len(sessions) == 2
    assert session1_id in sessions
    assert session2_id in sessions
    
    # Close session 1
    close_session(session1_id)
    
    # Verify session 1 was closed
    sessions = get_active_sessions()
    assert len(sessions) == 1
    assert session1_id not in sessions
    assert session2_id in sessions
    
    # Current session should still be session 2
    assert get_current_session_id() == session2_id


def test_close_current_session_switches_to_another(mock_st):
    """Test that closing current session switches to another session."""
    from gui.utils.session_manager import (
        create_new_session, close_session, get_active_sessions, get_current_session_id
    )
    
    # Create two sessions
    session1_id = create_new_session()
    session2_id = create_new_session()
    
    # Current session is session 2
    assert get_current_session_id() == session2_id
    
    # Close current session
    close_session(session2_id)
    
    # Should have switched to session 1
    assert get_current_session_id() == session1_id
    
    # Verify only session 1 remains
    sessions = get_active_sessions()
    assert len(sessions) == 1
    assert session1_id in sessions


def test_session_state_independence(mock_st):
    """Test that sessions maintain independent state."""
    from gui.utils.session_manager import (
        create_new_session, switch_session, get_current_session_id
    )
    
    # Create first session
    session1_id = create_new_session()
    mock_st.session_state.current_structure = "structure1"
    mock_st.session_state.workflow_config = {"ecutwfc": 50}
    mock_st.session_state.working_directory = "/tmp/dir1"
    
    # Create second session
    session2_id = create_new_session()
    
    # Session 2 should have clean state
    assert 'current_structure' not in mock_st.session_state or mock_st.session_state.current_structure is None
    assert mock_st.session_state.get('workflow_config', {}) == {}
    
    # Set different values in session 2
    mock_st.session_state.current_structure = "structure2"
    mock_st.session_state.workflow_config = {"ecutwfc": 80}
    mock_st.session_state.working_directory = "/tmp/dir2"
    
    # Switch back to session 1
    switch_session(session1_id)
    
    # Verify session 1's state was restored
    assert mock_st.session_state.current_structure == "structure1"
    assert mock_st.session_state.workflow_config == {"ecutwfc": 50}
    assert mock_st.session_state.working_directory == "/tmp/dir1"
    
    # Switch to session 2
    switch_session(session2_id)
    
    # Verify session 2's state
    assert mock_st.session_state.current_structure == "structure2"
    assert mock_st.session_state.workflow_config == {"ecutwfc": 80}
    assert mock_st.session_state.working_directory == "/tmp/dir2"


def test_get_current_session_id_creates_default(mock_st):
    """Test that get_current_session_id creates default session if none exists."""
    from gui.utils.session_manager import get_current_session_id
    
    # Should create session_1 automatically
    session_id = get_current_session_id()
    assert session_id == 'session_1'
    assert '_current_session_id' in mock_st.session_state
    assert mock_st.session_state._current_session_id == 'session_1'


def test_get_active_sessions_initializes_empty(mock_st):
    """Test that get_active_sessions initializes empty dict if needed."""
    from gui.utils.session_manager import get_active_sessions
    
    # Should initialize empty dict
    sessions = get_active_sessions()
    assert isinstance(sessions, dict)
    assert len(sessions) == 0
    assert '_active_sessions' in mock_st.session_state

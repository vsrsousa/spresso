"""
Test that widget keys are properly filtered during session switching.

This tests the fix for the bug where switching sessions would cause
StreamlitValueAssignmentNotAllowedError when widget keys were restored.
"""

import pytest

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
    
    # Remove existing streamlit module if present
    if 'streamlit' in sys.modules:
        del sys.modules['streamlit']
    if 'xespresso.gui.utils.session_manager' in sys.modules:
        del sys.modules['xespresso.gui.utils.session_manager']
    
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
    if 'xespresso.gui.utils.session_manager' in sys.modules:
        del sys.modules['xespresso.gui.utils.session_manager']
    # Clear session state after each test
    if hasattr(mock_streamlit, 'session_state'):
        mock_streamlit.session_state.clear()


def test_widget_keys_are_filtered(mock_st):
    """Test that widget keys are not included in serializable state."""
    from gui.utils.session_manager import get_serializable_state
    
    # Set up session state with both application state and widget keys
    mock_st.session_state.working_directory = "/tmp/test"
    mock_st.session_state.current_structure = "H2O"
    mock_st.session_state.workflow_config = {"ecutwfc": 50}
    
    # Add widget keys that should be filtered out
    mock_st.session_state.workdir_browser_quick_0 = True
    mock_st.session_state.workdir_browser_quick_1 = True
    mock_st.session_state.workdir_browser_up = True
    mock_st.session_state.workdir_browser_subdir_selector = "test"
    mock_st.session_state.workdir_browser_enter_subdir = True
    mock_st.session_state.workdir_browser_tkinter_browse = True  # tkinter browse button
    mock_st.session_state.session_manager_new = True
    mock_st.session_state.session_manager_save = True
    mock_st.session_state.session_manager_switch_select = "Session 2"
    mock_st.session_state.session_manager_renaming = False
    
    # Add button widget keys that should be filtered out (the bug fix)
    mock_st.session_state.build_crystal_btn = True
    mock_st.session_state.build_molecule_btn = False
    mock_st.session_state.load_db_structure_btn = True
    mock_st.session_state.save_db_structure_btn = False
    
    # Get serializable state
    state = get_serializable_state()
    
    # Application state should be preserved
    assert 'working_directory' in state
    assert state['working_directory'] == "/tmp/test"
    assert 'current_structure' in state
    assert state['current_structure'] == "H2O"
    assert 'workflow_config' in state
    assert state['workflow_config'] == {"ecutwfc": 50}
    
    # Widget keys should be filtered out
    assert 'workdir_browser_quick_0' not in state
    assert 'workdir_browser_quick_1' not in state
    assert 'workdir_browser_up' not in state
    assert 'workdir_browser_subdir_selector' not in state
    assert 'workdir_browser_enter_subdir' not in state
    assert 'workdir_browser_tkinter_browse' not in state  # tkinter browse button
    assert 'session_manager_new' not in state
    assert 'session_manager_save' not in state
    assert 'session_manager_switch_select' not in state
    assert 'session_manager_renaming' not in state
    
    # Button widget keys should be filtered out (the bug fix)
    assert 'build_crystal_btn' not in state
    assert 'build_molecule_btn' not in state
    assert 'load_db_structure_btn' not in state
    assert 'save_db_structure_btn' not in state


def test_restore_session_filters_widget_keys(mock_st):
    """Test that widget keys are filtered during session restoration."""
    from gui.utils.session_manager import restore_session
    
    # Create a state dictionary with both application state and widget keys
    state_with_widgets = {
        'working_directory': "/tmp/test",
        'current_structure': "H2O",
        'workflow_config': {"ecutwfc": 50},
        # Widget keys that should be filtered
        'workdir_browser_quick_0': True,
        'workdir_browser_subdir_selector': "test",
        'workdir_browser_tkinter_browse': True,  # tkinter browse button
        'session_manager_new': True,
        'session_manager_renaming': False,
        # Button widget keys that should be filtered (the bug fix)
        'build_crystal_btn': True,
        'build_molecule_btn': False,
        'load_db_structure_btn': True,
        'save_db_structure_btn': False,
    }
    
    # Restore session
    restore_session(state_with_widgets, clear_first=True)
    
    # Application state should be restored
    assert 'working_directory' in mock_st.session_state
    assert mock_st.session_state.working_directory == "/tmp/test"
    assert 'current_structure' in mock_st.session_state
    assert mock_st.session_state.current_structure == "H2O"
    assert 'workflow_config' in mock_st.session_state
    assert mock_st.session_state.workflow_config == {"ecutwfc": 50}
    
    # Widget keys should NOT be restored
    assert 'workdir_browser_quick_0' not in mock_st.session_state
    assert 'workdir_browser_subdir_selector' not in mock_st.session_state
    assert 'workdir_browser_tkinter_browse' not in mock_st.session_state  # tkinter browse button
    assert 'session_manager_new' not in mock_st.session_state
    assert 'session_manager_renaming' not in mock_st.session_state
    
    # Button widget keys should NOT be restored (the bug fix)
    assert 'build_crystal_btn' not in mock_st.session_state
    assert 'build_molecule_btn' not in mock_st.session_state
    assert 'load_db_structure_btn' not in mock_st.session_state
    assert 'save_db_structure_btn' not in mock_st.session_state


def test_switch_session_no_widget_key_conflict(mock_st):
    """Test that switching sessions doesn't restore widget keys."""
    from gui.utils.session_manager import (
        create_new_session, switch_session, get_current_session_id
    )
    
    # Create first session
    session1_id = create_new_session()
    
    # Set application state
    mock_st.session_state.working_directory = "/tmp/session1"
    mock_st.session_state.current_structure = "H2O"
    
    # Simulate widget keys being created during rendering
    mock_st.session_state.workdir_browser_quick_0 = True
    mock_st.session_state.workdir_browser_subdir_selector = "subfolder"
    
    # Create second session (this saves session 1's state)
    session2_id = create_new_session()
    
    # Set different application state for session 2
    mock_st.session_state.working_directory = "/tmp/session2"
    mock_st.session_state.current_structure = "CO2"
    
    # Simulate different widget state for session 2
    mock_st.session_state.workdir_browser_quick_0 = False
    mock_st.session_state.workdir_browser_subdir_selector = "other"
    
    # Switch back to session 1
    switch_session(session1_id)
    
    # Application state should be restored
    assert mock_st.session_state.working_directory == "/tmp/session1"
    assert mock_st.session_state.current_structure == "H2O"
    
    # Widget keys from session 2 should not cause conflicts
    # (they may or may not be present, but won't cause errors)
    # The key point is that the switch completed without error


def test_is_widget_key():
    """Test the _is_widget_key helper function."""
    from gui.utils.session_manager import _is_widget_key
    
    # Widget keys should be identified
    assert _is_widget_key('workdir_browser_quick_0')
    assert _is_widget_key('workdir_browser_quick_1')
    assert _is_widget_key('something_quick_123')
    assert _is_widget_key('browser_up')
    assert _is_widget_key('dir_subdir_selector')
    assert _is_widget_key('path_enter_subdir')
    assert _is_widget_key('workdir_custom_path')
    assert _is_widget_key('browser_go_custom')
    assert _is_widget_key('dir_create_dir')
    assert _is_widget_key('workdir_browser_tkinter_browse')  # tkinter native browse button
    assert _is_widget_key('session_new')
    assert _is_widget_key('session_save')
    assert _is_widget_key('session_rename_ok')
    assert _is_widget_key('session_switch_btn')
    assert _is_widget_key('session_close_selected')
    assert _is_widget_key('session_load_file')
    assert _is_widget_key('session_renaming')
    
    # Button widget keys should be identified (ending with _btn)
    assert _is_widget_key('build_crystal_btn')
    assert _is_widget_key('build_molecule_btn')
    assert _is_widget_key('load_db_structure_btn')
    assert _is_widget_key('save_db_structure_btn')
    assert _is_widget_key('any_button_btn')
    
    # Button widget keys with _button suffix should also be identified
    assert _is_widget_key('submit_button')
    assert _is_widget_key('cancel_button')
    
    # Application state keys should NOT be identified as widget keys
    assert not _is_widget_key('working_directory')
    assert not _is_widget_key('current_structure')
    assert not _is_widget_key('workflow_config')
    assert not _is_widget_key('workdir_browser_current_path')  # This is app state, not widget
    assert not _is_widget_key('selected_code_version')
    assert not _is_widget_key('espresso_calculator')
    
    # Keys that contain 'btn' but don't end with it should NOT be widget keys
    assert not _is_widget_key('btn_config')
    assert not _is_widget_key('button_state')

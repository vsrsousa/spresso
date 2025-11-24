"""
Tests for session management functionality.
"""

import pytest
import tempfile
import shutil
import os
import json
from pathlib import Path

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


@pytest.fixture
def temp_session_dir():
    """Create temporary directory for session files."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


def test_get_serializable_state(mock_st):
    """Test extraction of serializable state."""
    from gui.utils.session_manager import get_serializable_state
    
    # Set up session state
    mock_st.session_state['test_string'] = 'hello'
    mock_st.session_state['test_number'] = 42
    mock_st.session_state['test_list'] = [1, 2, 3]
    mock_st.session_state['test_dict'] = {'key': 'value'}
    mock_st.session_state['_internal'] = 'should_be_excluded'
    
    # Get serializable state
    state = get_serializable_state()
    
    # Verify results
    assert 'test_string' in state
    assert state['test_string'] == 'hello'
    assert 'test_number' in state
    assert state['test_number'] == 42
    assert 'test_list' in state
    assert state['test_list'] == [1, 2, 3]
    assert 'test_dict' in state
    assert state['test_dict'] == {'key': 'value'}
    assert '_internal' not in state  # Should be excluded


def test_save_and_load_session(mock_st, temp_session_dir):
    """Test saving and loading session."""
    from gui.utils.session_manager import save_session, load_session
    
    # Set up session state
    mock_st.session_state['machine_name'] = 'test_machine'
    mock_st.session_state['qe_version'] = '7.2'
    mock_st.session_state['test_config'] = {'ecutwfc': 50, 'kpts': [4, 4, 4]}
    
    # Save session with session name
    filepath = save_session(filename='test_session.json', session_dir=temp_session_dir, session_name='Test Session')
    
    # Verify file was created
    assert os.path.exists(filepath)
    assert filepath.endswith('test_session.json')
    
    # Load session
    loaded_state, session_name = load_session(filepath)
    
    # Verify loaded state
    assert 'machine_name' in loaded_state
    assert loaded_state['machine_name'] == 'test_machine'
    assert 'qe_version' in loaded_state
    assert loaded_state['qe_version'] == '7.2'
    assert 'test_config' in loaded_state
    assert loaded_state['test_config'] == {'ecutwfc': 50, 'kpts': [4, 4, 4]}
    # Verify session name
    assert session_name == 'Test Session'


def test_save_session_with_metadata(mock_st, temp_session_dir):
    """Test that saved session includes metadata."""
    from gui.utils.session_manager import save_session
    
    mock_st.session_state['test_key'] = 'test_value'
    
    # Save session with session name
    filepath = save_session(session_dir=temp_session_dir, session_name='My Session')
    
    # Load and verify metadata
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    assert 'metadata' in data
    assert 'saved_at' in data['metadata']
    assert 'version' in data['metadata']
    assert 'session_name' in data['metadata']
    assert data['metadata']['session_name'] == 'My Session'
    assert 'state' in data


def test_restore_session(mock_st):
    """Test restoring session state."""
    from gui.utils.session_manager import restore_session
    
    # Set initial state
    mock_st.session_state['old_key'] = 'old_value'
    
    # Restore new state
    new_state = {
        'new_key': 'new_value',
        'machine_name': 'cluster1',
    }
    restore_session(new_state, clear_first=True)
    
    # Verify state was restored
    assert 'new_key' in mock_st.session_state
    assert mock_st.session_state['new_key'] == 'new_value'
    assert 'machine_name' in mock_st.session_state
    assert mock_st.session_state['machine_name'] == 'cluster1'
    # Old key should be cleared
    assert 'old_key' not in mock_st.session_state


def test_reset_session(mock_st):
    """Test resetting session state."""
    from gui.utils.session_manager import reset_session
    
    # Set up state
    mock_st.session_state['key1'] = 'value1'
    mock_st.session_state['key2'] = 'value2'
    mock_st.session_state['key3'] = 'value3'
    
    # Reset but keep key2
    reset_session(keep_keys=['key2'])
    
    # Verify results
    assert 'key1' not in mock_st.session_state
    assert 'key2' in mock_st.session_state
    assert mock_st.session_state['key2'] == 'value2'
    assert 'key3' not in mock_st.session_state


def test_reset_session_all(mock_st):
    """Test resetting all session state."""
    from gui.utils.session_manager import reset_session
    
    # Set up state
    mock_st.session_state['key1'] = 'value1'
    mock_st.session_state['key2'] = 'value2'
    mock_st.session_state['_internal'] = 'internal_value'
    
    # Reset all
    reset_session()
    
    # Verify all keys (except internal) are cleared
    assert 'key1' not in mock_st.session_state
    assert 'key2' not in mock_st.session_state
    # Internal keys should be kept
    assert '_internal' in mock_st.session_state


def test_list_sessions(temp_session_dir):
    """Test listing available sessions."""
    from gui.utils.session_manager import list_sessions
    import time
    
    # Create some test session files
    session1 = {
        'metadata': {'saved_at': '2024-01-01T10:00:00'},
        'state': {'key': 'value1'}
    }
    session2 = {
        'metadata': {'saved_at': '2024-01-02T10:00:00'},
        'state': {'key': 'value2'}
    }
    
    with open(os.path.join(temp_session_dir, 'session1.json'), 'w') as f:
        json.dump(session1, f)
    
    time.sleep(0.1)  # Small delay to ensure different timestamps
    
    with open(os.path.join(temp_session_dir, 'session2.json'), 'w') as f:
        json.dump(session2, f)
    
    # List sessions
    sessions = list_sessions(session_dir=temp_session_dir)
    
    # Verify results
    assert len(sessions) == 2
    assert all('filename' in s for s in sessions)
    assert all('path' in s for s in sessions)
    assert all('saved_at' in s for s in sessions)
    # Should be sorted by saved_at (most recent first)
    assert sessions[0]['saved_at'] >= sessions[1]['saved_at']


def test_list_sessions_empty_dir(temp_session_dir):
    """Test listing sessions in empty directory."""
    from gui.utils.session_manager import list_sessions
    
    sessions = list_sessions(session_dir=temp_session_dir)
    assert sessions == []


def test_list_sessions_nonexistent_dir():
    """Test listing sessions in non-existent directory."""
    from gui.utils.session_manager import list_sessions
    
    sessions = list_sessions(session_dir='/nonexistent/directory')
    assert sessions == []


def test_save_session_auto_extension(mock_st, temp_session_dir):
    """Test that .json extension is added automatically."""
    from gui.utils.session_manager import save_session
    
    mock_st.session_state['test'] = 'value'
    
    # Save without .json extension
    filepath = save_session(filename='test_session', session_dir=temp_session_dir)
    
    # Verify .json was added
    assert filepath.endswith('.json')
    assert os.path.exists(filepath)


def test_load_session_nonexistent_file(temp_session_dir):
    """Test loading non-existent session file."""
    from gui.utils.session_manager import load_session
    
    with pytest.raises(IOError):
        load_session(os.path.join(temp_session_dir, 'nonexistent.json'))


def test_load_session_invalid_format(temp_session_dir):
    """Test loading session with invalid format."""
    from gui.utils.session_manager import load_session
    
    # Create invalid session file (missing 'state' key)
    invalid_file = os.path.join(temp_session_dir, 'invalid.json')
    with open(invalid_file, 'w') as f:
        json.dump({'metadata': {}}, f)
    
    with pytest.raises(ValueError):
        load_session(invalid_file)


def test_session_name_in_filename(mock_st, temp_session_dir):
    """Test that session name is used as filename."""
    from gui.utils.session_manager import save_session
    
    mock_st.session_state['test_key'] = 'test_value'
    
    # Save session with custom name
    filepath = save_session(session_dir=temp_session_dir, session_name='Al_scf')
    
    # Verify filename matches session name
    assert os.path.basename(filepath) == 'Al_scf.json'
    assert os.path.exists(filepath)


def test_load_session_preserves_name(mock_st, temp_session_dir):
    """Test that loading a session preserves the session name."""
    from gui.utils.session_manager import save_session, load_session
    
    mock_st.session_state['test_key'] = 'test_value'
    
    # Save with a specific session name
    filepath = save_session(session_dir=temp_session_dir, session_name='Al_scf')
    
    # Load the session
    state, session_name = load_session(filepath)
    
    # Verify session name is preserved
    assert session_name == 'Al_scf'
    assert 'test_key' in state
    assert state['test_key'] == 'test_value'


def test_save_load_structure(mock_st, temp_session_dir):
    """Test that ASE Atoms structures are properly saved and loaded."""
    try:
        from ase import Atoms
        from gui.utils.session_manager import save_session, load_session, restore_session
    except ImportError:
        pytest.skip("ASE not available")
    
    # Create a simple structure
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    mock_st.session_state['current_structure'] = atoms
    mock_st.session_state['test_key'] = 'test_value'
    
    # Save session
    filepath = save_session(session_dir=temp_session_dir, session_name='H2_test')
    
    # Load session
    state, session_name = load_session(filepath)
    
    # Verify session name
    assert session_name == 'H2_test'
    
    # Verify structure is in state with special format
    assert 'current_structure' in state
    assert isinstance(state['current_structure'], dict)
    assert state['current_structure']['__type__'] == 'ase.Atoms'
    assert '__data__' in state['current_structure']
    
    # Test restoration
    restore_session(state, clear_first=True)
    
    # Verify structure was properly restored
    assert 'current_structure' in mock_st.session_state
    restored_atoms = mock_st.session_state['current_structure']
    assert isinstance(restored_atoms, Atoms)
    assert restored_atoms.get_chemical_formula() == 'H2'
    assert len(restored_atoms) == 2

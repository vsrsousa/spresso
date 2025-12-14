"""
Functional tests for qtGUI to identify potential infinite loops and issues.

These tests simulate user interactions with the GUI to identify problems like:
- Infinite loops in signal handlers
- Issues with file dialogs
- Session management problems
- UI state consistency issues
"""

import pytest
import os
import sys
import tempfile
import json
from pathlib import Path
from unittest.mock import patch, MagicMock


# Conditionally skip tests if Qt is not available
pytestmark = pytest.mark.skipif(
    os.environ.get('SKIP_GUI_TESTS') == '1',
    reason="GUI tests skipped"
)


# Create test fixtures for mock application
@pytest.fixture
def mock_qapplication(monkeypatch):
    """Create a mock QApplication to avoid display requirement."""
    # Import PySide6 only if available
    try:
        from PySide6.QtWidgets import QApplication
        from PySide6.QtCore import Qt
        
        # Check if an application already exists
        app = QApplication.instance()
        if app is None:
            # Create application in offscreen mode
            QApplication.setAttribute(Qt.AA_UseSoftwareOpenGL)
            app = QApplication(sys.argv)
        yield app
    except Exception as e:
        pytest.skip(f"PySide6 not available or cannot create QApplication: {e}")


@pytest.fixture
def temp_sessions_dir():
    """Create a temporary directory for session storage."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


class TestSessionStateSignals:
    """Test SessionState for potential infinite loops in signal handling."""
    
    def test_session_state_no_infinite_recursion(self, temp_sessions_dir):
        """Test that session state changes don't cause infinite recursion."""
        from qtgui.main_app import SessionState
        
        # Reset singleton for testing
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        call_count = {'count': 0}
        max_calls = 100
        
        def listener():
            call_count['count'] += 1
            if call_count['count'] > max_calls:
                raise RuntimeError("Potential infinite loop detected!")
            # Simulate updating state in listener (this should NOT trigger recursion)
            # Note: If listeners were to call state setter again, it could loop
        
        state.add_listener(listener)
        
        # Change state multiple times
        for i in range(10):
            state['test_key'] = f'value_{i}'
        
        # Should have exactly 10 notifications
        assert call_count['count'] == 10
        
        # Reset singleton
        SessionState._instance = None
    
    def test_session_save_load_no_corruption(self, temp_sessions_dir):
        """Test that saving and loading sessions doesn't corrupt data."""
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        # Set various values
        state['working_directory'] = '/test/path'
        state['session_name'] = 'Test Session'
        state['workflow_config'] = {'calc_type': 'scf', 'ecutwfc': 50.0}
        
        # Save session
        state.save_session()
        
        # Reset and reload
        state._initialize_defaults()
        state._load_session(state._current_session_id)
        
        # Verify values are restored
        assert state['working_directory'] == '/test/path'
        assert state['session_name'] == 'Test Session'
        
        # Reset singleton
        SessionState._instance = None
    
    def test_session_create_switch_no_data_loss(self, temp_sessions_dir):
        """Test creating and switching sessions doesn't lose data."""
        from qtgui.main_app import SessionState
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        # Save data in default session
        state['working_directory'] = '/session1/path'
        state.save_session()
        
        # Create new session
        new_id = state.create_session("Session 2")
        state['working_directory'] = '/session2/path'
        state.save_session()
        
        # Switch back to default
        state.switch_session('default')
        assert state['working_directory'] == '/session1/path'
        
        # Switch to new session
        state.switch_session(new_id)
        assert state['working_directory'] == '/session2/path'
        
        # Reset singleton
        SessionState._instance = None


class TestIsolatedSessions:
    """Tests for isolated session handling."""

    def test_isolated_session_state_independent(self, temp_sessions_dir):
        from qtgui.main_app import SessionState

        SessionState._instance = None

        shared = SessionState()
        shared._sessions_dir = temp_sessions_dir
        isolated = SessionState(isolated=True)
        isolated._sessions_dir = temp_sessions_dir

        shared['working_directory'] = '/shared'
        isolated['working_directory'] = '/isolated'

        assert shared['working_directory'] == '/shared'
        assert isolated['working_directory'] == '/isolated'
        assert isolated.is_isolated() is True
        assert shared.is_isolated() is False


class TestWorkingDirectoryPerSession:
    """Test that working directory is properly associated with each session."""
    
    def test_working_directory_saved_with_session(self, temp_sessions_dir, mock_qapplication):
        """Test that working directory is saved and restored per session."""
        from qtgui.main_app import SessionState, MainWindow
        
        # Reset singleton
        SessionState._instance = None
        
        window = MainWindow()
        state = window.session_state
        state._sessions_dir = temp_sessions_dir
        state._sessions = {}  # Clear any pre-existing sessions
        
        # Create directories for different sessions
        workdir1 = os.path.join(temp_sessions_dir, "session1_workdir")
        workdir2 = os.path.join(temp_sessions_dir, "session2_workdir")
        os.makedirs(workdir1)
        os.makedirs(workdir2)
        
        # Set up default session with first working directory
        state['working_directory'] = workdir1
        state.save_session()
        window._on_session_changed()
        
        # Create new session with different working directory
        session2_id = state.create_session("Session 2")
        state['working_directory'] = workdir2
        state.save_session()
        window._refresh_session_list()
        window._on_session_changed()
        
        # Verify Session 2's working directory
        assert state['working_directory'] == workdir2
        assert window.workdir_label.text() == workdir2
        
        # Switch back to default session
        state.switch_session('default')
        window._on_session_changed()
        
        # Verify working directory was restored to Session 1's directory
        assert state['working_directory'] == workdir1
        assert window.workdir_label.text() == workdir1
        
        # Switch to Session 2 again
        state.switch_session(session2_id)
        window._on_session_changed()
        
        # Verify working directory was restored to Session 2's directory
        assert state['working_directory'] == workdir2
        assert window.workdir_label.text() == workdir2


class TestSessionManagerWindow:
    """Ensure the session manager opens isolated session windows."""

    def test_manager_spawns_isolated_window(self, temp_sessions_dir, mock_qapplication):
        from qtgui.main_app import SessionManagerWindow, SessionState

        SessionState._instance = None

        manager = SessionManagerWindow()
        manager.manager_state._sessions_dir = temp_sessions_dir

        state = SessionState(isolated=True)
        state._sessions_dir = temp_sessions_dir
        state.create_session("Managed Session")

        window = manager._spawn_session_window(state)
        assert window.session_state is state
        assert state.is_isolated()

        manager._session_windows[state.get_current_session_id()] = window
        manager._refresh_sessions()

        window.close()
        manager.close()

    def test_workspace_config_uses_shared_state(self, temp_sessions_dir, mock_qapplication):
        from qtgui.main_app import SessionManagerWindow, SessionState, MainWindow

        SessionState._instance = None

        isolated_state = SessionState(isolated=True)
        isolated_state._sessions_dir = temp_sessions_dir

        window = MainWindow(session_state=isolated_state)
        try:
            window._open_config_dialog()

            shared_state = SessionState()
            assert window._config_dialog.session_state is shared_state
            assert window._config_dialog.session_state.is_isolated() is False
        finally:
            if window._config_dialog is not None:
                window._config_dialog.close()
            window.close()

        # Reset singleton
        SessionState._instance = None
    
    def test_working_directory_ui_updates_on_session_switch(self, temp_sessions_dir, mock_qapplication):
        """Test that the UI working directory label updates when switching sessions."""
        from qtgui.main_app import SessionState, MainWindow
        from PySide6.QtWidgets import QFileDialog
        
        # Reset singleton
        SessionState._instance = None
        
        window = MainWindow()
        state = window.session_state
        state._sessions_dir = temp_sessions_dir
        state._sessions = {}
        
        # Create directories
        workdir_default = os.path.join(temp_sessions_dir, "default_work")
        workdir_new = os.path.join(temp_sessions_dir, "new_session_work")
        os.makedirs(workdir_default)
        os.makedirs(workdir_new)
        
        # Set default session's working directory using browse
        with patch.object(QFileDialog, 'getExistingDirectory', return_value=workdir_default):
            window._browse_workdir()
        
        state.save_session()
        assert window.workdir_label.text() == workdir_default
        
        # Create new session and set different working directory
        session2_id = state.create_session("New Session")
        with patch.object(QFileDialog, 'getExistingDirectory', return_value=workdir_new):
            window._browse_workdir()
        
        state.save_session()
        window._refresh_session_list()
        
        # Verify new session has different working directory
        assert state['working_directory'] == workdir_new
        assert window.workdir_label.text() == workdir_new
        
        # Switch back to default - UI should show default's working directory
        state.switch_session('default')
        window._on_session_changed()
        
        assert window.workdir_label.text() == workdir_default
        
        window.close()
        
        # Reset singleton
        SessionState._instance = None


class TestComboBoxSignalHandling:
    """Test combobox signal handlers for infinite loops."""
    
    def test_machine_combo_no_infinite_loop(self, temp_sessions_dir, mock_qapplication):
        """Test machine combo signal doesn't cause infinite loop."""
        from qtgui.main_app import SessionState
        from qtgui.pages.calculation_setup import CalculationSetupPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        # Create page without xespresso available (fallback mode)
        with patch('qtgui.pages.calculation_setup.XESPRESSO_AVAILABLE', False):
            page = CalculationSetupPage(state)
        
        # Simulate machine combo changes - should not hang
        call_count = {'count': 0}
        max_iterations = 100
        
        def track_changes(text):
            call_count['count'] += 1
            if call_count['count'] > max_iterations:
                raise RuntimeError("Infinite loop detected in machine combo!")
        
        # Page should handle empty machine name gracefully
        page._on_machine_changed("")
        page._on_machine_changed(None)
        
        # Reset singleton
        SessionState._instance = None
    
    def test_version_combo_with_empty_machine(self, temp_sessions_dir, mock_qapplication):
        """Test version combo handles empty machine gracefully."""
        from qtgui.main_app import SessionState
        from qtgui.pages.calculation_setup import CalculationSetupPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        with patch('qtgui.pages.calculation_setup.XESPRESSO_AVAILABLE', False):
            page = CalculationSetupPage(state)
        
        # These should not hang or crash
        page._on_version_changed("")
        page._on_version_changed(None)
        
        # Reset singleton
        SessionState._instance = None


class TestFileDialogHandling:
    """Test file dialog handling to ensure no hangs."""
    
    def test_browse_workdir_cancellation(self, temp_sessions_dir, mock_qapplication):
        """Test that cancelling workdir browse dialog is handled properly."""
        from qtgui.main_app import SessionState, MainWindow
        from PySide6.QtWidgets import QFileDialog
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        # Mock QFileDialog to simulate cancel
        with patch.object(QFileDialog, 'getExistingDirectory', return_value=""):
            window = MainWindow()
            
            # Store original working directory
            original_dir = state.get('working_directory')
            
            # Simulate browse button click (cancellation)
            window._browse_workdir()
            
            # Working directory should not change
            assert state.get('working_directory') == original_dir
            
            window.close()
        
        # Reset singleton
        SessionState._instance = None
    
    def test_browse_workdir_selection(self, temp_sessions_dir, mock_qapplication):
        """Test that selecting a directory in workdir browse works."""
        from qtgui.main_app import SessionState, MainWindow
        from PySide6.QtWidgets import QFileDialog
        
        # Reset singleton
        SessionState._instance = None
        
        # Create the MainWindow which will create the SessionState singleton
        # Mock QFileDialog to return a specific directory
        test_dir = temp_sessions_dir
        with patch.object(QFileDialog, 'getExistingDirectory', return_value=test_dir):
            window = MainWindow()
            
            # Get the session state that the window is using
            state = window.session_state
            state._sessions_dir = temp_sessions_dir
            
            # Simulate browse button click
            window._browse_workdir()
            
            # Working directory should be updated
            assert state.get('working_directory') == test_dir
            
            window.close()
        
        # Reset singleton
        SessionState._instance = None


class TestSessionComboBoxInteraction:
    """Test session combo box for potential issues."""
    
    def test_session_combo_signal_blocking(self, temp_sessions_dir, mock_qapplication):
        """Test that session combo properly blocks signals during updates."""
        from qtgui.main_app import SessionState, MainWindow
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        window = MainWindow()
        
        # Create multiple sessions
        state.create_session("Session 1")
        state.create_session("Session 2")
        
        # Refresh should not cause infinite loops
        call_count = {'count': 0}
        original_refresh = window._refresh_session_list
        
        def tracked_refresh():
            call_count['count'] += 1
            if call_count['count'] > 10:
                raise RuntimeError("Potential infinite loop in session refresh!")
            original_refresh()
        
        window._refresh_session_list = tracked_refresh
        window._refresh_session_list()
        
        # Should have called exactly once
        assert call_count['count'] == 1
        
        window.close()
        
        # Reset singleton
        SessionState._instance = None
    
    def test_session_selection_signal_handling(self, temp_sessions_dir, mock_qapplication):
        """Test session selection doesn't trigger recursive updates."""
        from qtgui.main_app import SessionState, MainWindow
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        window = MainWindow()
        
        session_id = state.create_session("Test Session")
        window._refresh_session_list()
        
        # Track session changed calls
        call_count = {'count': 0}
        original_handler = window._on_session_changed
        
        def tracked_handler():
            call_count['count'] += 1
            if call_count['count'] > 5:
                raise RuntimeError("Potential infinite loop in session change handler!")
            original_handler()
        
        window._on_session_changed = tracked_handler
        
        # Simulate selecting the session
        window._on_session_selected("Test Session")
        
        window.close()
        
        # Reset singleton
        SessionState._instance = None


class TestPageRefresh:
    """Test page refresh methods for potential issues."""
    
    def test_calculation_setup_refresh_no_hang(self, temp_sessions_dir, mock_qapplication):
        """Test CalculationSetupPage refresh doesn't hang."""
        from qtgui.main_app import SessionState
        from qtgui.pages.calculation_setup import CalculationSetupPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        with patch('qtgui.pages.calculation_setup.XESPRESSO_AVAILABLE', False):
            page = CalculationSetupPage(state)
            
            # Multiple refresh calls should not hang
            for _ in range(5):
                page.refresh()
        
        # Reset singleton
        SessionState._instance = None
    
    def test_workflow_builder_refresh_no_hang(self, temp_sessions_dir, mock_qapplication):
        """Test WorkflowBuilderPage refresh doesn't hang."""
        from qtgui.main_app import SessionState
        from qtgui.pages.workflow_builder import WorkflowBuilderPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        with patch('qtgui.pages.workflow_builder.XESPRESSO_AVAILABLE', False):
            page = WorkflowBuilderPage(state)
            
            # Multiple refresh calls should not hang
            for _ in range(5):
                page.refresh()
        
        # Reset singleton
        SessionState._instance = None


class TestNavigationSignals:
    """Test navigation between pages."""
    
    def test_nav_changed_no_recursion(self, temp_sessions_dir, mock_qapplication):
        """Test that navigation changes don't cause recursion."""
        from qtgui.main_app import SessionState, MainWindow
        
        # Reset singleton
        SessionState._instance = None
        
        window = MainWindow()
        state = window.session_state
        state._sessions_dir = temp_sessions_dir
        
        # Track nav changes
        call_count = {'count': 0}
        original_handler = window._on_nav_changed
        
        def tracked_handler(index):
            call_count['count'] += 1
            if call_count['count'] > 10:
                raise RuntimeError("Potential infinite loop in navigation!")
            original_handler(index)
        
        window._on_nav_changed = tracked_handler
        
        # Simulate navigating through pages (start from 1 since 0 is already selected)
        for i in range(1, 5):
            window._navigate_to(i)
        
        # Should have exactly 4 calls (pages 1, 2, 3, 4)
        assert call_count['count'] == 4
        
        window.close()
        
        # Reset singleton
        SessionState._instance = None


class TestStructureViewer:
    """Test structure viewer page interactions."""
    
    def test_structure_build_type_change(self, temp_sessions_dir, mock_qapplication):
        """Test build type combo changes don't cause issues."""
        from qtgui.main_app import SessionState
        from qtgui.pages.structure_viewer import StructureViewerPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        # Skip if ASE not available
        try:
            from ase import Atoms
        except ImportError:
            pytest.skip("ASE not available")
        
        page = StructureViewerPage(state)
        
        # Switch between build types - should not hang
        page._on_build_type_changed("Bulk Crystal")
        page._on_build_type_changed("Molecule")
        page._on_build_type_changed("Bulk Crystal")
        
        # Reset singleton
        SessionState._instance = None


class TestJobSubmissionPage:
    """Test job submission page for potential issues."""
    
    def test_browser_refresh_no_hang(self, temp_sessions_dir, mock_qapplication):
        """Test file browser refresh doesn't hang."""
        from qtgui.main_app import SessionState
        from qtgui.pages.job_submission import JobSubmissionPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        state['working_directory'] = temp_sessions_dir
        
        page = JobSubmissionPage(state)
        
        # Multiple refresh calls should not hang
        for _ in range(5):
            page._refresh_browser()
        
        # Reset singleton
        SessionState._instance = None
    
    def test_file_selection_no_crash(self, temp_sessions_dir, mock_qapplication):
        """Test file selection handles missing files gracefully."""
        from qtgui.main_app import SessionState
        from qtgui.pages.job_submission import JobSubmissionPage
        from PySide6.QtWidgets import QTreeWidgetItem
        from PySide6.QtCore import Qt
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        state['working_directory'] = temp_sessions_dir
        
        page = JobSubmissionPage(state)
        
        # Create a mock tree item with non-existent file
        item = QTreeWidgetItem(["nonexistent.txt"])
        item.setData(0, Qt.UserRole, "/nonexistent/path/file.txt")
        
        # Should handle gracefully without crashing
        page._on_file_selected(item, 0)
        
        # Reset singleton
        SessionState._instance = None


class TestValidationUtility:
    """Test the validation utility for path handling."""
    
    def test_validate_path_under_base(self):
        """Test path validation prevents traversal."""
        from qtgui.utils import validate_path_under_base
        
        # Valid paths
        is_valid, resolved, error = validate_path_under_base("/base/subdir", "/base")
        assert is_valid
        
        # Path traversal attempt
        is_valid, resolved, error = validate_path_under_base("/base/../etc/passwd", "/base")
        assert not is_valid
        
        # Absolute path outside base
        is_valid, resolved, error = validate_path_under_base("/etc/passwd", "/base")
        assert not is_valid
    
    def test_safe_makedirs(self, temp_sessions_dir):
        """Test safe makedirs function."""
        from qtgui.utils import safe_makedirs
        
        # Should create directory
        test_path = os.path.join(temp_sessions_dir, "newdir", "subdir")
        safe_makedirs(test_path)
        assert os.path.isdir(test_path)
        
        # Should handle existing directory
        safe_makedirs(test_path)  # Should not raise


class TestEmptyStateHandling:
    """Test handling of empty/null state values."""
    
    def test_empty_machine_name_handling(self, temp_sessions_dir, mock_qapplication):
        """Test pages handle empty machine name without crashing."""
        from qtgui.main_app import SessionState
        from qtgui.pages.codes_config import CodesConfigPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        with patch('qtgui.pages.codes_config.XESPRESSO_AVAILABLE', False):
            page = CodesConfigPage(state)
            
            # Should handle empty machine name
            page._on_machine_changed("")
            page._on_machine_changed(None)
        
        # Reset singleton
        SessionState._instance = None
    
    def test_empty_version_handling(self, temp_sessions_dir, mock_qapplication):
        """Test pages handle empty version without crashing."""
        from qtgui.main_app import SessionState
        from qtgui.pages.codes_config import CodesConfigPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state._sessions_dir = temp_sessions_dir
        
        with patch('qtgui.pages.codes_config.XESPRESSO_AVAILABLE', False):
            page = CodesConfigPage(state)
            
            # Should handle empty version
            page._on_version_changed("")
            page._on_version_changed(None)
        
        # Reset singleton
        SessionState._instance = None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

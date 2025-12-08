"""
Test file browser functionality for qtgui job submission page.

This test verifies that the file browser correctly displays:
1. All top-level directories in the working directory
2. Subdirectories containing calculation files
3. Calculation files themselves
"""

import os
import sys
import tempfile
import pytest
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestFileBrowser:
    """Test file browser visibility and functionality."""
    
    @pytest.fixture
    def test_workdir(self):
        """Create a test working directory structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create directory structure similar to user's problem
            workdir = Path(tmpdir) / "workdir"
            workdir.mkdir()
            
            # Top-level directories
            (workdir / "AlGdNi4" / "subdir").mkdir(parents=True)
            (workdir / "CuGd" / "calculation").mkdir(parents=True)
            (workdir / "scf" / "AlGdNi4").mkdir(parents=True)
            (workdir / "temp_AlGdNi4").mkdir()
            (workdir / "workflow" / "step1").mkdir(parents=True)
            
            # Add calculation files in subdirectories
            (workdir / "AlGdNi4" / "subdir" / "test.pwi").write_text("test input")
            (workdir / "CuGd" / "calculation" / "calc.pwi").write_text("test input")
            (workdir / "scf" / "AlGdNi4" / "scf.pwi").write_text("test input")
            (workdir / "workflow" / "step1" / "workflow.pwi").write_text("test input")
            
            # Add a file at root level
            (workdir / "nohup.out").write_text("nohup output")
            
            yield str(workdir)
    
    def test_all_top_level_directories_visible(self, test_workdir, mock_qapplication):
        """Test that all top-level directories are visible in the browser."""
        from qtgui.main_app import SessionState
        from qtgui.pages.job_submission import JobSubmissionPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state['working_directory'] = test_workdir
        
        page = JobSubmissionPage(state)
        page._refresh_browser()
        
        # Get all top-level items in the tree
        top_level_items = []
        for i in range(page.file_tree.topLevelItemCount()):
            item = page.file_tree.topLevelItem(i)
            top_level_items.append(item.text(0))
        
        # Verify all directories are visible
        assert "AlGdNi4/" in top_level_items, "AlGdNi4/ should be visible"
        assert "CuGd/" in top_level_items, "CuGd/ should be visible"
        assert "scf/" in top_level_items, "scf/ should be visible"
        assert "temp_AlGdNi4/" in top_level_items, "temp_AlGdNi4/ should be visible"
        assert "workflow/" in top_level_items, "workflow/ should be visible"
        
        # Reset singleton
        SessionState._instance = None
    
    def test_subdirectories_with_files_shown(self, test_workdir, mock_qapplication):
        """Test that subdirectories containing calculation files are shown."""
        from qtgui.main_app import SessionState
        from qtgui.pages.job_submission import JobSubmissionPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state['working_directory'] = test_workdir
        
        page = JobSubmissionPage(state)
        page._refresh_browser()
        
        # Find the AlGdNi4 directory item
        algdni4_item = None
        for i in range(page.file_tree.topLevelItemCount()):
            item = page.file_tree.topLevelItem(i)
            if item.text(0) == "AlGdNi4/":
                algdni4_item = item
                break
        
        assert algdni4_item is not None, "AlGdNi4/ directory should be in tree"
        assert algdni4_item.childCount() > 0, "AlGdNi4/ should have child items"
        
        # Check that subdir is shown as a child
        child_names = [algdni4_item.child(i).text(0) for i in range(algdni4_item.childCount())]
        assert "subdir" in child_names, "subdir should be shown as child of AlGdNi4/"
        
        # Reset singleton
        SessionState._instance = None
    
    def test_empty_directories_shown(self, test_workdir, mock_qapplication):
        """Test that empty directories (with no calculation files) are still shown."""
        from qtgui.main_app import SessionState
        from qtgui.pages.job_submission import JobSubmissionPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state['working_directory'] = test_workdir
        
        page = JobSubmissionPage(state)
        page._refresh_browser()
        
        # Get all top-level items
        top_level_items = []
        for i in range(page.file_tree.topLevelItemCount()):
            item = page.file_tree.topLevelItem(i)
            top_level_items.append(item.text(0))
        
        # temp_AlGdNi4 has no calculation files but should still be visible
        assert "temp_AlGdNi4/" in top_level_items, "temp_AlGdNi4/ (empty dir) should be visible"
        
        # Reset singleton
        SessionState._instance = None
    
    def test_calculation_files_shown(self, test_workdir, mock_qapplication):
        """Test that calculation files are shown under their parent directories."""
        from qtgui.main_app import SessionState
        from qtgui.pages.job_submission import JobSubmissionPage
        
        # Reset singleton
        SessionState._instance = None
        
        state = SessionState()
        state['working_directory'] = test_workdir
        
        page = JobSubmissionPage(state)
        page._refresh_browser()
        
        # Find the scf directory
        scf_item = None
        for i in range(page.file_tree.topLevelItemCount()):
            item = page.file_tree.topLevelItem(i)
            if item.text(0) == "scf/":
                scf_item = item
                break
        
        assert scf_item is not None, "scf/ directory should be in tree"
        assert scf_item.childCount() > 0, "scf/ should have children"
        
        # Find AlGdNi4 subdirectory under scf
        algdni4_child = None
        for i in range(scf_item.childCount()):
            child = scf_item.child(i)
            if child.text(0) == "AlGdNi4":
                algdni4_child = child
                break
        
        assert algdni4_child is not None, "AlGdNi4 should be under scf/"
        
        # Check that .pwi file is shown
        if algdni4_child.childCount() > 0:
            file_names = [algdni4_child.child(i).text(0) for i in range(algdni4_child.childCount())]
            assert any('.pwi' in name for name in file_names), "At least one .pwi file should be shown"
        
        # Reset singleton
        SessionState._instance = None


@pytest.fixture
def mock_qapplication():
    """Mock QApplication for testing."""
    from PySide6.QtWidgets import QApplication
    import sys
    
    # Check if QApplication already exists
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    
    yield app
    
    # Don't quit the application as it may be shared


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

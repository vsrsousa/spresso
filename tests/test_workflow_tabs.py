import pytest
import sys
from qtgui.main_app import SessionState, MainWindow


@pytest.fixture
def mock_qapplication(monkeypatch):
    try:
        from PySide6.QtWidgets import QApplication
        from PySide6.QtCore import Qt
        app = QApplication.instance()
        if app is None:
            QApplication.setAttribute(Qt.AA_UseSoftwareOpenGL)
            app = QApplication(sys.argv)
        yield app
    except Exception as e:
        pytest.skip(f"PySide6 not available or cannot create QApplication: {e}")


def test_independent_workflow_tabs(mock_qapplication):
    # Reset singleton
    SessionState._instance = None
    state = SessionState(isolated=True)

    w = MainWindow(session_state=state)

    # Launch the same workflow twice
    w.launch_workflow('SCF')
    w.launch_workflow('SCF')

    # Find the WorkflowTabsPage and assert only one SCF tab exists
    tabs_page = None
    for i in range(w.content_stack.count()):
        widget = w.content_stack.widget(i)
        if widget.__class__.__name__ == 'WorkflowTabsPage':
            tabs_page = widget
            break

    assert tabs_page is not None
    assert tabs_page.tabs.count() == 1
    assert tabs_page.tabs.tabText(0) == 'SCF'

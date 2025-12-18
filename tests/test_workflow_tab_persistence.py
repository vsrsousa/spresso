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


def test_workflow_tab_persistence(mock_qapplication):
    # Reset singleton
    SessionState._instance = None
    state = SessionState(isolated=True)

    w = MainWindow(session_state=state)

    # Launch SCF and set some values
    w.launch_workflow('SCF')

    # Find the WorkflowTabsPage
    tabs_page = None
    for i in range(w.content_stack.count()):
        widget = w.content_stack.widget(i)
        if widget.__class__.__name__ == 'WorkflowTabsPage':
            tabs_page = widget
            break

    assert tabs_page is not None
    page = tabs_page.tabs.widget(0)

    # Set values and explicitly save draft
    page.machine_edit.setText('medusa')
    page.pseudo_edit.setText('Fe=Fe.pseudo')
    page.proto_combo.setCurrentText('fast')
    page.magnetism_chk.setChecked(True)
    page._save_draft()

    # Close the tab
    tabs_page._on_tab_close_requested(0)

    # Re-open SCF
    w.launch_workflow('SCF')

    # Find new page and assert values restored
    tabs_page = None
    for i in range(w.content_stack.count()):
        widget = w.content_stack.widget(i)
        if widget.__class__.__name__ == 'WorkflowTabsPage':
            tabs_page = widget
            break

    assert tabs_page is not None
    page2 = tabs_page.tabs.widget(0)
    assert page2.machine_edit.text() == 'medusa'
    assert 'Fe=Fe.pseudo' in page2.pseudo_edit.text()
    assert page2.proto_combo.currentText() == 'fast'
    assert page2.magnetism_chk.isChecked() is True

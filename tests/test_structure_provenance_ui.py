import pytest
from ase import Atoms

from qtgui.pages.structure_viewer import StructureViewerPage


def test_structure_viewer_provenance_ui(monkeypatch, tmp_path):
    # Minimal session_state dict-like
    session_state = {}
    # Ensure a QApplication exists for Qt widgets in headless test environments
    app_created = False
    try:
        from PySide6.QtWidgets import QApplication
        if QApplication.instance() is None:
            _app = QApplication([])
            app_created = True
    except Exception:
        # If Qt isn't available or cannot create an application, skip GUI instantiation
        pytest.skip("Qt not available in this environment")

    page = StructureViewerPage(session_state, view_only=True)

    # Build a simple H2 molecule
    atoms = Atoms('H2', positions=[[0,0,0],[0,0,0.74]])

    # Set structure and ensure provenance label updates (no exceptions)
    page._set_structure(atoms, 'test')
    assert hasattr(page, 'provenance_count_label')
    # Text should be string (may be 'No provenance records' or count)
    assert isinstance(page.provenance_count_label.text(), str)

    if app_created:
        try:
            _app.quit()
        except Exception:
            pass

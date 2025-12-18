import os
import json
from pathlib import Path

from qtgui.pages.workflow_builder import WorkflowBuilderPage
from qtgui.session.state_manager import SessionState


def test_save_and_load_workflow(tmp_path, monkeypatch):
    # Ensure QApplication exists for widget instantiation
    try:
        from PySide6.QtWidgets import QApplication
        app = QApplication.instance() or QApplication([])
    except Exception:
        app = None

    ss = SessionState(isolated=True)
    ss.create_session('test')
    page = WorkflowBuilderPage(ss)
    # set some UI fields to known values
    page.ecutwfc_spin.setValue(60.0)
    page.occupations_combo.setCurrentText('fixed')
    page.workflow_combo.setCurrentText('SCF + Relaxation')

    # save to a temp workflows dir by monkeypatching _workflows_dir
    wf_dir = tmp_path / 'workflows'
    wf_dir.mkdir()
    monkeypatch.setattr(page, '_workflows_dir', lambda: str(wf_dir))

    page._save_workflow()
    # file should be created
    files = list(wf_dir.glob('*.json'))
    assert len(files) == 1
    data = json.loads(files[0].read_text())
    assert data['config']['ecutwfc'] == 60.0
    assert data['config']['occupations'] == 'fixed'
    assert len(data['calculations']) == 2

    # now test load
    page.workflow_combo.setCurrentText('Single SCF')
    page.summary_text.clear()
    page._apply_loaded_workflow(data)
    assert page.workflow_combo.currentText() == 'SCF + Relaxation'
    assert float(page.ecutwfc_spin.value()) == 60.0
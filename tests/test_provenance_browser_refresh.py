from PySide6.QtWidgets import QApplication


def test_config_change_resets_manager_browser():
    # Ensure QApplication exists
    app = QApplication.instance() or QApplication([])
    from qtgui.main_app import SessionManagerWindow
    from qtgui.dialogs import ConfigurationDialog

    mgr = SessionManagerWindow()
    # simulate opening browser
    try:
        pb = mgr._get_provenance_browser()
    except Exception:
        pb = None
    # open manager's configuration dialog so it's the connected instance
    mgr._open_config_dialog()
    cfg = mgr._config_dialog
    # ensure precondition: manager has provenance browser attribute
    assert hasattr(mgr, '_provenance_browser')
    # emit change and verify the manager's browser reference becomes None
    # emit the signal and verify the manager's browser reference becomes None
    try:
        cfg.configuration_changed.emit()
    except Exception:
        pass
    assert mgr._provenance_browser is None

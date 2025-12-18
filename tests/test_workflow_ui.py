def test_workflow_ui_imports():
    # Ensure UI placeholders import cleanly in headless/test environment
    from qtgui.pages.workflow_tasks.palette import TaskPalette
    from qtgui.pages.workflow_tasks.editor import TaskEditor
    from qtgui.pages.workflow_tasks.canvas import WorkflowCanvas

    # If Qt is available we must ensure a QApplication exists before
    # instantiating widgets; otherwise the module-level placeholders
    # set QWidget to 'object' and instantiation is a no-op.
    try:
        from PySide6.QtWidgets import QApplication
        app = QApplication.instance() or QApplication([])
    except Exception:
        app = None

    # Instantiate (should not raise)
    p = TaskPalette()
    e = TaskEditor()
    c = WorkflowCanvas()

    # Basic sanity checks
    assert p is not None
    assert e is not None
    assert c is not None

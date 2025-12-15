"""Toolbar helper."""
from PySide6.QtWidgets import QToolBar, QWidget
from PySide6.QtWidgets import QSizePolicy
from PySide6.QtGui import QAction


def setup_toolbar(main_window):
    toolbar = QToolBar("Main Toolbar")
    toolbar.setMovable(False)
    main_window.addToolBar(toolbar)
    # Minimal toolbar for session windows: keep only session-local actions.
    if getattr(main_window, 'is_session_window', False):
        # Provide quick access to configuration when in a session window
        config_action = QAction("⚙️ Configuration", main_window)
        config_action.setToolTip("Open configuration dialog")
        config_action.triggered.connect(main_window._open_config_dialog)
        toolbar.addAction(config_action)

        rename_action = QAction("✏️ Rename", main_window)
        rename_action.setToolTip("Rename current session")
        rename_action.triggered.connect(main_window._rename_session)
        toolbar.addAction(rename_action)

        # Toggle sidebar visibility (collapsible sidebar)
        toggle_action = QAction("⇔ Toggle Sidebar", main_window)
        toggle_action.setToolTip("Show/hide the left sidebar")
        toggle_action.triggered.connect(getattr(main_window, '_toggle_sidebar', lambda: None))
        toolbar.addAction(toggle_action)

        toolbar.addSeparator()

        job_monitor_action = QAction("🔍 Job Monitor", main_window)
        job_monitor_action.setToolTip("Track remote job submissions")
        job_monitor_action.triggered.connect(main_window._open_job_monitor)
        toolbar.addAction(job_monitor_action)
    else:
        # Full toolbar for manager/launcher windows
        config_action = QAction("⚙️ Configuration", main_window)
        config_action.setToolTip("Open configuration dialog")
        config_action.triggered.connect(main_window._open_config_dialog)
        toolbar.addAction(config_action)
        toolbar.addSeparator()

        # New and Load are manager-level actions and should appear here
        new_action = QAction("📋 New Session", main_window)
        new_action.setToolTip("Create a new session (opens a new window)")
        new_action.triggered.connect(main_window._new_session)
        toolbar.addAction(new_action)

        load_action = QAction("📂 Load Session", main_window)
        load_action.setToolTip("Load a saved session (opens a new window)")
        load_action.triggered.connect(main_window._load_session_file)
        toolbar.addAction(load_action)

        toolbar.addSeparator()

        rename_action = QAction("✏️ Rename", main_window)
        rename_action.setToolTip("Rename current session")
        rename_action.triggered.connect(main_window._rename_session)
        toolbar.addAction(rename_action)

        toggle_action = QAction("⇔ Toggle Sidebar", main_window)
        toggle_action.setToolTip("Show/hide the left sidebar")
        toggle_action.triggered.connect(getattr(main_window, '_toggle_sidebar', lambda: None))
        toolbar.addAction(toggle_action)

    toolbar.addSeparator()
    # Add Job Monitor for manager windows here; session windows already
    # include a job monitor action above to avoid duplication.
    if not getattr(main_window, 'is_session_window', False):
        job_monitor_action = QAction("🔍 Job Monitor", main_window)
        job_monitor_action.setToolTip("Track remote job submissions")
        job_monitor_action.triggered.connect(main_window._open_job_monitor)
        toolbar.addAction(job_monitor_action)

    spacer = QWidget()
    spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
    toolbar.addWidget(spacer)

    # Quit only present on non-session manager windows
    if not getattr(main_window, 'is_session_window', False):
        quit_action = QAction("🚪 Quit", main_window)
        quit_action.setToolTip("Quit the application")
        quit_action.triggered.connect(main_window.close)
        toolbar.addAction(quit_action)

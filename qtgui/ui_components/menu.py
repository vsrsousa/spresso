"""Menu bar helper."""
from qtpy.QtGui import QAction


def setup_menu(main_window):
    menubar = main_window.menuBar()

    file_menu = menubar.addMenu("&File")
    new_session_action = QAction("New Session", main_window)
    new_session_action.setShortcut("Ctrl+N")
    new_session_action.triggered.connect(main_window._new_session)
    file_menu.addAction(new_session_action)

    save_session_action = QAction("Save Session", main_window)
    save_session_action.setShortcut("Ctrl+S")
    save_session_action.triggered.connect(main_window._save_session)
    file_menu.addAction(save_session_action)

    file_menu.addSeparator()
    exit_action = QAction("Exit", main_window)
    exit_action.setShortcut("Ctrl+Q")
    exit_action.triggered.connect(main_window.close)
    file_menu.addAction(exit_action)

    edit_menu = menubar.addMenu("&Edit")
    config_action = QAction("Configuration...", main_window)
    config_action.setShortcut("Ctrl+,")
    config_action.triggered.connect(main_window._open_config_dialog)
    edit_menu.addAction(config_action)

    view_menu = menubar.addMenu("&View")
    nav_actions = [
        ("Calculation Setup", "Ctrl+1"),
        ("Workflow Builder", "Ctrl+2"),
        ("Job Submission", "Ctrl+3"),
        ("Results", "Ctrl+4"),
    ]
    for i, (name, shortcut) in enumerate(nav_actions):
        action = QAction(name, main_window)
        action.setShortcut(shortcut)
        action.triggered.connect(lambda checked, idx=i: main_window._navigate_to(idx))
        view_menu.addAction(action)

    help_menu = menubar.addMenu("&Help")
    about_action = QAction("About", main_window)
    about_action.triggered.connect(main_window._show_about)
    help_menu.addAction(about_action)

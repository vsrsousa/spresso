"""Sidebar UI helper."""
from PySide6.QtWidgets import (
    QGroupBox, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QListWidget, QListWidgetItem, QWidget, QComboBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont


def create_sidebar(main_window):
    """Create sidebar widget and populate navigation elements."""
    sidebar = QWidget()
    sidebar.setMaximumWidth(300)
    sidebar.setMinimumWidth(200)
    layout = QVBoxLayout(sidebar)
    layout.setContentsMargins(5, 5, 5, 5)

    # Title
    title_label = QLabel("⚛️ xespresso")
    title_font = QFont()
    title_font.setPointSize(16)
    title_font.setBold(True)
    title_label.setFont(title_font)
    title_label.setAlignment(Qt.AlignCenter)
    layout.addWidget(title_label)

    # (Session section is added below; avoid duplication)

    # Session section
    session_group = QGroupBox("📋 Session")
    session_layout = QVBoxLayout(session_group)

    # Manager windows keep the Active selector; session windows show static name/ID
    if not getattr(main_window, 'is_session_window', False):
        session_selector_layout = QHBoxLayout()
        session_selector_layout.addWidget(QLabel("Active:"))
        main_window.session_combo = QComboBox()
        main_window.session_combo.setToolTip("Active session (create new or load saved session)")
        main_window.session_combo.currentTextChanged.connect(main_window._on_session_selected)
        session_selector_layout.addWidget(main_window.session_combo, 1)
        session_layout.addLayout(session_selector_layout)

        main_window.session_name_label = QLabel("")
        main_window.session_name_label.setWordWrap(True)
        main_window._update_session_name_label()
        session_layout.addWidget(main_window.session_name_label)

        # Session buttons for manager: New / Load
        session_btn_layout1 = QHBoxLayout()
        new_session_btn = QPushButton("New")
        new_session_btn.clicked.connect(main_window._new_session)
        session_btn_layout1.addWidget(new_session_btn)
        load_session_btn = QPushButton("Load")
        load_session_btn.clicked.connect(main_window._load_session_dialog)
        session_btn_layout1.addWidget(load_session_btn)
        session_layout.addLayout(session_btn_layout1)

        session_btn_layout2 = QHBoxLayout()
        rename_session_btn = QPushButton("Rename")
        rename_session_btn.clicked.connect(main_window._rename_session)
        session_btn_layout2.addWidget(rename_session_btn)
        save_session_btn = QPushButton("Save")
        save_session_btn.clicked.connect(main_window._save_session)
        session_btn_layout2.addWidget(save_session_btn)
        session_layout.addLayout(session_btn_layout2)

        close_session_btn = QPushButton("Close")
        close_session_btn.clicked.connect(main_window._close_session)
        session_layout.addWidget(close_session_btn)
    else:
        # Session window: show name/ID label, Save, Rename, Edit, and Close (which closes window)
        main_window.session_name_label = QLabel("")
        main_window.session_name_label.setWordWrap(True)
        main_window._update_session_name_label()
        session_layout.addWidget(main_window.session_name_label)

        session_btn_layout = QHBoxLayout()
        rename_session_btn = QPushButton("Rename")
        rename_session_btn.clicked.connect(main_window._rename_session)
        session_btn_layout.addWidget(rename_session_btn)
        save_session_btn = QPushButton("Save")
        save_session_btn.clicked.connect(main_window._save_session)
        session_btn_layout.addWidget(save_session_btn)
        edit_session_btn = QPushButton("Edit")
        edit_session_btn.clicked.connect(getattr(main_window, '_toggle_edit_mode', lambda: None))
        session_btn_layout.addWidget(edit_session_btn)
        session_layout.addLayout(session_btn_layout)

        close_session_btn = QPushButton("Close")
        close_session_btn.clicked.connect(getattr(main_window, 'close', lambda: None))
        session_layout.addWidget(close_session_btn)

    # Add the session group to the layout so Session appears before Structure
    layout.addWidget(session_group)

    # Structure section (session windows only)
    if getattr(main_window, 'is_session_window', False):
        struct_group = QGroupBox("🧱 Structure")
        struct_layout = QVBoxLayout(struct_group)
        # structure_label will be updated by the session window via _update_structure_label
        main_window.structure_label = QLabel(main_window.session_state.get('structure_source', 'No structure selected'))
        main_window.structure_label.setWordWrap(True)
        struct_layout.addWidget(main_window.structure_label)

        struct_btns = QHBoxLayout()
        load_struct_btn = QPushButton("Load")
        load_struct_btn.setToolTip("Load structure from file for this session (once)")
        load_struct_btn.clicked.connect(getattr(main_window, '_choose_structure_file', lambda: None))
        struct_btns.addWidget(load_struct_btn)

        db_struct_btn = QPushButton("From DB")
        db_struct_btn.setToolTip("Choose structure from ASE database for this session (once)")
        db_struct_btn.clicked.connect(getattr(main_window, '_choose_structure_db', lambda: None))
        struct_btns.addWidget(db_struct_btn)

        build_struct_btn = QPushButton("Build")
        build_struct_btn.setToolTip("Build a structure (opens builder)")
        build_struct_btn.clicked.connect(getattr(main_window, '_build_structure', lambda: None))
        struct_btns.addWidget(build_struct_btn)

        view_btn = QPushButton("View")
        view_btn.setToolTip("View structure in the main area")
        view_btn.clicked.connect(getattr(main_window, '_view_structure', lambda: None))
        struct_btns.addWidget(view_btn)

        # Expose buttons so the session window can enable/disable them
        main_window.load_struct_btn = load_struct_btn
        main_window.db_struct_btn = db_struct_btn
        main_window.build_struct_btn = build_struct_btn

        # Disable buttons if a structure is already set or locked
        try:
            has_structure = main_window.session_state.get('current_structure') is not None
            locked = main_window.session_state.is_structure_locked()
            enabled = not (has_structure or locked)
            load_struct_btn.setEnabled(enabled)
            db_struct_btn.setEnabled(enabled)
            build_struct_btn.setEnabled(enabled)
        except Exception:
            pass

        struct_layout.addLayout(struct_btns)
        layout.addWidget(struct_group)

    # Configuration removed from sidebar; configuration is available in toolbar

    # Working directory
    workdir_group = QGroupBox("📁 Working Directory")
    workdir_layout = QVBoxLayout(workdir_group)
    main_window.workdir_label = QLabel(main_window.session_state.get('working_directory', '~'))
    main_window.workdir_label.setWordWrap(True)
    workdir_layout.addWidget(main_window.workdir_label)
    browse_btn = QPushButton("Browse...")
    browse_btn.clicked.connect(main_window._browse_workdir)
    workdir_layout.addWidget(browse_btn)
    layout.addWidget(workdir_group)

    # Workflow navigation
    workflow_group = QGroupBox("🔬 Workflows")
    workflow_layout = QVBoxLayout(workflow_group)
    main_window.nav_list = QListWidget()
    # Show independent workflow presets directly in the sidebar. Selecting
    # an item will launch the chosen workflow (no Calculation Setup page).
    presets = [
        "SCF",
        "Relax",
        "SCF+Relax",
        "Convergence Test",
        "Geometry Optimization",
        "NEB",
        "Post-Processing",
    ]
    for p in presets:
        main_window.nav_list.addItem(QListWidgetItem(p))
    # Connect navigation changes to the session window handler if available
    try:
        # Use currentRowChanged so programmatic changes trigger the handler
        main_window.nav_list.currentRowChanged.connect(getattr(main_window, '_on_nav_changed', lambda idx: None))
    except Exception:
        pass
    # When a preset is activated, call the session window launcher directly.
    def _on_preset_activated(item):
        try:
            name = item.text()
            getattr(main_window, 'launch_workflow', lambda n: None)(name)
        except Exception:
            pass

    main_window.nav_list.itemActivated.connect(_on_preset_activated)
    workflow_layout.addWidget(main_window.nav_list)
    layout.addWidget(workflow_group)

    layout.addStretch()

    about_label = QLabel(
        """
<b>About</b><br>
<b>xespresso GUI</b> - PySide6 interface for Quantum ESPRESSO calculations<br>
<br>
Version: 1.2.0<br>
<a href="https://github.com/vsrsousa/spresso">Documentation</a> | 
<a href="https://github.com/vsrsousa/spresso/issues">Report Issue</a>
"""
    )
    about_label.setTextFormat(Qt.RichText)
    about_label.setOpenExternalLinks(True)
    about_label.setWordWrap(True)
    layout.addWidget(about_label)

    return sidebar

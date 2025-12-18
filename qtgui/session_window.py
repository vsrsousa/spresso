"""
Session window (MainWindow) extracted from main_app.py.

This module contains the session workspace window so the session
manager and session window responsibilities are separated.
"""
import os
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QStackedWidget, QLabel, QFileDialog, QMessageBox, QSplitter,
    QInputDialog, QLineEdit
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction

from .session_state import SessionState
from .session import controllers as session_controllers
from .ui_components import create_sidebar, setup_menu, setup_toolbar, setup_statusbar
from .ui_components.window_focus import apply_focus_decorator


# Lazy import helper for page modules (copied here to avoid circular imports)
_page_modules = {}

def _get_page_class(name):
    """Lazy import page class to improve startup time."""
    if name not in _page_modules:
        if name == 'StructureViewerPage':
            from qtgui.pages.structure_viewer import StructureViewerPage
            _page_modules[name] = StructureViewerPage
        elif name == 'CalculationSetupPage':
            # For the workflow-first prototype replace the old CalculationSetup
            # page with the new WorkflowsPage. This is a deliberate replacement
            # (no back-compat shim) while we iterate on the workflow UX.
            from qtgui.pages.workflows_page import WorkflowsPage
            _page_modules[name] = WorkflowsPage
        elif name == 'WorkflowBuilderPage':
            from qtgui.pages.workflow_builder import WorkflowBuilderPage
            _page_modules[name] = WorkflowBuilderPage
        elif name == 'JobSubmissionPage':
            from qtgui.pages.job_submission import JobSubmissionPage
            _page_modules[name] = JobSubmissionPage
        elif name == 'ResultsPostprocessingPage':
            from qtgui.pages.results_postprocessing import ResultsPostprocessingPage
            _page_modules[name] = ResultsPostprocessingPage
        else:
            raise ValueError(f"Unknown page class: {name}")
    return _page_modules[name]


class MainWindow(QMainWindow):
    """Main application window for xespresso PySide6 GUI (session window)."""

    def __init__(self, session_state=None, parent=None):
        super().__init__(parent)
        # This class represents a session workspace window.
        # Mark early so toolbar/menu builders can adapt their contents.
        self.is_session_window = True
        self.setWindowTitle("⚛️ xespresso - Quantum ESPRESSO Configuration GUI")

        # Set window size as a proportion of the screen (80% width, 80% height)
        # with a reasonable minimum for usability
        screen = QApplication.primaryScreen().availableGeometry()
        width = max(int(screen.width() * 0.8), 800)
        height = max(int(screen.height() * 0.8), 600)
        self.resize(width, height)

        # Set minimum size to ensure usability on small screens
        self.setMinimumSize(800, 600)

        # Center the window on the screen
        self.move(
            (screen.width() - width) // 2,
            (screen.height() - height) // 2
        )

        # Initialize session state - allow isolated instances for independent windows
        self.session_state = session_state or SessionState(isolated=True)

        # Configuration dialog (created on demand)
        self._config_dialog = None
        self._open_prov_btn = None

        # Structure label and lock state updater
        self.structure_label = None

        # Job monitor dialog (created on demand, accessible without session)
        self._job_monitor = None

        # Guard to prevent recursive updates during session changes
        self._updating = False

        # Setup menus, toolbar, and statusbar via ui_components helpers
        setup_menu(self)
        setup_toolbar(self)
        setup_statusbar(self)
        # add per-session provenance toolbar button (if toolbar exists)
        try:
            self._add_provenance_toolbar_button()
        except Exception:
            pass

        # Build sidebar widget using helper and create pages/content stack
        self.sidebar = create_sidebar(self)

        # Central content stack for pages
        self.content_stack = QStackedWidget()

        # Use a splitter so the sidebar can be resized/collapsed
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.sidebar)
        splitter.addWidget(self.content_stack)
        splitter.setSizes([300, 800])

        central = QWidget()
        central_layout = QHBoxLayout(central)
        central_layout.addWidget(splitter)
        self.setCentralWidget(central)

        # Create pages (lazy imports)
        self._create_pages()
        # Structure viewer (created on demand when the user clicks View)
        self._structure_viewer = None
        # Update structure display (if any)
        try:
            self._update_structure_label()
        except Exception:
            pass
        # Apply visual focus decorator so active windows get a clear drop-shadow
        try:
            apply_focus_decorator(self)
        except Exception:
            pass

    def _update_structure_label(self):
        """Update structure display in the sidebar if present."""
        try:
            if hasattr(self, 'structure_label') and self.structure_label is not None:
                source = self.session_state.get('structure_source') or 'No structure selected'
                atoms = self.session_state.get('current_structure')
                formula = None
                name = source
                if atoms is not None:
                    try:
                        if hasattr(atoms, 'get_chemical_formula'):
                            formula = atoms.get_chemical_formula()
                        else:
                            # fallback: count symbols
                            from collections import Counter
                            syms = atoms.get_chemical_symbols()
                            counts = Counter(syms)
                            formula = ''.join(f"{el}{counts[el] if counts[el]>1 else ''}" for el in sorted(counts))
                    except Exception:
                        formula = None

                locked = self.session_state.is_structure_locked()
                parts = [f"Source: {source}"]
                if formula:
                    parts.append(f"Formula: {formula}")
                # show atom count if available
                try:
                    if atoms is not None:
                        atom_count = len(atoms)
                        parts.append(f"Atoms: {atom_count}")
                except Exception:
                    pass
                # do not display locked/editable marker in label
                text = "\n".join(parts)
                self.structure_label.setText(text)

                # Disable/enable structure buttons if present
                try:
                    enabled = not (atoms is not None or locked)
                    if hasattr(self, 'load_struct_btn'):
                        self.load_struct_btn.setEnabled(enabled)
                    if hasattr(self, 'db_struct_btn'):
                        self.db_struct_btn.setEnabled(enabled)
                    if hasattr(self, 'build_struct_btn'):
                        self.build_struct_btn.setEnabled(enabled)
                except Exception:
                    pass
        except Exception:
            pass

    def _toggle_sidebar(self):
        """Toggle visibility of the left sidebar."""
        try:
            if hasattr(self, 'sidebar') and self.sidebar is not None:
                self.sidebar.setVisible(not self.sidebar.isVisible())
        except Exception:
            pass

    def _create_pages(self):
        """Create all page widgets (workflow pages only) using lazy imports."""
        # Create pages in order matching navigation list (Structure viewer is created on demand)
        self.pages = [
            _get_page_class('CalculationSetupPage')(self.session_state),
            _get_page_class('WorkflowBuilderPage')(self.session_state),
            _get_page_class('JobSubmissionPage')(self.session_state),
            _get_page_class('ResultsPostprocessingPage')(self.session_state)
        ]

        # Store reference to JobSubmissionPage for setting job monitor later
        self._job_submission_page = self.pages[2]

        for page in self.pages:
            self.content_stack.addWidget(page)
        try:
            if hasattr(self, 'sidebar') and self.sidebar is not None:
                self.sidebar.setVisible(not self.sidebar.isVisible())
        except Exception:
            pass
        # Restore any previously open workflow tabs saved in the session state
        try:
            open_tabs = self.session_state.get('open_workflow_tabs') or []
            if open_tabs:
                from qtgui.pages.workflow_tabs_page import WorkflowTabsPage
                tabs_page = None
                for i in range(self.content_stack.count()):
                    w = self.content_stack.widget(i)
                    if isinstance(w, WorkflowTabsPage):
                        tabs_page = w
                        break
                if tabs_page is None:
                    tabs_page = WorkflowTabsPage(self.session_state)
                    self.content_stack.addWidget(tabs_page)
                for name in open_tabs:
                    try:
                        tabs_page.add_workflow_tab(name)
                    except Exception:
                        pass
        except Exception:
            pass

    def launch_workflow(self, workflow_name: str):
        """Open the workflow launcher for a named preset (called from sidebar).

        This deliberately uses the `WorkflowLauncherDialog` to allow the user
        to confirm parameters before starting a run. For automated tests the
        launcher runs in debug mode so no external binaries are invoked.
        """
        try:
            # Always create a new independent workflow tab for each launch.
            # Use a timestamp to ensure uniqueness so multiple instances can
            # coexist concurrently in the session area.
            from datetime import datetime
            ts = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            tab_name = f"workflow::{workflow_name}::{ts}"

            # Ensure there's a WorkflowTabsPage in the content stack and add a tab
            from qtgui.pages.workflow_tabs_page import WorkflowTabsPage

            # Find existing WorkflowTabsPage
            tabs_page = None
            for i in range(self.content_stack.count()):
                w = self.content_stack.widget(i)
                if isinstance(w, WorkflowTabsPage):
                    tabs_page = w
                    break

            if tabs_page is None:
                tabs_page = WorkflowTabsPage(self.session_state)
                self.content_stack.addWidget(tabs_page)

            # Add a new independent tab for this workflow
            tabs_page.add_workflow_tab(workflow_name)
            self.content_stack.setCurrentWidget(tabs_page)
        except Exception:
            return

    def _choose_structure_file(self):
        """Prompt user to choose a structure file and set it for this session (once)."""
        sessions_dir = os.path.expanduser("~")
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Structure File", sessions_dir, "Structure Files (*.cif *.xsd *.xyz *.vasp *.pdb);;All Files (*)")
        if not file_path:
            return False
        # Open the standalone Structure Manager dialog on the Upload tab so the
        # user can upload and save to the central DB. When the manager emits a
        # selected DB id we load it into this session.
        try:
            from qtgui.dialogs import StructureManagerDialog
            dlg = StructureManagerDialog(self, for_session=True)
            def _on_selected(db_id):
                try:
                    ok = self.session_state.set_structure_from_database(db_id)
                    if ok:
                        self._update_structure_label()
                        QMessageBox.information(self, "Structure Loaded", f"Structure loaded into session from DB ID {db_id}")
                        try:
                            dlg.close()
                        except Exception:
                            pass
                        # Bring session window to front
                        try:
                            self.show()
                            self.raise_()
                            self.activateWindow()
                        except Exception:
                            pass
                    else:
                        QMessageBox.warning(self, "Error", self.session_state.get_last_error_message() or "Could not load structure from DB")
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"Could not load structure into session: {e}")
            dlg.structure_selected.connect(_on_selected)
            dlg.show_upload_tab()
            return True
        except Exception:
            # Fallback: behave as before using non-GUI controller helpers
            db_id, loaded, msg = session_controllers.add_structure_file_and_load(self.session_state, file_path)
            if loaded:
                try:
                    self._update_structure_label()
                except Exception:
                    pass
                if db_id:
                    QMessageBox.information(self, "Structure Saved", f"Structure saved to DB with ID {db_id} and loaded into session.")
                return True
            else:
                QMessageBox.warning(self, "Error", msg or "Could not set structure for this session (maybe it is locked)")
                return False

    def _choose_structure_db(self):
        # Open the standalone Structure Manager dialog on the DB tab
        try:
            from qtgui.dialogs import StructureManagerDialog
            dlg = StructureManagerDialog(self, for_session=True)
            def _on_selected(db_id):
                try:
                    ok = self.session_state.set_structure_from_database(db_id)
                    if ok:
                        self._update_structure_label()
                        QMessageBox.information(self, "Structure Loaded", f"Structure loaded into session from DB ID {db_id}")
                        try:
                            dlg.close()
                        except Exception:
                            pass
                        try:
                            self.show()
                            self.raise_()
                            self.activateWindow()
                        except Exception:
                            pass
                    else:
                        QMessageBox.warning(self, "Error", self.session_state.get_last_error_message() or "Could not load structure from DB")
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"Could not load structure into session: {e}")
            dlg.structure_selected.connect(_on_selected)
            dlg.show_db_tab()
            return True
        except Exception:
            # fallback to previous simple DB selection if manager unavailable
            try:
                rows = session_controllers.list_db_rows(self.session_state.get('structure_db_path') or None)
            except Exception as e:
                QMessageBox.warning(self, "DB Error", f"Could not open structures DB: {e}")
                return False
            from PySide6.QtWidgets import QDialog, QVBoxLayout, QListWidget, QPushButton, QHBoxLayout
            dlg = QDialog(self)
            dlg.setWindowTitle("Select Structure from DB")
            dlg.setMinimumSize(400, 300)
            v = QVBoxLayout(dlg)
            listw = QListWidget()
            for row in rows:
                try:
                    atoms = row.toatoms()
                    formula = atoms.get_chemical_formula() if hasattr(atoms, 'get_chemical_formula') else ''
                except Exception:
                    formula = ''
                listw.addItem(f"{row.id}  {formula}  {row.key_value_pairs.get('name','') if row.key_value_pairs else ''}")
            v.addWidget(listw)
            btn_row = QHBoxLayout()
            ok_btn = QPushButton("OK")
            cancel_btn = QPushButton("Cancel")
            btn_row.addWidget(ok_btn)
            btn_row.addWidget(cancel_btn)
            v.addLayout(btn_row)
            def on_ok():
                sel = listw.currentItem()
                if not sel:
                    return
                text = sel.text().split()[0]
                try:
                    sid = int(text)
                except Exception:
                    sid = None
                if sid is not None:
                    ok2 = self.session_state.set_structure_from_database(sid)
                    if ok2:
                        try:
                            self._update_structure_label()
                        except Exception:
                            pass
                        dlg.accept()
                        return
                QMessageBox.warning(self, "Error", "Could not load selected structure")
            ok_btn.clicked.connect(on_ok)
            cancel_btn.clicked.connect(dlg.reject)
            res = dlg.exec()
            return res == QDialog.Accepted

    def _build_structure(self):
        """Open the standalone Structure Manager on the Build tab."""
        try:
            from qtgui.dialogs import StructureManagerDialog
            dlg = StructureManagerDialog(self, for_session=True)
            def _on_selected(db_id):
                try:
                    ok = self.session_state.set_structure_from_database(db_id)
                    if ok:
                        self._update_structure_label()
                        QMessageBox.information(self, "Structure Loaded", f"Structure loaded into session from DB ID {db_id}")
                        try:
                            dlg.close()
                        except Exception:
                            pass
                        try:
                            self.show()
                            self.raise_()
                            self.activateWindow()
                        except Exception:
                            pass
                    else:
                        QMessageBox.warning(self, "Error", self.session_state.get_last_error_message() or "Could not load structure from DB")
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"Could not load structure into session: {e}")
            dlg.structure_selected.connect(_on_selected)
            dlg.show_build_tab()
            return True
        except Exception:
            QMessageBox.information(self, "Build", "Structure builder not available.")
            return False

    def _view_structure(self):
        """Switch the main content to the Structure Viewer page if available."""
        try:
            # Create the StructureViewerPage on demand and show it
            if self._structure_viewer is None:
                from qtgui.pages.structure_viewer import StructureViewerPage
                # Create view-only structure viewer so only visualization is shown
                self._structure_viewer = StructureViewerPage(self.session_state, view_only=True)
                self.content_stack.addWidget(self._structure_viewer)
                # If a structure is already loaded in session, ensure viewer displays it
                atoms = self.session_state.get('current_structure')
                source = self.session_state.get('structure_source')
                try:
                    if atoms is not None:
                        self._structure_viewer._set_structure(atoms, source or '')
                except Exception:
                    pass
            # Switch to the dynamic structure viewer widget
            self.content_stack.setCurrentWidget(self._structure_viewer)
            return True
        except Exception:
            QMessageBox.information(self, "View", "Structure viewer not available.")
            return False

    def _setup_menu(self):
        """Setup the menu bar."""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("&File")

        new_session_action = QAction("New Session", self)
        new_session_action.setShortcut("Ctrl+N")
        new_session_action.triggered.connect(self._new_session)
        file_menu.addAction(new_session_action)

        save_session_action = QAction("Save Session", self)
        save_session_action.setShortcut("Ctrl+S")
        save_session_action.triggered.connect(self._save_session)
        file_menu.addAction(save_session_action)

        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Edit menu
        edit_menu = menubar.addMenu("&Edit")

        config_action = QAction("Configuration...", self)
        config_action.setShortcut("Ctrl+,")
        config_action.triggered.connect(self._open_config_dialog)
        edit_menu.addAction(config_action)

        # View menu
        view_menu = menubar.addMenu("&View")

        # Navigate menu items
        nav_actions = [
            ("Structure", "Ctrl+1"),
            ("Calculation Setup", "Ctrl+2"),
            ("Workflow Builder", "Ctrl+3"),
            ("Job Submission", "Ctrl+4"),
            ("Results", "Ctrl+5"),
        ]

        for i, (name, shortcut) in enumerate(nav_actions):
            action = QAction(name, self)
            action.setShortcut(shortcut)
            action.triggered.connect(lambda checked, idx=i: self._navigate_to(idx))
            view_menu.addAction(action)

        # Help menu
        help_menu = menubar.addMenu("&Help")

        about_action = QAction("About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

        # Structure Manager accessible from session windows as well
        try:
            struct_action = QAction("Structure Manager...", self)
            struct_action.triggered.connect(lambda: self._open_structure_manager())
            tools_menu = menubar.addMenu("&Tools")
            tools_menu.addAction(struct_action)
        except Exception:
            pass

    def _setup_toolbar(self):
        """Delegate toolbar setup to `ui_components.setup_toolbar`."""
        setup_toolbar(self)

    def _setup_statusbar(self):
        """Delegate statusbar setup to `ui_components.setup_statusbar`."""
        setup_statusbar(self)

    def _update_statusbar(self):
        """Update status bar with current session info."""
        session_name = self.session_state.get_session_name()
        self.statusbar.showMessage(f"Session: {session_name} | Ready")

    def _update_session_name_label(self):
        """Update the session name label in the sidebar."""
        session_name = self.session_state.get_session_name()
        session_id = self.session_state.get_current_session_id()
        try:
            # Some sidebars use QLabel named session_name_label
            self.session_name_label.setText(f"<b>{session_name}</b><br><small>ID: {session_id}</small>")
            self.session_name_label.setTextFormat(Qt.RichText)
        except Exception:
            pass

    def _refresh_session_list(self):
        """Populate the session combo with active sessions (session-local)."""
        # If a session combo exists (manager window), refresh it; otherwise noop
        if not hasattr(self, 'session_combo'):
            return
        try:
            self.session_combo.blockSignals(True)
            self.session_combo.clear()
            active = self.session_state.get_active_session_names()
            for name in active:
                self.session_combo.addItem(name, name)
            current = self.session_state.get_session_name()
            if current in active:
                idx = active.index(current)
                self.session_combo.setCurrentIndex(idx)
        finally:
            try:
                self.session_combo.blockSignals(False)
            except Exception:
                pass

    def _on_session_selected(self, text):
        """Handle session selection from combo box (no-op for isolated session windows)."""
        # Session windows have isolated state; switching active session via
        # the sidebar combo isn't supported here. Ignore selection changes.
        return

    def _toggle_edit_mode(self):
        """Toggle a simple edit mode flag and notify pages that support editing."""
        editing = getattr(self, '_editing', False)
        self._editing = not editing
        for page in getattr(self, 'pages', []):
            if hasattr(page, 'set_edit_mode'):
                try:
                    page.set_edit_mode(self._editing)
                except Exception:
                    pass
        try:
            self.statusbar.showMessage("Edit mode on" if self._editing else "Edit mode off")
        except Exception:
            pass

    def _load_session_dialog(self):
        """Open dialog to load a session file and open it in a new window."""
        sessions_dir = self.session_state.get_sessions_dir()
        os.makedirs(sessions_dir, exist_ok=True)
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Session File",
            sessions_dir,
            "Session Files (*.json);;All Files (*)",
        )
        if not file_path:
            return
        state = SessionState(isolated=True)
        if state.load_session_from_file(file_path):
            new_win = MainWindow(session_state=state)
            new_win.setWindowTitle(f"⚛️ xespresso - {state.get_session_name()}")
            new_win.show()
        else:
            QMessageBox.warning(self, "Error", f"Could not load session file:\n{file_path}")

    def _open_config_dialog(self):
        """Open the configuration dialog (uses shared configuration state)."""
        from qtgui.dialogs import ConfigurationDialog

        if self._config_dialog is None:
            self._config_dialog = ConfigurationDialog(self.session_state.get_config_state(), self)
            try:
                self._config_dialog.configuration_changed.connect(self._on_config_changed)
            except Exception:
                pass

        self._config_dialog.show()
        self._config_dialog.raise_()
        self._config_dialog.activateWindow()

    def _on_config_changed(self):
        """Handle configuration changes from the dialog by refreshing pages."""
        for page in getattr(self, 'pages', []):
            if hasattr(page, 'refresh'):
                try:
                    page.refresh()
                except Exception:
                    pass
        try:
            self.statusbar.showMessage("Configuration updated")
        except Exception:
            pass

    def _add_provenance_toolbar_button(self):
        try:
            from PySide6.QtWidgets import QPushButton
            btn = QPushButton("Provenance")
            btn.setToolTip("Open Provenance Browser for current structure")
            btn.clicked.connect(self._open_provenance_for_current)
            # add to existing toolbar if available
            try:
                tb = getattr(self, 'toolbar', None)
                if tb is not None and hasattr(tb, 'addWidget'):
                    tb.addWidget(btn)
            except Exception:
                pass
            self._open_prov_btn = btn
        except Exception:
            pass

    def _open_provenance_for_current(self):
        try:
            atoms = self.session_state.get('current_structure')
            manager = getattr(self, '_manager', None)
            if manager is None:
                # fallback: open a local provenance panel if manager not present
                try:
                    from qtgui.pages.workflow_tasks.provenance_panel import ProvenancePanel
                    win = None
                    if hasattr(self, 'centralWidget'):
                        win = self
                    panel = ProvenancePanel(parent=self)
                    if atoms is not None and hasattr(panel, 'update_for_atoms'):
                        panel.update_for_atoms(atoms)
                    # attempt to show as a dialog-like widget
                    try:
                        panel.show()
                    except Exception:
                        pass
                    return
                except Exception:
                    return

            # Use manager's provenance browser so it's shared globally
            try:
                pb = manager._get_provenance_browser()
                if pb is None:
                    return
                # If centralWidget is a ProvenancePanel, update it
                try:
                    panel = pb.centralWidget()
                    if atoms is not None and hasattr(panel, 'update_for_atoms'):
                        panel.update_for_atoms(atoms)
                except Exception:
                    pass
                manager._open_provenance_browser()
            except Exception:
                pass
        except Exception:
            pass

    def _rename_session(self):
        """Rename the current session and update manager if present."""
        current_name = self.session_state.get_session_name()
        name, ok = QInputDialog.getText(
            self,
            "Rename Session",
            "Enter a new name for the session:",
            QLineEdit.Normal,
            current_name,
        )
        if not ok or not name:
            return
        try:
            # Rename underlying session storage
            self.session_state.rename_session(name)
            # Notify manager to refresh list if available
            try:
                if hasattr(self, '_manager') and self._manager is not None:
                    try:
                        self._manager._refresh_sessions()
                    except Exception:
                        pass
            except Exception:
                pass
            try:
                self._update_session_name_label()
            except Exception:
                pass
            try:
                self.statusbar.showMessage(f"Session renamed to: {name}")
            except Exception:
                pass
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not rename session:\n{e}")

    def _on_nav_changed(self, index):
        """Handle navigation change."""
        self.content_stack.setCurrentIndex(index)

        # Update status bar (structure viewer is shown on demand via View)
        page_names = [
            "Calculation Setup",
            "Workflow Builder",
            "Job Submission",
            "Results & Post-Processing"
        ]
        if 0 <= index < len(page_names):
            session_name = self.session_state.get_session_name()
            self.statusbar.showMessage(f"Session: {session_name} | Page: {page_names[index]}")
        # Refresh the target page if it provides a refresh method
        try:
            if hasattr(self, 'pages') and 0 <= index < len(self.pages):
                page = self.pages[index]
                if hasattr(page, 'refresh'):
                    try:
                        page.refresh()
                    except Exception:
                        pass
        except Exception:
            pass

    def _on_session_changed(self):
        """Handle UI updates when the current session changes.

        Tests call this hook directly to force the UI to refresh labels
        after programmatic session switches. Keep implementation minimal
        and robust for headless test runs.
        """
        try:
            # Update session name label and working directory display
            try:
                self._update_session_name_label()
            except Exception:
                pass
            try:
                wd = self.session_state.get('working_directory', os.path.expanduser("~"))
                if hasattr(self, 'workdir_label') and self.workdir_label is not None:
                    self.workdir_label.setText(wd)
            except Exception:
                pass

            # Refresh pages to reflect new session state (if they implement refresh)
            for page in getattr(self, 'pages', []):
                if hasattr(page, 'refresh'):
                    try:
                        page.refresh()
                    except Exception:
                        pass
        except Exception:
            pass

    def _navigate_to(self, index):
        """Navigate to a specific page."""
        self.nav_list.setCurrentRow(index)

    def _browse_workdir(self):
        """Open directory browser for working directory."""
        current_dir = self.session_state.get('working_directory', os.path.expanduser("~"))
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Working Directory",
            current_dir,
            QFileDialog.ShowDirsOnly
        )
        if directory:
            self.session_state['working_directory'] = directory
            try:
                self.workdir_label.setText(directory)
            except Exception:
                pass
            self.statusbar.showMessage(f"Working directory set to: {directory}")

    def _new_session(self):
        """Create a new session in a new window (sessions are independent)."""
        name, ok = QInputDialog.getText(
            self,
            "New Session",
            "Enter a name for the new session:",
            QLineEdit.Normal,
            f"Session {len(self.session_state.list_sessions()) + 1}"
        )

        if ok and name:
            # Create isolated session state and open it in a new MainWindow
            # Prevent overwriting an existing saved session silently
            existing = None
            for sid, meta in self.session_state._sessions.items():
                if isinstance(meta.get('name'), str) and meta.get('name').lower() == name.lower():
                    existing = sid
                    break
            if existing is not None:
                resp = QMessageBox.question(self, "Overwrite Session?", f"A saved session named '{name}' already exists. Overwrite it?", QMessageBox.Yes | QMessageBox.No)
                if resp != QMessageBox.Yes:
                    return

            state = SessionState(isolated=True)
            try:
                state.create_session(name)
            except Exception as exc:
                QMessageBox.warning(self, "Error", f"Could not create session:\n{exc}")
                return

            new_window = MainWindow(session_state=state)
            new_window.setWindowTitle(f"⚛️ xespresso - {name}")
            new_window.show()
            self.statusbar.showMessage(f"Opened new session window: {name}")

    def _save_session(self):
        """Save the current session.

        First collects current state from all pages, then saves to disk.
        """
        self._updating = True
        try:
            for page in self.pages:
                if hasattr(page, 'save_state'):
                    try:
                        page.save_state()
                    except Exception as e:
                        print(f"Warning: Could not save state from page: {e}")

            self.session_state.save_session()
            self.statusbar.showMessage("Session saved")
        finally:
            self._updating = False

    def _close_session(self):
        """Close the current session and clear all data."""
        reply = QMessageBox.question(
            self,
            "Close Session",
            "Close the current session and clear all data?\n\n"
            "Make sure to save your session first if you want to keep it.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            self.session_state.clear_state()
            try:
                self.session_combo.blockSignals(True)
                self.session_combo.clear()
                self.session_combo.blockSignals(False)
            except Exception:
                pass

            from PySide6.QtWidgets import QApplication
            QApplication.processEvents()

            for page in self.pages:
                if hasattr(page, 'refresh'):
                    try:
                        page.refresh()
                        QApplication.processEvents()
                    except Exception as e:
                        print(f"Warning: Could not refresh page: {e}")

            try:
                self._update_session_name_label()
                self.workdir_label.setText(self.session_state.get('working_directory', '~'))
            except Exception:
                pass

            self.statusbar.showMessage("Session closed")

    def _show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About xespresso GUI",
            """<h2>xespresso GUI</h2>
<p><b>Version:</b> 1.2.0</p>
<p>PySide6 interface for Quantum ESPRESSO calculations.</p>
"""
        )

    def _get_job_monitor(self):
        """
        Get or create the Job Monitor dialog instance.
        """
        if self._job_monitor is None:
            from qtgui.dialogs.job_monitor_dialog import JobMonitorDialog
            xespresso_dir = os.path.expanduser("~/.xespresso")
            self._job_monitor = JobMonitorDialog(config_dir=xespresso_dir, parent=self)

            if hasattr(self, '_job_submission_page'):
                self._job_submission_page.set_job_monitor(self._job_monitor)

        return self._job_monitor

    def _open_job_monitor(self):
        job_monitor = self._get_job_monitor()
        job_monitor.show()
        job_monitor.raise_()
        job_monitor.activateWindow()

    def _open_structure_manager(self):
        try:
            from qtgui.dialogs import StructureManagerDialog
            if not hasattr(self, '_structure_manager') or self._structure_manager is None:
                self._structure_manager = StructureManagerDialog(self)
            self._structure_manager.show_db_tab()
            self._structure_manager.raise_()
            self._structure_manager.activateWindow()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not open Structure Manager:\n{e}")

    def closeEvent(self, event):
        """Handle window close event."""
        self.session_state.save_session()
        if self._config_dialog is not None:
            self._config_dialog.close()

        try:
            if hasattr(self, '_manager') and hasattr(self, '_manager_session_id'):
                try:
                    self._manager._on_session_window_closed(self._manager_session_id)
                except Exception:
                    pass
        except Exception:
            pass

        event.accept()

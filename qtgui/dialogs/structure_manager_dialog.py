"""
Structure Manager Dialog

Wraps the existing `StructureViewerPage` into a non-modal dialog that acts as
an independent structure manager (Upload / Build / DB / View). The dialog
emits `structure_selected(int)` when the user chooses a structure to use in a
session (i.e. selects or saves a structure to the ASE DB and clicks "Use in session").
"""
from qtpy.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QPushButton
from qtpy.QtCore import Signal

from qtgui.pages.structure_viewer import StructureViewerPage
from qtgui.session_state import SessionState


class StructureManagerDialog(QDialog):
    structure_selected = Signal(int)

    def __init__(self, parent=None, for_session=False):
        super().__init__(parent)
        self.setWindowTitle("Structure Manager")
        self.setModal(False)
        self.setMinimumSize(900, 700)

        # Use the shared session state instance for structure management
        shared_state = SessionState()

        layout = QVBoxLayout(self)
        self.page = StructureViewerPage(shared_state)
        # remember if this dialog was opened from a session button; when True,
        # actions that 'load' a structure from the DB should behave like
        # 'Use in Session' — emit `structure_selected` and close the dialog.
        self._for_session = bool(for_session)

        # Wrap the page's DB-load method so that when used from a session
        # the selected DB entry will be emitted and the dialog closed.
        try:
            orig_load = getattr(self.page, '_load_from_database', None)
            if callable(orig_load):
                def _wrapped_load_from_database(*args, **kwargs):
                    res = orig_load(*args, **kwargs)
                    try:
                        if self._for_session:
                            src = self.page.session_state.get('structure_source', '')
                            import re
                            m = re.match(r'^Database:\s*ID\s*(\d+)$', src)
                            if m:
                                sid = int(m.group(1))
                                # emit and close
                                try:
                                    self.structure_selected.emit(sid)
                                except Exception:
                                    pass
                                try:
                                    self.close()
                                except Exception:
                                    pass
                    except Exception:
                        pass
                    return res
                setattr(self.page, '_load_from_database', _wrapped_load_from_database)
        except Exception:
            pass
        layout.addWidget(self.page)

        # Action row
        row = QHBoxLayout()
        row.addStretch()
        self.use_btn = QPushButton("Use in Session")
        self.use_btn.clicked.connect(self._use_in_session)
        row.addWidget(self.use_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        row.addWidget(close_btn)

        layout.addLayout(row)

    def show_upload_tab(self):
        try:
            self.page.tabs.setCurrentIndex(0)
            self.show()
            self.raise_()
            self.activateWindow()
        except Exception:
            pass

    def show_build_tab(self):
        try:
            # Build tab is index 1 in StructureViewerPage
            self.page.tabs.setCurrentIndex(1)
            self.show()
            self.raise_()
            self.activateWindow()
        except Exception:
            pass

    def show_db_tab(self):
        try:
            # Database tab is index 2
            self.page.tabs.setCurrentIndex(2)
            # refresh list
            try:
                self.page._refresh_db_list()
            except Exception:
                pass
            self.show()
            self.raise_()
            self.activateWindow()
        except Exception:
            pass

    def show_view_tab(self):
        try:
            # View tab is index 3
            self.page.tabs.setCurrentIndex(3)
            self.show()
            self.raise_()
            self.activateWindow()
        except Exception:
            pass

    def _use_in_session(self):
        """Determine the currently loaded structure and emit its DB id.

        Behavior:
        - If the page indicates the structure was loaded from DB (structure_source), parse and emit that id.
        - Else, attempt to save the current structure to the DB (using the upload or build save helpers)
          and emit the new id when available.
        """
        src = self.page.session_state.get('structure_source', '')
        # Case: already Database: ID N
        import re
        m = re.match(r'^Database:\s*ID\s*(\d+)$', src)
        if m:
            try:
                sid = int(m.group(1))
                self.structure_selected.emit(sid)
                return
            except Exception:
                pass

        # Otherwise try to save current structure using available save widgets
        # Prefer the upload tab widgets if present
        row_id = None
        try:
            # Try upload save widgets
            if hasattr(self.page, 'upload_db_path_edit'):
                row_id = self.page._save_structure_to_db(
                    getattr(self.page, 'upload_db_path_edit'),
                    getattr(self.page, 'upload_name_edit'),
                    getattr(self.page, 'upload_tags_edit')
                )
        except Exception:
            row_id = None

        if not row_id:
            try:
                # Try build save widgets
                if hasattr(self.page, 'build_db_path_edit'):
                    row_id = self.page._save_structure_to_db(
                        getattr(self.page, 'build_db_path_edit'),
                        getattr(self.page, 'build_name_edit'),
                        getattr(self.page, 'build_tags_edit')
                    )
            except Exception:
                row_id = None

        if row_id:
            try:
                self.structure_selected.emit(int(row_id))
                return
            except Exception:
                pass

        # As a last-resort, if current_structure exists but couldn't be saved, notify user
        try:
            if self.page.session_state.get('current_structure') is not None:
                from qtpy.QtWidgets import QMessageBox
                QMessageBox.information(self, "Info", "Structure is loaded in the manager but could not be saved to DB. Use the ASE Database tab to save it manually.")
                return
        except Exception:
            pass

        # Nothing to do
        from qtpy.QtWidgets import QMessageBox
        QMessageBox.warning(self, "No Structure", "No structure available to use in session.")

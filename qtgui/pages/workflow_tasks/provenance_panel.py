from __future__ import annotations

from pathlib import Path
from typing import List

try:
    from qtpy.QtWidgets import QWidget, QVBoxLayout, QListWidget, QPushButton, QTextEdit, QFileDialog
except Exception:
    # If Qt isn't available in this environment, provide a noop placeholder
    QWidget = object

from ase import Atoms
from xespresso.provenance import ProvenanceDB, ProvenanceRecord


class ProvenancePanel(QWidget):
    """Simple panel to query and display provenance records for a structure.

    This is a lightweight Qt widget showing a list of records and a details
    view. It expects to be embedded in the main application and that the
    caller provides the current `Atoms` object when requesting updates.
    """

    def __init__(self, parent=None, db_path: str | None = None):
        super().__init__() if hasattr(self.__class__, "__mro__") else None
        # Prefer the default cached DB (which respects env var). Allow
        # explicit override via db_path parameter.
        if db_path:
            self.db = ProvenanceDB.get_default(path=db_path)
        else:
            self.db = ProvenanceDB.get_default()
        self._build_ui()

    def _build_ui(self):
        try:
            self.layout = QVBoxLayout()
            self.list_widget = QListWidget()
            self.details = QTextEdit()
            self.details.setReadOnly(True)
            self.refresh_btn = QPushButton("Refresh for current structure")
            self.open_db_btn = QPushButton("Open DB Entry")
            self.open_db_btn.setToolTip("Open the ASE DB entry linked to the selected provenance record")
            self.open_db_btn.clicked.connect(self._on_open_db_entry)
            self.save_btn = QPushButton("Save selected snapshot as...")

            self.layout.addWidget(self.list_widget)
            self.layout.addWidget(self.details)
            self.layout.addWidget(self.refresh_btn)
            self.layout.addWidget(self.open_db_btn)
            self.layout.addWidget(self.save_btn)
            self.setLayout(self.layout)

            self.refresh_btn.clicked.connect(self._on_refresh)
            self.save_btn.clicked.connect(self._on_save_snapshot)
            self.list_widget.itemSelectionChanged.connect(self._on_select)
            self._last_atoms = None
            self._last_records: List[ProvenanceRecord] = []
        except Exception:
            pass

    def update_for_atoms(self, atoms: Atoms):
        """Query the provenance DB for `atoms` and populate the list."""
        try:
            # If the default DB path may have changed externally, refresh
            # the cached default to ensure we query the right file. This is
            # inexpensive and safe.
            self.db = ProvenanceDB.get_default()
            records = self.db.get_records(atoms)
            self.list_widget.clear()
            self._last_atoms = atoms
            self._last_records = records
            for r in records:
                display = f"{r.timestamp} - {r.calc_type} - {r.label} - E={r.energy}"
                self.list_widget.addItem(display)
            self.details.clear()
        except Exception:
            pass

    def _on_refresh(self):
        # Placeholder: the caller should call `update_for_atoms(atoms)`
        pass

    def _on_select(self):
        try:
            idx = self.list_widget.currentRow()
            if idx < 0 or idx >= len(self._last_records):
                self.details.clear()
                return
            rec = self._last_records[idx]
            txt = []
            txt.append(f"Timestamp: {rec.timestamp}")
            txt.append(f"Type: {rec.calc_type}")
            txt.append(f"Label: {rec.label}")
            txt.append(f"Energy: {rec.energy}")
            txt.append(f"Fermi: {rec.efermi}")
            txt.append(f"Params: {rec.params}")
            txt.append(f"Meta: {rec.meta}")
            # If linked ASE DB id is present, show it and enable Open DB button
            db_id = None
            try:
                db_id = rec.meta.get('db_id') if rec.meta else None
            except Exception:
                db_id = None
            if db_id:
                txt.append(f"Linked ASE DB id: {db_id}")
                try:
                    self.open_db_btn.setEnabled(True)
                except Exception:
                    pass
            else:
                try:
                    self.open_db_btn.setEnabled(False)
                except Exception:
                    pass
            self.details.setPlainText("\n".join(txt))
        except Exception:
            pass

    def _on_save_snapshot(self):
        try:
            idx = self.list_widget.currentRow()
            if idx < 0:
                return
            path, _ = QFileDialog.getSaveFileName(self, "Save snapshot", "", "Trajectory (*.traj);;All Files (*)")
            if not path:
                return
            if idx >= len(self._last_records):
                return
            rec = self._last_records[idx]
            snap = rec.meta.get("snapshot") if rec and rec.meta else None
            if not snap:
                return
            snap_path = Path(snap)
            if not snap_path.exists():
                return
            # copy file to target
            import shutil

            shutil.copy(str(snap_path), path)
        except Exception:
            pass

    def _on_open_db_entry(self):
        """Open the ASE DB entry linked to the selected provenance record.

        This is a lightweight helper: it opens the ASE DB file (if configured in
        `ProvenanceDB` path or prompts the user) and displays the row info in a
        message box. It's best-effort and headless-safe.
        """
        try:
            idx = self.list_widget.currentRow()
            if idx < 0 or idx >= len(self._last_records):
                return
            rec = self._last_records[idx]
            db_id = rec.meta.get('db_id') if rec and rec.meta else None
            if not db_id:
                QMessageBox.warning(self, "DB Entry", "No linked ASE DB id for this record")
                return

            # Prefer to open via parent MainWindow's Structure Viewer if available
            parent = self.parent()
            mw = None
            while parent is not None:
                if hasattr(parent, '_view_structure') and hasattr(parent, 'session_state'):
                    mw = parent
                    break
                parent = getattr(parent, 'parent', None)() if callable(getattr(parent, 'parent', None)) else getattr(parent, 'parent', None)

            # Determine candidate DB file: try session_state stored path first
            candidate_db = None
            try:
                if mw is not None:
                    candidate_db = mw.session_state.get('structure_db_path')
            except Exception:
                candidate_db = None

            # If we have a main window and a DB path, open the structure viewer and select
            if mw is not None and candidate_db:
                try:
                    mw._view_structure()
                    sv = getattr(mw, '_structure_viewer', None)
                    if sv is not None and hasattr(sv, 'open_db_and_select'):
                        ok = sv.open_db_and_select(candidate_db, int(db_id))
                        if ok:
                            return
                except Exception:
                    pass

            # Fallback: open ASE DB file chooser and show row info (previous behavior)
            try:
                from ase.db import connect as ase_db_connect
            except Exception:
                QMessageBox.warning(self, "DB Entry", "ASE DB support not available in this environment")
                return

            # Try best-effort guess from provenance DB dir
            db_file = None
            try:
                pdb_path = getattr(self.db, 'path', None)
                if pdb_path:
                    pdir = Path(pdb_path).parent
                    candidates = list(pdir.glob('*.db'))
                    if len(candidates) == 1:
                        db_file = str(candidates[0])
            except Exception:
                db_file = None

            if not db_file:
                db_file, _ = QFileDialog.getOpenFileName(self, "Select ASE DB file", os.path.expanduser("~"), "ASE DB files (*.db);;All files (*)")
                if not db_file:
                    return

            try:
                db = ase_db_connect(db_file)
                try:
                    row = db.get(id=int(db_id))
                except Exception:
                    QMessageBox.warning(self, "DB Entry", f"Could not find row id {db_id} in {db_file}")
                    return

                info = [f"DB: {db_file}", f"ID: {db_id}", f"Formula: {row.formula}", f"Atoms: {row.natoms}"]
                if getattr(row, 'key_value_pairs', None):
                    info.append(f"Meta: {row.key_value_pairs}")
                QMessageBox.information(self, "DB Entry", "\n".join(info))
            except Exception as e:
                QMessageBox.warning(self, "DB Entry", f"Could not open ASE DB: {e}")
                return
        except Exception:
            pass

from qtpy.QtWidgets import QWidget, QHBoxLayout, QLabel, QLineEdit, QComboBox
from qtpy.QtCore import Signal


class CalculationHeaderWidget(QWidget):
    """Lightweight header widget for workflow tabs.

    Provides top-level inputs for machine selection, pseudopotentials and
    a run label. This is intentionally simple and does not import the
    heavy `CalculationSetupPage`.
    """

    changed = Signal()

    def __init__(self, session_state, preset_name=None, parent=None):
        super().__init__(parent)
        self.session = session_state
        self.preset_name = preset_name
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        # Machine selector (editable combo so it can display current selection)
        layout.addWidget(QLabel("Machine:"))
        self.machine_combo = QComboBox()
        self.machine_combo.setEditable(True)
        cur = None
        try:
            cur = self.session.get('current_machine_name')
        except Exception:
            cur = None
        if cur:
            self.machine_combo.addItem(cur)
            self.machine_combo.setCurrentText(cur)
        layout.addWidget(self.machine_combo)

        # Pseudopotentials quick entry
        layout.addWidget(QLabel("Pseudopotentials:"))
        self.pseudo_edit = QLineEdit()
        self.pseudo_edit.setPlaceholderText("e.g. Fe=Fe.UPF,O=O.UPF")
        layout.addWidget(self.pseudo_edit)

        # Label
        layout.addWidget(QLabel("Label:"))
        self.label_edit = QLineEdit()
        if self.preset_name:
            self.label_edit.setPlaceholderText(f"{self.preset_name} run label (optional)")
        else:
            self.label_edit.setPlaceholderText("Run label (optional)")
        layout.addWidget(self.label_edit)

        # Connect simple change signals
        try:
            self.machine_combo.currentTextChanged.connect(self.changed)
            self.pseudo_edit.editingFinished.connect(self.changed)
            self.label_edit.editingFinished.connect(self.changed)
        except Exception:
            pass

    def get_header_config(self):
        """Return a dict with machine_name, pseudopotentials and label."""
        cfg = {}
        try:
            m = self.machine_combo.currentText().strip()
            if m:
                cfg['machine_name'] = m
        except Exception:
            pass
        try:
            s = self.pseudo_edit.text().strip()
            if s:
                parts = [p.strip() for p in s.split(',') if p.strip()]
                pp = {}
                for p in parts:
                    if '=' in p:
                        k, v = p.split('=', 1)
                        pp[k.strip()] = v.strip()
                if pp:
                    cfg['pseudopotentials'] = pp
        except Exception:
            pass
        try:
            lbl = self.label_edit.text().strip()
            if lbl:
                cfg['label'] = lbl
        except Exception:
            pass
        return cfg

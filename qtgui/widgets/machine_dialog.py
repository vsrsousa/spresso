from qtpy.QtWidgets import QDialog, QVBoxLayout, QFormLayout, QComboBox, QLabel, QDialogButtonBox

# Optional xespresso helpers
try:
    from xespresso.machines.config.loader import list_machines, load_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR
    XESPRESSO_AVAILABLE = True
except Exception:
    list_machines = lambda *a, **k: []
    load_machine = lambda *a, **k: None
    load_codes_config = lambda *a, **k: None
    DEFAULT_CONFIG_PATH = None
    DEFAULT_MACHINES_DIR = None
    DEFAULT_CODES_DIR = None
    XESPRESSO_AVAILABLE = False


class MachineDialog(QDialog):
    """Minimal dialog to select machine, QE version and code.

    Returns selected values via `get_config()` after accepted.
    """

    def __init__(self, session=None, parent=None):
        super().__init__(parent)
        self.session = session
        self.setWindowTitle('Select Machine')
        self._init_ui()
        try:
            self._load_machines()
        except Exception:
            pass

    def _init_ui(self):
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.machine_combo = QComboBox(); self.machine_combo.setEditable(True)
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        self.machine_info = QLabel('')
        self.version_combo = QComboBox()
        self.version_combo.currentTextChanged.connect(self._on_version_changed)
        self.code_combo = QComboBox()
        form.addRow('Machine:', self.machine_combo)
        form.addRow(self.machine_info)
        form.addRow('QE Version:', self.version_combo)
        form.addRow('Code:', self.code_combo)
        layout.addLayout(form)
        # Keep OK/Cancel buttons; WorkflowInstancePage will show this dialog
        # modelessly and listen for the finished(int) signal to apply changes.
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _load_machines(self):
        if not XESPRESSO_AVAILABLE:
            return
        try:
            machines = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            for m in machines:
                self.machine_combo.addItem(m)
            # try to populate codes if session has current_machine_name
            cur = None
            try:
                cur = self.session.get('current_machine_name')
            except Exception:
                cur = None
            if cur:
                idx = self.machine_combo.findText(cur)
                if idx >= 0:
                    self.machine_combo.setCurrentIndex(idx)
            # trigger a change to populate codes
            try:
                if cur:
                    self._on_machine_changed(cur)
            except Exception:
                pass
        except Exception:
            pass

    def get_config(self):
        return {
            'machine_name': self.machine_combo.currentText(),
            'version': self.version_combo.currentText(),
            'code': self.code_combo.currentText(),
        }

    def _on_machine_changed(self, machine_name):
        if not machine_name or not XESPRESSO_AVAILABLE:
            return
        try:
            machine = load_machine(DEFAULT_CONFIG_PATH, machine_name, DEFAULT_MACHINES_DIR, return_object=True)
            info = f"Type: {getattr(machine, 'execution', '')}"
            if getattr(machine, 'scheduler', None):
                info += f", Scheduler: {machine.scheduler}"
            try:
                self.machine_info.setText(info)
            except Exception:
                pass
            # load codes for this machine if available
            try:
                codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, verbose=False)
                self.version_combo.clear()
                self.code_combo.clear()
                if codes and getattr(codes, 'has_any_codes', None) and codes.has_any_codes():
                    versions = codes.list_versions()
                    if versions:
                        for version in versions:
                            self.version_combo.addItem(version)
                    else:
                        all_codes = codes.get_all_codes()
                        for code_name in all_codes.keys():
                            self.code_combo.addItem(code_name)
            except Exception:
                pass
        except Exception:
            pass

    def _on_version_changed(self, version):
        if not version or not XESPRESSO_AVAILABLE:
            return
        try:
            machine_name = self.machine_combo.currentText()
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, version=version, verbose=False)
            self.code_combo.clear()
            if codes:
                for code_name in codes.get_all_codes(version=version).keys():
                    self.code_combo.addItem(code_name)
        except Exception:
            pass

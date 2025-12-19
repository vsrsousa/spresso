import json
import logging
from typing import Dict, Iterable

from qtpy.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QComboBox,
    QPushButton, QGroupBox, QFormLayout, QDoubleSpinBox, QPlainTextEdit, QDialog
)
from qtpy.QtCore import Signal

logger = logging.getLogger(__name__)

# Try to import optional helpers; if missing provide safe fallbacks.
try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except Exception:
    PSEUDO_SELECTOR_AVAILABLE = False

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


class FallbackPseudoSelector(QWidget):
    """Simple text-based pseudopotential editor used when the real selector is missing.

    API compatibility with the real selector:
    - `changed` Signal
    - `get_pseudopotentials()` -> Dict[str, str]
    - `set_elements(iterable_of_symbols)`
    - `set_pseudopotentials(mapping)`
    """

    changed = Signal()

    def __init__(self, session_state=None, parent=None):
        super().__init__(parent)
        self.session = session_state
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel('Pseudopotentials (one per line: Element: filename)'))
        self.text = QPlainTextEdit()
        self.text.textChanged.connect(self.changed)
        layout.addWidget(self.text)

    def get_pseudopotentials(self) -> Dict[str, str]:
        out = {}
        for line in self.text.toPlainText().splitlines():
            line = line.strip()
            if not line:
                continue
            if ':' in line:
                el, fn = line.split(':', 1)
                out[el.strip()] = fn.strip()
        return out

    def set_elements(self, elements: Iterable[str]):
        # No-op for fallback, but keep a hint in the text area if empty.
        if not self.text.toPlainText().strip():
            hint = '\n'.join(f"{el}: " for el in sorted(set(elements)))
            self.text.setPlainText(hint)

    def set_pseudopotentials(self, mapping: Dict[str, str]):
        try:
            s = '\n'.join(f"{k}: {v}" for k, v in sorted(mapping.items()))
            self.text.setPlainText(s)
        except Exception:
            logger.debug('Failed to set pseudopotentials in fallback', exc_info=True)


class CalculationConfigWidget(QWidget):
    """Configuration widget used inside workflow tabs.

    Provides a safe, import-friendly UI for selecting machine, code, basic
    QE parameters and pseudopotentials. Magnetism/Hubbard controls are omitted
    per current agreement.
    """

    changed = Signal()

    def __init__(self, session_state, preset_name: str = None, parent=None):
        super().__init__(parent)
        self.session = session_state
        self.preset_name = preset_name
        self._init_ui()
        try:
            self._load_machines()
        except Exception:
            logger.debug('Could not load machines during init', exc_info=True)

    def _init_ui(self):
        main = QVBoxLayout(self)

        # Machine / version / code group
        env = QGroupBox('Machine')
        env_layout = QFormLayout(env)
        self.machine_combo = QComboBox()
        self.machine_combo.setEditable(True)
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        # Machine row: include short info label and a Configure button
        mrow = QHBoxLayout()
        mrow.addWidget(self.machine_combo, 1)
        self.machine_info_label = QLabel("")
        self.machine_info_label.setStyleSheet("color: #555; font-size: 11px;")
        mrow.addWidget(self.machine_info_label)
        cfg_btn = QPushButton("Configure...")
        cfg_btn.setFixedWidth(110)
        cfg_btn.clicked.connect(self._open_machine_config_dialog)
        mrow.addWidget(cfg_btn)
        env_layout.addRow('Machine:', mrow)

        self.version_combo = QComboBox()
        self.version_combo.currentTextChanged.connect(self._on_version_changed)
        env_layout.addRow('QE Version:', self.version_combo)

        self.code_combo = QComboBox()
        self.code_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        env_layout.addRow('Code:', self.code_combo)

        main.addWidget(env)

        # Pseudopotentials selector (optional)
        pseudo_group = QGroupBox('Pseudopotentials')
        pseudo_layout = QVBoxLayout(pseudo_group)
        if PSEUDO_SELECTOR_AVAILABLE:
            try:
                self.pseudo_selector = PseudopotentialsSelectorWidget(self.session)
                self.pseudo_selector.changed.connect(self.changed)
                pseudo_layout.addWidget(self.pseudo_selector)
            except Exception:
                logger.debug('Failed to init real pseudo selector', exc_info=True)
                self.pseudo_selector = FallbackPseudoSelector(self.session)
                self.pseudo_selector.changed.connect(self.changed)
                pseudo_layout.addWidget(self.pseudo_selector)
        else:
            self.pseudo_selector = FallbackPseudoSelector(self.session)
            self.pseudo_selector.changed.connect(self.changed)
            pseudo_layout.addWidget(self.pseudo_selector)
        main.addWidget(pseudo_group)

        # Protocol + label
        row = QHBoxLayout()
        row.addWidget(QLabel('Protocol:'))
        self.protocol_combo = QComboBox()
        # use PRESETS if available, otherwise fallbacks
        try:
            from xespresso.workflow.simple_workflow import PRESETS
            self.protocol_combo.addItems(sorted(PRESETS.keys()))
        except Exception:
            self.protocol_combo.addItems(['fast', 'moderate', 'accurate'])
        self.protocol_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        row.addWidget(self.protocol_combo)
        row.addWidget(QLabel('Label:'))
        self.label_edit = QLineEdit()
        self.label_edit.editingFinished.connect(self.changed)
        row.addWidget(self.label_edit)
        main.addLayout(row)

        # Basic parameters
        params = QGroupBox('Basic Parameters')
        pform = QFormLayout(params)
        self.ecutwfc_spin = QDoubleSpinBox(); self.ecutwfc_spin.setRange(10.0, 200.0); self.ecutwfc_spin.setValue(50.0)
        pform.addRow('Ecutwfc:', self.ecutwfc_spin)
        self.ecutrho_spin = QDoubleSpinBox(); self.ecutrho_spin.setRange(40.0, 1600.0); self.ecutrho_spin.setValue(400.0)
        pform.addRow('Ecutrho:', self.ecutrho_spin)
        self.kspacing_spin = QDoubleSpinBox(); self.kspacing_spin.setRange(0.01, 5.0); self.kspacing_spin.setDecimals(3); self.kspacing_spin.setValue(0.3)
        pform.addRow('K-spacing (Å⁻¹):', self.kspacing_spin)
        self.conv_thr_spin = QDoubleSpinBox(); self.conv_thr_spin.setRange(1e-12, 1e-2); self.conv_thr_spin.setDecimals(12); self.conv_thr_spin.setValue(1e-8)
        pform.addRow('conv_thr:', self.conv_thr_spin)
        self.forc_conv_thr_spin = QDoubleSpinBox(); self.forc_conv_thr_spin.setRange(1e-8, 1e-1); self.forc_conv_thr_spin.setDecimals(8); self.forc_conv_thr_spin.setValue(1e-3)
        pform.addRow('forc_conv_thr:', self.forc_conv_thr_spin)
        main.addWidget(params)

        # Wire changed signals
        try:
            self.ecutwfc_spin.valueChanged.connect(self.changed)
            self.ecutrho_spin.valueChanged.connect(self.changed)
            self.kspacing_spin.valueChanged.connect(self.changed)
            self.conv_thr_spin.valueChanged.connect(self.changed)
            self.forc_conv_thr_spin.valueChanged.connect(self.changed)
        except Exception:
            pass

    def _load_machines(self):
        if not XESPRESSO_AVAILABLE:
            return
        try:
            machines = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            for m in machines:
                self.machine_combo.addItem(m)
            cur = self.session.get('current_machine_name')
            if cur:
                idx = self.machine_combo.findText(cur)
                if idx >= 0:
                    self.machine_combo.setCurrentIndex(idx)
        except Exception:
            logger.debug('Failed to load machines', exc_info=True)

    def _on_machine_changed(self, machine_name):
        # update versions and codes when the machine changes
        if not machine_name or not XESPRESSO_AVAILABLE:
            return
        try:
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, verbose=False)
            self.version_combo.blockSignals(True)
            self.code_combo.blockSignals(True)
            try:
                self.version_combo.clear()
                self.code_combo.clear()
                if codes and getattr(codes, 'list_versions', None):
                    for v in codes.list_versions():
                        self.version_combo.addItem(v)
                else:
                    all_codes = codes.get_all_codes() if codes else {}
                    for cname in all_codes.keys():
                        self.code_combo.addItem(cname)
            finally:
                self.version_combo.blockSignals(False)
                self.code_combo.blockSignals(False)
            self.changed.emit()
        except Exception:
            logger.debug('Failed to load codes for machine', exc_info=True)

    def _on_version_changed(self, version):
        if not version or not XESPRESSO_AVAILABLE:
            return
        try:
            machine_name = self.machine_combo.currentText()
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, version=version, verbose=False)
            self.code_combo.blockSignals(True)
            try:
                self.code_combo.clear()
                if codes and getattr(codes, 'get_all_codes', None):
                    for code_name in codes.get_all_codes(version=version).keys():
                        self.code_combo.addItem(code_name)
            finally:
                self.code_combo.blockSignals(False)
            self.changed.emit()
        except Exception:
            logger.debug('Failed to update codes for version', exc_info=True)

    def update_for_structure(self, atoms):
        if atoms is None:
            return
        try:
            elements = set(atoms.get_chemical_symbols())
            if getattr(self, 'pseudo_selector', None) is not None:
                try:
                    self.pseudo_selector.set_elements(elements)
                except Exception:
                    logger.debug('Failed to update pseudo selector elements', exc_info=True)
            # populate label if empty
            try:
                if (not self.label_edit.text().strip()) and self.preset_name:
                    formula = None
                    try:
                        formula = atoms.get_chemical_formula()
                    except Exception:
                        pass
                    if formula:
                        self.label_edit.setText(f"{self.preset_name.lower()}/{formula}")
            except Exception:
                pass
        except Exception:
            logger.debug('Failed during update_for_structure', exc_info=True)

    def get_config(self):
        cfg = {}
        cfg['machine_name'] = self.machine_combo.currentText()
        cfg['qe_version'] = self.version_combo.currentText()
        cfg['selected_code'] = self.code_combo.currentText()
        cfg['label'] = self.label_edit.text().strip()
        cfg['ecutwfc'] = float(self.ecutwfc_spin.value())
        cfg['ecutrho'] = float(self.ecutrho_spin.value())
        try:
            cfg['kspacing'] = float(self.kspacing_spin.value())
        except Exception:
            pass
        try:
            cfg['conv_thr'] = float(self.conv_thr_spin.value())
        except Exception:
            pass
        try:
            cfg['forc_conv_thr'] = float(self.forc_conv_thr_spin.value())
        except Exception:
            pass
        try:
            cfg['pseudopotentials'] = self.pseudo_selector.get_pseudopotentials() if self.pseudo_selector is not None else {}
        except Exception:
            cfg['pseudopotentials'] = {}
        cfg['calc_type'] = cfg.get('calc_type', 'scf')
        return cfg

    def restore_from_dict(self, cfg: dict):
        if not cfg or not isinstance(cfg, dict):
            return
        try:
            if cfg.get('machine_name'):
                idx = self.machine_combo.findText(cfg['machine_name'])
                if idx >= 0:
                    self.machine_combo.setCurrentIndex(idx)
        except Exception:
            pass
        try:
            if cfg.get('qe_version'):
                idx = self.version_combo.findText(cfg['qe_version'])
                if idx >= 0:
                    self.version_combo.setCurrentIndex(idx)
        except Exception:
            pass
        try:
            if cfg.get('selected_code'):
                idx = self.code_combo.findText(cfg['selected_code'])
                if idx >= 0:
                    self.code_combo.setCurrentIndex(idx)
        except Exception:
            pass
        try:
            if cfg.get('label'):
                self.label_edit.setText(cfg['label'])
        except Exception:
            pass
        try:
            if 'ecutwfc' in cfg:
                self.ecutwfc_spin.setValue(float(cfg['ecutwfc']))
        except Exception:
            pass
        try:
            if 'ecutrho' in cfg:
                self.ecutrho_spin.setValue(float(cfg['ecutrho']))
        except Exception:
            pass
        try:
            if cfg.get('pseudopotentials') and getattr(self, 'pseudo_selector', None) is not None:
                try:
                    self.pseudo_selector.set_pseudopotentials(cfg['pseudopotentials'])
                except Exception:
                    pass
        except Exception:
            pass
import logging
from qtpy.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QComboBox,
    QPushButton, QGroupBox, QFormLayout, QDoubleSpinBox, QSpinBox,
    QCheckBox, QPlainTextEdit
)
from qtpy.QtCore import Signal

logger = logging.getLogger(__name__)

# Predefined magnetic moments and typical Hubbard U values (keep in sync with CalculationSetupPage)
PREDEFINED_MAGNETIC_MOMENTS = {
    'Fe': 2.2, 'Co': 1.7, 'Ni': 0.6, 'Mn': 5.0, 'Cr': 3.0,
    'V': 2.0, 'Ti': 1.0, 'Gd': 7.0, 'Nd': 3.0
}

TYPICAL_HUBBARD_U = {
    'Fe': 4.0, 'Co': 3.5, 'Ni': 3.0, 'Mn': 4.0, 'Cr': 3.5,
    'V': 3.0, 'Ti': 2.5, 'Cu': 4.0, 'Zn': 4.0,
    'Gd': 6.0, 'Nd': 5.0, 'Ce': 5.0, 'O': 0.0
}

try:
    from xespresso.machines.config.loader import list_machines, load_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR
    XESPRESSO_AVAILABLE = True
except Exception:
    list_machines = lambda *args, **kwargs: []
    load_machine = lambda *args, **kwargs: None
    load_codes_config = lambda *args, **kwargs: None
    XESPRESSO_AVAILABLE = False

try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except Exception:
    PSEUDO_SELECTOR_AVAILABLE = False

try:
    from xespresso.workflow.simple_workflow import PRESETS
    WORKFLOW_PRESETS_AVAILABLE = True
except Exception:
    PRESETS = {}
    WORKFLOW_PRESETS_AVAILABLE = False


class CalculationConfigWidget(QWidget):
    """Recreates the subset of CalculationSetupPage needed inside workflow tabs.

    Exposes `get_config()` and `restore_from_dict(cfg)`.
    """

    changed = Signal()

    def __init__(self, session_state, preset_name: str = None, parent=None):
        super().__init__(parent)
        self.session = session_state
        self.preset_name = preset_name
        self._init_ui()
        self._load_machines()

    def _init_ui(self):
        main = QVBoxLayout(self)

        # Top row: machine, pseudopotentials, QE version, code
        env = QGroupBox("Machine")
        env_layout = QFormLayout(env)
        self.machine_combo = QComboBox()
        self.machine_combo.setEditable(True)
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        env_layout.addRow("Machine:", self.machine_combo)

        # (Pseudopotentials selector will be placed in its own group below)
import logging
from qtpy.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit, QComboBox,
    QGroupBox, QFormLayout, QDoubleSpinBox
)
from qtpy.QtCore import Signal

logger = logging.getLogger(__name__)

# Defer heavy xespresso / pseudo selector imports: use guards so importing
# this widget during tests doesn't trigger optional dependency failures.
try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except Exception:
    PSEUDO_SELECTOR_AVAILABLE = False

try:
    from xespresso.machines.config.loader import list_machines, load_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR
    XESPRESSO_AVAILABLE = True
except Exception:
    list_machines = lambda *a, **k: []
    load_machine = lambda *a, **k: None
    load_codes_config = lambda *a, **k: None
    XESPRESSO_AVAILABLE = False


class CalculationConfigWidget(QWidget):
    """Configuration widget used inside workflow tabs (magnetism/hubbard omitted).

    This provides the public API expected by `WorkflowInstancePage`:
    - constructor(session_state, preset_name)
    - `get_config()` -> dict
    - `restore_from_dict(cfg)`
    - `update_for_structure(atoms)`
    """

    changed = Signal()

    def __init__(self, session_state, preset_name: str = None, parent=None):
        super().__init__(parent)
        self.session = session_state
        self.preset_name = preset_name
        self._init_ui()
        # try to load machines if available
        try:
            self._load_machines()
        except Exception:
            logger.debug('Could not load machines during init', exc_info=True)

    def _init_ui(self):
        main = QVBoxLayout(self)

        # Machine / version / code group
        env = QGroupBox('Machine')
        env_layout = QFormLayout(env)
        self.machine_combo = QComboBox()
        self.machine_combo.setEditable(True)
        self.machine_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        env_layout.addRow('Machine:', self.machine_combo)

        self.version_combo = QComboBox()
        self.version_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        env_layout.addRow('QE Version:', self.version_combo)

        self.code_combo = QComboBox()
        self.code_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        env_layout.addRow('Code:', self.code_combo)

        main.addWidget(env)

        # Pseudopotentials selector (optional)
        pseudo_group = QGroupBox('Pseudopotentials')
        pseudo_layout = QVBoxLayout(pseudo_group)
        if PSEUDO_SELECTOR_AVAILABLE:
            try:
                self.pseudo_selector = PseudopotentialsSelectorWidget(self.session)
                self.pseudo_selector.changed.connect(self.changed)
                pseudo_layout.addWidget(self.pseudo_selector)
            except Exception:
                # Fall back to the simple text-based selector if the real
                # widget failed to initialize.
                try:
                    self.pseudo_selector = FallbackPseudoSelector(self.session)
                    self.pseudo_selector.changed.connect(self.changed)
                    pseudo_layout.addWidget(self.pseudo_selector)
                except Exception:
                    self.pseudo_selector = None
        else:
            # Use the fallback selector when the richer widget isn't available
            try:
                self.pseudo_selector = FallbackPseudoSelector(self.session)
                self.pseudo_selector.changed.connect(self.changed)
                pseudo_layout.addWidget(self.pseudo_selector)
            except Exception:
                self.pseudo_selector = None
        main.addWidget(pseudo_group)

        # Protocol + label
        row = QHBoxLayout()
        row.addWidget(QLabel('Protocol:'))
        self.protocol_combo = QComboBox()
        self.protocol_combo.addItems(['fast', 'moderate', 'accurate'])
        self.protocol_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        row.addWidget(self.protocol_combo)
        row.addWidget(QLabel('Label:'))
        self.label_edit = QLineEdit()
        self.label_edit.editingFinished.connect(self.changed)
        row.addWidget(self.label_edit)
        main.addLayout(row)

        # Basic parameters
        params = QGroupBox('Basic Parameters')
        pform = QFormLayout(params)
        self.ecutwfc_spin = QDoubleSpinBox(); self.ecutwfc_spin.setRange(10.0, 200.0); self.ecutwfc_spin.setValue(50.0)
        pform.addRow('Ecutwfc:', self.ecutwfc_spin)
        self.ecutrho_spin = QDoubleSpinBox(); self.ecutrho_spin.setRange(40.0, 1600.0); self.ecutrho_spin.setValue(400.0)
        pform.addRow('Ecutrho:', self.ecutrho_spin)
        self.kspacing_spin = QDoubleSpinBox(); self.kspacing_spin.setRange(0.01, 5.0); self.kspacing_spin.setDecimals(3); self.kspacing_spin.setValue(0.3)
        pform.addRow('K-spacing (Å⁻¹):', self.kspacing_spin)
        self.conv_thr_spin = QDoubleSpinBox(); self.conv_thr_spin.setRange(1e-12, 1e-2); self.conv_thr_spin.setDecimals(12); self.conv_thr_spin.setValue(1e-8)
        pform.addRow('conv_thr:', self.conv_thr_spin)
        self.forc_conv_thr_spin = QDoubleSpinBox(); self.forc_conv_thr_spin.setRange(1e-8, 1e-1); self.forc_conv_thr_spin.setDecimals(8); self.forc_conv_thr_spin.setValue(1e-3)
        pform.addRow('forc_conv_thr:', self.forc_conv_thr_spin)
        main.addWidget(params)

        # Wire changed signals
        try:
            self.ecutwfc_spin.valueChanged.connect(self.changed)
            self.ecutrho_spin.valueChanged.connect(self.changed)
            self.kspacing_spin.valueChanged.connect(self.changed)
            self.conv_thr_spin.valueChanged.connect(self.changed)
            self.forc_conv_thr_spin.valueChanged.connect(self.changed)
        except Exception:
            pass

    def _load_machines(self):
        if not XESPRESSO_AVAILABLE:
            return
        try:
            machines = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            for m in machines:
                self.machine_combo.addItem(m)
            cur = self.session.get('current_machine_name')
            if cur:
                idx = self.machine_combo.findText(cur)
                if idx >= 0:
                    self.machine_combo.setCurrentIndex(idx)
        except Exception:
            logger.debug('Failed to load machines', exc_info=True)

    def _on_machine_changed(self, machine_name):
        if not machine_name or not XESPRESSO_AVAILABLE:
            return
        try:
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, verbose=False)
            self.version_combo.clear(); self.code_combo.clear()
            if codes and getattr(codes, 'versions', None):
                for v in codes.list_versions():
                    self.version_combo.addItem(v)
            else:
                all_codes = codes.get_all_codes() if codes else {}
                for cname in all_codes.keys():
                    self.code_combo.addItem(cname)
        except Exception:
            logger.debug('Failed to load codes for machine', exc_info=True)

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
            logger.debug('Failed to update codes for version', exc_info=True)

    def _open_machine_config_dialog(self):
        """Open the MachineConfigPage in a dialog so users can create/edit machines."""
        try:
            from qtgui.pages.machine_config import MachineConfigPage
        except Exception:
            return
        dlg = QDialog(self)
        dlg.setWindowTitle("Machine Configuration")
        layout = QVBoxLayout(dlg)
        page = MachineConfigPage(self.session)
        layout.addWidget(page)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.accept)
        layout.addWidget(close_btn)
        dlg.exec_()
        # reload machine list/info after possible edits
        try:
            self._load_machines()
        except Exception:
            logger.debug('Failed to reload machines after dialog', exc_info=True)

    def update_for_structure(self, atoms):
        if atoms is None:
            return
        try:
            elements = set(atoms.get_chemical_symbols())
            if getattr(self, 'pseudo_selector', None) is not None:
                try:
                    self.pseudo_selector.set_elements(elements)
                except Exception:
                    logger.debug('Failed to update pseudo selector elements', exc_info=True)
            # populate label if empty
            try:
                if (not self.label_edit.text().strip()) and self.preset_name:
                    formula = None
                    try:
                        formula = atoms.get_chemical_formula()
                    except Exception:
                        pass
                    if formula:
                        self.label_edit.setText(f"{self.preset_name.lower()}/{formula}")
            except Exception:
                pass
        except Exception:
            logger.debug('Failed during update_for_structure', exc_info=True)

    def get_config(self):
        cfg = {}
        cfg['machine_name'] = self.machine_combo.currentText()
        cfg['qe_version'] = self.version_combo.currentText()
        cfg['selected_code'] = self.code_combo.currentText()
        cfg['label'] = self.label_edit.text().strip()
        cfg['ecutwfc'] = float(self.ecutwfc_spin.value())
        cfg['ecutrho'] = float(self.ecutrho_spin.value())
        try:
            cfg['kspacing'] = float(self.kspacing_spin.value())
        except Exception:
            pass
        try:
            cfg['conv_thr'] = float(self.conv_thr_spin.value())
        except Exception:
            pass
        try:
            cfg['forc_conv_thr'] = float(self.forc_conv_thr_spin.value())
        except Exception:
            pass
        if getattr(self, 'pseudo_selector', None) is not None:
            try:
                cfg['pseudopotentials'] = self.pseudo_selector.get_pseudopotentials()
            except Exception:
                cfg['pseudopotentials'] = {}
        else:
            cfg['pseudopotentials'] = {}
        cfg['calc_type'] = cfg.get('calc_type', 'scf')
        return cfg

    def restore_from_dict(self, cfg: dict):
        if not cfg or not isinstance(cfg, dict):
            return
        try:
            if cfg.get('machine_name'):
                idx = self.machine_combo.findText(cfg['machine_name'])
                if idx >= 0:
                    self.machine_combo.setCurrentIndex(idx)
        except Exception:
            pass
        try:
            if cfg.get('qe_version'):
                idx = self.version_combo.findText(cfg['qe_version'])
                if idx >= 0:
                    self.version_combo.setCurrentIndex(idx)
        except Exception:
            pass
        try:
            if cfg.get('selected_code'):
                idx = self.code_combo.findText(cfg['selected_code'])
                if idx >= 0:
                    self.code_combo.setCurrentIndex(idx)
        except Exception:
            pass
        try:
            if cfg.get('label'):
                self.label_edit.setText(cfg['label'])
        except Exception:
            pass
        try:
            if 'ecutwfc' in cfg:
                self.ecutwfc_spin.setValue(float(cfg['ecutwfc']))
        except Exception:
            pass
        try:
            if 'ecutrho' in cfg:
                self.ecutrho_spin.setValue(float(cfg['ecutrho']))
        except Exception:
            pass
        try:
            if cfg.get('pseudopotentials') and getattr(self, 'pseudo_selector', None) is not None:
                try:
                    self.pseudo_selector.set_pseudopotentials(cfg['pseudopotentials'])
                except Exception:
                    pass
        except Exception:
            pass

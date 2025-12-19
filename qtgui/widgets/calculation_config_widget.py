"""
CalculationConfigWidget

Clean single-definition implementation mirroring pseudopotential UI from
`qtgui/pages/calculation_setup.py`. Ensures a single code path and adds
debug logging to help diagnose missing selector instances.
"""

import logging
from typing import Dict

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QLineEdit, QComboBox,
    QGroupBox, QFormLayout, QDoubleSpinBox, QHBoxLayout, QPushButton,
    QDialog
)
from qtpy.QtCore import Signal

logger = logging.getLogger(__name__)

# Optional selector import
try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except Exception:
    PSEUDO_SELECTOR_AVAILABLE = False

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


class CalculationConfigWidget(QWidget):
    """Configuration widget used inside workflow tabs.

    Mirrors `CalculationSetupPage` behavior for pseudopotentials: prefer
    `PseudopotentialsSelectorWidget` when available; otherwise provide
    per-element manual inputs.
    """

    changed = Signal()

    def __init__(self, session_state: Dict, preset_name: str = None, parent=None):
        super().__init__(parent)
        self.session = session_state or {}
        self.preset_name = preset_name
        self.pseudo_edits = {}
        # keep attribute present so later code can inspect it
        self.pseudo_selector = None
        self._init_ui()
        try:
            self._load_machines()
        except Exception:
            logger.debug('Could not load machines during init', exc_info=True)

    def _init_ui(self):
        main = QVBoxLayout(self)

        # Machine / version / code group (minimal)
        env = QGroupBox('Machine')
        env_layout = QFormLayout(env)
        self.machine_combo = QComboBox(); self.machine_combo.setEditable(True)
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        env_layout.addRow('Machine:', self.machine_combo)
        self.machine_info_label = QLabel("")
        env_layout.addRow(self.machine_info_label)
        self.version_combo = QComboBox(); self.version_combo.currentTextChanged.connect(self._on_version_changed)
        env_layout.addRow('QE Version:', self.version_combo)
        self.code_combo = QComboBox(); self.code_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        env_layout.addRow('Code:', self.code_combo)
        main.addWidget(env)

        # Pseudopotentials group
        pseudo_group = QGroupBox('🔬 Pseudopotentials')
        pseudo_layout = QVBoxLayout(pseudo_group)

        if PSEUDO_SELECTOR_AVAILABLE:
            # instantiate selector exactly as in CalculationSetupPage
            try:
                self.pseudo_selector = PseudopotentialsSelectorWidget(self.session)
                try:
                    self.pseudo_selector.changed.connect(self.changed)
                except Exception:
                    pass
                pseudo_layout.addWidget(self.pseudo_selector)
            except Exception:
                pass
                self.pseudo_selector = None

            # safeguard: ensure reference exists after parenting
            try:
                if getattr(self, 'pseudo_selector', None) is None:
                    found = self.findChild(PseudopotentialsSelectorWidget)
                    if found is not None:
                        self.pseudo_selector = found
            except Exception:
                pass

        else:
            # manual per-element inputs
            self.pseudo_selector = None
            self.pseudo_info_label = QLabel('Load a structure to configure pseudopotentials')
            pseudo_layout.addWidget(self.pseudo_info_label)
            self.pseudo_container = QWidget()
            self.pseudo_container_layout = QFormLayout(self.pseudo_container)
            pseudo_layout.addWidget(self.pseudo_container)

        main.addWidget(pseudo_group)

        # Protocol + label
        row = QHBoxLayout()
        row.addWidget(QLabel('Protocol:'))
        self.protocol_combo = QComboBox()
        try:
            from xespresso.workflow.simple_workflow import PRESETS
            self.protocol_combo.addItems(sorted(PRESETS.keys()))
        except Exception:
            self.protocol_combo.addItems(['fast', 'moderate', 'accurate'])
        self.protocol_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        row.addWidget(self.protocol_combo)
        row.addWidget(QLabel('Label:'))
        self.label_edit = QLineEdit(); self.label_edit.editingFinished.connect(self.changed)
        row.addWidget(self.label_edit)
        main.addLayout(row)

        # Basic parameters
        params = QGroupBox('Basic Parameters')
        pform = QFormLayout(params)
        self.ecutwfc_spin = QDoubleSpinBox(); self.ecutwfc_spin.setRange(10.0, 200.0); self.ecutwfc_spin.setValue(50.0)
        pform.addRow('Ecutwfc:', self.ecutwfc_spin)
        main.addWidget(params)

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
                try:
                    self._on_machine_changed(cur)
                except Exception:
                    logger.debug('Failed to initialize machine codes', exc_info=True)
        except Exception:
            logger.debug('Failed to load machines', exc_info=True)

    def _on_machine_changed(self, machine_name):
        if not machine_name or not XESPRESSO_AVAILABLE:
            return
        try:
            machine = load_machine(DEFAULT_CONFIG_PATH, machine_name, DEFAULT_MACHINES_DIR, return_object=True)
            try:
                self.session['calc_machine'] = machine
                self.session['selected_machine'] = machine_name
            except Exception:
                pass
            info = f"Type: {getattr(machine, 'execution', '')}"
            if getattr(machine, 'scheduler', None):
                info += f", Scheduler: {machine.scheduler}"
            self.machine_info_label.setText(info)
            try:
                self._load_codes(machine_name)
            except Exception:
                logger.debug('Failed to load codes for machine', exc_info=True)
        except Exception:
            logger.debug('Failed to handle machine change', exc_info=True)

    def _load_codes(self, machine_name):
        if not XESPRESSO_AVAILABLE:
            return
        self.version_combo.blockSignals(True)
        self.code_combo.blockSignals(True)
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
            logger.debug('Failed to load codes', exc_info=True)
        finally:
            try:
                self.version_combo.blockSignals(False)
                self.code_combo.blockSignals(False)
            except Exception:
                pass
        if self.version_combo.count() > 0:
            try:
                self._on_version_changed(self.version_combo.currentText())
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
            logger.debug('Failed to update codes for version', exc_info=True)

    def update_for_structure(self, atoms):
        """
        Update element-dependent inputs based on the provided `atoms` object.

        This mirrors the behavior in CalculationSetupPage: when a structure
        is available, compute its unique element symbols and populate the
        pseudopotential selector (or manual inputs) accordingly.
        """
        elements = set()
        if atoms is None:
            elements = set()
        else:
            try:
                # ASE Atoms provides get_chemical_symbols / get_chemical_formula
                symbols = []
                try:
                    symbols = atoms.get_chemical_symbols()
                except Exception:
                    # fallback: some objects expose .symbols or .get_elements
                    symbols = getattr(atoms, 'symbols', []) or getattr(atoms, 'get_elements', lambda: [])()
                elements = set(symbols)
            except Exception:
                elements = set()

        if self.pseudo_selector is not None:
            try:
                self.pseudo_selector.set_elements(elements)
            except Exception:
                logger.debug('Failed to set elements on pseudo_selector', exc_info=True)
        else:
            # manual inputs path
            try:
                self.pseudo_info_label.setText('Configure pseudopotentials for: ' + ', '.join(sorted(elements)) if elements else 'Load a structure to configure pseudopotentials')
                # rebuild manual inputs
                while self.pseudo_container_layout.count():
                    item = self.pseudo_container_layout.takeAt(0)
                    if item.widget():
                        item.widget().deleteLater()
                self.pseudo_edits = {}
                for el in sorted(elements):
                    le = QLineEdit()
                    le.setPlaceholderText(f"e.g., {el}.UPF")
                    self.pseudo_edits[el] = le
                    self.pseudo_container_layout.addRow(f"{el}:", le)
            except Exception:
                logger.debug('Failed to update manual pseudo inputs', exc_info=True)

    def get_config(self) -> Dict:
        """
        Collect current header/config values from the form.

        Returns a dict suitable for merging into a workflow configuration.
        """
        out = {}
        try:
            out['machine_name'] = self.machine_combo.currentText()
        except Exception:
            out['machine_name'] = None
        try:
            out['code'] = self.code_combo.currentText()
        except Exception:
            out['code'] = None
        try:
            out['protocol'] = self.protocol_combo.currentText()
        except Exception:
            out['protocol'] = None
        try:
            out['label'] = self.label_edit.text().strip() or None
        except Exception:
            out['label'] = None
        try:
            out['ecutwfc'] = float(self.ecutwfc_spin.value())
        except Exception:
            out['ecutwfc'] = None

        # pseudopotentials mapping
        try:
            if self.pseudo_selector is not None and hasattr(self.pseudo_selector, 'get_pseudopotentials'):
                out['pseudopotentials'] = self.pseudo_selector.get_pseudopotentials() or {}
            else:
                out['pseudopotentials'] = {k: v.text().strip() for k, v in self.pseudo_edits.items() if v.text().strip()}
        except Exception:
            out['pseudopotentials'] = {}

        return out

    def restore_from_dict(self, cfg: Dict):
        """
        Restore UI state from a config dictionary saved in session state.
        """
        if not cfg:
            return
        # restore begins
        try:
            if 'machine_name' in cfg and cfg.get('machine_name'):
                try:
                    self.machine_combo.setCurrentText(cfg.get('machine_name'))
                except Exception:
                    pass
        except Exception:
            pass
        try:
            if 'code' in cfg and cfg.get('code'):
                try:
                    self.code_combo.setCurrentText(cfg.get('code'))
                except Exception:
                    pass
        except Exception:
            pass

        try:
            if 'protocol' in cfg and cfg.get('protocol'):
                try:
                    self.protocol_combo.setCurrentText(cfg.get('protocol'))
                except Exception:
                    pass
        except Exception:
            pass

        try:
            if 'label' in cfg and cfg.get('label'):
                try:
                    self.label_edit.setText(cfg.get('label'))
                except Exception:
                    pass
        except Exception:
            pass

        # restore pseudopotentials mapping if present
        try:
            mapping = cfg.get('pseudopotentials') or {}
            if mapping:
                if self.pseudo_selector is not None and hasattr(self.pseudo_selector, 'set_pseudopotentials'):
                    try:
                        # ensure elements are set first
                        elems = set(mapping.keys())
                        try:
                            self.pseudo_selector.set_elements(elems)
                            # Ensure inputs are rebuilt by explicitly invoking
                            # the config change handler in case internal paths
                            # did not trigger UI rebuild earlier.
                            try:
                                self.pseudo_selector._on_config_changed(getattr(self.pseudo_selector, 'config_combo').currentText())
                            except Exception:
                                pass
                        except Exception:
                            pass
                        try:
                            self.pseudo_selector.set_pseudopotentials(mapping)
                        except Exception:
                            pass
                    except Exception:
                        logger.debug('Failed to restore pseudopotentials via selector', exc_info=True)
                else:
                    # manual path
                    for k, v in mapping.items():
                        if k in self.pseudo_edits:
                            try:
                                self.pseudo_edits[k].setText(v)
                            except Exception:
                                pass
        except Exception:
            logger.debug('Failed to restore pseudopotentials mapping', exc_info=True)

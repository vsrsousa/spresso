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
    QFormLayout, QDoubleSpinBox, QHBoxLayout, QPushButton,
    QDialog, QWidget as QtWidget, QGroupBox, QSpinBox,
    QRadioButton
)
from qtpy.QtCore import Signal, Qt

logger = logging.getLogger(__name__)

from qtgui.widgets.curtain_widget import CurtainWidget

# Use xespresso.tools.kpts_from_spacing for k-mesh computation
# k-mesh computation will import xespresso.tools.kpts_from_spacing lazily

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

        # Machine / version / code group (minimal) inside a curtain
        env_curtain = CurtainWidget('Execution Environment')
        env_widget = QtWidget()
        env_layout = QFormLayout(env_widget)
        self.machine_combo = QComboBox(); self.machine_combo.setEditable(True)
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        env_layout.addRow('Machine:', self.machine_combo)
        self.machine_info_label = QLabel("")
        env_layout.addRow(self.machine_info_label)
        self.version_combo = QComboBox(); self.version_combo.currentTextChanged.connect(self._on_version_changed)
        env_layout.addRow('QE Version:', self.version_combo)
        self.code_combo = QComboBox(); self.code_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        env_layout.addRow('Code:', self.code_combo)
        env_curtain.content_layout.addWidget(env_widget)
        main.addWidget(env_curtain)

        # Pseudopotentials group inside a curtain
        pseudo_curtain = CurtainWidget('Pseudopotentials')
        pseudo_widget = QtWidget()
        pseudo_layout = QVBoxLayout(pseudo_widget)

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
            self.pseudo_container = QtWidget()
            self.pseudo_container_layout = QFormLayout(self.pseudo_container)
            pseudo_layout.addWidget(self.pseudo_container)

        pseudo_curtain.content_layout.addWidget(pseudo_widget)
        main.addWidget(pseudo_curtain)

        # Label (separate row) - show above protocol per user request
        label_row = QHBoxLayout()
        label_row.addWidget(QLabel('Label:'))
        self.label_edit = QLineEdit(); self.label_edit.editingFinished.connect(self.changed)
        label_row.addWidget(self.label_edit)
        main.addLayout(label_row)

        # Protocol
        proto_row = QHBoxLayout()
        proto_row.addWidget(QLabel('Protocol:'))
        self.protocol_combo = QComboBox()
        try:
            from xespresso.workflow.simple_workflow import PRESETS
            self._PRESETS = PRESETS
            # Preserve meaningful ordering: fast, moderate, accurate
            desired = ['fast', 'moderate', 'accurate']
            ordered = [k for k in desired if k in PRESETS]
            # Append any additional presets after the desired ones (keep their original order)
            ordered += [k for k in PRESETS.keys() if k not in ordered]
            if ordered:
                self.protocol_combo.addItems(ordered)
            else:
                self.protocol_combo.addItems(list(PRESETS.keys()))
        except Exception:
            self._PRESETS = {}
            # fallback to the precision order
            self.protocol_combo.addItems(['fast', 'moderate', 'accurate'])
        self.protocol_combo.currentTextChanged.connect(self._on_protocol_changed)
        proto_row.addWidget(self.protocol_combo)
        main.addLayout(proto_row)

        # Basic parameters in a curtain
        params_curtain = CurtainWidget('Basic Parameters')
        params_widget = QtWidget()
        pform = QFormLayout(params_widget)
        # Compact layout: place related controls side-by-side to save space
        self.ecutwfc_spin = QDoubleSpinBox(); self.ecutwfc_spin.setRange(10.0, 200.0); self.ecutwfc_spin.setValue(50.0)
        self.ecutwfc_spin.valueChanged.connect(self._on_ecutwfc_changed)
        self.ecutwfc_spin.valueChanged.connect(lambda v: self.changed.emit())
        self.ecutrho_spin = QDoubleSpinBox(); self.ecutrho_spin.setRange(40.0, 1600.0); self.ecutrho_spin.setValue(400.0)
        self.ecutrho_spin.valueChanged.connect(lambda v: self.changed.emit())

        # Row 1: Ecutwfc | Ecutrho | Convergence
        self.conv_thr_edit = QLineEdit('1.0e-8')
        self.conv_thr_edit.editingFinished.connect(lambda: self.changed.emit())

        row1 = QtWidget(); row1_l = QHBoxLayout(row1); row1_l.setContentsMargins(0, 0, 0, 0); row1_l.setSpacing(6)
        lbl_ecutw = QLabel('Ecutwfc:'); lbl_ecutw.setFixedWidth(110)
        lbl_ecutrho = QLabel('Ecutrho:'); lbl_ecutrho.setFixedWidth(110)
        lbl_conv = QLabel('Convergence:'); lbl_conv.setFixedWidth(110)
        row1_l.addWidget(lbl_ecutw); row1_l.addWidget(self.ecutwfc_spin)
        row1_l.addWidget(lbl_ecutrho); row1_l.addWidget(self.ecutrho_spin)
        row1_l.addWidget(lbl_conv); row1_l.addWidget(self.conv_thr_edit)
        pform.addRow(row1)

        # Row 2: Occupations | Smearing | Degauss
        self.occupations_combo = QComboBox()
        self.occupations_combo.addItems(['smearing', 'fixed', 'tetrahedra'])
        self.occupations_combo.currentTextChanged.connect(self._on_occupations_changed)
        self.occupations_combo.currentTextChanged.connect(lambda v: self.changed.emit())

        self.smearing_combo = QComboBox()
        self.smearing_combo.addItems(['marzari-vanderbilt', 'gaussian', 'methfessel-paxton', 'fermi-dirac'])
        self.smearing_combo.currentTextChanged.connect(lambda v: self.changed.emit())
        self.degauss_spin = QDoubleSpinBox(); self.degauss_spin.setRange(0.001, 0.1); self.degauss_spin.setValue(0.02)
        self.degauss_spin.setSingleStep(0.005); self.degauss_spin.setDecimals(4)
        self.degauss_spin.valueChanged.connect(lambda v: self.changed.emit())

        row2 = QtWidget(); row2_l = QHBoxLayout(row2); row2_l.setContentsMargins(0, 0, 0, 0); row2_l.setSpacing(6)
        lbl_occ = QLabel('Occupations:'); lbl_occ.setFixedWidth(110)
        lbl_smear = QLabel('Smearing:'); lbl_smear.setFixedWidth(110)
        lbl_degauss = QLabel('Degauss:'); lbl_degauss.setFixedWidth(110)
        row2_l.addWidget(lbl_occ); row2_l.addWidget(self.occupations_combo)
        row2_l.addWidget(lbl_smear); row2_l.addWidget(self.smearing_combo)
        row2_l.addWidget(lbl_degauss); row2_l.addWidget(self.degauss_spin)
        pform.addRow(row2)

        # Kpoints mode and values compact
        # K-spacing: keep only k-spacing control and show computed k-mesh
        kpts_widget = QtWidget(); kpts_layout = QHBoxLayout(kpts_widget); kpts_layout.setContentsMargins(0, 0, 0, 0)
        # Only K-spacing control (explicit mode not supported here)
        self.kspacing_spin = QDoubleSpinBox(); self.kspacing_spin.setRange(0.05, 2.0); self.kspacing_spin.setValue(0.3)
        self.kspacing_spin.setSingleStep(0.01)
        self.kspacing_spin.valueChanged.connect(lambda v: (self._update_kmesh_display(), self.changed.emit()))
        kpts_layout.addWidget(QLabel('K-spacing:'))
        kpts_layout.addWidget(self.kspacing_spin)

        # computed k-mesh label (only output, no k1/k2/k3 controls)
        self.kmesh_label = QLabel('k-mesh: -')
        kpts_layout.addWidget(self.kmesh_label)

        # No explicit k1/k2/k3 controls — only show computed k-mesh

        pform.addRow(kpts_widget)
        params_curtain.content_layout.addWidget(params_widget)
        main.addWidget(params_curtain)

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

    def _on_protocol_changed(self, proto: str):
        """Populate basic parameters from selected protocol preset."""
        try:
            preset = (self._PRESETS or {}).get(proto, None)
            if preset is None:
                return
            if 'ecutwfc' in preset:
                try:
                    self.ecutwfc_spin.setValue(float(preset.get('ecutwfc')))
                except Exception:
                    pass
            if 'ecutrho' in preset and hasattr(self, 'ecutrho_spin'):
                try:
                    self.ecutrho_spin.setValue(float(preset.get('ecutrho')))
                except Exception:
                    pass
            if 'conv_thr' in preset and hasattr(self, 'conv_thr_edit'):
                try:
                    self.conv_thr_edit.setText(str(preset.get('conv_thr')))
                except Exception:
                    pass
            # additional mappings from CalculationSetupPage PRESETS
            if 'occupations' in preset and hasattr(self, 'occupations_combo'):
                try:
                    occ = preset.get('occupations')
                    if occ in ['smearing', 'fixed', 'tetrahedra']:
                        self.occupations_combo.setCurrentText(occ)
                except Exception:
                    pass
            if 'smearing' in preset and hasattr(self, 'smearing_combo'):
                try:
                    s = preset.get('smearing')
                    if s:
                        self.smearing_combo.setCurrentText(s)
                except Exception:
                    pass
            if 'degauss' in preset and hasattr(self, 'degauss_spin'):
                try:
                    self.degauss_spin.setValue(float(preset.get('degauss')))
                except Exception:
                    pass
            if 'kspacing' in preset and hasattr(self, 'kspacing_spin'):
                try:
                    self.kspacing_spin.setValue(float(preset.get('kspacing')))
                except Exception:
                    pass
            if 'kpoints' in preset and isinstance(preset.get('kpoints'), (list, tuple)):
                try:
                    k = preset.get('kpoints')
                    if len(k) >= 3 and hasattr(self, 'k1_spin'):
                        self.k1_spin.setValue(int(k[0])); self.k2_spin.setValue(int(k[1])); self.k3_spin.setValue(int(k[2]))
                except Exception:
                    pass
            if 'mixing_beta' in preset and hasattr(self, 'mixing_beta'):
                try:
                    self.mixing_beta.setValue(float(preset.get('mixing_beta')))
                except Exception:
                    pass
            if 'electron_maxstep' in preset and hasattr(self, 'electron_maxstep'):
                try:
                    self.electron_maxstep.setValue(int(preset.get('electron_maxstep')))
                except Exception:
                    pass
            # emit changed so draft saves
            try:
                self.changed.emit()
            except Exception:
                pass
        except Exception:
            pass

    def _on_ecutwfc_changed(self, value: float):
        """Update ecutrho when ecutwfc changes (same rule as CalculationSetupPage)."""
        try:
            if hasattr(self, 'ecutrho_spin'):
                self.ecutrho_spin.setValue(value * 8)
        except Exception:
            pass

    def _on_occupations_changed(self, occupations: str):
        try:
            if hasattr(self, 'smearing_combo') and hasattr(self, 'degauss_spin') and hasattr(self, 'smearing_combo'):
                # show/hide smearing-specific controls
                visible = (occupations == 'smearing')
                try:
                    self.smearing_combo.setVisible(visible)
                    self.degauss_spin.setVisible(visible)
                except Exception:
                    pass
        except Exception:
            pass

    def _on_kpts_mode_changed(self, checked: bool):
        # No-op: explicit mode not supported in this simplified UI
        return

    def _update_kmesh_display(self):
        """Compute and display the k-mesh corresponding to the current k-spacing.

        Uses `kpts_from_spacing` (xespresso.tools) if available and requires a
        structure present in `self._last_atoms` or `self.session['current_structure']`.
        If unavailable, shows '-' placeholder.
        """
        try:
            kspacing = float(getattr(self, 'kspacing_spin', None).value()) if getattr(self, 'kspacing_spin', None) is not None else None
        except Exception:
            kspacing = None
        atoms = getattr(self, '_last_atoms', None) or self.session.get('current_structure')
        if kspacing is None or atoms is None:
            try:
                if hasattr(self, 'kmesh_label'):
                    self.kmesh_label.setText('k-mesh: -')
            except Exception:
                pass
            return

        # compute grid using xespresso.kpts_from_spacing (no fallbacks)
        try:
            from xespresso import kpts_from_spacing
            # call the helper as intended by xespresso; it returns a tuple
            kmesh = kpts_from_spacing(atoms, kspacing)
        except Exception:
            kmesh = None

        try:
            if kmesh is None:
                self.kmesh_label.setText('k-mesh: -')
            else:
                # format as returned by helper (typically a tuple)
                self.kmesh_label.setText(f'k-mesh: {kmesh}')
        except Exception:
            pass

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

        # Auto-generate label similar to CalculationSetupPage: formula/protocol
        try:
            if atoms is not None:
                formula = None
                try:
                    formula = atoms.get_chemical_formula()
                except Exception:
                    try:
                        # fallback to chemical symbols join
                        syms = getattr(atoms, 'get_chemical_symbols', None)
                        if callable(syms):
                            formula = ''.join(sorted(set(syms())))
                    except Exception:
                        formula = None
                if formula:
                    current_label = (self.label_edit.text() or '').strip()
                    # if empty or looks like a previous auto-generated label, replace
                    proto = (self.protocol_combo.currentText() or '').strip() or 'scf'
                    auto = f"{formula}/{proto}"
                    if not current_label or '/' in current_label:
                        try:
                            self.label_edit.setText(auto)
                            try:
                                self.changed.emit()
                            except Exception:
                                pass
                        except Exception:
                            pass
        except Exception:
            pass

        # remember last atoms for k-mesh computation
        try:
            self._last_atoms = atoms
            # update computed k-mesh display whenever structure changes
            try:
                self._update_kmesh_display()
            except Exception:
                pass
        except Exception:
            pass

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

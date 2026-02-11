"""Single, clean CalculationWindow implementation.

This file was overwritten to remove duplicated and corrupted sections and
to ensure a single `CalculationWindow` class provides the requested
tabs. It intentionally keeps logic defensive to work in minimal test
environments.
"""
import os
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QFormLayout,
    QComboBox, QLineEdit, QTextEdit, QPushButton, QCheckBox,
    QTabWidget, QApplication, QMessageBox, QGroupBox
)
from qtpy.QtCore import Qt

try:
    from xespresso.workflow.simple_workflow import PRESETS
except Exception:
    PRESETS = {}

try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except Exception:
    PSEUDO_SELECTOR_AVAILABLE = False

try:
    from xespresso.machines.config.loader import list_machines, load_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
except Exception:
    list_machines = lambda *a, **k: []
    load_machine = lambda *a, **k: None
    DEFAULT_CONFIG_PATH = None
    DEFAULT_MACHINES_DIR = None

try:
    from xespresso.codes.manager import load_codes_config
    CODES_AVAILABLE = True
except Exception:
    load_codes_config = lambda *a, **k: None
    CODES_AVAILABLE = False

try:
    from xespresso.utils.auth import test_ssh_connection
except Exception:
    def test_ssh_connection(username, host, key_path=None, port=22):
        import subprocess
        key_path = os.path.expanduser(key_path) if key_path else None
        cmd = ["ssh", "-p", str(port), "-o", "PasswordAuthentication=no", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5"]
        if key_path:
            cmd += ["-i", key_path]
        cmd += [f"{username}@{host}", "echo 'Connection successful'"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except Exception:
            return False


class CalculationWindow(QWidget):
    """Modeless calculation configuration window with required tabs."""

    def __init__(self, calc_name: str, session_state=None, parent=None):
        super().__init__(parent)
        self.calc_name = calc_name
        self.session_state = session_state or {}
        self.setWindowTitle(f"Calculation: {calc_name}")
        self.resize(900, 700)

        # Determine calculation type from calc_name
        self.calculation_type = self._get_calculation_type_from_name(calc_name)

        main_layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        self._build_machine_tab()
        self._build_pseudopotentials_tab()
        self._build_basic_tab()
        self._build_magnetism_tab()
        self._build_hubbard_tab()
        self._build_prepare_tab()
        self._build_submit_tab()

        # Listen for structure changes to update pseudopotentials
        if self.session_state and hasattr(self.session_state, 'add_listener'):
            self.session_state.add_listener(self._update_pseudopotentials_for_structure)

    def _get_calculation_type_from_name(self, calc_name: str) -> str:
        """Map calculation window name to Quantum ESPRESSO calculation type."""
        name_lower = calc_name.lower()
        
        # Map button names to calculation types
        if 'scf' in name_lower:
            return 'scf'
        elif 'relax' in name_lower:
            return 'relax'
        elif 'geometry optimization' in name_lower or 'vc-relax' in name_lower:
            return 'vc-relax'
        elif 'md' in name_lower or 'molecular dynamics' in name_lower:
            return 'md'
        elif 'neb' in name_lower:
            return 'neb'
        else:
            # Default to SCF for unknown types
            return 'scf'

    def _build_machine_tab(self):
        w = QWidget()
        form = QFormLayout(w)

        self.machine_combo = QComboBox()
        try:
            for m in list_machines() or []:
                self.machine_combo.addItem(str(m))
            cur = (self.session_state.get('current_machine_name') if self.session_state else None)
            if cur and self.machine_combo.findText(cur) >= 0:
                self.machine_combo.setCurrentText(cur)
        except Exception:
            pass
        form.addRow('Machine:', self.machine_combo)

        self.version_combo = QComboBox()
        self.codes_combo = QComboBox()
        row = QWidget(); row_l = QHBoxLayout(row); row_l.setContentsMargins(0,0,0,0)
        row_l.addWidget(QLabel('Version:')); row_l.addWidget(self.version_combo)
        row_l.addWidget(QLabel('Code:')); row_l.addWidget(self.codes_combo)
        form.addRow(row)

        self.test_conn_btn = QPushButton('Test Connection')
        self.results_label = QLabel('')
        self.results_label.setWordWrap(True)
        form.addRow(self.test_conn_btn)
        form.addRow(self.results_label)

        try:
            self.test_conn_btn.clicked.connect(self._machine_test_connection)
        except Exception:
            pass

        try:
            self.machine_combo.currentTextChanged.connect(lambda name: self._load_codes_for_machine(name))
        except Exception:
            pass

        # Load codes for the initially selected machine
        try:
            current_machine = self.machine_combo.currentText()
            if current_machine:
                self._load_codes_for_machine(current_machine)
        except Exception:
            pass

        self.tabs.addTab(w, 'Machine')

    def _machine_test_connection(self):
        """Test the connection to the machine."""
        mn = self.machine_combo.currentText() or ''
        if not mn and self.session_state:
            mn = self.session_state.get('current_machine_name') or ''

        machine_obj = None
        try:
            machine_obj = load_machine(DEFAULT_CONFIG_PATH, mn, DEFAULT_MACHINES_DIR, return_object=True)
        except Exception:
            try:
                machine_obj = load_machine(mn)
            except Exception:
                machine_obj = None

        if machine_obj is None:
            self.results_label.setText("❌ Failed to load machine configuration")
            self.results_label.setStyleSheet("color: red;")
            return

        execution = getattr(machine_obj, 'execution', 'remote')
        
        if execution == "local":
            workdir = getattr(machine_obj, 'workdir', './calculations')
            user = os.environ.get('USER', 'unknown')
            self.results_label.setText(
                f"✅ Local machine - connection OK\n"
                f"Working directory: {workdir}\n"
                f"Current user: {user}"
            )
            self.results_label.setStyleSheet("color: green;")
        else:
            # Test remote connection
            host = getattr(machine_obj, 'host', None)
            username = getattr(machine_obj, 'username', None)
            port = getattr(machine_obj, 'port', 22)
            auth = getattr(machine_obj, 'auth', None)
            ssh_key = None
            if isinstance(auth, dict):
                ssh_key = auth.get('ssh_key') or auth.get('key')
            
            if not host or not username:
                QMessageBox.warning(self, "Warning", "Please enter host and username")
                return
            
            key_path = os.path.expanduser(ssh_key) if ssh_key else None
            
            if not key_path or not os.path.isfile(key_path):
                self.results_label.setText(f"❌ SSH key not found: {key_path}\n💡 Check the SSH key path")
                self.results_label.setStyleSheet("color: red;")
                return
            
            self.results_label.setText("Testing SSH connection...")
            self.results_label.setStyleSheet("color: blue;")
            
            # Process events to update UI
            QApplication.processEvents()
            
            try:
                success = test_ssh_connection(username, host, key_path, port)
                
                if success:
                    self.results_label.setText(
                        f"✅ SSH connection successful!\n"
                        f"Connected to: {username}@{host}:{port}"
                    )
                    self.results_label.setStyleSheet("color: green;")
                else:
                    self.results_label.setText(
                        "❌ SSH connection failed.\n"
                        "💡 Check your credentials and SSH key configuration"
                    )
                    self.results_label.setStyleSheet("color: red;")
            except Exception as e:
                self.results_label.setText(f"❌ Test failed: {e}")
                self.results_label.setStyleSheet("color: red;")

    def _build_pseudopotentials_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        
        if PSEUDO_SELECTOR_AVAILABLE:
            self.pseudo_selector = PseudopotentialsSelectorWidget(self.session_state)
            self.pseudo_selector.changed.connect(self._on_pseudopotentials_changed)
            layout.addWidget(self.pseudo_selector)
            
            # Update pseudopotentials for current structure
            self._update_pseudopotentials_for_structure()
        else:
            # Fallback to simple text editor
            form = QFormLayout()
            self.pseudo_editor = QTextEdit()
            form.addRow('Pseudopotentials:', self.pseudo_editor)
            layout.addLayout(form)
        
        self.tabs.addTab(w, 'Pseudopotentials')

    def _update_pseudopotentials_for_structure(self):
        """Update pseudopotentials selector based on current structure."""
        if not PSEUDO_SELECTOR_AVAILABLE:
            return
            
        try:
            atoms = self.session_state.get('current_structure')
            if atoms is not None:
                # Get unique elements from structure
                symbols = atoms.get_chemical_symbols()
                elements = set(symbols)
                self.pseudo_selector.set_elements(elements)
            else:
                self.pseudo_selector.set_elements(set())
        except Exception as e:
            # If anything fails, just clear elements
            try:
                self.pseudo_selector.set_elements(set())
            except Exception:
                pass

    def _on_pseudopotentials_changed(self):
        """Handle pseudopotentials configuration changes."""
        # Could add validation or other logic here if needed
        pass

    def _build_basic_tab(self):
        w = QWidget()
        form = QFormLayout(w)
        
        # Protocol selection
        self.protocol_combo = QComboBox()
        try:
            self.protocol_combo.addItems(list(PRESETS.keys()))
            self.protocol_combo.currentTextChanged.connect(self._on_protocol_changed)
            # Set fast as the default selection
            fast_index = self.protocol_combo.findText('fast')
            if fast_index >= 0:
                self.protocol_combo.setCurrentIndex(fast_index)
        except Exception:
            pass
        form.addRow('Protocol:', self.protocol_combo)
        
        # Energy cutoffs
        self.ecutwfc_edit = QLineEdit('50')
        form.addRow('ecutwfc (Ry):', self.ecutwfc_edit)
        
        self.ecutrho_edit = QLineEdit('400')
        form.addRow('ecutrho (Ry):', self.ecutrho_edit)
        
        # Convergence and mixing
        self.conv_thr_edit = QLineEdit('1.0e-8')
        form.addRow('conv_thr:', self.conv_thr_edit)
        
        self.mixing_beta_edit = QLineEdit('0.5')
        form.addRow('mixing_beta:', self.mixing_beta_edit)
        
        # K-points and SCF
        self.kspacing_edit = QLineEdit('0.3')
        form.addRow('kspacing (Å⁻¹):', self.kspacing_edit)
        
        self.electron_maxstep_edit = QLineEdit('200')
        form.addRow('electron_maxstep:', self.electron_maxstep_edit)
        
        # Additional SCF parameters
        self.verbosity_combo = QComboBox()
        self.verbosity_combo.addItems(['low', 'high'])
        self.verbosity_combo.setCurrentText('low')
        form.addRow('verbosity:', self.verbosity_combo)
        
        self.restart_mode_combo = QComboBox()
        self.restart_mode_combo.addItems(['from_scratch', 'restart'])
        self.restart_mode_combo.setCurrentText('from_scratch')
        form.addRow('restart_mode:', self.restart_mode_combo)
        
        self.disk_io_combo = QComboBox()
        self.disk_io_combo.addItems(['low', 'medium', 'high', 'none'])
        self.disk_io_combo.setCurrentText('low')
        form.addRow('disk_io:', self.disk_io_combo)
        
        self.mixing_mode_combo = QComboBox()
        self.mixing_mode_combo.addItems(['plain', 'TF', 'local-TF'])
        self.mixing_mode_combo.setCurrentText('plain')
        form.addRow('mixing_mode:', self.mixing_mode_combo)
        
        # Forces and stress calculation
        self.calc_forces_check = QCheckBox('Calculate forces')
        self.calc_forces_check.setChecked(False)  # Disabled by default for speed
        form.addRow('', self.calc_forces_check)
        
        self.calc_stress_check = QCheckBox('Calculate stress')
        self.calc_stress_check.setChecked(False)  # Disabled by default for speed
        form.addRow('', self.calc_stress_check)
        
        # Relaxation parameters
        self.ion_dynamics_combo = QComboBox()
        self.ion_dynamics_combo.addItems(['bfgs', 'damp', 'verlet', 'langevin', 'none'])
        self.ion_dynamics_combo.setCurrentText('bfgs')
        form.addRow('ion_dynamics:', self.ion_dynamics_combo)
        
        self.cell_dynamics_combo = QComboBox()
        self.cell_dynamics_combo.addItems(['none', 'bfgs', 'damp-pr', 'damp-w'])
        self.cell_dynamics_combo.setCurrentText('none')
        form.addRow('cell_dynamics:', self.cell_dynamics_combo)
        
        self.cell_dofree_combo = QComboBox()
        self.cell_dofree_combo.addItems(['all', 'shape', 'volume', 'x', 'y', 'z', 'xy', 'xz', 'yz'])
        self.cell_dofree_combo.setCurrentText('all')
        form.addRow('cell_dofree:', self.cell_dofree_combo)
        
        self.press_edit = QLineEdit('0.0')
        form.addRow('press (kbar):', self.press_edit)
        
        self.press_conv_test_edit = QLineEdit('0.5')
        form.addRow('press_conv_test (kbar):', self.press_conv_test_edit)
        
        # Initialize with fast preset (default)
        try:
            self._on_protocol_changed('fast')
        except Exception:
            pass
        
        # Initialize calculation type visibility based on the determined type
        try:
            self._on_calculation_type_changed(self.calculation_type)
        except Exception:
            pass
        
        self.tabs.addTab(w, 'Basic Parameters')

    def _build_magnetism_tab(self):
        w = QWidget()
        form = QFormLayout(w)
        self.magnetism_chk = QCheckBox('Enable magnetism')
        form.addRow(self.magnetism_chk)
        self.tabs.addTab(w, 'Magnetism')

    def _build_hubbard_tab(self):
        w = QWidget()
        form = QFormLayout(w)
        self.hubbard_text = QTextEdit()
        self.hubbard_text.setPlaceholderText('Enter Hubbard U values e.g. Fe:5.3')
        form.addRow('Hubbard U:', self.hubbard_text)
        self.tabs.addTab(w, 'Hubbard')

    def _build_prepare_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        box = QGroupBox('Preparation Actions')
        bl = QVBoxLayout(box)
        bl.addWidget(QLabel('Actions: generate input, validate settings, dry-run'))
        layout.addWidget(box)
        self.tabs.addTab(w, 'Prepare')

    def _build_submit_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        info = QLabel('Job submission configuration and controls')
        layout.addWidget(info)
        row = QHBoxLayout()
        self.submit_btn = QPushButton('Submit Job')
        self.close_btn = QPushButton('Close')
        row.addWidget(self.submit_btn)
        row.addWidget(self.close_btn)
        layout.addLayout(row)
        try:
            self.close_btn.clicked.connect(self.close)
            self.submit_btn.clicked.connect(self._on_submit)
        except Exception:
            pass
        self.tabs.addTab(w, 'Submit')

    def _load_codes_for_machine(self, machine_name: str):
        try:
            self.codes_combo.clear()
            try:
                self.version_combo.clear()
            except Exception:
                pass

            if not CODES_AVAILABLE:
                try:
                    self.version_combo.setVisible(False)
                except Exception:
                    pass
                return

            cfg = load_codes_config(machine_name)
            if not cfg:
                try:
                    self.version_combo.setVisible(False)
                except Exception:
                    pass
                return

            try:
                versions = cfg.list_versions() if hasattr(cfg, 'list_versions') else (getattr(cfg, 'versions', {}) or {}).keys()
                versions = list(versions) if versions is not None else []
            except Exception:
                versions = []

            if versions:
                try:
                    self.version_combo.setVisible(True)
                except Exception:
                    pass
                self.version_combo.blockSignals(True)
                self.version_combo.clear()
                for v in versions:
                    try:
                        self.version_combo.addItem(str(v))
                    except Exception:
                        pass
                try:
                    default_v = getattr(cfg, 'qe_version', None)
                    if default_v:
                        idx = self.version_combo.findText(str(default_v))
                        if idx >= 0:
                            self.version_combo.setCurrentIndex(idx)
                except Exception:
                    pass
                self.version_combo.blockSignals(False)
                try:
                    sel_v = self.version_combo.currentText() or (versions[0] if versions else None)
                    vcodes = cfg.get_all_codes(version=sel_v) if hasattr(cfg, 'get_all_codes') else {}
                    self.codes_combo.clear()
                    for name in (vcodes.keys() if isinstance(vcodes, dict) else list(vcodes)):
                        try:
                            self.codes_combo.addItem(str(name), (sel_v, name))
                        except Exception:
                            self.codes_combo.addItem(str(name))
                except Exception:
                    pass
                try:
                    self.version_combo.currentTextChanged.disconnect()
                except Exception:
                    pass
                try:
                    self.version_combo.currentTextChanged.connect(lambda v: self._on_version_changed(v, cfg))
                except Exception:
                    pass
            else:
                try:
                    self.version_combo.setVisible(False)
                except Exception:
                    pass
                all_codes = cfg.get_all_codes() if hasattr(cfg, 'get_all_codes') else getattr(cfg, 'codes', {})
                self.codes_combo.clear()
                if isinstance(all_codes, dict):
                    for k in all_codes.keys():
                        try:
                            self.codes_combo.addItem(str(k), (None, k))
                        except Exception:
                            self.codes_combo.addItem(str(k))
                else:
                    for k in list(all_codes):
                        try:
                            self.codes_combo.addItem(str(k), (None, k))
                        except Exception:
                            self.codes_combo.addItem(str(k))
        except Exception:
            pass

    def _on_version_changed(self, version, cfg):
        try:
            self.codes_combo.clear()
            if not cfg:
                return
            try:
                vcodes = cfg.get_all_codes(version=version) if hasattr(cfg, 'get_all_codes') else {}
                for name in (vcodes.keys() if isinstance(vcodes, dict) else list(vcodes)):
                    try:
                        self.codes_combo.addItem(str(name), (version, name))
                    except Exception:
                        self.codes_combo.addItem(str(name))
            except Exception:
                pass
        except Exception:
            pass

    def _on_submit(self):
        try:
            params = {
                'calc_name': self.calc_name,
                'machine': self.machine_combo.currentText() if hasattr(self, 'machine_combo') else None,
            }
            try:
                raw = self.pseudo_editor.toPlainText()
                mapping = {}
                for part in raw.split(','):
                    part = part.strip()
                    if not part:
                        continue
                    if '=' in part:
                        k, v = part.split('=', 1)
                        mapping[k.strip()] = v.strip()
                    elif ':' in part:
                        k, v = part.split(':', 1)
                        mapping[k.strip()] = v.strip()
                params['pseudopotentials'] = mapping
            except Exception:
                pass
            try:
                if hasattr(self.session_state, 'update'):
                    self.session_state.update({'_last_submission': params})
                else:
                    self.session_state['_last_submission'] = params
            except Exception:
                pass
        except Exception:
            pass

    def _on_protocol_changed(self, proto: str):
        try:
            preset = PRESETS.get(proto) if PRESETS else None
            if not preset:
                return
            
            # Update all parameter fields based on preset
            if hasattr(self, 'ecutwfc_edit') and 'ecutwfc' in preset:
                self.ecutwfc_edit.setText(str(preset['ecutwfc']))
            
            if hasattr(self, 'ecutrho_edit') and 'ecutrho' in preset:
                self.ecutrho_edit.setText(str(preset['ecutrho']))
            
            if hasattr(self, 'conv_thr_edit') and 'conv_thr' in preset:
                self.conv_thr_edit.setText(f"{preset['conv_thr']:.0e}")
            
            if hasattr(self, 'mixing_beta_edit') and 'mixing_beta' in preset:
                self.mixing_beta_edit.setText(str(preset['mixing_beta']))
            
            if hasattr(self, 'kspacing_edit') and 'kspacing' in preset:
                self.kspacing_edit.setText(str(preset['kspacing']))
            
            if hasattr(self, 'electron_maxstep_edit') and 'electron_maxstep' in preset:
                self.electron_maxstep_edit.setText(str(preset['electron_maxstep']))
                
        except Exception:
            pass

    def _on_calculation_type_changed(self, calc_type: str):
        """Show/hide parameters based on calculation type."""
        try:
            # Relaxation parameters are only relevant for relax, vc-relax, and md calculations
            show_relaxation = calc_type in ['relax', 'vc-relax', 'md']
            
            # Hide/show relaxation parameter widgets
            relaxation_widgets = [
                self.ion_dynamics_combo, self.cell_dynamics_combo, 
                self.cell_dofree_combo, self.press_edit, self.press_conv_test_edit
            ]
            
            for widget in relaxation_widgets:
                if hasattr(widget, 'setVisible'):
                    widget.setVisible(show_relaxation)
                    # Also hide the label by finding the parent layout item
                    try:
                        # Get the form layout and find the row containing this widget
                        form_layout = widget.parent().layout()
                        if isinstance(form_layout, QFormLayout):
                            # Find the row index for this widget
                            for i in range(form_layout.rowCount()):
                                field_item = form_layout.itemAt(i, QFormLayout.FieldRole)
                                if field_item and field_item.widget() == widget:
                                    label_item = form_layout.itemAt(i, QFormLayout.LabelRole)
                                    if label_item and label_item.widget():
                                        label_item.widget().setVisible(show_relaxation)
                                    break
                    except Exception:
                        pass
                        
        except Exception:
            pass

    def get_pseudopotentials(self):
        """Get the current pseudopotentials configuration as a dict."""
        if PSEUDO_SELECTOR_AVAILABLE and hasattr(self, 'pseudo_selector'):
            try:
                return self.pseudo_selector.get_pseudopotentials()
            except Exception:
                pass
        
        # Fallback to manual text input
        if hasattr(self, 'pseudo_editor'):
            try:
                text = self.pseudo_editor.toPlainText().strip()
                if text:
                    # Parse simple format like "Fe=Fe.UPF\nO=O.UPF"
                    pseudo_dict = {}
                    for line in text.split('\n'):
                        line = line.strip()
                        if '=' in line:
                            element, pseudo = line.split('=', 1)
                            pseudo_dict[element.strip()] = pseudo.strip()
                    return pseudo_dict
            except Exception:
                pass
        
        return {}

    def get_basic_parameters(self):
        """Get the current basic calculation parameters as a dict."""
        params = {}
        
        try:
            # Protocol
            if hasattr(self, 'protocol_combo'):
                protocol = self.protocol_combo.currentText()
                if protocol:
                    params['protocol'] = protocol
            
            # Calculation type (determined from window name)
            if hasattr(self, 'calculation_type') and self.calculation_type:
                params['calculation'] = self.calculation_type
            
            # Energy cutoffs
            if hasattr(self, 'ecutwfc_edit') and self.ecutwfc_edit.text():
                params['ecutwfc'] = float(self.ecutwfc_edit.text())
            
            if hasattr(self, 'ecutrho_edit') and self.ecutrho_edit.text():
                params['ecutrho'] = float(self.ecutrho_edit.text())
            
            # Convergence and mixing
            if hasattr(self, 'conv_thr_edit') and self.conv_thr_edit.text():
                params['conv_thr'] = float(self.conv_thr_edit.text())
            
            if hasattr(self, 'mixing_beta_edit') and self.mixing_beta_edit.text():
                params['mixing_beta'] = float(self.mixing_beta_edit.text())
            
            # K-points and SCF
            if hasattr(self, 'kspacing_edit') and self.kspacing_edit.text():
                params['kspacing'] = float(self.kspacing_edit.text())
            
            if hasattr(self, 'electron_maxstep_edit') and self.electron_maxstep_edit.text():
                params['electron_maxstep'] = int(self.electron_maxstep_edit.text())
            
            # Additional SCF parameters
            if hasattr(self, 'verbosity_combo'):
                params['verbosity'] = self.verbosity_combo.currentText()
            
            if hasattr(self, 'restart_mode_combo'):
                params['restart_mode'] = self.restart_mode_combo.currentText()
            
            if hasattr(self, 'disk_io_combo'):
                params['disk_io'] = self.disk_io_combo.currentText()
            
            if hasattr(self, 'mixing_mode_combo'):
                params['mixing_mode'] = self.mixing_mode_combo.currentText()
            
            # Forces and stress
            if hasattr(self, 'calc_forces_check'):
                params['tprnfor'] = self.calc_forces_check.isChecked()
            
            if hasattr(self, 'calc_stress_check'):
                params['tstress'] = self.calc_stress_check.isChecked()
            
            # Relaxation parameters
            if hasattr(self, 'ion_dynamics_combo'):
                ion_dyn = self.ion_dynamics_combo.currentText()
                if ion_dyn != 'none':  # Only add if not none
                    params['ion_dynamics'] = ion_dyn
            
            if hasattr(self, 'cell_dynamics_combo'):
                cell_dyn = self.cell_dynamics_combo.currentText()
                if cell_dyn != 'none':  # Only add if not none
                    params['cell_dynamics'] = cell_dyn
            
            if hasattr(self, 'cell_dofree_combo'):
                cell_free = self.cell_dofree_combo.currentText()
                params['cell_dofree'] = cell_free
            
            if hasattr(self, 'press_edit') and self.press_edit.text():
                press = float(self.press_edit.text())
                if press != 0.0:  # Only add if not zero
                    params['press'] = press
            
            if hasattr(self, 'press_conv_test_edit') and self.press_conv_test_edit.text():
                params['press_conv_test'] = float(self.press_conv_test_edit.text())
                
        except Exception:
            pass
        
        return params
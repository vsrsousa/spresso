"""Single, clean CalculationWindow implementation.

This file was overwritten to remove duplicated and corrupted sections and
to ensure a single `CalculationWindow` class provides the requested
tabs. It intentionally keeps logic defensive to work in minimal test
environments.
"""
import json
import os
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QFormLayout,
    QComboBox, QLineEdit, QTextEdit, QPushButton, QCheckBox,
    QTabWidget, QApplication, QMessageBox, QGroupBox, QScrollArea
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


try:
    import importlib.util
    spec = importlib.util.spec_from_file_location('pseudo_orbitals', os.path.join(os.path.dirname(__file__), '..', '..', 'xespresso', 'tools', 'pseudo_orbitals.py'))
    pseudo_orbitals = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pseudo_orbitals)
    parse_pseudopotential_orbitals = pseudo_orbitals.parse_pseudopotential_orbitals
    PSEUDO_ORBITALS_AVAILABLE = True
except Exception:
    PSEUDO_ORBITALS_AVAILABLE = False


class CalculationWindow(QWidget):
    """Modeless calculation configuration window with required tabs."""

    def __init__(self, calc_name: str, session_state=None, parent=None):
        super().__init__(parent)
        self.calc_name = calc_name
        self.session_state = session_state or {}
        self.setWindowTitle(f"Calculation: {calc_name}")
        self.resize(900, 700)

        # Initialize dictionaries for dynamic inputs
        self.magnetic_edits = {}
        self.magnetic_checkboxes = {}  # Track which elements have magnetism enabled

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

        # Listen for structure changes to update pseudopotentials and magnetic inputs
        if self.session_state and hasattr(self.session_state, 'add_listener'):
            self.session_state.add_listener(self._update_pseudo_inputs_for_structure)
            self.session_state.add_listener(self._update_magnetic_inputs_for_structure)

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
        layout.setSpacing(5)
        layout.setContentsMargins(5, 5, 5, 5)

        # Pseudopotentials Configuration
        self.pseudo_group = QGroupBox("📁 Pseudopotentials")
        self.pseudo_group.setCheckable(True)
        self.pseudo_group.setChecked(True)  # Enabled by default
        pseudo_layout = QVBoxLayout(self.pseudo_group)

        # Configuration selector
        config_layout = QHBoxLayout()
        config_layout.addWidget(QLabel("Config:"))
        self.pseudo_config_combo = QComboBox()
        self.pseudo_config_combo.currentTextChanged.connect(self._on_pseudo_config_changed)
        config_layout.addWidget(self.pseudo_config_combo)

        refresh_btn = QPushButton("↻")
        refresh_btn.setMaximumWidth(25)
        refresh_btn.clicked.connect(self._load_pseudo_configs)
        config_layout.addWidget(refresh_btn)
        pseudo_layout.addLayout(config_layout)

        # Status label
        self.pseudo_status_label = QLabel("")
        self.pseudo_status_label.setWordWrap(True)
        pseudo_layout.addWidget(self.pseudo_status_label)

        # Per-element pseudopotential inputs
        self.pseudo_container = QWidget()
        self.pseudo_form_layout = QFormLayout(self.pseudo_container)
        pseudo_layout.addWidget(self.pseudo_container)

        # Show/hide pseudo controls when the group checkbox is toggled
        self.pseudo_group.toggled.connect(lambda checked: self._on_pseudo_group_toggled(checked))

        layout.addWidget(self.pseudo_group)
        layout.addStretch()

        self.tabs.addTab(w, 'Pseudopotentials')

        # Initialize
        self.pseudo_edits = {}
        if PSEUDO_SELECTOR_AVAILABLE:
            self._load_pseudo_configs()
        self._update_pseudo_inputs_for_structure()

    def _load_pseudo_configs(self):
        """Load pseudopotential configurations."""
        if not PSEUDO_SELECTOR_AVAILABLE:
            return

        self.pseudo_config_combo.blockSignals(True)
        self.pseudo_config_combo.clear()
        self.pseudo_config_combo.addItem("Manual", None)

        try:
            from xespresso.pseudopotentials import PseudopotentialsManager
            configs = PseudopotentialsManager.list_configs()
            for config_name in configs:
                self.pseudo_config_combo.addItem(config_name, config_name)

            # Check for default
            if PseudopotentialsManager.has_default_config():
                self.pseudo_config_combo.insertItem(1, "Default", "default")
                self.pseudo_config_combo.setCurrentIndex(1)

        except Exception:
            pass

        self.pseudo_config_combo.blockSignals(False)
        self._on_pseudo_config_changed(self.pseudo_config_combo.currentText())

    def _on_pseudo_config_changed(self, config_name):
        """Handle pseudopotential config change."""
        if not PSEUDO_SELECTOR_AVAILABLE:
            return

        try:
            from xespresso.pseudopotentials import PseudopotentialsManager

            if config_name == "Manual":
                # Manual mode - clear any loaded pseudopotentials
                self.pseudo_status_label.setText("Manual configuration - enter pseudopotentials manually")
                self.pseudo_status_label.setStyleSheet("")
            elif config_name == "Default":
                # Load default config
                config = PseudopotentialsManager.load_config("default")
                if config:
                    self._load_pseudopotentials_from_config(config)
                    self.pseudo_status_label.setText(f"Loaded default config: {config.description or config.name}")
                    self.pseudo_status_label.setStyleSheet("color: green;")
                else:
                    self.pseudo_status_label.setText("❌ Default config not found")
                    self.pseudo_status_label.setStyleSheet("color: red;")
            else:
                # Load specific config
                config = PseudopotentialsManager.load_config(config_name)
                if config:
                    self._load_pseudopotentials_from_config(config)
                    desc = config.description or f"{config.library} {config.version}" if config.library else config.name
                    self.pseudo_status_label.setText(f"Loaded: {desc}")
                    self.pseudo_status_label.setStyleSheet("color: green;")
                else:
                    self.pseudo_status_label.setText(f"❌ Config '{config_name}' not found")
                    self.pseudo_status_label.setStyleSheet("color: red;")

        except Exception as e:
            self.pseudo_status_label.setText(f"❌ Error loading config: {e}")
            self.pseudo_status_label.setStyleSheet("color: red;")

        # Update inputs regardless
        self._update_pseudo_inputs_for_structure()

    def _load_pseudopotentials_from_config(self, config):
        """Load pseudopotentials from config into session state."""
        if not self.session_state:
            return

        pseudopotentials = {}
        for element, pseudo in config.pseudopotentials.items():
            pseudopotentials[element] = pseudo.filename

        self.session_state['pseudopotentials'] = pseudopotentials

    def _on_pseudo_group_toggled(self, checked):
        """Show/hide pseudopotential controls when the group is toggled."""
        self.pseudo_container.setVisible(checked)
        if checked:
            self._update_pseudo_inputs_for_structure()

    def _update_pseudo_inputs_for_structure(self):
        """Update pseudopotential inputs based on current structure."""
        # Clear existing inputs
        while self.pseudo_form_layout.count():
            item = self.pseudo_form_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        self.pseudo_edits = {}

        # Get current structure
        atoms = None
        if self.session_state:
            atoms = self.session_state.get('current_structure')

        if atoms:
            elements = sorted(set(atoms.get_chemical_symbols()))

            # Get existing pseudopotentials from session state
            existing_pseudos = {}
            if self.session_state:
                existing_pseudos = self.session_state.get('pseudopotentials', {})

            # Create inputs for each element
            for element in elements:
                edit = QLineEdit()
                edit.setPlaceholderText(f"{element}.UPF")
                edit.textChanged.connect(self._on_pseudopotentials_changed)

                # Pre-fill with existing pseudopotential if available
                if element in existing_pseudos:
                    edit.setText(existing_pseudos[element])

                self.pseudo_edits[element] = edit
                self.pseudo_form_layout.addRow(f"{element}:", edit)
        else:
            # No structure loaded
            info_label = QLabel("Load a structure to configure pseudopotentials")
            info_label.setStyleSheet("font-style: italic; color: gray;")
            self.pseudo_form_layout.addRow(info_label)

    def _on_pseudopotentials_changed(self):
        """Handle pseudopotential input changes."""
        pseudopotentials = {}
        for element, edit in self.pseudo_edits.items():
            pseudo = edit.text().strip()
            if pseudo:
                pseudopotentials[element] = pseudo

        # Update session state
        if self.session_state:
            self.session_state['pseudopotentials'] = pseudopotentials

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
        layout = QVBoxLayout(w)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # Magnetic Configuration (optional)
        self.magnetic_group = QGroupBox("🧲 Magnetic Configuration (Optional)")
        self.magnetic_group.setCheckable(True)
        self.magnetic_group.setChecked(False)
        magnetic_layout = QVBoxLayout(self.magnetic_group)
        
        magnetic_info = QLabel("""
<p><b>Configure magnetic properties for spin-polarized calculations.</b></p>
<p>Select which elements should be magnetic and specify their magnetic moments.</p>
<p>• Single value: all atoms of that element get the same magnetization</p>
<p>• Multiple values (comma-separated): creates different magnetic species (e.g., 1.0,-1.0 for AFM)</p>
<p>• Auto-expand cell: automatically creates supercell for complex configurations</p>
""")
        magnetic_info.setTextFormat(Qt.RichText)
        magnetic_info.setWordWrap(True)
        magnetic_layout.addWidget(magnetic_info)
        
        # Preset selector (wrapped in a widget so we can show/hide it with the group)
        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Preset:"))
        self.magnetic_preset_combo = QComboBox()
        self.magnetic_preset_combo.addItems(["Custom", "Ferromagnetic", "Antiferromagnetic"])
        self.magnetic_preset_combo.currentTextChanged.connect(self._on_magnetic_preset_changed)
        preset_layout.addWidget(self.magnetic_preset_combo)
        self.magnetic_preset_widget = QWidget()
        self.magnetic_preset_widget.setLayout(preset_layout)
        self.magnetic_preset_widget.setVisible(False)
        magnetic_layout.addWidget(self.magnetic_preset_widget)
        
        # Container for per-element magnetic inputs (hidden until enabled)
        self.magnetic_container = QWidget()
        self.magnetic_container_layout = QFormLayout(self.magnetic_container)
        self.magnetic_container.setVisible(False)
        magnetic_layout.addWidget(self.magnetic_container)
        
        # Expand cell option
        self.expand_cell_check = QCheckBox("Auto-expand cell if more magnetic moments specified than atoms exist")
        self.expand_cell_check.setToolTip("Automatically create supercell to accommodate complex magnetic configurations")
        self.expand_cell_check.setVisible(False)
        magnetic_layout.addWidget(self.expand_cell_check)

        # Show/hide magnetic controls when the group checkbox is toggled
        self.magnetic_group.toggled.connect(lambda checked: self._on_magnetic_group_toggled(checked))
        
        layout.addWidget(self.magnetic_group)
        layout.addStretch()
        
        self.tabs.addTab(w, 'Magnetism')

    def _update_magnetic_inputs(self, elements):
        """Update magnetic input fields for structure elements."""
        # Clear existing inputs
        while self.magnetic_container_layout.count():
            item = self.magnetic_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.magnetic_edits = {}
        self.magnetic_checkboxes = {}
        
        for element in sorted(elements):
            # Create a container widget for each element with checkbox and magnetic moments input
            container = QWidget()
            hlayout = QHBoxLayout(container)
            hlayout.setContentsMargins(0, 0, 0, 0)
            
            # Checkbox to enable/disable magnetism for this element
            checkbox = QCheckBox()
            checkbox.setChecked(False)  # Default to unchecked
            checkbox.setToolTip(f"Enable magnetic configuration for {element}")
            checkbox.stateChanged.connect(lambda state, elem=element: self._on_magnetic_element_toggled(elem, state))
            hlayout.addWidget(checkbox)
            
            # Magnetic moments input (comma-separated values)
            edit = QLineEdit()
            edit.setPlaceholderText("e.g., 1.0 or 1.0,-1.0 or 1.0,1.0,-1.0,-1.0")
            edit.setToolTip(f"Magnetic moments for {element} (comma-separated)")
            edit.setEnabled(False)  # Disabled by default
            hlayout.addWidget(edit)
            
            self.magnetic_checkboxes[element] = checkbox
            self.magnetic_edits[element] = edit
            self.magnetic_container_layout.addRow(f"{element}:", container)

    def _on_magnetic_preset_changed(self, preset):
        """Handle magnetic preset selection."""
        if preset == "Custom":
            return
        
        # Clear all current settings
        for element in self.magnetic_checkboxes:
            self.magnetic_checkboxes[element].setChecked(False)
            self.magnetic_edits[element].setText("")
        
        # Apply preset
        for element in self.magnetic_checkboxes:
            if element in ['Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Ti', 'Gd', 'Nd', 'Ce']:
                self.magnetic_checkboxes[element].setChecked(True)
                if preset == "Ferromagnetic":
                    mag_val = {'Fe': 2.2, 'Co': 1.7, 'Ni': 0.6, 'Mn': 5.0, 'Cr': 3.0, 'V': 2.0, 'Ti': 1.0, 'Gd': 7.0, 'Nd': 3.0, 'Ce': 5.0}.get(element, 1.0)
                    self.magnetic_edits[element].setText(f"{mag_val:.2f}")
                elif preset == "Antiferromagnetic":
                    # For AFM, we need two opposite values
                    mag_val = {'Fe': 2.2, 'Co': 1.7, 'Ni': 0.6, 'Mn': 5.0, 'Cr': 3.0, 'V': 2.0, 'Ti': 1.0, 'Gd': 7.0, 'Nd': 3.0, 'Ce': 5.0}.get(element, 1.0)
                    self.magnetic_edits[element].setText(f"{mag_val:.2f},{-mag_val:.2f}")

    def _on_magnetic_element_toggled(self, element, state):
        """Enable/disable magnetic moment input when element checkbox is toggled."""
        if element in self.magnetic_edits:
            self.magnetic_edits[element].setEnabled(state == 2)  # Qt.CheckState.Checked
    
    def _on_magnetic_group_toggled(self, checked):
        """Show or hide magnetic controls when the magnetic group is toggled."""
        # Show preset selector and per-element inputs directly under the group
        try:
            self.magnetic_preset_widget.setVisible(checked)
        except Exception:
            pass
        try:
            self.magnetic_container.setVisible(checked)
        except Exception:
            pass
        try:
            self.expand_cell_check.setVisible(checked)
        except Exception:
            pass
        
        # Update magnetic inputs when enabled
        if checked:
            self._update_magnetic_inputs_for_structure()

    def _update_magnetic_inputs_for_structure(self):
        """Update magnetic inputs based on current structure."""
        try:
            atoms = None
            if self.session_state:
                atoms = self.session_state.get('current_structure')
            
            if atoms:
                elements = set(atoms.get_chemical_symbols())
                self._update_magnetic_inputs(elements)
            else:
                # No structure loaded
                self._update_magnetic_inputs([])
        except Exception:
            pass

    def _on_magnetism_toggled(self, state):
        """Show/hide magnetic configuration options when magnetism is enabled/disabled."""
        enabled = state == 2  # Qt.CheckState.Checked
        self._show_magnetic_config(enabled)
        if enabled:
            self._update_atom_magnetizations()

    def _build_hubbard_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(5)
        layout.setContentsMargins(5, 5, 5, 5)

        # Hubbard Configuration
        self.hubbard_group = QGroupBox("🧲 Hubbard Parameters (DFT+U)")
        self.hubbard_group.setCheckable(True)
        self.hubbard_group.setChecked(False)
        hubbard_layout = QVBoxLayout(self.hubbard_group)

        # Format selection and projector
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("Format:"))
        self.hubbard_format_combo = QComboBox()
        self.hubbard_format_combo.addItems(["Auto", "Old (QE < 7.0)", "New (QE >= 7.0)"])
        self.hubbard_format_combo.setToolTip("QE version determines format: Auto detects based on orbital specifications, Old uses SYSTEM namelist, New uses HUBBARD card")
        format_layout.addWidget(self.hubbard_format_combo)

        format_layout.addWidget(QLabel("Projector:"))
        self.hubbard_projector_combo = QComboBox()
        self.hubbard_projector_combo.addItems(["atomic", "ortho-atomic", "norm-atomic", "wf", "pseudo"])
        self.hubbard_projector_combo.setCurrentText("ortho-atomic")
        self.hubbard_projector_combo.setToolTip("Projector type for new format Hubbard calculations (ortho-atomic recommended by QE)")
        format_layout.addWidget(self.hubbard_projector_combo)
        hubbard_layout.addLayout(format_layout)

        # Status label
        self.hubbard_status_label = QLabel("")
        self.hubbard_status_label.setWordWrap(True)
        hubbard_layout.addWidget(self.hubbard_status_label)

        # U Parameters section
        u_group = QGroupBox("U Parameters (On-site)")
        u_layout = QVBoxLayout(u_group)

        u_info = QLabel("Configure Hubbard U values for each element/orbital combination.\n"
                       "• Auto format: Automatically chooses based on orbital specifications\n"
                       "• Old format: U values per element (e.g., Fe: 4.3)\n"
                       "• New format: U values per element-orbital (e.g., Fe-3d: 4.3)")
        u_info.setWordWrap(True)
        u_layout.addWidget(u_info)

        self.hubbard_u_container = QWidget()
        self.hubbard_u_form_layout = QFormLayout(self.hubbard_u_container)
        u_layout.addWidget(self.hubbard_u_container)

        # Add U parameter button
        add_u_layout = QHBoxLayout()
        self.hubbard_u_element_combo = QComboBox()
        self.hubbard_u_element_combo.setEditable(True)
        self.hubbard_u_element_combo.setPlaceholderText("Element or Element-Orbital")
        add_u_layout.addWidget(self.hubbard_u_element_combo)

        self.hubbard_u_value_edit = QLineEdit()
        self.hubbard_u_value_edit.setPlaceholderText("U value (eV)")
        add_u_layout.addWidget(self.hubbard_u_value_edit)

        add_u_btn = QPushButton("Add U")
        add_u_btn.clicked.connect(self._add_hubbard_u_parameter)
        add_u_layout.addWidget(add_u_btn)
        u_layout.addLayout(add_u_layout)

        hubbard_layout.addWidget(u_group)

        # V Parameters section
        v_group = QGroupBox("V Parameters (Inter-site)")
        v_layout = QVBoxLayout(v_group)

        v_info = QLabel("Configure Hubbard V values for inter-site interactions.\n"
                       "• Old format: V(na,nb,k) parameters\n"
                       "• New format: V between specific element-orbital pairs")
        v_info.setWordWrap(True)
        v_layout.addWidget(v_info)

        self.hubbard_v_container = QWidget()
        self.hubbard_v_form_layout = QFormLayout(self.hubbard_v_container)
        v_layout.addWidget(self.hubbard_v_container)

        # Add V parameter button
        add_v_layout = QHBoxLayout()
        self.hubbard_v_spec1_combo = QComboBox()
        self.hubbard_v_spec1_combo.setEditable(True)
        self.hubbard_v_spec1_combo.setPlaceholderText("Species1-Orbital")
        add_v_layout.addWidget(self.hubbard_v_spec1_combo)

        self.hubbard_v_spec2_combo = QComboBox()
        self.hubbard_v_spec2_combo.setEditable(True)
        self.hubbard_v_spec2_combo.setPlaceholderText("Species2-Orbital")
        add_v_layout.addWidget(self.hubbard_v_spec2_combo)

        self.hubbard_v_value_edit = QLineEdit()
        self.hubbard_v_value_edit.setPlaceholderText("V value (eV)")
        add_v_layout.addWidget(self.hubbard_v_value_edit)

        add_v_btn = QPushButton("Add V")
        add_v_btn.clicked.connect(self._add_hubbard_v_parameter)
        add_v_layout.addWidget(add_v_btn)
        v_layout.addLayout(add_v_layout)

        hubbard_layout.addWidget(v_group)

        # Advanced parameters section (collapsible)
        self.hubbard_advanced_group = QGroupBox("Advanced Parameters (J, α, β)")
        self.hubbard_advanced_group.setCheckable(True)
        self.hubbard_advanced_group.setChecked(False)
        advanced_layout = QVBoxLayout(self.hubbard_advanced_group)

        self.hubbard_advanced_container = QWidget()
        self.hubbard_advanced_form_layout = QFormLayout(self.hubbard_advanced_container)
        advanced_layout.addWidget(self.hubbard_advanced_container)

        # Add advanced parameter controls
        add_advanced_layout = QHBoxLayout()
        self.hubbard_advanced_type_combo = QComboBox()
        self.hubbard_advanced_type_combo.addItems(["J", "α (alpha)", "β (beta)"])
        add_advanced_layout.addWidget(self.hubbard_advanced_type_combo)

        self.hubbard_advanced_element_combo = QComboBox()
        self.hubbard_advanced_element_combo.setEditable(True)
        self.hubbard_advanced_element_combo.setPlaceholderText("Element")
        add_advanced_layout.addWidget(self.hubbard_advanced_element_combo)

        self.hubbard_advanced_value_edit = QLineEdit()
        self.hubbard_advanced_value_edit.setPlaceholderText("Value (eV)")
        add_advanced_layout.addWidget(self.hubbard_advanced_value_edit)

        add_advanced_btn = QPushButton("Add")
        add_advanced_btn.clicked.connect(self._add_hubbard_advanced_parameter)
        add_advanced_layout.addWidget(add_advanced_btn)
        advanced_layout.addLayout(add_advanced_layout)

        hubbard_layout.addWidget(self.hubbard_advanced_group)

        # Show/hide hubbard controls when the group checkbox is toggled
        self.hubbard_group.toggled.connect(lambda checked: self._on_hubbard_group_toggled(checked))

        layout.addWidget(self.hubbard_group)
        layout.addStretch()

        self.tabs.addTab(w, 'Hubbard')

        # Initialize
        self.hubbard_u_edits = {}
        self.hubbard_v_edits = {}
        self.hubbard_advanced_edits = {}
        self._update_hubbard_inputs_for_structure()

    def _on_hubbard_group_toggled(self, checked):
        """Show/hide Hubbard controls when the group is toggled."""
        # Update inputs when enabled
        if checked:
            self._update_hubbard_inputs_for_structure()

    def _update_hubbard_inputs_for_structure(self):
        """Update Hubbard inputs based on current structure."""
        try:
            atoms = None
            if self.session_state:
                atoms = self.session_state.get('current_structure')

            if atoms:
                elements = sorted(set(atoms.get_chemical_symbols()))

                # Update element combos
                self.hubbard_u_element_combo.clear()
                self.hubbard_v_spec1_combo.clear()
                self.hubbard_v_spec2_combo.clear()
                self.hubbard_advanced_element_combo.clear()

                # Get available orbitals based on selected pseudopotentials
                available_orbitals = self._get_available_hubbard_orbitals()

                # Add elements and element-orbital combinations
                for element in elements:
                    self.hubbard_u_element_combo.addItem(element)
                    self.hubbard_v_spec1_combo.addItem(element)
                    self.hubbard_v_spec2_combo.addItem(element)
                    self.hubbard_advanced_element_combo.addItem(element)

                    # Add orbital combinations based on pseudopotential analysis
                    element_orbitals = available_orbitals.get(element, ['3d'])  # Default fallback
                    for orbital in element_orbitals:
                        orbital_combo = f"{element}-{orbital}"
                        self.hubbard_u_element_combo.addItem(orbital_combo)
                        self.hubbard_v_spec1_combo.addItem(orbital_combo)
                        self.hubbard_v_spec2_combo.addItem(orbital_combo)

                # Get existing Hubbard parameters from session state
                existing_hubbard = {}
                if self.session_state:
                    existing_hubbard = self.session_state.get('hubbard', {})

                # Pre-populate U parameters
                self._clear_hubbard_u_inputs()
                if 'u' in existing_hubbard:
                    for key, value in existing_hubbard['u'].items():
                        self._add_hubbard_u_input(key, value)

                # Pre-populate V parameters
                self._clear_hubbard_v_inputs()
                if 'v' in existing_hubbard:
                    for v_param in existing_hubbard['v']:
                        if isinstance(v_param, dict):
                            spec1 = f"{v_param.get('species1', '')}-{v_param.get('orbital1', '')}"
                            spec2 = f"{v_param.get('species2', '')}-{v_param.get('orbital2', '')}"
                            value = v_param.get('value', 0.0)
                            self._add_hubbard_v_input(spec1, spec2, value)

                # Update status
                self.hubbard_status_label.setText(f"Configured for elements: {', '.join(elements)}")
                self.hubbard_status_label.setStyleSheet("color: green;")
            else:
                # No structure loaded
                self._clear_hubbard_u_inputs()
                self._clear_hubbard_v_inputs()
                self._clear_hubbard_advanced_inputs()
                self.hubbard_status_label.setText("Load a structure to configure Hubbard parameters")
                self.hubbard_status_label.setStyleSheet("")

        except Exception as e:
            self.hubbard_status_label.setText(f"Error updating Hubbard inputs: {e}")
            self.hubbard_status_label.setStyleSheet("color: red;")

    def _clear_hubbard_u_inputs(self):
        """Clear all U parameter inputs."""
        while self.hubbard_u_form_layout.count():
            item = self.hubbard_u_form_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.hubbard_u_edits = {}

    def _clear_hubbard_v_inputs(self):
        """Clear all V parameter inputs."""
        while self.hubbard_v_form_layout.count():
            item = self.hubbard_v_form_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.hubbard_v_edits = {}

    def _clear_hubbard_advanced_inputs(self):
        """Clear all advanced parameter inputs."""
        while self.hubbard_advanced_form_layout.count():
            item = self.hubbard_advanced_form_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.hubbard_advanced_edits = {}

    def _add_hubbard_u_parameter(self):
        """Add a U parameter from the input fields."""
        element_orbital = self.hubbard_u_element_combo.currentText().strip()
        value_text = self.hubbard_u_value_edit.text().strip()

        if not element_orbital or not value_text:
            return

        try:
            value = float(value_text)
            self._add_hubbard_u_input(element_orbital, value)
            self.hubbard_u_value_edit.clear()
            self._update_hubbard_session_state()
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "U value must be a number")

    def _add_hubbard_u_input(self, element_orbital, value):
        """Add a U parameter input field."""
        edit = QLineEdit(str(value))
        edit.setMaximumWidth(80)
        edit.textChanged.connect(self._update_hubbard_session_state)

        remove_btn = QPushButton("×")
        remove_btn.setMaximumWidth(25)
        remove_btn.clicked.connect(lambda: self._remove_hubbard_u_parameter(element_orbital))

        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit)
        layout.addWidget(remove_btn)

        self.hubbard_u_form_layout.addRow(f"{element_orbital}:", container)
        self.hubbard_u_edits[element_orbital] = edit

    def _remove_hubbard_u_parameter(self, element_orbital):
        """Remove a U parameter."""
        if element_orbital in self.hubbard_u_edits:
            # Find and remove the row
            for i in range(self.hubbard_u_form_layout.rowCount()):
                label_item = self.hubbard_u_form_layout.itemAt(i, QFormLayout.LabelRole)
                if label_item and label_item.widget():
                    label_text = label_item.widget().text()
                    if label_text.startswith(f"{element_orbital}:"):
                        # Remove the row
                        self.hubbard_u_form_layout.removeRow(i)
                        break

            del self.hubbard_u_edits[element_orbital]
            self._update_hubbard_session_state()

    def _add_hubbard_v_parameter(self):
        """Add a V parameter from the input fields."""
        spec1 = self.hubbard_v_spec1_combo.currentText().strip()
        spec2 = self.hubbard_v_spec2_combo.currentText().strip()
        value_text = self.hubbard_v_value_edit.text().strip()

        if not spec1 or not spec2 or not value_text:
            return

        try:
            value = float(value_text)
            self._add_hubbard_v_input(spec1, spec2, value)
            self.hubbard_v_value_edit.clear()
            self._update_hubbard_session_state()
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "V value must be a number")

    def _add_hubbard_v_input(self, spec1, spec2, value):
        """Add a V parameter input field."""
        edit = QLineEdit(str(value))
        edit.setMaximumWidth(80)
        edit.textChanged.connect(self._update_hubbard_session_state)

        remove_btn = QPushButton("×")
        remove_btn.setMaximumWidth(25)
        remove_btn.clicked.connect(lambda: self._remove_hubbard_v_parameter(spec1, spec2))

        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit)
        layout.addWidget(remove_btn)

        label = f"{spec1} ↔ {spec2}:"
        self.hubbard_v_form_layout.addRow(label, container)
        self.hubbard_v_edits[(spec1, spec2)] = edit

    def _remove_hubbard_v_parameter(self, spec1, spec2):
        """Remove a V parameter."""
        key = (spec1, spec2)
        if key in self.hubbard_v_edits:
            # Find and remove the row
            for i in range(self.hubbard_v_form_layout.rowCount()):
                label_item = self.hubbard_v_form_layout.itemAt(i, QFormLayout.LabelRole)
                if label_item and label_item.widget():
                    label_text = label_item.widget().text()
                    if f"{spec1} ↔ {spec2}:" in label_text:
                        # Remove the row
                        self.hubbard_v_form_layout.removeRow(i)
                        break

            del self.hubbard_v_edits[key]
            self._update_hubbard_session_state()

    def _add_hubbard_advanced_parameter(self):
        """Add an advanced parameter (J, alpha, beta)."""
        param_type = self.hubbard_advanced_type_combo.currentText()
        element = self.hubbard_advanced_element_combo.currentText().strip()
        value_text = self.hubbard_advanced_value_edit.text().strip()

        if not element or not value_text:
            return

        try:
            value = float(value_text)

            # Map display names to parameter keys
            type_map = {"J": "j", "α (alpha)": "alpha", "β (beta)": "beta"}
            param_key = type_map.get(param_type, param_type.lower())

            self._add_hubbard_advanced_input(param_key, element, value)
            self.hubbard_advanced_value_edit.clear()
            self._update_hubbard_session_state()
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Parameter value must be a number")

    def _add_hubbard_advanced_input(self, param_type, element, value):
        """Add an advanced parameter input field."""
        edit = QLineEdit(str(value))
        edit.setMaximumWidth(80)
        edit.textChanged.connect(self._update_hubbard_session_state)

        remove_btn = QPushButton("×")
        remove_btn.setMaximumWidth(25)
        remove_btn.clicked.connect(lambda: self._remove_hubbard_advanced_parameter(param_type, element))

        container = QWidget()
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit)
        layout.addWidget(remove_btn)

        # Display name mapping
        display_map = {"j": "J", "alpha": "α", "beta": "β"}
        display_type = display_map.get(param_type, param_type.upper())

        label = f"{display_type} ({element}):"
        self.hubbard_advanced_form_layout.addRow(label, container)
        self.hubbard_advanced_edits[(param_type, element)] = edit

    def _remove_hubbard_advanced_parameter(self, param_type, element):
        """Remove an advanced parameter."""
        key = (param_type, element)
        if key in self.hubbard_advanced_edits:
            # Find and remove the row
            display_map = {"j": "J", "alpha": "α", "beta": "β"}
            display_type = display_map.get(param_type, param_type.upper())

            for i in range(self.hubbard_advanced_form_layout.rowCount()):
                label_item = self.hubbard_advanced_form_layout.itemAt(i, QFormLayout.LabelRole)
                if label_item and label_item.widget():
                    label_text = label_item.widget().text()
                    if f"{display_type} ({element}):" in label_text:
                        # Remove the row
                        self.hubbard_advanced_form_layout.removeRow(i)
                        break

            del self.hubbard_advanced_edits[key]
            self._update_hubbard_session_state()

    def _update_hubbard_session_state(self):
        """Update session state with current Hubbard parameters."""
        if not self.session_state:
            return

        hubbard_config = {}

        # Format selection
        format_map = {"Auto": None, "Old (QE < 7.0)": "old", "New (QE >= 7.0)": "new"}
        selected_format = self.hubbard_format_combo.currentText()
        if format_map[selected_format] == "old":
            hubbard_config["use_new_format"] = False
        elif format_map[selected_format] == "new":
            hubbard_config["use_new_format"] = True

        # Projector
        hubbard_config["projector"] = self.hubbard_projector_combo.currentText()

        # U parameters
        if self.hubbard_u_edits:
            hubbard_config["u"] = {}
            for element_orbital, edit in self.hubbard_u_edits.items():
                try:
                    value = float(edit.text().strip())
                    hubbard_config["u"][element_orbital] = value
                except ValueError:
                    pass

        # V parameters
        if self.hubbard_v_edits:
            hubbard_config["v"] = []
            for (spec1, spec2), edit in self.hubbard_v_edits.items():
                try:
                    value = float(edit.text().strip())
                    # Parse species and orbitals
                    parts1 = spec1.split('-', 1)
                    parts2 = spec2.split('-', 1)

                    v_param = {
                        "species1": parts1[0],
                        "orbital1": parts1[1] if len(parts1) > 1 else "",
                        "species2": parts2[0],
                        "orbital2": parts2[1] if len(parts2) > 1 else "",
                        "i": 1,
                        "j": 1,
                        "value": value
                    }
                    hubbard_config["v"].append(v_param)
                except (ValueError, IndexError):
                    pass

        # Advanced parameters
        for (param_type, element), edit in self.hubbard_advanced_edits.items():
            try:
                value = float(edit.text().strip())
                if param_type not in hubbard_config:
                    hubbard_config[param_type] = {}
                hubbard_config[param_type][element] = value
            except ValueError:
                pass

        self.session_state['hubbard'] = hubbard_config

    def _get_available_hubbard_orbitals(self):
        """Get available Hubbard orbitals based on selected pseudopotentials."""
        if not PSEUDO_ORBITALS_AVAILABLE or not self.session_state:
            return {}

        pseudopotentials = self.session_state.get('pseudopotentials', {})
        if not pseudopotentials:
            return {}

        available_orbitals = {}

        # Try to find pseudopotential files in common locations
        # This is more generic than hardcoded paths
        for element, pseudo_filename in pseudopotentials.items():
            if not pseudo_filename:
                continue

            # Try to find the pseudopotential file
            pseudo_path = None
            
            # Check if it's an absolute path
            if os.path.isfile(pseudo_filename):
                pseudo_path = pseudo_filename
            else:
                # Try relative to current working directory
                cwd_path = os.path.join(os.getcwd(), pseudo_filename)
                if os.path.isfile(cwd_path):
                    pseudo_path = cwd_path
                else:
                    # Try in common pseudopotential directories
                    # This could be extended with environment variables or config
                    common_dirs = [
                        'pseudo',
                        'pseudopotentials', 
                        'data/pseudo',
                        'examples/datas/pseudo'
                    ]
                    for common_dir in common_dirs:
                        candidate_path = os.path.join(os.getcwd(), common_dir, pseudo_filename)
                        if os.path.isfile(candidate_path):
                            pseudo_path = candidate_path
                            break

            if pseudo_path:
                try:
                    # Parse the pseudopotential file to extract orbitals
                    orbital_info = parse_pseudopotential_orbitals(pseudo_path)
                    orbitals = orbital_info.get('orbitals', [])

                    # Convert to the format expected by Hubbard (e.g., '3d', '4f')
                    # Filter out s and p orbitals as they're typically not used for Hubbard U
                    hubbard_orbitals = [orb for orb in orbitals if orb[-1] in ['d', 'f']]

                    if hubbard_orbitals:
                        available_orbitals[element] = hubbard_orbitals
                    else:
                        # Fallback to JSON file data
                        available_orbitals[element] = self._get_fallback_orbitals_for_element(element)

                except Exception as e:
                    # If parsing fails, use fallback from JSON
                    available_orbitals[element] = self._get_fallback_orbitals_for_element(element)
            else:
                # If pseudopotential file not found, use fallback from JSON
                available_orbitals[element] = self._get_fallback_orbitals_for_element(element)

        return available_orbitals

    def _get_fallback_orbitals_for_element(self, element):
        """Get fallback orbital suggestions from JSON file."""
        try:
            # Load the hubbard orbitals JSON file
            json_path = os.path.join(os.path.dirname(__file__), '..', '..', 'xespresso', 'data', 'hubbard_orbitals.json')
            with open(json_path, 'r') as f:
                orbital_data = json.load(f)
            
            # Filter to only d and f orbitals for Hubbard U
            all_orbitals = orbital_data.get("orbitals", {}).get(element, ['3d'])
            hubbard_orbitals = [orb for orb in all_orbitals if orb[-1] in ['d', 'f']]
            return hubbard_orbitals if hubbard_orbitals else ['3d']
        except Exception:
            # Ultimate fallback if JSON loading fails
            return ['3d']

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

    def get_magnetism_config(self):
        """Get the current magnetism configuration as a dict."""
        config = {}
        
        try:
            # Check if magnetism is enabled
            if not hasattr(self, 'magnetic_group') or not self.magnetic_group.isChecked():
                return config
            
            # Get atoms from session state
            atoms = None
            if self.session_state:
                atoms = self.session_state.get('current_structure')
            
            if not atoms:
                return config
            
            # Build magnetic configuration from element-based inputs
            magnetic_config = {}
            expand_cell = hasattr(self, 'expand_cell_check') and self.expand_cell_check.isChecked()
            
            for element, checkbox in self.magnetic_checkboxes.items():
                if checkbox.isChecked():
                    edit = self.magnetic_edits.get(element)
                    if edit and edit.text().strip():
                        try:
                            # Parse comma-separated values
                            values = [float(x.strip()) for x in edit.text().split(',') if x.strip()]
                            if values:
                                magnetic_config[element] = values
                        except ValueError:
                            # Skip invalid values
                            continue
            
            if magnetic_config:
                config['magnetic_config'] = magnetic_config
                config['expand_cell'] = expand_cell
                
                # Generate the actual magnetic moments using xespresso's setup_magnetic_config
                try:
                    from xespresso import setup_magnetic_config
                    mag_config = setup_magnetic_config(atoms, magnetic_config, expand_cell=expand_cell)
                    config.update(mag_config)
                except Exception:
                    # Fallback: create basic config
                    config['input_ntyp'] = {'starting_magnetization': {}}
                    species_counter = {}
                    for element, values in magnetic_config.items():
                        if not values:
                            continue
                        # Use first value for each element (simplified fallback)
                        mag = values[0]
                        if element not in species_counter:
                            species_counter[element] = 0
                        species_counter[element] += 1
                        species_name = f"{element}{species_counter[element]}" if species_counter[element] > 1 else element
                        config['input_ntyp']['starting_magnetization'][species_name] = mag
                
        except Exception:
            pass
        
        return config
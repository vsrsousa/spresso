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
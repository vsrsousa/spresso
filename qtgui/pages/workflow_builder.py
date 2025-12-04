"""
Workflow Builder Page for xespresso PySide6 GUI.

This page creates multi-step workflows using the GUIWorkflow class.
"""

import os

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QDoubleSpinBox,
    QSpinBox, QCheckBox, QTextEdit, QRadioButton
)
from PySide6.QtCore import Qt

try:
    from xespresso.machines.config.loader import (
        list_machines, load_machine,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False

try:
    from ase import Atoms
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

# Import the pseudopotentials selector widget
try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except ImportError:
    PSEUDO_SELECTOR_AVAILABLE = False

# Import GUIWorkflow for creating workflow objects
try:
    from gui.workflows.base import GUIWorkflow
    WORKFLOW_AVAILABLE = True
except ImportError:
    WORKFLOW_AVAILABLE = False


class WorkflowBuilderPage(QWidget):
    """Workflow builder page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self._loading = False  # Guard to prevent infinite loops
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the user interface."""
        main_layout = QVBoxLayout(self)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Header
        header_label = QLabel("<h2>🔄 Workflow Builder</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p>Build multi-step workflows using <b>workflow modules</b> that orchestrate calculations.</p>
<p><b>Modular Design:</b></p>
<ul>
<li>Workflow modules coordinate multiple calculations</li>
<li>Each calculation step uses calculation modules to prepare objects</li>
<li>Job submission executes the prepared workflow steps</li>
</ul>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        # Structure Status
        self.structure_status = QLabel("")
        self.structure_status.setWordWrap(True)
        scroll_layout.addWidget(self.structure_status)
        
        # Execution Environment
        env_group = QGroupBox("🖥️ Execution Environment")
        env_layout = QFormLayout(env_group)
        
        self.machine_combo = QComboBox()
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        env_layout.addRow("Machine:", self.machine_combo)
        
        self.machine_info_label = QLabel("")
        env_layout.addRow(self.machine_info_label)
        
        self.version_combo = QComboBox()
        self.version_combo.currentTextChanged.connect(self._on_version_changed)
        env_layout.addRow("QE Version:", self.version_combo)
        
        self.code_combo = QComboBox()
        env_layout.addRow("Code:", self.code_combo)
        
        scroll_layout.addWidget(env_group)
        
        # Workflow Type
        workflow_group = QGroupBox("🔧 Workflow Type")
        workflow_layout = QFormLayout(workflow_group)
        
        self.workflow_combo = QComboBox()
        self.workflow_combo.addItems([
            "Single SCF",
            "SCF + Relaxation",
            "SCF + Relax + SCF (on relaxed structure)",
            "Custom Multi-Step"
        ])
        self.workflow_combo.currentTextChanged.connect(self._on_workflow_changed)
        workflow_layout.addRow("Select Workflow:", self.workflow_combo)
        
        scroll_layout.addWidget(workflow_group)
        
        # Common Parameters
        params_group = QGroupBox("⚙️ Common Parameters")
        params_layout = QFormLayout(params_group)
        
        self.ecutwfc_spin = QDoubleSpinBox()
        self.ecutwfc_spin.setRange(10.0, 200.0)
        self.ecutwfc_spin.setValue(50.0)
        self.ecutwfc_spin.setSingleStep(5.0)
        self.ecutwfc_spin.setSuffix(" Ry")
        params_layout.addRow("Energy Cutoff:", self.ecutwfc_spin)
        
        self.ecutrho_spin = QDoubleSpinBox()
        self.ecutrho_spin.setRange(40.0, 1600.0)
        self.ecutrho_spin.setValue(400.0)
        self.ecutrho_spin.setSingleStep(20.0)
        self.ecutrho_spin.setSuffix(" Ry")
        params_layout.addRow("Charge Density Cutoff:", self.ecutrho_spin)
        
        self.occupations_combo = QComboBox()
        self.occupations_combo.addItems(["smearing", "fixed", "tetrahedra"])
        params_layout.addRow("Occupations:", self.occupations_combo)
        
        scroll_layout.addWidget(params_group)
        
        # K-points
        kpts_group = QGroupBox("🔷 K-points")
        kpts_layout = QVBoxLayout(kpts_group)
        
        kpts_mode_layout = QHBoxLayout()
        self.kspacing_radio = QRadioButton("K-spacing")
        self.kspacing_radio.setChecked(True)
        self.kspacing_radio.toggled.connect(self._on_kpts_mode_changed)
        kpts_mode_layout.addWidget(self.kspacing_radio)
        
        self.explicit_radio = QRadioButton("Explicit Grid")
        kpts_mode_layout.addWidget(self.explicit_radio)
        kpts_layout.addLayout(kpts_mode_layout)
        
        # K-spacing widget
        self.kspacing_widget = QWidget()
        kspacing_inner = QHBoxLayout(self.kspacing_widget)
        kspacing_inner.setContentsMargins(0, 0, 0, 0)
        kspacing_inner.addWidget(QLabel("K-spacing (Å⁻¹):"))
        self.kspacing_spin = QDoubleSpinBox()
        self.kspacing_spin.setRange(0.1, 1.0)
        self.kspacing_spin.setValue(0.3)
        self.kspacing_spin.setSingleStep(0.05)
        kspacing_inner.addWidget(self.kspacing_spin)
        kpts_layout.addWidget(self.kspacing_widget)
        
        # Explicit k-points widget
        self.explicit_widget = QWidget()
        explicit_inner = QHBoxLayout(self.explicit_widget)
        explicit_inner.setContentsMargins(0, 0, 0, 0)
        
        explicit_inner.addWidget(QLabel("k₁:"))
        self.k1_spin = QSpinBox()
        self.k1_spin.setRange(1, 20)
        self.k1_spin.setValue(4)
        explicit_inner.addWidget(self.k1_spin)
        
        explicit_inner.addWidget(QLabel("k₂:"))
        self.k2_spin = QSpinBox()
        self.k2_spin.setRange(1, 20)
        self.k2_spin.setValue(4)
        explicit_inner.addWidget(self.k2_spin)
        
        explicit_inner.addWidget(QLabel("k₃:"))
        self.k3_spin = QSpinBox()
        self.k3_spin.setRange(1, 20)
        self.k3_spin.setValue(4)
        explicit_inner.addWidget(self.k3_spin)
        
        self.explicit_widget.setVisible(False)
        kpts_layout.addWidget(self.explicit_widget)
        
        scroll_layout.addWidget(kpts_group)
        
        # Pseudopotentials
        pseudo_group = QGroupBox("🔬 Pseudopotentials")
        pseudo_layout = QVBoxLayout(pseudo_group)
        
        if PSEUDO_SELECTOR_AVAILABLE:
            # Use the advanced selector widget
            self.pseudo_selector = PseudopotentialsSelectorWidget(self.session_state)
            pseudo_layout.addWidget(self.pseudo_selector)
        else:
            # Fallback to simple manual inputs
            self.pseudo_info_label = QLabel("Load a structure to configure pseudopotentials")
            pseudo_layout.addWidget(self.pseudo_info_label)
            
            self.pseudo_container = QWidget()
            self.pseudo_container_layout = QFormLayout(self.pseudo_container)
            pseudo_layout.addWidget(self.pseudo_container)
        
        scroll_layout.addWidget(pseudo_group)
        
        # Magnetic Configuration (optional)
        self.magnetic_group = QGroupBox("🧲 Magnetic Configuration (Optional)")
        self.magnetic_group.setCheckable(True)
        self.magnetic_group.setChecked(False)
        magnetic_layout = QVBoxLayout(self.magnetic_group)
        
        magnetic_info = QLabel("""
<p><b>Configure magnetic properties for spin-polarized calculations.</b></p>
<p>Specify starting magnetization for each element. Values typically range from -1 to 1.</p>
""")
        magnetic_info.setTextFormat(Qt.RichText)
        magnetic_info.setWordWrap(True)
        magnetic_layout.addWidget(magnetic_info)
        
        # Preset selector
        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Preset:"))
        self.magnetic_preset_combo = QComboBox()
        self.magnetic_preset_combo.addItems(["Custom", "Ferromagnetic", "Antiferromagnetic"])
        self.magnetic_preset_combo.currentTextChanged.connect(self._on_magnetic_preset_changed)
        preset_layout.addWidget(self.magnetic_preset_combo)
        magnetic_layout.addLayout(preset_layout)
        
        # Container for per-element magnetic inputs
        self.magnetic_container = QWidget()
        self.magnetic_container_layout = QFormLayout(self.magnetic_container)
        magnetic_layout.addWidget(self.magnetic_container)
        
        # Store magnetic inputs
        self.magnetic_edits = {}
        
        scroll_layout.addWidget(self.magnetic_group)
        
        # Hubbard (DFT+U) Configuration (optional)
        self.hubbard_group = QGroupBox("⚛️ Hubbard (DFT+U) Configuration (Optional)")
        self.hubbard_group.setCheckable(True)
        self.hubbard_group.setChecked(False)
        hubbard_layout = QVBoxLayout(self.hubbard_group)
        
        hubbard_info = QLabel("""
<p><b>Configure Hubbard U corrections for strongly correlated systems.</b></p>
<p>Specify U values (in eV) for transition metals and rare earths.</p>
""")
        hubbard_info.setTextFormat(Qt.RichText)
        hubbard_info.setWordWrap(True)
        hubbard_layout.addWidget(hubbard_info)
        
        # Hubbard format selector
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("Format:"))
        self.hubbard_format_combo = QComboBox()
        self.hubbard_format_combo.addItems(["New (QE 7.x)", "Old (QE 6.x)"])
        format_layout.addWidget(self.hubbard_format_combo)
        hubbard_layout.addLayout(format_layout)
        
        # Container for per-element Hubbard inputs
        self.hubbard_container = QWidget()
        self.hubbard_container_layout = QFormLayout(self.hubbard_container)
        hubbard_layout.addWidget(self.hubbard_container)
        
        # Store Hubbard inputs
        self.hubbard_u_edits = {}
        self.hubbard_orbital_edits = {}
        
        scroll_layout.addWidget(self.hubbard_group)
        
        # Workflow-specific options
        self.relax_group = QGroupBox("🔄 Relaxation Options")
        relax_layout = QFormLayout(self.relax_group)
        
        self.relax_type_combo = QComboBox()
        self.relax_type_combo.addItems(["relax", "vc-relax"])
        relax_layout.addRow("Relaxation Type:", self.relax_type_combo)
        
        self.forc_conv_edit = QLineEdit("1.0e-3")
        relax_layout.addRow("Force Convergence (Ry/bohr):", self.forc_conv_edit)
        
        self.relax_group.setVisible(False)
        scroll_layout.addWidget(self.relax_group)
        
        # Build Workflow Button
        build_group = QGroupBox("✨ Build Workflow")
        build_layout = QVBoxLayout(build_group)
        
        build_info = QLabel("""
<p>Click below to build the workflow using <b>workflow modules</b>.</p>
""")
        build_info.setTextFormat(Qt.RichText)
        build_layout.addWidget(build_info)
        
        build_btn = QPushButton("🔄 Build Workflow")
        build_btn.clicked.connect(self._build_workflow)
        build_layout.addWidget(build_btn)
        
        scroll_layout.addWidget(build_group)
        
        # Results area
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        scroll_layout.addWidget(self.results_label)
        
        # Workflow summary
        summary_group = QGroupBox("📊 Current Workflow")
        summary_layout = QVBoxLayout(summary_group)
        
        self.summary_text = QTextEdit()
        self.summary_text.setReadOnly(True)
        self.summary_text.setMaximumHeight(150)
        summary_layout.addWidget(self.summary_text)
        
        scroll_layout.addWidget(summary_group)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
        
        # Load initial data
        self._load_machines()
        self._update_structure_status()
    
    def _load_machines(self):
        """Load available machines."""
        if not XESPRESSO_AVAILABLE:
            return
        
        self.machine_combo.blockSignals(True)
        try:
            machines = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            for machine in machines:
                self.machine_combo.addItem(machine)
        except Exception:
            pass
        finally:
            self.machine_combo.blockSignals(False)
        # Manually trigger the handler for the first item if any
        if self.machine_combo.count() > 0:
            self._on_machine_changed(self.machine_combo.currentText())
    
    def _on_machine_changed(self, machine_name):
        """Handle machine selection change."""
        if not machine_name or not XESPRESSO_AVAILABLE:
            return
        
        # Prevent recursive updates
        if self._loading:
            return
        
        try:
            self._loading = True
            machine = load_machine(DEFAULT_CONFIG_PATH, machine_name, DEFAULT_MACHINES_DIR, return_object=True)
            self.session_state['workflow_machine'] = machine
            self.session_state['selected_machine'] = machine_name
            
            info = f"Type: {machine.execution}"
            if machine.scheduler:
                info += f", Scheduler: {machine.scheduler}"
            self.machine_info_label.setText(info)
            
            self._load_codes(machine_name)
        except Exception as e:
            self.machine_info_label.setText(f"Error: {e}")
        finally:
            self._loading = False
    
    def _load_codes(self, machine_name):
        """Load codes for a machine."""
        if not XESPRESSO_AVAILABLE:
            return
        
        self.version_combo.blockSignals(True)
        self.code_combo.blockSignals(True)
        try:
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, verbose=False)
            
            self.version_combo.clear()
            self.code_combo.clear()
            
            if codes and codes.has_any_codes():
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
        finally:
            self.version_combo.blockSignals(False)
            self.code_combo.blockSignals(False)
        # Manually trigger the handler for the first item if any
        if self.version_combo.count() > 0:
            self._on_version_changed(self.version_combo.currentText())
    
    def _on_version_changed(self, version):
        """Handle version selection change."""
        if not version or not XESPRESSO_AVAILABLE:
            return
        
        machine_name = self.machine_combo.currentText()
        if not machine_name:
            return
        
        self.code_combo.blockSignals(True)
        try:
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, version=version, verbose=False)
            
            self.code_combo.clear()
            if codes:
                version_codes = codes.get_all_codes(version=version)
                for code_name in version_codes.keys():
                    self.code_combo.addItem(code_name)
                
                idx = self.code_combo.findText('pw')
                if idx >= 0:
                    self.code_combo.setCurrentIndex(idx)
        except Exception:
            pass
        finally:
            self.code_combo.blockSignals(False)
    
    def _on_workflow_changed(self, workflow_type):
        """Handle workflow type change."""
        show_relax = workflow_type in ["SCF + Relaxation", "SCF + Relax + SCF (on relaxed structure)"]
        self.relax_group.setVisible(show_relax)
    
    def _on_kpts_mode_changed(self, checked):
        """Handle k-points mode change."""
        self.kspacing_widget.setVisible(checked)
        self.explicit_widget.setVisible(not checked)
    
    def _update_structure_status(self):
        """Update structure status display."""
        atoms = self.session_state.get('current_structure')
        
        if atoms is not None and ASE_AVAILABLE:
            if isinstance(atoms, Atoms):
                formula = atoms.get_chemical_formula()
                natoms = len(atoms)
                self.structure_status.setText(f"✅ Structure loaded: {formula} ({natoms} atoms)")
                self.structure_status.setStyleSheet("color: green;")
                
                self._update_pseudo_inputs(atoms)
            else:
                self.structure_status.setText("⚠️ Invalid structure type")
                self.structure_status.setStyleSheet("color: orange;")
        else:
            self.structure_status.setText("⚠️ No structure loaded. Please load a structure first.")
            self.structure_status.setStyleSheet("color: orange;")
    
    def _update_pseudo_inputs(self, atoms):
        """Update pseudopotential input fields."""
        elements = set(atoms.get_chemical_symbols())
        
        if PSEUDO_SELECTOR_AVAILABLE and hasattr(self, 'pseudo_selector'):
            # Use the advanced selector widget
            self.pseudo_selector.set_elements(elements)
        else:
            # Fallback to simple manual inputs
            while self.pseudo_container_layout.count():
                item = self.pseudo_container_layout.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
            
            self.pseudo_edits = {}
            
            self.pseudo_info_label.setText(f"Configure pseudopotentials for: {', '.join(sorted(elements))}")
            
            for element in sorted(elements):
                edit = QLineEdit()
                edit.setPlaceholderText(f"e.g., {element}.UPF")
                self.pseudo_edits[element] = edit
                self.pseudo_container_layout.addRow(f"{element}:", edit)
        
        # Update magnetic inputs
        self._update_magnetic_inputs(elements)
        
        # Update Hubbard inputs
        self._update_hubbard_inputs(elements)
    
    def _update_magnetic_inputs(self, elements):
        """Update magnetic configuration input fields for the given elements."""
        # Clear existing inputs
        while self.magnetic_container_layout.count():
            item = self.magnetic_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.magnetic_edits = {}
        
        # Predefined magnetic moments for common elements
        predefined_moments = {
            'Fe': 2.2, 'Co': 1.7, 'Ni': 0.6, 'Mn': 5.0, 'Cr': 3.0,
            'V': 2.0, 'Ti': 1.0, 'Gd': 7.0, 'Nd': 3.0, 'Sm': 5.0
        }
        
        for element in sorted(elements):
            edit = QDoubleSpinBox()
            edit.setRange(-10.0, 10.0)
            edit.setSingleStep(0.1)
            edit.setDecimals(2)
            # Set default value based on predefined moments or 0
            default_val = predefined_moments.get(element, 0.0)
            edit.setValue(default_val)
            self.magnetic_edits[element] = edit
            self.magnetic_container_layout.addRow(f"{element}:", edit)
    
    def _update_hubbard_inputs(self, elements):
        """Update Hubbard U input fields for the given elements."""
        # Clear existing inputs
        while self.hubbard_container_layout.count():
            item = self.hubbard_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.hubbard_u_edits = {}
        self.hubbard_orbital_edits = {}
        
        # Common Hubbard U values and orbitals
        common_hubbard = {
            'Fe': (4.0, '3d'), 'Co': (3.5, '3d'), 'Ni': (3.0, '3d'),
            'Mn': (4.0, '3d'), 'Cr': (3.5, '3d'), 'V': (3.0, '3d'),
            'Ti': (2.5, '3d'), 'Cu': (4.0, '3d'), 'Zn': (3.0, '3d'),
            'Ce': (4.5, '4f'), 'Nd': (5.0, '4f'), 'Gd': (6.0, '4f')
        }
        
        for element in sorted(elements):
            # Create a horizontal layout for each element
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            
            # U value spinbox
            u_edit = QDoubleSpinBox()
            u_edit.setRange(0.0, 15.0)
            u_edit.setSingleStep(0.5)
            u_edit.setDecimals(1)
            u_edit.setSuffix(" eV")
            default_u, default_orbital = common_hubbard.get(element, (0.0, '3d'))
            u_edit.setValue(default_u)
            row_layout.addWidget(QLabel("U:"))
            row_layout.addWidget(u_edit)
            
            # Orbital selector
            orbital_edit = QComboBox()
            orbital_edit.addItems(['3d', '4d', '5d', '4f', '5f', '2p', '3p'])
            orbital_edit.setCurrentText(default_orbital)
            row_layout.addWidget(QLabel("Orbital:"))
            row_layout.addWidget(orbital_edit)
            
            self.hubbard_u_edits[element] = u_edit
            self.hubbard_orbital_edits[element] = orbital_edit
            self.hubbard_container_layout.addRow(f"{element}:", row_widget)
    
    def _on_magnetic_preset_changed(self, preset):
        """Handle magnetic preset selection change."""
        if preset == "Custom":
            return
        
        # Predefined magnetic moments
        predefined_moments = {
            'Fe': 2.2, 'Co': 1.7, 'Ni': 0.6, 'Mn': 5.0, 'Cr': 3.0,
            'V': 2.0, 'Ti': 1.0, 'Gd': 7.0, 'Nd': 3.0, 'Sm': 5.0
        }
        
        for element, edit in self.magnetic_edits.items():
            base_moment = predefined_moments.get(element, 0.0)
            if preset == "Ferromagnetic":
                edit.setValue(base_moment)
            elif preset == "Antiferromagnetic":
                # For AFM, we'd need to alternate signs - simplified here
                edit.setValue(base_moment)
    
    def _get_config(self):
        """Get the current configuration as a dictionary."""
        config = {}
        
        config['workflow_type'] = self.workflow_combo.currentText()
        config['ecutwfc'] = self.ecutwfc_spin.value()
        config['ecutrho'] = self.ecutrho_spin.value()
        config['occupations'] = self.occupations_combo.currentText()
        
        # K-points
        if self.kspacing_radio.isChecked():
            config['kspacing'] = self.kspacing_spin.value()
        else:
            config['kpts'] = (self.k1_spin.value(), self.k2_spin.value(), self.k3_spin.value())
        
        # Relaxation options
        if self.relax_group.isVisible():
            config['relax_type'] = self.relax_type_combo.currentText()
            config['forc_conv_thr'] = float(self.forc_conv_edit.text())
        
        # Pseudopotentials
        if PSEUDO_SELECTOR_AVAILABLE and hasattr(self, 'pseudo_selector'):
            config['pseudopotentials'] = self.pseudo_selector.get_pseudopotentials()
        else:
            config['pseudopotentials'] = {}
            for element, edit in getattr(self, 'pseudo_edits', {}).items():
                pseudo = edit.text().strip()
                if pseudo:
                    config['pseudopotentials'][element] = pseudo
        
        # Magnetic configuration (optional)
        config['enable_magnetism'] = self.magnetic_group.isChecked()
        if config['enable_magnetism']:
            config['magnetic_config'] = {}
            for element, edit in getattr(self, 'magnetic_edits', {}).items():
                config['magnetic_config'][element] = [edit.value()]
        
        # Hubbard (DFT+U) configuration (optional)
        config['enable_hubbard'] = self.hubbard_group.isChecked()
        if config['enable_hubbard']:
            config['hubbard_format'] = 'new' if 'New' in self.hubbard_format_combo.currentText() else 'old'
            config['hubbard_u'] = {}
            config['hubbard_orbitals'] = {}
            for element in getattr(self, 'hubbard_u_edits', {}):
                u_value = self.hubbard_u_edits[element].value()
                if u_value > 0:
                    config['hubbard_u'][element] = u_value
                    config['hubbard_orbitals'][element] = self.hubbard_orbital_edits[element].currentText()
        
        # Machine and Queue configuration
        config['machine_name'] = self.machine_combo.currentText()
        config['qe_version'] = self.version_combo.currentText()
        config['selected_code'] = self.code_combo.currentText()
        
        # Get the machine object and convert to queue configuration
        # This ensures the scheduler settings (slurm/direct, local/remote) are passed to workflow
        machine = self.session_state.get('calc_machine')
        if machine:
            try:
                # Machine's to_queue() method returns the complete queue configuration
                # including execution mode, scheduler type, and all settings
                config['queue'] = machine.to_queue()
            except (AttributeError, TypeError) as e:
                import logging
                logging.warning(f"Could not convert machine to queue: {e}")
                # Don't set a default queue - let workflow handle it
        
        return config
    
    def _build_workflow(self):
        """Build the workflow using GUIWorkflow class to create calculator objects."""
        atoms = self.session_state.get('current_structure')
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded. Please load a structure first.")
            return
        
        config = self._get_config()
        
        # Validate pseudopotentials
        if not config.get('pseudopotentials'):
            QMessageBox.warning(self, "Warning", "Please specify pseudopotentials for all elements")
            return
        
        # Validate machine selection
        machine = self.session_state.get('calc_machine')
        if not machine:
            QMessageBox.warning(
                self, 
                "No Machine Selected",
                "Please select a machine configuration before building the workflow.\n\n"
                "A machine defines where and how the calculations will run (local/remote, scheduler type, resources, etc.).\n\n"
                "Go to Machine Configuration page to create or select a machine."
            )
            return
        
        # Check if GUIWorkflow is available
        if not WORKFLOW_AVAILABLE:
            QMessageBox.critical(
                self,
                "Workflow Module Not Available",
                "The gui.workflows.base module is not available.\n\n"
                "Workflow functionality requires the GUIWorkflow class."
            )
            return
        
        try:
            # Get workflow type
            workflow_type = config['workflow_type']
            
            # Create base label for workflow
            base_label = config.get('label', f"workflow/{atoms.get_chemical_formula()}")
            
            # Create GUIWorkflow object
            workflow = GUIWorkflow(atoms, config, base_label)
            
            # Add calculation steps based on workflow type
            if workflow_type == "Single SCF":
                # Single SCF calculation
                scf_config = self._create_step_config(config, 'scf')
                workflow.add_calculation("scf", scf_config)
                steps = ["scf"]
                
            elif workflow_type == "SCF + Relaxation":
                # SCF followed by relaxation
                scf_config = self._create_step_config(config, 'scf')
                workflow.add_calculation("scf", scf_config)
                
                relax_config = self._create_step_config(config, config.get('relax_type', 'relax'))
                workflow.add_calculation("relax", relax_config)
                steps = ["scf", config.get('relax_type', 'relax')]
                
            elif workflow_type == "SCF + Relax + SCF (on relaxed structure)":
                # Initial SCF
                scf_initial_config = self._create_step_config(config, 'scf')
                workflow.add_calculation("scf_initial", scf_initial_config)
                
                # Relaxation
                relax_config = self._create_step_config(config, config.get('relax_type', 'relax'))
                workflow.add_calculation("relax", relax_config)
                
                # Final SCF on relaxed structure
                scf_final_config = self._create_step_config(config, 'scf')
                workflow.add_calculation("scf_final", scf_final_config)
                steps = ["scf_initial", config.get('relax_type', 'relax'), "scf_final"]
                
            else:
                # Custom workflow - just store config
                steps = ["custom"]
            
            # Store workflow object in session state
            self.session_state['gui_workflow'] = workflow
            self.session_state['workflow_config'] = config
            
            # Update summary display
            self.summary_text.setText(f"""
Workflow Type: {workflow_type}
Steps: {len(steps)}
  - {chr(10).join([f'  {i+1}. {s}' for i, s in enumerate(steps)])}

Energy Cutoff: {config['ecutwfc']} Ry
Machine: {config.get('machine_name', 'Not selected')}
QE Version: {config.get('qe_version', 'Not selected')}

Calculators created: {len(workflow.calculations)}
""")
            
            self.results_label.setText(f"""
✅ <b>Workflow built successfully!</b>

The GUIWorkflow object has been created with {len(workflow.calculations)} calculation step(s).
Each step has its own Espresso calculator object ready to use.

Go to <b>Job Submission</b> page to execute the workflow steps.
""")
            self.results_label.setStyleSheet("color: green;")
            self.results_label.setTextFormat(Qt.RichText)
            
            QMessageBox.information(
                self,
                "Workflow Built",
                f"Workflow '{workflow_type}' has been built with {len(workflow.calculations)} calculator(s).\n\n"
                f"Each step has a prepared Espresso calculator object.\n\n"
                "Go to Job Submission page to execute the workflow."
            )
            
        except Exception as e:
            import traceback
            error_details = traceback.format_exc()
            QMessageBox.critical(
                self,
                "Workflow Build Failed",
                f"Failed to build workflow: {e}\n\nSee console for details."
            )
            print(f"Error building workflow:\n{error_details}")
    
    def _create_step_config(self, base_config, calc_type):
        """
        Create configuration for a specific calculation step.
        
        Args:
            base_config: Base workflow configuration
            calc_type: Type of calculation (scf, relax, vc-relax, etc.)
            
        Returns:
            dict: Configuration for this step
        """
        step_config = base_config.copy()
        step_config['calc_type'] = calc_type
        
        # Add relaxation-specific parameters if needed
        if calc_type in ['relax', 'vc-relax']:
            step_config['forc_conv_thr'] = base_config.get('forc_conv_thr', 1.0e-3)
        
        return step_config
    
    def refresh(self):
        """Refresh the page."""
        # Use loading guard to prevent infinite loops
        if self._loading:
            return
        self._loading = True
        try:
            self._load_machines()
            self._update_structure_status()
        finally:
            self._loading = False
    
    def save_state(self):
        """Save current page state to session state.
        
        This is called before the session is saved to disk to ensure
        all current UI values are captured in the session state.
        """
        config = self._get_config()
        self.session_state['workflow_config'] = config

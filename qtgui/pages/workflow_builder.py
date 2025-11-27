"""
Workflow Builder Page for xespresso PyQt GUI.

This page creates multi-step workflows.
"""

import os

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QDoubleSpinBox,
    QSpinBox, QCheckBox, QTextEdit, QRadioButton
)
from PyQt5.QtCore import Qt

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


class WorkflowBuilderPage(QWidget):
    """Workflow builder page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
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
        
        self.pseudo_info_label = QLabel("Load a structure to configure pseudopotentials")
        pseudo_layout.addWidget(self.pseudo_info_label)
        
        self.pseudo_container = QWidget()
        self.pseudo_container_layout = QFormLayout(self.pseudo_container)
        pseudo_layout.addWidget(self.pseudo_container)
        
        scroll_layout.addWidget(pseudo_group)
        
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
        
        try:
            machines = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            for machine in machines:
                self.machine_combo.addItem(machine)
        except Exception as e:
            pass
    
    def _on_machine_changed(self, machine_name):
        """Handle machine selection change."""
        if not machine_name or not XESPRESSO_AVAILABLE:
            return
        
        try:
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
    
    def _load_codes(self, machine_name):
        """Load codes for a machine."""
        if not XESPRESSO_AVAILABLE:
            return
        
        try:
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR)
            
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
        except Exception as e:
            pass
    
    def _on_version_changed(self, version):
        """Handle version selection change."""
        if not version or not XESPRESSO_AVAILABLE:
            return
        
        machine_name = self.machine_combo.currentText()
        if not machine_name:
            return
        
        try:
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, version=version)
            
            self.code_combo.clear()
            if codes:
                version_codes = codes.get_all_codes(version=version)
                for code_name in version_codes.keys():
                    self.code_combo.addItem(code_name)
                
                idx = self.code_combo.findText('pw')
                if idx >= 0:
                    self.code_combo.setCurrentIndex(idx)
        except Exception as e:
            pass
    
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
        while self.pseudo_container_layout.count():
            item = self.pseudo_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        elements = set(atoms.get_chemical_symbols())
        self.pseudo_edits = {}
        
        self.pseudo_info_label.setText(f"Configure pseudopotentials for: {', '.join(sorted(elements))}")
        
        for element in sorted(elements):
            edit = QLineEdit()
            edit.setPlaceholderText(f"e.g., {element}.UPF")
            self.pseudo_edits[element] = edit
            self.pseudo_container_layout.addRow(f"{element}:", edit)
    
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
        config['pseudopotentials'] = {}
        for element, edit in getattr(self, 'pseudo_edits', {}).items():
            pseudo = edit.text().strip()
            if pseudo:
                config['pseudopotentials'][element] = pseudo
        
        # Machine
        config['machine_name'] = self.machine_combo.currentText()
        config['qe_version'] = self.version_combo.currentText()
        config['selected_code'] = self.code_combo.currentText()
        
        return config
    
    def _build_workflow(self):
        """Build the workflow."""
        atoms = self.session_state.get('current_structure')
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded. Please load a structure first.")
            return
        
        config = self._get_config()
        
        # Validate pseudopotentials
        if not config.get('pseudopotentials'):
            QMessageBox.warning(self, "Warning", "Please specify pseudopotentials for all elements")
            return
        
        # Store config in session state
        self.session_state['workflow_config'] = config
        
        # Build workflow summary
        workflow_type = config['workflow_type']
        
        steps = []
        if workflow_type == "Single SCF":
            steps = ["scf"]
        elif workflow_type == "SCF + Relaxation":
            steps = ["scf", config.get('relax_type', 'relax')]
        elif workflow_type == "SCF + Relax + SCF (on relaxed structure)":
            steps = ["scf_initial", config.get('relax_type', 'relax'), "scf_final"]
        else:
            steps = ["custom"]
        
        self.summary_text.setText(f"""
Workflow Type: {workflow_type}
Steps: {len(steps)}
  - {chr(10).join([f'  {i+1}. {s}' for i, s in enumerate(steps)])}

Energy Cutoff: {config['ecutwfc']} Ry
Machine: {config.get('machine_name', 'Not selected')}
QE Version: {config.get('qe_version', 'Not selected')}
""")
        
        self.results_label.setText("""
✅ <b>Workflow built successfully!</b>

The workflow configuration has been stored in session state.
Go to <b>Job Submission</b> page to execute the workflow steps.
""")
        self.results_label.setStyleSheet("color: green;")
        self.results_label.setTextFormat(Qt.RichText)
        
        QMessageBox.information(
            self,
            "Workflow Built",
            f"Workflow '{workflow_type}' has been built with {len(steps)} step(s).\n\nGo to Job Submission page to execute."
        )
    
    def refresh(self):
        """Refresh the page."""
        self._load_machines()
        self._update_structure_status()

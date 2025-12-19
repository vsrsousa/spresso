"""
Workflow Builder Page for xespresso PySide6 GUI.

This page creates multi-step workflows using the GUIWorkflow class.
"""

import os

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QDoubleSpinBox,
    QSpinBox, QCheckBox, QTextEdit, QRadioButton
)
from qtpy.QtCore import Qt

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
    from qtgui.workflows.base import GUIWorkflow
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
        
        # Mode Status (shows if calculation setup is active)
        self.mode_status = QLabel("")
        self.mode_status.setWordWrap(True)
        scroll_layout.addWidget(self.mode_status)
        
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

        # Save / Load workflow buttons
        save_btn = QPushButton("💾 Save Workflow")
        save_btn.clicked.connect(self._save_workflow)
        build_layout.addWidget(save_btn)

        load_btn = QPushButton("📂 Load Workflow")
        load_btn.clicked.connect(self._load_workflow)
        build_layout.addWidget(load_btn)
        
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

    # -- Persistence helpers -------------------------------------------------
    def _workflows_dir(self):
        # prefer session working directory under .xespresso/workflows
        base = os.path.expanduser("~/.xespresso/workflows")
        try:
            os.makedirs(base, exist_ok=True)
        except Exception:
            pass
        return base

    def _save_workflow(self):
        # Build a lightweight descriptor from current UI selections
        wf = self._build_workflow_descriptor(return_obj=True)
        if wf is None:
            QMessageBox.warning(self, "Save Workflow", "No workflow to save")
            return
        # Ask user for filename
        path = os.path.join(self._workflows_dir(), f"{wf.get('base_label','workflow')}.json")
        try:
            self._save_workflow_to_path(path)
            # Avoid showing modal dialogs during automated tests (pytest)
            if not os.environ.get("PYTEST_CURRENT_TEST"):
                QMessageBox.information(self, "Save Workflow", f"Saved workflow to: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Save Workflow", f"Could not save workflow: {e}")

    def _load_workflow(self):
        # Let user choose a workflow file from workflows dir
        start = self._workflows_dir()
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Workflow", start, "Workflow Files (*.json);;All Files (*)")
        if not file_path:
            return
        try:
            data = self._load_workflow_from_path(file_path)
            self._apply_loaded_workflow(data)
            QMessageBox.information(self, "Load Workflow", f"Loaded workflow from: {file_path}")
        except Exception as e:
            QMessageBox.warning(self, "Load Workflow", f"Could not load workflow: {e}")

    # Programmatic helpers -------------------------------------------------
    def _save_workflow_to_path(self, path: str):
        """Save current workflow descriptor to a given path (programmatic helper)."""
        wf = self._build_workflow_descriptor(return_obj=True)
        if wf is None:
            raise RuntimeError("No workflow to save")
        import json

        dpath = os.path.dirname(path)
        if dpath:
            os.makedirs(dpath, exist_ok=True)
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(wf, f, indent=2)
        return path

    def _load_workflow_from_path(self, path: str):
        """Load workflow descriptor from path and return the data."""
        import json

        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        # apply to UI for convenience
        try:
            self._apply_loaded_workflow(data)
        except Exception:
            pass
        return data

    def _build_workflow_descriptor(self, return_obj=False):
        # Build a minimal serializable workflow descriptor instead of GUIWorkflow
        base_label = f"wf_{self.session_state.get('session_name','session')}"
        cfg = {
            'ecutwfc': float(self.ecutwfc_spin.value()),
            'ecutrho': float(self.ecutrho_spin.value()),
            'occupations': str(self.occupations_combo.currentText()),
        }
        wf = {'base_label': base_label, 'config': cfg, 'calculations': []}
        # simple mapping based on workflow_combo
        choice = self.workflow_combo.currentText()
        if choice == 'Single SCF':
            wf['calculations'].append({'name': 'scf', 'config': cfg})
        elif choice == 'SCF + Relaxation':
            wf['calculations'].append({'name': 'scf', 'config': cfg})
            wf['calculations'].append({'name': 'relax', 'config': cfg})
        elif choice == 'SCF + Relax + SCF (on relaxed structure)':
            wf['calculations'].append({'name': 'scf_init', 'config': cfg})
            wf['calculations'].append({'name': 'relax', 'config': cfg})
            wf['calculations'].append({'name': 'scf_relaxed', 'config': cfg})
        else:
            # custom: no steps until user defines them
            pass

        # update summary text
        try:
            import json

            self.summary_text.setPlainText(json.dumps(wf, indent=2))
        except Exception:
            self.summary_text.setPlainText(str(wf))

        if return_obj:
            return wf
        return None

    def _apply_loaded_workflow(self, data):
        try:
            # apply basic config
            cfg = data.get('config', {})
            if 'ecutwfc' in cfg:
                try:
                    self.ecutwfc_spin.setValue(float(cfg.get('ecutwfc')))
                except Exception:
                    pass
            if 'occupations' in cfg:
                try:
                    idx = self.occupations_combo.findText(cfg.get('occupations'))
                    if idx >= 0:
                        self.occupations_combo.setCurrentIndex(idx)
                except Exception:
                    pass
            # set workflow selection based on number of calculations
            steps = data.get('calculations', [])
            if len(steps) == 1:
                self.workflow_combo.setCurrentText('Single SCF')
            elif len(steps) == 2:
                self.workflow_combo.setCurrentText('SCF + Relaxation')
            elif len(steps) == 3:
                self.workflow_combo.setCurrentText('SCF + Relax + SCF (on relaxed structure)')
            # update summary
            try:
                import json

                self.summary_text.setPlainText(json.dumps(data, indent=2))
            except Exception:
                self.summary_text.setPlainText(str(data))
        except Exception:
            pass
    
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
                self.structure_status.setStyleSheet("color: #d97706;")
        else:
            self.structure_status.setText("⚠️ No structure loaded. Please load a structure first.")
            self.structure_status.setStyleSheet("color: #d97706;")
    
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
        
        # Check if calculation setup mode is active
        workflow_mode = self.session_state.get('workflow_mode')
        if workflow_mode == 'single':
            reply = QMessageBox.question(
                self,
                "Switch to Multi-Step Workflow Mode?",
                "You previously prepared a single calculation.\n\n"
                "Using Workflow Builder will switch to multi-step workflow mode and replace the calculation configuration.\n\n"
                "Do you want to continue?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        # Check if GUIWorkflow is available
        if not WORKFLOW_AVAILABLE:
            QMessageBox.critical(
                self,
                "Workflow Module Not Available",
                "The qtgui.workflows.base module is not available.\n\n"
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
            # Mark that we're using multi-step workflow mode
            self.session_state['workflow_mode'] = 'multi'
            
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
            # Restore UI state from saved configuration
            self._restore_config_to_ui()
            # Update mode indicator
            self._update_mode_indicator()
        finally:
            self._loading = False
    
    def _restore_config_to_ui(self):
        """Restore UI state from saved workflow_config in session state.
        
        This ensures that when a session is loaded, the magnetic and hubbard
        checkboxes and other UI elements reflect the saved configuration,
        allowing continued editing after session save/reload.
        """
        config = self.session_state.get('workflow_config')
        if not config:
            return
        
        # Restore machine selection from config (check both current_machine_name and machine_name)
        # Session state may have current_machine_name at the top level
        machine_name = self.session_state.get('current_machine_name') or config.get('machine_name')
        if machine_name:
            idx = self.machine_combo.findText(machine_name)
            if idx >= 0:
                self.machine_combo.blockSignals(True)
                self.machine_combo.setCurrentIndex(idx)
                self.machine_combo.blockSignals(False)
                # Manually trigger machine changed to load codes
                self._on_machine_changed(machine_name)
        
        # Restore QE version selection
        qe_version = config.get('qe_version')
        if qe_version:
            idx = self.version_combo.findText(qe_version)
            if idx >= 0:
                self.version_combo.blockSignals(True)
                self.version_combo.setCurrentIndex(idx)
                self.version_combo.blockSignals(False)
                # Manually trigger version changed to load codes for this version
                self._on_version_changed(qe_version)
        
        # Restore selected code
        selected_code = config.get('selected_code')
        if selected_code:
            idx = self.code_combo.findText(selected_code)
            if idx >= 0:
                self.code_combo.setCurrentIndex(idx)
        
        # Restore magnetic configuration checkbox state
        if config.get('enable_magnetism'):
            self.magnetic_group.setChecked(True)
            # Restore magnetic values if elements exist
            if config.get('magnetic_config'):
                for element, mag_value in config['magnetic_config'].items():
                    if element in self.magnetic_edits:
                        # mag_value is stored as a list in the config
                        value = mag_value[0] if isinstance(mag_value, list) else mag_value
                        self.magnetic_edits[element].setValue(value)
        else:
            # Explicitly uncheck if the config says magnetism is disabled
            self.magnetic_group.setChecked(False)
        
        # Restore Hubbard configuration checkbox state
        if config.get('enable_hubbard'):
            self.hubbard_group.setChecked(True)
            # Restore Hubbard U values if elements exist
            if config.get('hubbard_u'):
                for element, u_value in config['hubbard_u'].items():
                    if element in self.hubbard_u_edits:
                        self.hubbard_u_edits[element].setValue(u_value)
            # Restore orbital selections
            if config.get('hubbard_orbitals'):
                for element, orbital in config['hubbard_orbitals'].items():
                    if element in self.hubbard_orbital_edits:
                        idx = self.hubbard_orbital_edits[element].findText(orbital)
                        if idx >= 0:
                            self.hubbard_orbital_edits[element].setCurrentIndex(idx)
            # Restore Hubbard format
            if config.get('hubbard_format'):
                format_str = "New (QE 7.x)" if config['hubbard_format'] == 'new' else "Old (QE 6.x)"
                idx = self.hubbard_format_combo.findText(format_str)
                if idx >= 0:
                    self.hubbard_format_combo.setCurrentIndex(idx)
        else:
            # Explicitly uncheck if the config says Hubbard is disabled
            self.hubbard_group.setChecked(False)
    
    def save_state(self):
        """Save current page state to session state.
        
        This is called before the session is saved to disk to ensure
        all current UI values are captured in the session state.
        
        Only save if we're in multi-step workflow mode or if mode isn't set yet.
        """
        workflow_mode = self.session_state.get('workflow_mode')
        
        # Only save if we're in multi mode, or mode is not set (for backward compatibility)
        if workflow_mode == 'single':
            # Don't overwrite calculation setup's config
            return
        
        config = self._get_config()
        existing_config = self.session_state.get('workflow_config')
        
        # Merge with existing config to preserve settings from other pages
        # (e.g., Calculation Setup page might have set magnetic/Hubbard settings)
        merged_config = self._merge_configs(config, existing_config)
        if merged_config is not None:
            self.session_state['workflow_config'] = merged_config
            # If mode isn't set, set it to multi since we're saving from here
            if workflow_mode is None:
                self.session_state['workflow_mode'] = 'multi'
    
    def _merge_configs(self, config, existing_config):
        """Merge current config with existing config to avoid overwriting other pages' settings.
        
        Args:
            config (dict): New configuration from current UI state
            existing_config (dict or None): Existing workflow configuration in session state
            
        Returns:
            dict: Merged configuration to save
            
        Decision logic:
            - If no existing config: Save current config
            - If existing config exists: Merge current config with existing, preferring
              existing values for magnetic/hubbard if they were set there
        """
        # If no existing config, save current one
        if not existing_config:
            return config
        
        # Start with a copy of the current config
        merged = config.copy()
        
        # If existing config has magnetic/hubbard settings that are enabled,
        # preserve them unless current config also has them enabled
        if existing_config.get('enable_magnetism') and not config.get('enable_magnetism'):
            # Preserve existing magnetic config if current page doesn't have it enabled
            merged['enable_magnetism'] = existing_config['enable_magnetism']
            merged['magnetic_config'] = existing_config.get('magnetic_config', {})
        
        if existing_config.get('enable_hubbard') and not config.get('enable_hubbard'):
            # Preserve existing Hubbard config if current page doesn't have it enabled
            merged['enable_hubbard'] = existing_config['enable_hubbard']
            merged['hubbard_u'] = existing_config.get('hubbard_u', {})
            merged['hubbard_orbitals'] = existing_config.get('hubbard_orbitals', {})
            merged['hubbard_format'] = existing_config.get('hubbard_format', 'new')
        
        return merged
    
    def _update_mode_indicator(self):
        """Update the mode status indicator to show which workflow mode is active."""
        workflow_mode = self.session_state.get('workflow_mode')
        
        if workflow_mode == 'single':
            self.mode_status.setText(
                "ℹ️ <b>Single Calculation Mode Active</b> - "
                "Currently using Calculation Setup. To use this page, click 'Build Workflow' and confirm the mode switch."
            )
            self.mode_status.setStyleSheet("color: #ff9800; font-weight: bold; background-color: #fff3e0; padding: 10px; border-radius: 5px;")
            self.mode_status.setTextFormat(Qt.RichText)
        elif workflow_mode == 'multi':
            self.mode_status.setText("✅ <b>Multi-Step Workflow Mode Active</b>")
            self.mode_status.setStyleSheet("color: green; font-weight: bold;")
            self.mode_status.setTextFormat(Qt.RichText)
        else:
            # No mode set yet
            self.mode_status.setText("")
            self.mode_status.setStyleSheet("")


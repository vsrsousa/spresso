"""
Calculation Setup Page for xespresso PySide6 GUI.

This page is responsible for configuring calculations.
"""

import os

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QDoubleSpinBox,
    QSpinBox, QCheckBox, QTextEdit, QRadioButton, QButtonGroup
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QPalette, QColor

from qtgui.utils import ASE_ESPRESSO_COMMAND_TEMPLATE

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


# Predefined magnetic moments for common elements
PREDEFINED_MAGNETIC_MOMENTS = {
    'Fe': 2.2, 'Co': 1.7, 'Ni': 0.6, 'Mn': 5.0, 'Cr': 3.0,
    'V': 2.0, 'Ti': 1.0, 'Gd': 7.0, 'Nd': 3.0
}

# Typical Hubbard U values for common elements
TYPICAL_HUBBARD_U = {
    'Fe': 4.0, 'Co': 3.5, 'Ni': 3.0, 'Mn': 4.0, 'Cr': 3.5,
    'V': 3.0, 'Ti': 2.5, 'Cu': 4.0, 'Zn': 4.0,
    'Gd': 6.0, 'Nd': 5.0, 'Ce': 5.0, 'O': 0.0
}


class CalculationSetupPage(QWidget):
    """Calculation setup page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self._loading = False  # Guard to prevent infinite loops
        # Initialize dictionaries for dynamic inputs
        self.magnetic_edits = {}
        self.hubbard_edits = {}
        self.hubbard_orbital_edits = {}
        self.hubbard_checkboxes = {}  # Track which elements have Hubbard enabled
        self._setup_ui()
    
    def _get_combobox_stylesheet(self):
        """Get the stylesheet for comboboxes with proper contrast."""
        return """
            QComboBox {
                padding: 2px;
            }
            QComboBox::drop-down {
                border: 0px;
            }
            QComboBox QAbstractItemView {
                selection-background-color: #B3D9FF;
                selection-color: black;
                background-color: #F5F5F5;
                color: black;
            }
            QComboBox QAbstractItemView::item:hover {
                background-color: #D6EBFF;
                color: black;
            }
        """
    
    def _apply_combobox_styling(self):
        """Apply custom styling to comboboxes for better visibility."""
        combobox_style = self._get_combobox_stylesheet()
        
        # Apply to all comboboxes in the page
        for widget in self.findChildren(QComboBox):
            widget.setStyleSheet(combobox_style)
    
    def _setup_ui(self):
        """Setup the user interface."""
        main_layout = QVBoxLayout(self)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Header
        header_label = QLabel("<h2>📊 Calculation Setup</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p>Configure your calculation parameters. This page uses <b>calculation modules</b>
to prepare atoms and Espresso calculator objects following xespresso's design patterns.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        # Structure Status
        self.structure_status = QLabel("")
        self.structure_status.setWordWrap(True)
        scroll_layout.addWidget(self.structure_status)
        
        # Mode Status (shows if workflow builder is active)
        self.mode_status = QLabel("")
        self.mode_status.setWordWrap(True)
        scroll_layout.addWidget(self.mode_status)
        
        if not XESPRESSO_AVAILABLE:
            error_label = QLabel("❌ xespresso modules not available.")
            error_label.setStyleSheet("color: red; font-weight: bold;")
            scroll_layout.addWidget(error_label)
        
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
        
        # Calculation Type
        calc_group = QGroupBox("⚙️ Calculation Type")
        calc_layout = QFormLayout(calc_group)
        
        self.calc_type_combo = QComboBox()
        self.calc_type_combo.addItems(["scf", "relax", "vc-relax"])
        calc_layout.addRow("Calculation:", self.calc_type_combo)
        
        self.label_edit = QLineEdit()
        self.label_edit.setPlaceholderText("e.g., Al/scf")
        calc_layout.addRow("Label (subfolder name):", self.label_edit)
        
        scroll_layout.addWidget(calc_group)
        
        # Basic Parameters
        params_group = QGroupBox("🔧 Basic Parameters")
        params_layout = QFormLayout(params_group)
        
        self.ecutwfc_spin = QDoubleSpinBox()
        self.ecutwfc_spin.setRange(10.0, 200.0)
        self.ecutwfc_spin.setValue(50.0)
        self.ecutwfc_spin.setSingleStep(5.0)
        self.ecutwfc_spin.setSuffix(" Ry")
        self.ecutwfc_spin.valueChanged.connect(self._on_ecutwfc_changed)
        params_layout.addRow("Energy Cutoff:", self.ecutwfc_spin)
        
        self.ecutrho_spin = QDoubleSpinBox()
        self.ecutrho_spin.setRange(40.0, 1600.0)
        self.ecutrho_spin.setValue(400.0)
        self.ecutrho_spin.setSingleStep(20.0)
        self.ecutrho_spin.setSuffix(" Ry")
        params_layout.addRow("Charge Density Cutoff:", self.ecutrho_spin)
        
        self.occupations_combo = QComboBox()
        self.occupations_combo.addItems(["smearing", "fixed", "tetrahedra"])
        self.occupations_combo.currentTextChanged.connect(self._on_occupations_changed)
        params_layout.addRow("Occupations:", self.occupations_combo)
        
        self.conv_thr_edit = QLineEdit("1.0e-8")
        params_layout.addRow("Convergence Threshold:", self.conv_thr_edit)
        
        scroll_layout.addWidget(params_group)
        
        # Smearing Parameters
        self.smearing_group = QGroupBox("📊 Smearing Parameters")
        smearing_layout = QFormLayout(self.smearing_group)
        
        self.smearing_combo = QComboBox()
        self.smearing_combo.addItems(["marzari-vanderbilt", "gaussian", "methfessel-paxton", "fermi-dirac"])
        smearing_layout.addRow("Smearing Type:", self.smearing_combo)
        
        self.degauss_spin = QDoubleSpinBox()
        self.degauss_spin.setRange(0.001, 0.1)
        self.degauss_spin.setValue(0.02)
        self.degauss_spin.setSingleStep(0.005)
        self.degauss_spin.setDecimals(4)
        self.degauss_spin.setSuffix(" Ry")
        smearing_layout.addRow("Degauss:", self.degauss_spin)
        
        scroll_layout.addWidget(self.smearing_group)
        
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
        
        scroll_layout.addWidget(self.magnetic_group)
        
        # Hubbard (DFT+U) Configuration (optional)
        self.hubbard_group = QGroupBox("⚛️ Hubbard (DFT+U) Configuration (Optional)")
        self.hubbard_group.setCheckable(True)
        self.hubbard_group.setChecked(False)
        hubbard_layout = QVBoxLayout(self.hubbard_group)
        
        hubbard_info = QLabel("""
<p><b>Configure Hubbard U corrections for strongly correlated systems.</b></p>
<p>Common for transition metals (Fe, Mn, Co, Ni) and rare earths. U values typically range from 2-8 eV.</p>
""")
        hubbard_info.setTextFormat(Qt.RichText)
        hubbard_info.setWordWrap(True)
        hubbard_layout.addWidget(hubbard_info)
        
        # Hubbard format selector
        format_layout = QHBoxLayout()
        format_layout.addWidget(QLabel("Format:"))
        self.hubbard_format_combo = QComboBox()
        self.hubbard_format_combo.addItems(["New (QE >= 7.0)", "Old (QE < 7.0)"])
        format_layout.addWidget(self.hubbard_format_combo)
        hubbard_layout.addLayout(format_layout)
        
        # Container for per-element Hubbard U inputs
        self.hubbard_container = QWidget()
        self.hubbard_container_layout = QFormLayout(self.hubbard_container)
        hubbard_layout.addWidget(self.hubbard_container)
        
        scroll_layout.addWidget(self.hubbard_group)
        
        # Resources Configuration
        resources_group = QGroupBox("⚙️ Resources Configuration")
        resources_layout = QVBoxLayout(resources_group)
        
        self.adjust_resources_check = QCheckBox("Adjust Resources")
        self.adjust_resources_check.stateChanged.connect(self._on_adjust_resources_changed)
        resources_layout.addWidget(self.adjust_resources_check)
        
        self.resources_widget = QWidget()
        resources_form = QFormLayout(self.resources_widget)
        
        self.nprocs_spin = QSpinBox()
        self.nprocs_spin.setRange(1, 256)
        self.nprocs_spin.setValue(1)
        resources_form.addRow("Number of Processors:", self.nprocs_spin)
        
        self.nodes_spin = QSpinBox()
        self.nodes_spin.setRange(1, 1000)
        self.nodes_spin.setValue(1)
        resources_form.addRow("Nodes:", self.nodes_spin)
        
        self.ntasks_spin = QSpinBox()
        self.ntasks_spin.setRange(1, 256)
        self.ntasks_spin.setValue(16)
        resources_form.addRow("Tasks per Node:", self.ntasks_spin)
        
        self.time_edit = QLineEdit("02:00:00")
        resources_form.addRow("Time Limit:", self.time_edit)
        
        self.partition_edit = QLineEdit()
        self.partition_edit.setPlaceholderText("compute")
        resources_form.addRow("Partition/Queue:", self.partition_edit)
        
        self.resources_widget.setVisible(False)
        resources_layout.addWidget(self.resources_widget)
        
        scroll_layout.addWidget(resources_group)
        
        # Prepare Calculation Button
        prepare_group = QGroupBox("✨ Prepare Calculation")
        prepare_layout = QVBoxLayout(prepare_group)
        
        prepare_info = QLabel("""
<p>Click below to prepare atoms and Espresso calculator objects using the <b>calculation module</b>.</p>
""")
        prepare_info.setTextFormat(Qt.RichText)
        prepare_layout.addWidget(prepare_info)
        
        prepare_btn = QPushButton("🔧 Prepare Calculation")
        prepare_btn.clicked.connect(self._prepare_calculation)
        prepare_layout.addWidget(prepare_btn)
        
        scroll_layout.addWidget(prepare_group)
        
        # Results area
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        scroll_layout.addWidget(self.results_label)
        
        # Configuration summary
        config_group = QGroupBox("📝 Current Configuration")
        config_layout = QVBoxLayout(config_group)
        
        self.config_text = QTextEdit()
        self.config_text.setReadOnly(True)
        self.config_text.setMaximumHeight(150)
        config_layout.addWidget(self.config_text)
        
        scroll_layout.addWidget(config_group)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
        
        # Load initial data
        self._load_machines()
        self._update_structure_status()
        
        # Apply styling after all widgets are created
        self._apply_combobox_styling()
    
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
        except Exception as e:
            self.results_label.setText(f"⚠️ Could not load machines: {e}")
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
            self.session_state['calc_machine'] = machine
            self.session_state['selected_machine'] = machine_name
            
            info = f"Type: {machine.execution}"
            if machine.scheduler:
                info += f", Scheduler: {machine.scheduler}"
            self.machine_info_label.setText(info)
            
            # Load codes for this machine
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
                    # No versions, load codes directly
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
                
                # Try to select 'pw' by default
                idx = self.code_combo.findText('pw')
                if idx >= 0:
                    self.code_combo.setCurrentIndex(idx)
        except Exception:
            pass
        finally:
            self.code_combo.blockSignals(False)
    
    def _on_ecutwfc_changed(self, value):
        """Update ecutrho when ecutwfc changes."""
        self.ecutrho_spin.setValue(value * 8)
    
    def _on_occupations_changed(self, occupations):
        """Handle occupations change."""
        self.smearing_group.setVisible(occupations == "smearing")
    
    def _on_kpts_mode_changed(self, checked):
        """Handle k-points mode change."""
        self.kspacing_widget.setVisible(checked)
        self.explicit_widget.setVisible(not checked)
    
    def _on_adjust_resources_changed(self, state):
        """Handle adjust resources checkbox change."""
        self.resources_widget.setVisible(state == Qt.Checked)
    
    def _update_structure_status(self):
        """Update structure status display."""
        atoms = self.session_state.get('current_structure')
        
        if atoms is not None and ASE_AVAILABLE:
            if isinstance(atoms, Atoms):
                formula = atoms.get_chemical_formula()
                natoms = len(atoms)
                self.structure_status.setText(f"✅ Structure loaded: {formula} ({natoms} atoms)")
                self.structure_status.setStyleSheet("color: green;")
                
                # Update label default (structure/calculation_type format)
                calc_type = self.calc_type_combo.currentText()
                self.label_edit.setText(f"{formula}/{calc_type}")
                
                # Update pseudopotential inputs
                self._update_pseudo_inputs(atoms)
            else:
                self.structure_status.setText("⚠️ Invalid structure type")
                self.structure_status.setStyleSheet("color: orange;")
        else:
            self.structure_status.setText("⚠️ No structure loaded. Please load a structure first.")
            self.structure_status.setStyleSheet("color: orange;")
    
    def _update_pseudo_inputs(self, atoms):
        """Update pseudopotential input fields for structure elements."""
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
        """Update magnetic input fields for structure elements."""
        # Clear existing inputs
        while self.magnetic_container_layout.count():
            item = self.magnetic_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.magnetic_edits = {}
        
        for element in sorted(elements):
            spin = QDoubleSpinBox()
            spin.setRange(-10.0, 10.0)
            spin.setSingleStep(0.1)
            spin.setDecimals(2)
            default_val = PREDEFINED_MAGNETIC_MOMENTS.get(element, 0.0)
            spin.setValue(default_val)
            self.magnetic_edits[element] = spin
            self.magnetic_container_layout.addRow(f"{element}:", spin)
    
    def _update_hubbard_inputs(self, elements):
        """Update Hubbard U input fields for structure elements."""
        # Clear existing inputs
        while self.hubbard_container_layout.count():
            item = self.hubbard_container_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.hubbard_edits = {}
        self.hubbard_orbital_edits = {}
        self.hubbard_checkboxes = {}
        
        for element in sorted(elements):
            # Create a container widget for each element with checkbox, U value and orbital
            container = QWidget()
            hlayout = QHBoxLayout(container)
            hlayout.setContentsMargins(0, 0, 0, 0)
            
            # Checkbox to enable/disable Hubbard for this element
            checkbox = QCheckBox()
            checkbox.setChecked(False)  # Default to unchecked
            checkbox.setToolTip(f"Enable Hubbard U correction for {element}")
            hlayout.addWidget(checkbox)
            
            # U value spinbox
            spin = QDoubleSpinBox()
            spin.setRange(0.0, 20.0)
            spin.setSingleStep(0.5)
            spin.setDecimals(1)
            spin.setSuffix(" eV")
            default_val = TYPICAL_HUBBARD_U.get(element, 0.0)
            spin.setValue(default_val)
            hlayout.addWidget(spin)
            
            # Orbital selector (for new format)
            orbital_combo = QComboBox()
            # Common orbitals for the element
            orbitals = self._get_common_orbitals_for_element(element)
            orbital_combo.addItems(orbitals)
            # Apply styling to this combobox
            orbital_combo.setStyleSheet(self._get_combobox_stylesheet())
            hlayout.addWidget(QLabel("Orbital:"))
            hlayout.addWidget(orbital_combo)
            
            self.hubbard_checkboxes[element] = checkbox
            self.hubbard_edits[element] = spin
            self.hubbard_orbital_edits[element] = orbital_combo
            self.hubbard_container_layout.addRow(f"{element}:", container)
    
    def _get_common_orbitals_for_element(self, element):
        """Get common orbital choices for an element.
        
        Returns a list of orbital strings appropriate for the element
        for use in DFT+U calculations.
        """
        # Mapping of elements to their common Hubbard orbitals
        orbital_map = {
            # 3d transition metals
            'Ti': ['3d', '4s', '3p'],
            'V': ['3d', '4s', '3p'],
            'Cr': ['3d', '4s', '3p'],
            'Mn': ['3d', '4s', '3p'],
            'Fe': ['3d', '4s', '3p'],
            'Co': ['3d', '4s', '3p'],
            'Ni': ['3d', '4s', '3p'],
            'Cu': ['3d', '4s', '3p'],
            'Zn': ['3d', '4s', '3p'],
            # 4d transition metals
            'Zr': ['4d', '5s', '4p'],
            'Nb': ['4d', '5s', '4p'],
            'Mo': ['4d', '5s', '4p'],
            'Tc': ['4d', '5s', '4p'],
            'Ru': ['4d', '5s', '4p'],
            'Rh': ['4d', '5s', '4p'],
            'Pd': ['4d', '5s', '4p'],
            'Ag': ['4d', '5s', '4p'],
            'Cd': ['4d', '5s', '4p'],
            # 4f rare earths
            'Ce': ['4f', '5d', '5p'],
            'Pr': ['4f', '5d', '5p'],
            'Nd': ['4f', '5d', '5p'],
            'Pm': ['4f', '5d', '5p'],
            'Sm': ['4f', '5d', '5p'],
            'Eu': ['4f', '5d', '5p'],
            'Gd': ['4f', '5d', '5p'],
            'Tb': ['4f', '5d', '5p'],
            'Dy': ['4f', '5d', '5p'],
            'Ho': ['4f', '5d', '5p'],
            'Er': ['4f', '5d', '5p'],
            'Tm': ['4f', '5d', '5p'],
            'Yb': ['4f', '5d', '5p'],
            'Lu': ['4f', '5d', '5p'],
            # p-block elements
            'O': ['2p', '2s'],
            'S': ['3p', '3s'],
            'Se': ['4p', '4s'],
            'N': ['2p', '2s'],
            'P': ['3p', '3s'],
            'As': ['4p', '4s'],
        }
        
        return orbital_map.get(element, ['3d', '4d', '4f', '2p', '3p'])
    
    def _on_magnetic_preset_changed(self, preset):
        """Handle magnetic preset selection."""
        if preset == "Custom":
            return
        
        for element, spin in self.magnetic_edits.items():
            if element in PREDEFINED_MAGNETIC_MOMENTS:
                if preset == "Ferromagnetic":
                    spin.setValue(PREDEFINED_MAGNETIC_MOMENTS[element])
                elif preset == "Antiferromagnetic":
                    spin.setValue(-PREDEFINED_MAGNETIC_MOMENTS[element])
            else:
                spin.setValue(0.0)
    
    def _load_modules_from_codes(self, config):
        """Load modules from codes configuration.
        
        Helper method to load modules for a specific machine and QE version
        from the codes configuration JSON file. This is used to automatically
        include the correct module load commands in generated job scripts.
        
        Args:
            config (dict): Configuration dictionary that must contain:
                - machine_name (str): Name of the machine configuration
                - qe_version (str): Quantum ESPRESSO version string (e.g., '7.2')
                
        Returns:
            list or None: List of module names to load (e.g., ['quantum-espresso/7.2']),
                         or None if modules not found or config is incomplete
                         
        Note:
            This method gracefully handles missing xespresso, missing config files,
            or incomplete config dictionaries by returning None. Modules are optional
            for local execution but typically required for HPC environments.
        """
        try:
            if not XESPRESSO_AVAILABLE:
                return None
            
            machine_name = config.get('machine_name')
            qe_version = config.get('qe_version')
            
            if not machine_name or not qe_version:
                return None
            
            codes = load_codes_config(machine_name, DEFAULT_CODES_DIR, verbose=False)
            
            if codes and codes.versions and qe_version in codes.versions:
                version_config = codes.versions[qe_version]
                return version_config.get('modules')
        except Exception:
            # Silently fail - modules are optional
            pass
        
        return None
    
    def _get_config(self):
        """Get the current configuration as a dictionary."""
        config = {}
        
        config['calc_type'] = self.calc_type_combo.currentText()
        config['label'] = self.label_edit.text()
        config['ecutwfc'] = self.ecutwfc_spin.value()
        config['ecutrho'] = self.ecutrho_spin.value()
        config['occupations'] = self.occupations_combo.currentText()
        config['conv_thr'] = float(self.conv_thr_edit.text())
        
        if config['occupations'] == 'smearing':
            config['smearing'] = self.smearing_combo.currentText()
            config['degauss'] = self.degauss_spin.value()
        
        # K-points
        if self.kspacing_radio.isChecked():
            config['kspacing'] = self.kspacing_spin.value()
        else:
            config['kpts'] = (self.k1_spin.value(), self.k2_spin.value(), self.k3_spin.value())
        
        # Pseudopotentials
        if PSEUDO_SELECTOR_AVAILABLE and hasattr(self, 'pseudo_selector'):
            config['pseudopotentials'] = self.pseudo_selector.get_pseudopotentials()
        else:
            config['pseudopotentials'] = {}
            for element, edit in getattr(self, 'pseudo_edits', {}).items():
                pseudo = edit.text().strip()
                if pseudo:
                    config['pseudopotentials'][element] = pseudo
        
        # Machine and Queue configuration
        config['machine_name'] = self.machine_combo.currentText()
        config['qe_version'] = self.version_combo.currentText()
        config['selected_code'] = self.code_combo.currentText()
        
        # Get the machine object and convert to queue configuration
        # This ensures the scheduler settings (slurm/direct, local/remote) are passed to Espresso
        machine = self.session_state.get('calc_machine')
        if machine:
            try:
                # Machine's to_queue() method returns the complete queue configuration
                # including execution mode, scheduler type, and all settings
                config['queue'] = machine.to_queue()
            except (AttributeError, TypeError) as e:
                import logging
                logging.warning(f"Could not convert machine to queue: {e}")
                # Don't set a default queue - let prepare_calculation_from_gui handle it
        
        # Load modules from codes configuration
        modules = self._load_modules_from_codes(config)
        if modules:
            config['modules'] = modules
        
        # Resources
        if self.adjust_resources_check.isChecked():
            config['adjust_resources'] = True
            config['nprocs'] = self.nprocs_spin.value()
            config['resources'] = {
                'nodes': self.nodes_spin.value(),
                'ntasks-per-node': self.ntasks_spin.value(),
                'time': self.time_edit.text(),
                'partition': self.partition_edit.text()
            }
        
        # Magnetic configuration
        if self.magnetic_group.isChecked():
            config['enable_magnetism'] = True
            config['magnetic_config'] = {}
            for element, spin in self.magnetic_edits.items():
                value = spin.value()
                if value != 0.0:
                    config['magnetic_config'][element] = [value]
        else:
            config['enable_magnetism'] = False
            config['magnetic_config'] = {}
        
        # Hubbard configuration
        if self.hubbard_group.isChecked():
            config['enable_hubbard'] = True
            config['hubbard_format'] = 'new' if 'New' in self.hubbard_format_combo.currentText() else 'old'
            config['hubbard_u'] = {}
            config['hubbard_orbitals'] = {}
            
            for element, checkbox in self.hubbard_checkboxes.items():
                # Only save Hubbard U if checkbox is checked
                if checkbox.isChecked():
                    spin = self.hubbard_edits[element]
                    value = spin.value()
                    config['hubbard_u'][element] = value
                    # Save orbital information if available (for new format)
                    if element in self.hubbard_orbital_edits:
                        orbital = self.hubbard_orbital_edits[element].currentText()
                        config['hubbard_orbitals'][element] = orbital
        else:
            config['enable_hubbard'] = False
            config['hubbard_u'] = {}
        
        # Always use crystal coordinates for atomic positions
        config['crystal_coordinates'] = True
        
        return config
    
    def _prepare_calculation(self):
        """Prepare the calculation - creates Espresso calculator and stores in session state."""
        atoms = self.session_state.get('current_structure')
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded. Please load a structure first.")
            return
        
        # Check if workflow builder mode is active
        workflow_mode = self.session_state.get('workflow_mode')
        if workflow_mode == 'multi':
            reply = QMessageBox.question(
                self,
                "Switch to Single Calculation Mode?",
                "You previously built a multi-step workflow.\n\n"
                "Using Calculation Setup will switch to single calculation mode and replace the workflow configuration.\n\n"
                "Do you want to continue?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
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
                "Please select a machine configuration before preparing the calculation.\n\n"
                "A machine defines where and how the calculation will run (local/remote, scheduler type, resources, etc.).\n\n"
                "Go to Machine Configuration page to create or select a machine."
            )
            return
        
        # Set up environment for Espresso calculator
        # 1. Set ASE_ESPRESSO_COMMAND environment variable (like Streamlit GUI does)
        # This tells xespresso's Espresso calculator how to run QE executables
        import os
        os.environ['ASE_ESPRESSO_COMMAND'] = ASE_ESPRESSO_COMMAND_TEMPLATE
        
        # 2. Load modules from codes configuration if available
        modules = self._load_modules_from_codes(config)
        if modules:
            config['modules'] = modules
        
        # Store config in session state
        self.session_state['workflow_config'] = config
        # Mark that we're using single calculation mode
        self.session_state['workflow_mode'] = 'single'
        
        # Store config in session state for use in job_submission
        # Note: Calculation preparation is now handled in job_submission page
        import json
        self.config_text.setText(json.dumps(config, indent=2))
        
        self.results_label.setText("""
✅ <b>Configuration saved successfully!</b>

The configuration has been stored and is ready for use.

Go to <b>Job Submission</b> page to prepare and run the calculation.
""")
        self.results_label.setStyleSheet("color: green;")
        self.results_label.setTextFormat(Qt.RichText)
        
        QMessageBox.information(
            self, 
            "Configuration Saved",
            "Configuration has been saved.\n\nGo to Job Submission page to prepare and run the calculation."
        )
    
    def refresh(self):
        """Refresh the page."""
        # Use loading guard to prevent infinite loops
        if self._loading:
            return
        self._loading = True
        try:
            self._load_machines()
            self._update_structure_status()
            # Restore UI state from saved configuration if exists
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
                        # mag_value might be a list from setup_magnetic_config
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
                    if element in self.hubbard_edits:
                        self.hubbard_edits[element].setValue(u_value)
                    # Check the checkbox for elements that have Hubbard U configured
                    if element in self.hubbard_checkboxes:
                        self.hubbard_checkboxes[element].setChecked(True)
            # Restore orbital selections
            if config.get('hubbard_orbitals'):
                for element, orbital in config['hubbard_orbitals'].items():
                    if element in self.hubbard_orbital_edits:
                        idx = self.hubbard_orbital_edits[element].findText(orbital)
                        if idx >= 0:
                            self.hubbard_orbital_edits[element].setCurrentIndex(idx)
            # Restore Hubbard format
            if config.get('hubbard_format'):
                format_str = "New (QE >= 7.0)" if config['hubbard_format'] == 'new' else "Old (QE < 7.0)"
                idx = self.hubbard_format_combo.findText(format_str)
                if idx >= 0:
                    self.hubbard_format_combo.setCurrentIndex(idx)
        else:
            # Explicitly uncheck if the config says Hubbard is disabled
            self.hubbard_group.setChecked(False)
        
        # Restore other configuration values
        if config.get('calc_type'):
            idx = self.calc_type_combo.findText(config['calc_type'])
            if idx >= 0:
                self.calc_type_combo.setCurrentIndex(idx)
        
        if config.get('label'):
            self.label_edit.setText(config['label'])
        
        if 'ecutwfc' in config:
            self.ecutwfc_spin.setValue(config['ecutwfc'])
        
        if 'ecutrho' in config:
            self.ecutrho_spin.setValue(config['ecutrho'])
        
        if config.get('occupations'):
            idx = self.occupations_combo.findText(config['occupations'])
            if idx >= 0:
                self.occupations_combo.setCurrentIndex(idx)
        
        if 'conv_thr' in config:
            self.conv_thr_edit.setText(str(config['conv_thr']))
        
        # Restore smearing parameters
        if config.get('smearing'):
            idx = self.smearing_combo.findText(config['smearing'])
            if idx >= 0:
                self.smearing_combo.setCurrentIndex(idx)
        
        if 'degauss' in config:
            self.degauss_spin.setValue(config['degauss'])
        
        # Restore k-points
        if 'kspacing' in config:
            self.kspacing_radio.setChecked(True)
            self.kspacing_spin.setValue(config['kspacing'])
        elif 'kpts' in config:
            self.explicit_radio.setChecked(True)
            kpts = config['kpts']
            if isinstance(kpts, (list, tuple)) and len(kpts) == 3:
                self.k1_spin.setValue(kpts[0])
                self.k2_spin.setValue(kpts[1])
                self.k3_spin.setValue(kpts[2])
        
        # Restore resources configuration
        if config.get('adjust_resources'):
            self.adjust_resources_check.setChecked(True)
            if 'nprocs' in config:
                self.nprocs_spin.setValue(config['nprocs'])
            resources = config.get('resources', {})
            if 'nodes' in resources:
                self.nodes_spin.setValue(resources['nodes'])
            if 'ntasks-per-node' in resources:
                self.ntasks_spin.setValue(resources['ntasks-per-node'])
            if 'time' in resources:
                self.time_edit.setText(resources['time'])
            if 'partition' in resources:
                self.partition_edit.setText(resources['partition'])
        
        # Restore pseudopotentials
        if config.get('pseudopotentials'):
            if PSEUDO_SELECTOR_AVAILABLE and hasattr(self, 'pseudo_selector'):
                # Use the advanced selector widget
                self.pseudo_selector.set_pseudopotentials(config['pseudopotentials'])
            else:
                # Fallback to simple manual inputs
                for element, pseudo_file in config['pseudopotentials'].items():
                    if element in getattr(self, 'pseudo_edits', {}):
                        self.pseudo_edits[element].setText(pseudo_file)
    
    def _merge_configs(self, config, existing_config):
        """Merge current config with existing config and return the result.
        
        A configuration is considered valid if it has pseudopotentials defined,
        which indicates the user has completed the "Prepare Calculation" step.
        
        Args:
            config (dict): New configuration from current UI state
            existing_config (dict or None): Existing workflow configuration in session state
            
        Returns:
            dict or None: Merged configuration to save, or None to skip saving
            
        Decision logic:
            - If new config has pseudopotentials: Save it (it's a valid prepared config)
            - If no existing config: Save new config even if incomplete (preserve UI state)
            - If existing config has pseudopotentials but new config doesn't:
              Merge the new config (especially magnetic/hubbard settings) with existing pseudopotentials
            - Otherwise: Return None to skip saving (don't overwrite valid with incomplete)
        """
        # Always save if new config has pseudopotentials (indicates a prepared config)
        if config.get('pseudopotentials'):
            return config
        
        # If no existing config, save current state even if incomplete
        # to preserve other fields like calc_type, ecutwfc, kpts, etc.
        if not existing_config:
            return config
        
        # If existing config has pseudopotentials but new one doesn't,
        # merge important fields (magnetic/hubbard settings, basic params) into existing
        # This allows users to edit magnetic/hubbard settings after session reload
        # without having to re-prepare the calculation
        if existing_config.get('pseudopotentials') and not config.get('pseudopotentials'):
            # Create a new merged config by copying the new config
            merged_config = config.copy()
            # Preserve pseudopotentials from existing config
            merged_config['pseudopotentials'] = existing_config['pseudopotentials']
            # Return the merged config
            return merged_config
        
        # Otherwise, return None to skip saving (don't overwrite valid with incomplete)
        return None
    
    def save_state(self):
        """Save current page state to session state.
        
        This is called before the session is saved to disk to ensure
        all current UI values are captured in the session state.
        
        Only save if we're in single calculation mode or if mode isn't set yet.
        """
        workflow_mode = self.session_state.get('workflow_mode')
        
        # Only save if we're in single mode, or mode is not set (for backward compatibility)
        if workflow_mode == 'multi':
            # Don't overwrite workflow builder's config
            return
        
        config = self._get_config()
        existing_config = self.session_state.get('workflow_config')
        
        # Merge configs if needed and get the final config to save
        merged_config = self._merge_configs(config, existing_config)
        if merged_config is not None:
            self.session_state['workflow_config'] = merged_config
            # If mode isn't set, set it to single since we're saving from here
            if workflow_mode is None:
                self.session_state['workflow_mode'] = 'single'
    
    def _update_mode_indicator(self):
        """Update the mode status indicator to show which workflow mode is active."""
        workflow_mode = self.session_state.get('workflow_mode')
        
        if workflow_mode == 'multi':
            self.mode_status.setText(
                "ℹ️ <b>Multi-Step Workflow Mode Active</b> - "
                "Currently using Workflow Builder. To use this page, click 'Prepare Calculation' and confirm the mode switch."
            )
            self.mode_status.setStyleSheet("color: #ff9800; font-weight: bold; background-color: #fff3e0; padding: 10px; border-radius: 5px;")
            self.mode_status.setTextFormat(Qt.RichText)
        elif workflow_mode == 'single':
            self.mode_status.setText("✅ <b>Single Calculation Mode Active</b>")
            self.mode_status.setStyleSheet("color: green; font-weight: bold;")
            self.mode_status.setTextFormat(Qt.RichText)
        else:
            # No mode set yet
            self.mode_status.setText("")
            self.mode_status.setStyleSheet("")


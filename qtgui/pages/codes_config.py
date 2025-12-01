"""
Codes Configuration Page for xespresso PySide6 GUI.

This module handles the Quantum ESPRESSO codes configuration interface,
allowing users to:
- Auto-detect QE executables on machines
- Configure multiple versions
- Save and load code configurations
"""

import os
import traceback

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QTextEdit, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTableWidget, QTableWidgetItem,
    QHeaderView, QApplication
)
from PySide6.QtCore import Qt

try:
    from xespresso.machines.config.loader import (
        list_machines,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    from xespresso.codes.manager import (
        detect_qe_codes, load_codes_config, CodesManager,
        DEFAULT_CODES_DIR
    )
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False


class CodesConfigPage(QWidget):
    """Codes configuration page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self.detected_codes = None
        self._setup_ui()
        self._load_machines_list()
    
    def _setup_ui(self):
        """Setup the user interface."""
        main_layout = QVBoxLayout(self)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Header
        header_label = QLabel("<h2>⚙️ Quantum ESPRESSO Codes Configuration</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p><b>Configure and auto-detect</b> Quantum ESPRESSO executable paths.</p>
<p>This page is for <b>configuration only</b> - once codes are saved, you can select versions
in the Calculation Setup or Workflow Builder pages.</p>
<p>Auto-detection is supported for both local and remote systems.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        info_label = QLabel("""
<p>💡 <b>Configuration vs. Selection:</b></p>
<ul>
<li><b>Configure</b> codes here (auto-detect and save)</li>
<li><b>Select</b> code versions in Calculation Setup or Workflow Builder</li>
<li>Multiple QE versions can coexist for the same machine</li>
<li>Configurations are saved to <code>~/.xespresso/codes/</code></li>
</ul>
""")
        info_label.setTextFormat(Qt.RichText)
        info_label.setWordWrap(True)
        info_label.setStyleSheet("background-color: #e3f2fd; padding: 10px; border-radius: 5px;")
        scroll_layout.addWidget(info_label)
        
        if not XESPRESSO_AVAILABLE:
            error_label = QLabel("❌ xespresso modules not available. Cannot configure codes.")
            error_label.setStyleSheet("color: red; font-weight: bold;")
            scroll_layout.addWidget(error_label)
            scroll_area.setWidget(scroll_widget)
            main_layout.addWidget(scroll_area)
            return
        
        # Machine Selection
        machine_group = QGroupBox("Select Machine")
        machine_layout = QHBoxLayout(machine_group)
        
        machine_layout.addWidget(QLabel("Machine:"))
        self.machine_combo = QComboBox()
        self.machine_combo.currentTextChanged.connect(self._on_machine_changed)
        machine_layout.addWidget(self.machine_combo, 1)
        
        scroll_layout.addWidget(machine_group)
        
        # Auto-Detection Section
        detect_group = QGroupBox("Auto-Detect Codes")
        detect_layout = QFormLayout(detect_group)
        
        self.qe_prefix_edit = QLineEdit()
        self.qe_prefix_edit.setPlaceholderText("e.g., /opt/qe-7.2/bin")
        detect_layout.addRow("QE Installation Prefix (optional):", self.qe_prefix_edit)
        
        self.qe_version_edit = QLineEdit()
        self.qe_version_edit.setPlaceholderText("e.g., 7.2, 7.1, 6.8")
        detect_layout.addRow("QE Version (optional but recommended):", self.qe_version_edit)
        
        self.label_edit = QLineEdit()
        self.label_edit.setPlaceholderText("e.g., 'production', 'dev', 'test'")
        detect_layout.addRow("Label (optional):", self.label_edit)
        
        self.modules_edit = QTextEdit()
        self.modules_edit.setPlaceholderText("Version-specific modules (e.g., 'qe/7.2' or 'quantum_espresso-7.4.1')")
        self.modules_edit.setMaximumHeight(80)
        detect_layout.addRow("Modules to Load (optional):", self.modules_edit)
        
        self.search_paths_edit = QTextEdit()
        self.search_paths_edit.setPlaceholderText("Additional directories to search for executables")
        self.search_paths_edit.setMaximumHeight(80)
        detect_layout.addRow("Additional Search Paths (optional):", self.search_paths_edit)
        
        tip_label = QLabel("""
<p>💡 <b>Tip: Explicit Version</b></p>
<p>Auto-detection may pick up compiler versions. It's recommended to specify the QE version explicitly!</p>
""")
        tip_label.setTextFormat(Qt.RichText)
        tip_label.setStyleSheet("background-color: #fff3e0; padding: 10px; border-radius: 5px;")
        detect_layout.addRow(tip_label)
        
        detect_btn = QPushButton("🔍 Auto-Detect Codes")
        detect_btn.clicked.connect(self._detect_codes)
        detect_layout.addRow(detect_btn)
        
        scroll_layout.addWidget(detect_group)
        
        # Detected Codes Section
        self.detected_group = QGroupBox("Detected Codes")
        detected_layout = QVBoxLayout(self.detected_group)
        
        self.codes_table = QTableWidget()
        self.codes_table.setColumnCount(4)
        self.codes_table.setHorizontalHeaderLabels(["Code", "Path", "Version", "Label"])
        self.codes_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        detected_layout.addWidget(self.codes_table)
        
        save_info = QLabel("""
<p>💾 <b>Saving Codes:</b></p>
<ul>
<li>Configuration will be saved to <code>{machine}.json</code></li>
<li>If <b>version/label</b> is specified, codes are stored in the <code>versions</code> structure</li>
<li>Multiple QE versions with different labels can coexist</li>
</ul>
""")
        save_info.setTextFormat(Qt.RichText)
        save_info.setWordWrap(True)
        detected_layout.addWidget(save_info)
        
        save_btn = QPushButton("💾 Save Codes Configuration")
        save_btn.clicked.connect(self._save_codes)
        detected_layout.addWidget(save_btn)
        
        self.detected_group.setVisible(False)
        scroll_layout.addWidget(self.detected_group)
        
        # Existing Configuration Section
        existing_group = QGroupBox("Existing Codes Configuration")
        existing_layout = QVBoxLayout(existing_group)
        
        self.existing_status_label = QLabel("")
        existing_layout.addWidget(self.existing_status_label)
        
        version_select_layout = QHBoxLayout()
        version_select_layout.addWidget(QLabel("Select QE Version:"))
        self.version_combo = QComboBox()
        self.version_combo.currentTextChanged.connect(self._on_version_changed)
        version_select_layout.addWidget(self.version_combo, 1)
        
        load_version_btn = QPushButton("Load Version")
        load_version_btn.clicked.connect(self._load_version)
        version_select_layout.addWidget(load_version_btn)
        existing_layout.addLayout(version_select_layout)
        
        self.existing_codes_table = QTableWidget()
        self.existing_codes_table.setColumnCount(3)
        self.existing_codes_table.setHorizontalHeaderLabels(["Code", "Path", "Version"])
        self.existing_codes_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        existing_layout.addWidget(self.existing_codes_table)
        
        scroll_layout.addWidget(existing_group)
        
        # Results area
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        scroll_layout.addWidget(self.results_label)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
    
    def _load_machines_list(self):
        """Load the list of available machines."""
        if not XESPRESSO_AVAILABLE:
            return
        
        try:
            machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.blockSignals(True)
            self.machine_combo.clear()
            
            if machines_list:
                for machine in machines_list:
                    self.machine_combo.addItem(machine)
            else:
                self.results_label.setText("⚠️ No machines configured. Please configure a machine first.")
                self.results_label.setStyleSheet("color: orange;")
            self.machine_combo.blockSignals(False)
            # Manually trigger the handler for the first item if any
            if self.machine_combo.count() > 0:
                self._on_machine_changed(self.machine_combo.currentText())
        except Exception as e:
            self.machine_combo.blockSignals(False)
            self.results_label.setText(f"⚠️ Could not load machines: {e}")
            self.results_label.setStyleSheet("color: orange;")
    
    def _on_machine_changed(self, machine_name):
        """Handle machine selection change."""
        if machine_name:
            self._load_existing_codes(machine_name)
    
    def _detect_codes(self):
        """Auto-detect Quantum ESPRESSO codes."""
        if not XESPRESSO_AVAILABLE:
            return
        
        machine_name = self.machine_combo.currentText()
        if not machine_name:
            QMessageBox.warning(self, "Warning", "Please select a machine first")
            return
        
        try:
            self.results_label.setText("Detecting Quantum ESPRESSO codes...")
            self.results_label.setStyleSheet("color: blue;")
            
            # Process events to update UI
            QApplication.processEvents()
            
            modules_text = self.modules_edit.toPlainText()
            modules = [m.strip() for m in modules_text.split("\n") if m.strip()] if modules_text else None
            
            search_paths_text = self.search_paths_edit.toPlainText()
            search_paths = [p.strip() for p in search_paths_text.split("\n") if p.strip()] if search_paths_text else None
            
            qe_prefix = self.qe_prefix_edit.text().strip() or None
            qe_version = self.qe_version_edit.text().strip() or None
            label = self.label_edit.text().strip() or None
            
            codes_config = detect_qe_codes(
                machine_name=machine_name,
                qe_prefix=qe_prefix,
                search_paths=search_paths,
                modules=modules,
                auto_load_machine=True,
                qe_version=qe_version,
                label=label
            )
            
            if codes_config and codes_config.has_any_codes():
                all_codes = codes_config.get_all_codes()
                self.results_label.setText(f"✅ Detected {len(all_codes)} codes!")
                self.results_label.setStyleSheet("color: green;")
                
                self.detected_codes = codes_config
                self._display_detected_codes(all_codes, codes_config.label)
            else:
                self.results_label.setText("⚠️ No codes detected. Check paths and modules.")
                self.results_label.setStyleSheet("color: orange;")
                self.detected_codes = None
                self.detected_group.setVisible(False)
                
        except Exception as e:
            self.results_label.setText(f"❌ Error detecting codes: {e}")
            self.results_label.setStyleSheet("color: red;")
            self.detected_codes = None
    
    def _display_detected_codes(self, codes, label):
        """Display detected codes in table."""
        self.codes_table.setRowCount(len(codes))
        
        for i, (name, code) in enumerate(codes.items()):
            self.codes_table.setItem(i, 0, QTableWidgetItem(name))
            self.codes_table.setItem(i, 1, QTableWidgetItem(code.path))
            self.codes_table.setItem(i, 2, QTableWidgetItem(code.version or "Unknown"))
            self.codes_table.setItem(i, 3, QTableWidgetItem(label or "default"))
        
        self.detected_group.setVisible(True)
    
    def _save_codes(self):
        """Save the detected codes configuration."""
        if not XESPRESSO_AVAILABLE or not self.detected_codes:
            QMessageBox.warning(self, "Warning", "No codes to save. Please detect codes first.")
            return
        
        try:
            filepath = CodesManager.save_config(
                self.detected_codes,
                output_dir=DEFAULT_CODES_DIR,
                overwrite=False,
                merge=True,
                interactive=False
            )
            
            self.results_label.setText(f"✅ Codes saved to: {filepath}")
            self.results_label.setStyleSheet("color: green;")
            
            if self.detected_codes.qe_version:
                version_msg = f"📦 Added/updated QE version: {self.detected_codes.qe_version}"
                if self.detected_codes.label:
                    version_msg += f" with label {self.detected_codes.label}"
                self.results_label.setText(self.results_label.text() + f"\n{version_msg}")
            
            # Clear detected codes and refresh existing
            self.detected_codes = None
            self.detected_group.setVisible(False)
            self._load_existing_codes(self.machine_combo.currentText())
            
            self.session_state['current_codes'] = self.detected_codes
            
        except Exception as e:
            self.results_label.setText(f"❌ Error saving codes: {e}")
            self.results_label.setStyleSheet("color: red;")
    
    def _load_existing_codes(self, machine_name):
        """Load existing codes configuration for a machine."""
        if not XESPRESSO_AVAILABLE or not machine_name:
            return
        
        try:
            existing_codes = load_codes_config(machine_name, DEFAULT_CODES_DIR)
            
            if existing_codes:
                self.existing_status_label.setText(f"✅ Loaded existing configuration for '{machine_name}'")
                self.existing_status_label.setStyleSheet("color: green;")
                
                # Check for versions
                if existing_codes.versions:
                    available_versions = existing_codes.list_versions()
                    self.version_combo.blockSignals(True)
                    self.version_combo.clear()
                    for version in available_versions:
                        self.version_combo.addItem(version)
                    self.version_combo.blockSignals(False)
                else:
                    # No version structure, show all codes
                    self.version_combo.blockSignals(True)
                    self.version_combo.clear()
                    self.version_combo.blockSignals(False)
                    self._display_existing_codes(existing_codes.codes)
                    self.session_state['current_codes'] = existing_codes
            else:
                self.existing_status_label.setText(f"ℹ️ No codes configuration found for this machine.")
                self.existing_status_label.setStyleSheet("color: blue;")
                self.version_combo.blockSignals(True)
                self.version_combo.clear()
                self.version_combo.blockSignals(False)
                self.existing_codes_table.setRowCount(0)
                
        except Exception as e:
            self.existing_status_label.setText(f"⚠️ Could not load codes configuration: {e}")
            self.existing_status_label.setStyleSheet("color: orange;")
    
    def _on_version_changed(self, version):
        """Handle version selection change."""
        pass  # Loading is done via button
    
    def _load_version(self):
        """Load codes for selected version."""
        machine_name = self.machine_combo.currentText()
        selected_version = self.version_combo.currentText()
        
        if not machine_name or not selected_version:
            return
        
        try:
            version_config = load_codes_config(
                machine_name,
                DEFAULT_CODES_DIR,
                version=selected_version
            )
            
            if version_config:
                version_codes = version_config.get_all_codes(version=selected_version)
                self._display_existing_codes(version_codes)
                
                self.session_state['current_codes'] = version_config
                self.session_state['selected_qe_version'] = selected_version
                
                self.results_label.setText(f"✅ Loaded QE {selected_version} configuration!")
                self.results_label.setStyleSheet("color: green;")
            else:
                self.results_label.setText(f"⚠️ Could not load QE {selected_version} configuration.")
                self.results_label.setStyleSheet("color: orange;")
                
        except Exception as e:
            self.results_label.setText(f"❌ Error loading version: {e}")
            self.results_label.setStyleSheet("color: red;")
    
    def _display_existing_codes(self, codes):
        """Display existing codes in table."""
        self.existing_codes_table.setRowCount(len(codes))
        
        for i, (name, code) in enumerate(codes.items()):
            self.existing_codes_table.setItem(i, 0, QTableWidgetItem(name))
            self.existing_codes_table.setItem(i, 1, QTableWidgetItem(code.path))
            self.existing_codes_table.setItem(i, 2, QTableWidgetItem(code.version or "Unknown"))
    
    def refresh(self):
        """Refresh the page."""
        self._load_machines_list()

"""
Pseudopotentials Configuration Page for xespresso PySide6 GUI.

This module handles the pseudopotential configuration interface,
allowing users to:
- Auto-detect pseudopotential files on machines
- Configure multiple pseudopotential libraries
- Save and load pseudopotential configurations
"""

import os

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTableWidget, QTableWidgetItem,
    QHeaderView, QFileDialog, QTextEdit, QCheckBox, QListWidget
)
from qtpy.QtCore import Qt

try:
    from xespresso.machines.config.loader import (
        list_machines, load_machine,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    from xespresso.pseudopotentials.manager import (
        create_pseudopotentials_config,
        load_pseudopotentials_config,
        PseudopotentialsManager,
        DEFAULT_PSEUDOPOTENTIALS_DIR
    )
    from xespresso.machines.machine import Machine
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False
    DEFAULT_PSEUDOPOTENTIALS_DIR = os.path.expanduser("~/.xespresso/pseudopotentials")


class PseudopotentialsConfigPage(QWidget):
    """Pseudopotentials configuration page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self.detected_pseudos = None
        self.current_pseudos = None
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
        header_label = QLabel("<h2>⚙️ Pseudopotentials Configuration</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p><b>Configure and auto-detect</b> pseudopotential files for Quantum ESPRESSO calculations.</p>
<p>Pseudopotentials are primarily stored on your <b>local computer</b> after downloading them.
During calculations, xespresso automatically copies the needed pseudopotentials to remote machines.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        # Info box
        info_label = QLabel("""
<p>💡 <b>How Pseudopotentials Work:</b></p>
<ul>
<li><b>Download</b> pseudopotentials to your local computer (e.g., from SSSP, PSLibrary, Pseudo Dojo)</li>
<li><b>Configure</b> them here by pointing to the directory containing .UPF files</li>
<li><b>Machine selection is optional</b> - only needed if pseudopotentials are stored on a remote machine</li>
<li>During calculations, xespresso automatically copies pseudopotentials to the remote machine</li>
<li>Configurations are saved to <code>~/.xespresso/pseudopotentials/</code></li>
</ul>
""")
        info_label.setTextFormat(Qt.RichText)
        info_label.setWordWrap(True)
        info_label.setOpenExternalLinks(True)
        info_label.setStyleSheet("background-color: #e7f3ff; color: #1a1a1a; padding: 10px; border-radius: 5px;")
        scroll_layout.addWidget(info_label)
        
        if not XESPRESSO_AVAILABLE:
            error_label = QLabel("❌ xespresso modules not available. Cannot configure pseudopotentials.")
            error_label.setStyleSheet("color: red; font-weight: bold;")
            scroll_layout.addWidget(error_label)
            scroll_area.setWidget(scroll_widget)
            main_layout.addWidget(scroll_area)
            return
        
        # Remote machine option
        remote_group = QGroupBox("Machine Selection (Optional)")
        remote_layout = QVBoxLayout(remote_group)
        
        self.use_remote_check = QCheckBox("Pseudopotentials are stored on a remote machine")
        self.use_remote_check.setToolTip("Check this only if pseudopotentials are on a remote machine, not your local computer")
        self.use_remote_check.stateChanged.connect(self._on_remote_changed)
        remote_layout.addWidget(self.use_remote_check)
        
        self.machine_row = QWidget()
        machine_layout = QHBoxLayout(self.machine_row)
        machine_layout.setContentsMargins(0, 0, 0, 0)
        machine_layout.addWidget(QLabel("Remote Machine:"))
        self.machine_combo = QComboBox()
        self.machine_combo.setToolTip("Choose the remote machine where pseudopotentials are stored")
        machine_layout.addWidget(self.machine_combo, 1)
        self.machine_row.setVisible(False)
        remote_layout.addWidget(self.machine_row)
        
        scroll_layout.addWidget(remote_group)
        
        # Auto-detection section
        detect_group = QGroupBox("🔍 Auto-Detect Pseudopotentials")
        detect_layout = QFormLayout(detect_group)
        
        self.config_name_edit = QLineEdit()
        self.config_name_edit.setPlaceholderText("e.g., SSSP_efficiency, pbe_standard, my_pseudos")
        detect_layout.addRow("Configuration Name *:", self.config_name_edit)
        
        base_path_layout = QHBoxLayout()
        self.base_path_edit = QLineEdit()
        self.base_path_edit.setPlaceholderText("e.g., /home/user/pseudopotentials/SSSP, ~/pseudo/pbe")
        base_path_layout.addWidget(self.base_path_edit)
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_pseudo_dir)
        base_path_layout.addWidget(browse_btn)
        detect_layout.addRow("Pseudopotentials Directory *:", base_path_layout)
        
        self.functional_edit = QLineEdit()
        self.functional_edit.setPlaceholderText("e.g., PBE, LDA, PBEsol")
        detect_layout.addRow("Functional (optional):", self.functional_edit)
        
        self.library_edit = QLineEdit()
        self.library_edit.setPlaceholderText("e.g., SSSP, PSLibrary, Pseudo Dojo")
        detect_layout.addRow("Library Name (optional):", self.library_edit)
        
        self.version_edit = QLineEdit()
        self.version_edit.setPlaceholderText("e.g., 1.1.2, 1.0.0, 0.4")
        detect_layout.addRow("Library Version (optional):", self.version_edit)
        
        self.description_edit = QTextEdit()
        self.description_edit.setPlaceholderText("e.g., SSSP 1.1.2 PBE efficiency set for production calculations")
        self.description_edit.setMaximumHeight(60)
        detect_layout.addRow("Description (optional):", self.description_edit)
        
        self.recursive_check = QCheckBox("Search subdirectories recursively")
        self.recursive_check.setChecked(True)
        self.recursive_check.setToolTip("Search for .UPF files in subdirectories")
        detect_layout.addRow(self.recursive_check)
        
        tips_label = QLabel("""
<p><b>💡 Tips:</b></p>
<ul>
<li>The directory should contain .UPF or .upf pseudopotential files</li>
<li>Files are typically named like <code>Fe.pbe-spn-kjpaw_psl.0.2.1.UPF</code></li>
<li>Auto-detection will extract element symbols from filenames</li>
<li><b>Local pseudopotentials</b> are recommended (stored on your computer)</li>
</ul>
""")
        tips_label.setTextFormat(Qt.RichText)
        tips_label.setWordWrap(True)
        detect_layout.addRow(tips_label)
        
        detect_btn = QPushButton("🔍 Auto-Detect Pseudopotentials")
        detect_btn.clicked.connect(self._detect_pseudopotentials)
        detect_layout.addRow(detect_btn)
        
        scroll_layout.addWidget(detect_group)
        
        # Detection results group
        self.results_group = QGroupBox("Detected Pseudopotentials")
        results_layout = QVBoxLayout(self.results_group)
        
        self.detection_status = QLabel("")
        results_layout.addWidget(self.detection_status)
        
        self.detection_info = QLabel("")
        self.detection_info.setWordWrap(True)
        results_layout.addWidget(self.detection_info)
        
        self.detected_table = QTableWidget()
        self.detected_table.setColumnCount(4)
        self.detected_table.setHorizontalHeaderLabels(["Element", "Filename", "Type", "Functional"])
        self.detected_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.detected_table.setMaximumHeight(200)
        results_layout.addWidget(self.detected_table)
        
        save_info = QLabel("""
<p><b>💾 Saving Pseudopotentials:</b></p>
<ul>
<li>Configuration will be saved to <code>~/.xespresso/pseudopotentials/{name}.json</code></li>
<li>Pseudopotential files remain in their original location (local or remote)</li>
<li>During calculations, xespresso automatically copies needed files to remote machines</li>
</ul>
""")
        save_info.setTextFormat(Qt.RichText)
        save_info.setWordWrap(True)
        results_layout.addWidget(save_info)
        
        save_btn = QPushButton("💾 Save Pseudopotentials Configuration")
        save_btn.clicked.connect(self._save_pseudopotentials)
        results_layout.addWidget(save_btn)
        
        self.save_status = QLabel("")
        self.save_status.setWordWrap(True)
        results_layout.addWidget(self.save_status)
        
        self.results_group.setVisible(False)
        scroll_layout.addWidget(self.results_group)
        
        # Existing configurations section
        existing_group = QGroupBox("📋 Existing Pseudopotential Configurations")
        existing_layout = QVBoxLayout(existing_group)
        
        self.existing_status = QLabel("")
        existing_layout.addWidget(self.existing_status)
        
        config_select_layout = QHBoxLayout()
        config_select_layout.addWidget(QLabel("Select Configuration:"))
        self.config_combo = QComboBox()
        self.config_combo.setToolTip("Choose a pseudopotential configuration to view or use")
        config_select_layout.addWidget(self.config_combo, 1)
        
        load_btn = QPushButton("Load Configuration")
        load_btn.clicked.connect(self._load_configuration)
        config_select_layout.addWidget(load_btn)
        existing_layout.addLayout(config_select_layout)
        
        # Loaded configuration details
        self.loaded_group = QWidget()
        loaded_layout = QVBoxLayout(self.loaded_group)
        loaded_layout.setContentsMargins(0, 0, 0, 0)
        
        self.loaded_info = QLabel("")
        self.loaded_info.setWordWrap(True)
        self.loaded_info.setTextFormat(Qt.RichText)
        loaded_layout.addWidget(self.loaded_info)
        
        self.loaded_table = QTableWidget()
        self.loaded_table.setColumnCount(4)
        self.loaded_table.setHorizontalHeaderLabels(["Element", "Filename", "Type", "Path"])
        self.loaded_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.loaded_table.setMaximumHeight(200)
        loaded_layout.addWidget(self.loaded_table)
        
        # Action buttons
        action_layout = QHBoxLayout()
        
        self.set_default_btn = QPushButton("⭐ Set as Default")
        self.set_default_btn.setToolTip("Set this configuration as the default for calculations")
        self.set_default_btn.clicked.connect(self._set_default)
        action_layout.addWidget(self.set_default_btn)
        
        self.clear_default_btn = QPushButton("🚫 Clear Default")
        self.clear_default_btn.setToolTip("Remove the current default configuration")
        self.clear_default_btn.clicked.connect(self._clear_default)
        action_layout.addWidget(self.clear_default_btn)
        
        self.delete_btn = QPushButton("🗑️ Delete Configuration")
        self.delete_btn.clicked.connect(self._delete_configuration)
        action_layout.addWidget(self.delete_btn)
        
        loaded_layout.addLayout(action_layout)
        
        self.action_status = QLabel("")
        self.action_status.setWordWrap(True)
        loaded_layout.addWidget(self.action_status)
        
        self.loaded_group.setVisible(False)
        existing_layout.addWidget(self.loaded_group)
        
        scroll_layout.addWidget(existing_group)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
        
        # Load initial data
        self._load_machines()
        self._load_existing_configs()
    
    def _on_remote_changed(self, state):
        """Handle remote checkbox state change."""
        self.machine_row.setVisible(state == Qt.Checked)
    
    def _load_machines(self):
        """Load available machines."""
        if not XESPRESSO_AVAILABLE:
            return
        
        try:
            machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            if machines_list:
                for machine in machines_list:
                    self.machine_combo.addItem(machine)
        except Exception:
            # Machine loading is optional - user can still configure local pseudopotentials
            self.machine_combo.clear()
    
    def _load_existing_configs(self):
        """Load existing pseudopotential configurations."""
        if not XESPRESSO_AVAILABLE:
            self.existing_status.setText("⚠️ xespresso modules not available")
            self.existing_status.setStyleSheet("color: #d97706;")
            return
        
        try:
            existing_configs = PseudopotentialsManager.list_configs(DEFAULT_PSEUDOPOTENTIALS_DIR)
            
            self.config_combo.clear()
            if existing_configs:
                self.existing_status.setText(f"✅ Found {len(existing_configs)} saved configuration(s)")
                self.existing_status.setStyleSheet("color: green;")
                for config in existing_configs:
                    self.config_combo.addItem(config)
            else:
                self.existing_status.setText("ℹ️ No pseudopotential configurations saved yet. Create one above!")
                self.existing_status.setStyleSheet("color: blue;")
        except Exception as e:
            self.existing_status.setText(f"⚠️ Could not load configurations: {e}")
            self.existing_status.setStyleSheet("color: #d97706;")
    
    def _browse_pseudo_dir(self):
        """Browse for pseudopotential directory."""
        current_dir = self.base_path_edit.text() or os.path.expanduser("~")
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Pseudopotentials Directory",
            current_dir,
            QFileDialog.ShowDirsOnly
        )
        if directory:
            self.base_path_edit.setText(directory)
    
    def _detect_pseudopotentials(self):
        """Auto-detect pseudopotentials."""
        config_name = self.config_name_edit.text().strip()
        base_path = self.base_path_edit.text().strip()
        
        if not config_name or not base_path:
            QMessageBox.warning(self, "Warning", "Configuration name and base path are required!")
            return
        
        try:
            # Determine if using remote machine
            ssh_conn = None
            selected_machine = None
            
            if self.use_remote_check.isChecked() and self.machine_combo.currentText():
                selected_machine = self.machine_combo.currentText()
                try:
                    machine_obj = load_machine(
                        DEFAULT_CONFIG_PATH,
                        selected_machine,
                        DEFAULT_MACHINES_DIR,
                        return_object=True
                    )
                    
                    if machine_obj and isinstance(machine_obj, Machine) and machine_obj.is_remote:
                        ssh_conn = {
                            'host': machine_obj.host,
                            'username': machine_obj.username,
                            'port': getattr(machine_obj, 'port', 22)
                        }
                except Exception as e:
                    self.detection_status.setText(f"⚠️ Could not load machine SSH details: {e}")
                    self.detection_status.setStyleSheet("color: #d97706;")
            
            # Get optional fields
            functional = self.functional_edit.text().strip() or None
            library = self.library_edit.text().strip() or None
            version = self.version_edit.text().strip() or None
            description = self.description_edit.toPlainText().strip() or None
            recursive = self.recursive_check.isChecked()
            
            # Detect pseudopotentials
            pseudo_config = create_pseudopotentials_config(
                name=config_name,
                base_path=base_path,
                machine_name=selected_machine if self.use_remote_check.isChecked() else None,
                ssh_connection=ssh_conn,
                recursive=recursive,
                description=description,
                functional=functional,
                library=library,
                version=version,
                save=False  # Don't save yet, let user confirm
            )
            
            if pseudo_config and pseudo_config.pseudopotentials:
                self.detected_pseudos = pseudo_config
                self._show_detected_pseudos()
            else:
                self.detected_pseudos = None
                self.results_group.setVisible(False)
                QMessageBox.warning(self, "Warning", "No pseudopotentials detected. Check the path and try again.")
                
        except Exception as e:
            self.detection_status.setText(f"⚠️ Error detecting pseudopotentials: {e}")
            self.detection_status.setStyleSheet("color: red;")
            self.results_group.setVisible(True)
            
            # Show user-friendly error message (avoid exposing full traceback)
            QMessageBox.critical(self, "Error", f"Error detecting pseudopotentials:\n{e}")
    
    def _show_detected_pseudos(self):
        """Display detected pseudopotentials."""
        if not self.detected_pseudos:
            return
        
        pseudo_config = self.detected_pseudos
        
        self.results_group.setVisible(True)
        self.detection_status.setText(f"✅ Detected {len(pseudo_config.pseudopotentials)} pseudopotentials!")
        self.detection_status.setStyleSheet("color: green;")
        
        # Build info text
        info_parts = [f"<b>Configuration:</b> {pseudo_config.name}"]
        if pseudo_config.machine_name:
            info_parts.append(f"<b>Remote Machine:</b> {pseudo_config.machine_name}")
        else:
            info_parts.append("<b>Location:</b> Local (your computer)")
        info_parts.append(f"<b>Base Path:</b> {pseudo_config.base_path}")
        if pseudo_config.functional:
            info_parts.append(f"<b>Functional:</b> {pseudo_config.functional}")
        if pseudo_config.library:
            lib_str = f"<b>Library:</b> {pseudo_config.library}"
            if pseudo_config.version:
                lib_str += f" v{pseudo_config.version}"
            info_parts.append(lib_str)
        if pseudo_config.description:
            info_parts.append(f"<b>Description:</b> {pseudo_config.description}")
        
        self.detection_info.setText("<br>".join(info_parts))
        
        # Populate table
        elements = sorted(pseudo_config.pseudopotentials.keys())
        self.detected_table.setRowCount(len(elements))
        
        for i, element in enumerate(elements):
            pseudo = pseudo_config.get_pseudopotential(element)
            self.detected_table.setItem(i, 0, QTableWidgetItem(element))
            self.detected_table.setItem(i, 1, QTableWidgetItem(pseudo.filename))
            self.detected_table.setItem(i, 2, QTableWidgetItem(pseudo.type or "Unknown"))
            self.detected_table.setItem(i, 3, QTableWidgetItem(pseudo.functional or "Unknown"))
    
    def _save_pseudopotentials(self):
        """Save the detected pseudopotentials configuration."""
        if not self.detected_pseudos:
            QMessageBox.warning(self, "Warning", "No pseudopotentials to save. Run detection first.")
            return
        
        try:
            filepath = PseudopotentialsManager.save_config(
                self.detected_pseudos,
                output_dir=DEFAULT_PSEUDOPOTENTIALS_DIR,
                overwrite=True  # Allow overwrite since user explicitly clicked save
            )
            
            location = f"on remote machine '{self.detected_pseudos.machine_name}'" if self.detected_pseudos.machine_name else "locally"
            
            self.save_status.setText(f"✅ Pseudopotentials saved to: {filepath}")
            self.save_status.setStyleSheet("color: green;")
            
            QMessageBox.information(
                self, "Success",
                f"Pseudopotentials configuration saved!\n\n"
                f"Configuration: {self.detected_pseudos.name}\n"
                f"Location: {location}\n"
                f"Elements: {len(self.detected_pseudos.pseudopotentials)}\n"
                f"Saved to: {filepath}"
            )
            
            # Clear detected and reload configs
            self.detected_pseudos = None
            self.results_group.setVisible(False)
            self._load_existing_configs()
            
        except Exception as e:
            self.save_status.setText(f"❌ Error saving: {e}")
            self.save_status.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error saving pseudopotentials:\n{e}")
    
    def _load_configuration(self):
        """Load a selected configuration."""
        selected_config = self.config_combo.currentText()
        if not selected_config:
            QMessageBox.warning(self, "Warning", "Please select a configuration")
            return
        
        try:
            loaded_config = load_pseudopotentials_config(
                selected_config,
                DEFAULT_PSEUDOPOTENTIALS_DIR
            )
            
            if loaded_config:
                self.current_pseudos = loaded_config
                self._show_loaded_config(selected_config)
            else:
                QMessageBox.warning(self, "Warning", f"Failed to load configuration: {selected_config}")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading configuration:\n{e}")
    
    def _show_loaded_config(self, config_name):
        """Display a loaded configuration."""
        if not self.current_pseudos:
            return
        
        config = self.current_pseudos
        self.loaded_group.setVisible(True)
        
        # Build info text
        info_parts = [f"<h3>Configuration: {config.name}</h3>"]
        if config.machine_name:
            info_parts.append(f"<b>Remote Machine:</b> {config.machine_name}")
        else:
            info_parts.append("<b>Location:</b> Local (your computer)")
        info_parts.append(f"<b>Base Path:</b> {config.base_path}")
        if config.functional:
            info_parts.append(f"<b>Functional:</b> {config.functional}")
        if config.library:
            lib_str = f"<b>Library:</b> {config.library}"
            if config.version:
                lib_str += f" v{config.version}"
            info_parts.append(lib_str)
        info_parts.append(f"<b>Elements:</b> {len(config.pseudopotentials)}")
        if config.description:
            info_parts.append(f"<b>Description:</b> {config.description}")
        
        # Check if this is the default
        is_default = config_name == "default"
        has_default = PseudopotentialsManager.has_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
        
        if is_default:
            info_parts.append("<p><b>⭐ This is the default configuration</b></p>")
        
        self.loaded_info.setText("<br>".join(info_parts))
        
        # Populate table
        elements = config.list_elements()
        self.loaded_table.setRowCount(len(elements))
        
        for i, element in enumerate(elements):
            pseudo = config.get_pseudopotential(element)
            self.loaded_table.setItem(i, 0, QTableWidgetItem(element))
            self.loaded_table.setItem(i, 1, QTableWidgetItem(pseudo.filename))
            self.loaded_table.setItem(i, 2, QTableWidgetItem(pseudo.type or "-"))
            self.loaded_table.setItem(i, 3, QTableWidgetItem(pseudo.path))
        
        # Update button visibility
        self.set_default_btn.setVisible(not is_default)
        self.clear_default_btn.setVisible(has_default and not is_default)
        self.delete_btn.setVisible(not is_default)
        
        if is_default:
            self.action_status.setText("ℹ️ Cannot delete the default configuration. Clear it first.")
            self.action_status.setStyleSheet("color: blue;")
        else:
            self.action_status.setText("")
    
    def _set_default(self):
        """Set the current configuration as default."""
        selected_config = self.config_combo.currentText()
        if not selected_config:
            return
        
        try:
            PseudopotentialsManager.set_default_config(
                selected_config,
                DEFAULT_PSEUDOPOTENTIALS_DIR
            )
            self.action_status.setText(f"✅ Set '{selected_config}' as default configuration")
            self.action_status.setStyleSheet("color: green;")
            self._load_existing_configs()
            
            QMessageBox.information(
                self, "Success",
                f"'{selected_config}' is now the default configuration.\n\n"
                "This configuration will be automatically selected in Calculation Setup and Workflow Builder."
            )
        except Exception as e:
            self.action_status.setText(f"❌ Error setting default: {e}")
            self.action_status.setStyleSheet("color: red;")
    
    def _clear_default(self):
        """Clear the default configuration."""
        try:
            PseudopotentialsManager.clear_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
            self.action_status.setText("✅ Default configuration cleared")
            self.action_status.setStyleSheet("color: green;")
            self._load_existing_configs()
        except Exception as e:
            self.action_status.setText(f"❌ Error clearing default: {e}")
            self.action_status.setStyleSheet("color: red;")
    
    def _delete_configuration(self):
        """Delete the selected configuration."""
        selected_config = self.config_combo.currentText()
        if not selected_config:
            return
        
        if selected_config == "default":
            QMessageBox.warning(
                self, "Warning",
                "Cannot delete the default configuration. Clear it first using the 'Clear Default' button."
            )
            return
        
        reply = QMessageBox.question(
            self, "Confirm Delete",
            f"Are you sure you want to delete the configuration '{selected_config}'?\n\n"
            "This will only delete the configuration file, not the pseudopotential files themselves.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            try:
                PseudopotentialsManager.delete_config(
                    selected_config,
                    DEFAULT_PSEUDOPOTENTIALS_DIR
                )
                self.action_status.setText(f"✅ Deleted configuration: {selected_config}")
                self.action_status.setStyleSheet("color: green;")
                self.loaded_group.setVisible(False)
                self.current_pseudos = None
                self._load_existing_configs()
            except Exception as e:
                self.action_status.setText(f"❌ Error deleting: {e}")
                self.action_status.setStyleSheet("color: red;")
    
    def refresh(self):
        """Refresh the page."""
        self._load_machines()
        self._load_existing_configs()

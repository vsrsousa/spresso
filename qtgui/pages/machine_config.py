"""
Machine Configuration Page for xespresso PySide6 GUI.

This module handles the machine configuration interface, allowing users to:
- Create and edit machine configurations
- Test connections to remote machines
- Save configurations for later use
"""

import os
import traceback

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QSpinBox, QTextEdit, QCheckBox, QPushButton,
    QGroupBox, QFormLayout, QMessageBox, QScrollArea,
    QFrame, QTableWidget, QTableWidgetItem, QHeaderView,
    QApplication
)
from PySide6.QtCore import Qt

try:
    from xespresso.machines.machine import Machine
    from xespresso.machines.config.loader import (
        load_machine, save_machine, list_machines,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False

# SSH connection testing is optional (requires paramiko)
try:
    from xespresso.utils.auth import test_ssh_connection
    SSH_TEST_AVAILABLE = True
except ImportError:
    SSH_TEST_AVAILABLE = False
    def test_ssh_connection(username, host, key_path=None, port=22):
        """Fallback SSH connection test when paramiko is not available."""
        import subprocess
        
        key_path = os.path.expanduser(key_path) if key_path else None
        cmd = ["ssh", "-p", str(port), "-o", "PasswordAuthentication=no", 
               "-o", "BatchMode=yes", "-o", "ConnectTimeout=5"]
        if key_path:
            cmd += ["-i", key_path]
        cmd += [f"{username}@{host}", "echo 'Connection successful'"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False


class MachineConfigPage(QWidget):
    """Machine configuration page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self._loading = False  # Guard to prevent infinite loops
        self._setup_ui()
        self._load_machines_list()
    
    def _setup_ui(self):
        """Setup the user interface."""
        # Main layout with scroll area
        main_layout = QVBoxLayout(self)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Header
        header_label = QLabel("<h2>⚙️ Machine Configuration</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p><b>Create and configure</b> computational machines/clusters for your calculations.</p>
<p>This page is for <b>configuration only</b> - once machines are saved, you can select them 
in the Calculation Setup or Workflow Builder pages.</p>
<p>Supports both local and remote (SSH) execution environments.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        info_label = QLabel("""
<p>💡 <b>Configuration vs. Selection:</b></p>
<ul>
<li><b>Configure</b> machines here (one-time setup or updates)</li>
<li><b>Select</b> configured machines in Calculation Setup or Workflow Builder</li>
<li>Configurations are saved to <code>~/.xespresso/machines/</code></li>
</ul>
""")
        info_label.setTextFormat(Qt.RichText)
        info_label.setWordWrap(True)
        info_label.setStyleSheet("background-color: #e3f2fd; color: #1a1a1a; padding: 10px; border-radius: 5px;")
        scroll_layout.addWidget(info_label)
        
        if not XESPRESSO_AVAILABLE:
            error_label = QLabel("❌ xespresso modules not available. Cannot configure machines.")
            error_label.setStyleSheet("color: red; font-weight: bold;")
            scroll_layout.addWidget(error_label)
            scroll_area.setWidget(scroll_widget)
            main_layout.addWidget(scroll_area)
            return
        
        # Existing Machines Section
        machines_group = QGroupBox("Existing Machines")
        machines_layout = QVBoxLayout(machines_group)
        
        self.machines_status_label = QLabel("")
        machines_layout.addWidget(self.machines_status_label)
        
        select_layout = QHBoxLayout()
        select_layout.addWidget(QLabel("Select a machine to edit:"))
        self.machine_combo = QComboBox()
        self.machine_combo.currentTextChanged.connect(self._on_machine_selected)
        select_layout.addWidget(self.machine_combo, 1)
        machines_layout.addLayout(select_layout)
        
        scroll_layout.addWidget(machines_group)
        
        # Machine Configuration Form
        config_group = QGroupBox("Machine Configuration Form")
        config_layout = QFormLayout(config_group)
        
        # Basic settings
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("Unique identifier for this machine")
        config_layout.addRow("Machine Name:", self.name_edit)
        
        self.execution_combo = QComboBox()
        self.execution_combo.addItems(["local", "remote"])
        self.execution_combo.currentTextChanged.connect(self._on_execution_changed)
        config_layout.addRow("Execution Mode:", self.execution_combo)
        
        self.scheduler_combo = QComboBox()
        self.scheduler_combo.addItems(["direct", "slurm", "pbs", "sge"])
        self.scheduler_combo.currentTextChanged.connect(self._on_scheduler_changed)
        config_layout.addRow("Scheduler Type:", self.scheduler_combo)
        
        self.workdir_edit = QLineEdit("./calculations")
        self.workdir_edit.setPlaceholderText("Directory for calculation files")
        config_layout.addRow("Working Directory:", self.workdir_edit)
        
        self.nprocs_spin = QSpinBox()
        self.nprocs_spin.setRange(1, 10000)
        self.nprocs_spin.setValue(1)
        config_layout.addRow("Number of Processors:", self.nprocs_spin)
        
        self.launcher_edit = QLineEdit("mpirun -np {nprocs}")
        self.launcher_edit.setPlaceholderText("MPI launch command template")
        config_layout.addRow("MPI Launcher:", self.launcher_edit)
        
        scroll_layout.addWidget(config_group)
        
        # Remote Connection Settings
        self.remote_group = QGroupBox("Remote Connection Settings")
        remote_layout = QFormLayout(self.remote_group)
        
        self.host_edit = QLineEdit()
        self.host_edit.setPlaceholderText("Remote hostname or IP")
        remote_layout.addRow("Host:", self.host_edit)
        
        self.username_edit = QLineEdit()
        self.username_edit.setPlaceholderText("SSH username")
        remote_layout.addRow("Username:", self.username_edit)
        
        self.port_spin = QSpinBox()
        self.port_spin.setRange(1, 65535)
        self.port_spin.setValue(22)
        remote_layout.addRow("SSH Port:", self.port_spin)
        
        self.ssh_key_edit = QLineEdit("~/.ssh/id_rsa")
        self.ssh_key_edit.setPlaceholderText("Path to SSH private key")
        remote_layout.addRow("SSH Key Path:", self.ssh_key_edit)
        
        self.remote_group.setVisible(False)
        scroll_layout.addWidget(self.remote_group)
        
        # Environment Modules
        modules_group = QGroupBox("Environment Modules")
        modules_layout = QVBoxLayout(modules_group)
        
        self.use_modules_check = QCheckBox("Use Environment Modules")
        self.use_modules_check.stateChanged.connect(self._on_modules_changed)
        modules_layout.addWidget(self.use_modules_check)
        
        self.modules_edit = QTextEdit()
        self.modules_edit.setPlaceholderText("Environment modules to load (one per line)")
        self.modules_edit.setMaximumHeight(100)
        self.modules_edit.setVisible(False)
        modules_layout.addWidget(self.modules_edit)
        
        scroll_layout.addWidget(modules_group)
        
        # Advanced Settings
        advanced_group = QGroupBox("Advanced Settings")
        advanced_layout = QFormLayout(advanced_group)
        
        self.prepend_edit = QTextEdit()
        self.prepend_edit.setPlaceholderText("Commands to run before calculation")
        self.prepend_edit.setMaximumHeight(80)
        advanced_layout.addRow("Prepend Commands:", self.prepend_edit)
        
        self.postpend_edit = QTextEdit()
        self.postpend_edit.setPlaceholderText("Commands to run after calculation")
        self.postpend_edit.setMaximumHeight(80)
        advanced_layout.addRow("Postpend Commands:", self.postpend_edit)
        
        self.env_setup_edit = QLineEdit()
        self.env_setup_edit.setPlaceholderText("Shell commands to setup environment")
        advanced_layout.addRow("Environment Setup:", self.env_setup_edit)
        
        scroll_layout.addWidget(advanced_group)
        
        # Scheduler Resources
        self.scheduler_group = QGroupBox("Scheduler Resources")
        scheduler_layout = QFormLayout(self.scheduler_group)
        
        self.nodes_spin = QSpinBox()
        self.nodes_spin.setRange(1, 10000)
        self.nodes_spin.setValue(1)
        scheduler_layout.addRow("Nodes:", self.nodes_spin)
        
        self.ntasks_spin = QSpinBox()
        self.ntasks_spin.setRange(1, 256)
        self.ntasks_spin.setValue(20)
        scheduler_layout.addRow("Tasks per Node:", self.ntasks_spin)
        
        self.time_edit = QLineEdit("24:00:00")
        scheduler_layout.addRow("Wall Time:", self.time_edit)
        
        self.partition_edit = QLineEdit()
        self.partition_edit.setPlaceholderText("Partition/Queue name")
        scheduler_layout.addRow("Partition/Queue:", self.partition_edit)
        
        self.scheduler_group.setVisible(False)
        scroll_layout.addWidget(self.scheduler_group)
        
        # Buttons
        buttons_layout = QHBoxLayout()
        
        save_btn = QPushButton("💾 Save Machine Configuration")
        save_btn.clicked.connect(self._save_machine)
        buttons_layout.addWidget(save_btn)
        
        test_btn = QPushButton("🔍 Test Connection")
        test_btn.clicked.connect(self._test_connection)
        buttons_layout.addWidget(test_btn)
        
        scroll_layout.addLayout(buttons_layout)
        
        # Results area
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        scroll_layout.addWidget(self.results_label)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
    
    def _load_machines_list(self):
        """Load the list of existing machines."""
        if not XESPRESSO_AVAILABLE:
            return
        
        self.machine_combo.blockSignals(True)
        try:
            machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            self.machine_combo.clear()
            self.machine_combo.addItem("[Create New Machine]")
            
            if machines_list:
                self.machines_status_label.setText(
                    f"✅ {len(machines_list)} machine(s) configured: {', '.join(machines_list)}"
                )
                self.machines_status_label.setStyleSheet("color: green;")
                for machine in machines_list:
                    self.machine_combo.addItem(machine)
            else:
                self.machines_status_label.setText(
                    "ℹ️ No machines configured yet. Create your first machine below."
                )
                self.machines_status_label.setStyleSheet("color: blue;")
        except Exception as e:
            self.machines_status_label.setText(f"⚠️ Could not load machines list: {e}")
            self.machines_status_label.setStyleSheet("color: #d97706;")
        finally:
            self.machine_combo.blockSignals(False)
    
    def _on_machine_selected(self, machine_name):
        """Handle machine selection."""
        if machine_name == "[Create New Machine]" or not machine_name:
            self._clear_form()
            return
        
        # Prevent recursive updates
        if self._loading:
            return
        
        try:
            self._loading = True
            machine = load_machine(DEFAULT_CONFIG_PATH, machine_name, DEFAULT_MACHINES_DIR, return_object=True)
            self._populate_form(machine)
            self.session_state['current_machine'] = machine
            self.session_state['current_machine_name'] = machine_name
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error loading machine: {e}")
        finally:
            self._loading = False
    
    def _populate_form(self, machine):
        """Populate the form with machine data."""
        self.name_edit.setText(machine.name)
        
        idx = self.execution_combo.findText(machine.execution)
        if idx >= 0:
            self.execution_combo.setCurrentIndex(idx)
        
        idx = self.scheduler_combo.findText(machine.scheduler)
        if idx >= 0:
            self.scheduler_combo.setCurrentIndex(idx)
        
        self.workdir_edit.setText(machine.workdir)
        self.nprocs_spin.setValue(machine.nprocs)
        self.launcher_edit.setText(machine.launcher)
        
        # Remote settings
        if machine.is_remote:
            self.host_edit.setText(machine.host or "")
            self.username_edit.setText(machine.username or "")
            self.port_spin.setValue(machine.port or 22)
            if hasattr(machine, 'auth') and machine.auth:
                self.ssh_key_edit.setText(machine.auth.get("ssh_key", "~/.ssh/id_rsa"))
        
        # Modules
        use_modules = getattr(machine, 'use_modules', False)
        self.use_modules_check.setChecked(use_modules)
        if hasattr(machine, 'modules') and machine.modules:
            self.modules_edit.setText("\n".join(machine.modules))
        
        # Advanced
        if hasattr(machine, 'prepend') and machine.prepend:
            if isinstance(machine.prepend, list):
                self.prepend_edit.setText("\n".join(machine.prepend))
            else:
                self.prepend_edit.setText(machine.prepend)
        
        if hasattr(machine, 'postpend') and machine.postpend:
            if isinstance(machine.postpend, list):
                self.postpend_edit.setText("\n".join(machine.postpend))
            else:
                self.postpend_edit.setText(machine.postpend)
        
        if hasattr(machine, 'env_setup'):
            self.env_setup_edit.setText(machine.env_setup or "")
        
        # Scheduler resources
        if hasattr(machine, 'resources') and machine.resources:
            res = machine.resources
            self.nodes_spin.setValue(res.get('nodes', 1))
            self.ntasks_spin.setValue(res.get('ntasks-per-node', 20))
            self.time_edit.setText(res.get('time', '24:00:00'))
            self.partition_edit.setText(res.get('partition', ''))
    
    def _clear_form(self):
        """Clear the configuration form."""
        self.name_edit.clear()
        self.execution_combo.setCurrentIndex(0)
        self.scheduler_combo.setCurrentIndex(0)
        self.workdir_edit.setText("./calculations")
        self.nprocs_spin.setValue(1)
        self.launcher_edit.setText("mpirun -np {nprocs}")
        self.host_edit.clear()
        self.username_edit.clear()
        self.port_spin.setValue(22)
        self.ssh_key_edit.setText("~/.ssh/id_rsa")
        self.use_modules_check.setChecked(False)
        self.modules_edit.clear()
        self.prepend_edit.clear()
        self.postpend_edit.clear()
        self.env_setup_edit.clear()
        self.nodes_spin.setValue(1)
        self.ntasks_spin.setValue(20)
        self.time_edit.setText("24:00:00")
        self.partition_edit.clear()
    
    def _on_execution_changed(self, execution):
        """Handle execution mode change."""
        self.remote_group.setVisible(execution == "remote")
    
    def _on_scheduler_changed(self, scheduler):
        """Handle scheduler type change."""
        self.scheduler_group.setVisible(scheduler != "direct")
    
    def _on_modules_changed(self, state):
        """Handle modules checkbox change."""
        self.modules_edit.setVisible(state == Qt.Checked)
    
    def _save_machine(self):
        """Save the machine configuration."""
        if not XESPRESSO_AVAILABLE:
            QMessageBox.critical(self, "Error", "xespresso modules not available")
            return
        
        try:
            machine_name = self.name_edit.text().strip()
            if not machine_name:
                QMessageBox.warning(self, "Warning", "Please enter a machine name")
                return
            
            # Build machine config
            machine_config = {
                "name": machine_name,
                "execution": self.execution_combo.currentText(),
                "scheduler": self.scheduler_combo.currentText(),
                "workdir": self.workdir_edit.text(),
                "nprocs": self.nprocs_spin.value(),
                "launcher": self.launcher_edit.text(),
                "use_modules": self.use_modules_check.isChecked(),
            }
            
            if self.use_modules_check.isChecked():
                modules_text = self.modules_edit.toPlainText()
                machine_config["modules"] = [m.strip() for m in modules_text.split("\n") if m.strip()]
            
            prepend_text = self.prepend_edit.toPlainText()
            if prepend_text:
                machine_config["prepend"] = [p.strip() for p in prepend_text.split("\n") if p.strip()]
            
            postpend_text = self.postpend_edit.toPlainText()
            if postpend_text:
                machine_config["postpend"] = [p.strip() for p in postpend_text.split("\n") if p.strip()]
            
            env_setup = self.env_setup_edit.text()
            if env_setup:
                machine_config["env_setup"] = env_setup
            
            if self.execution_combo.currentText() == "remote":
                machine_config["host"] = self.host_edit.text()
                machine_config["username"] = self.username_edit.text()
                machine_config["port"] = self.port_spin.value()
                machine_config["auth"] = {
                    "method": "key",
                    "ssh_key": self.ssh_key_edit.text()
                }
            
            if self.scheduler_combo.currentText() != "direct":
                machine_config["resources"] = {
                    "nodes": self.nodes_spin.value(),
                    "ntasks-per-node": self.ntasks_spin.value(),
                    "time": self.time_edit.text(),
                }
                partition = self.partition_edit.text()
                if partition:
                    machine_config["resources"]["partition"] = partition
            
            # Create Machine object
            new_machine = Machine(**machine_config)
            
            # Save machine
            save_machine(new_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            
            self.results_label.setText(
                f"✅ Machine '{machine_name}' saved successfully!\n"
                f"💾 Configuration saved to: ~/.xespresso/machines/{machine_name}.json"
            )
            self.results_label.setStyleSheet("color: green;")
            
            self.session_state['current_machine'] = new_machine
            self.session_state['current_machine_name'] = machine_name
            
            # Refresh machines list
            self._load_machines_list()
            
        except Exception as e:
            self.results_label.setText(f"❌ Error saving machine: {e}")
            self.results_label.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error saving machine:\n{traceback.format_exc()}")
    
    def _test_connection(self):
        """Test the connection to the machine."""
        execution = self.execution_combo.currentText()
        
        if execution == "local":
            workdir = self.workdir_edit.text()
            user = os.environ.get('USER', 'unknown')
            self.results_label.setText(
                f"✅ Local machine - connection OK\n"
                f"Working directory: {workdir}\n"
                f"Current user: {user}"
            )
            self.results_label.setStyleSheet("color: green;")
        else:
            # Test remote connection
            host = self.host_edit.text()
            username = self.username_edit.text()
            port = self.port_spin.value()
            ssh_key = self.ssh_key_edit.text()
            
            if not host or not username:
                QMessageBox.warning(self, "Warning", "Please enter host and username")
                return
            
            key_path = os.path.expanduser(ssh_key)
            
            if not os.path.isfile(key_path):
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
    
    def refresh(self):
        """Refresh the page."""
        # Use loading guard to prevent infinite loops
        if self._loading:
            return
        self._loading = True
        try:
            self._load_machines_list()
            self._clear_form()
        finally:
            self._loading = False

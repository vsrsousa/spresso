"""
Pseudopotentials selector widget for the Qt GUI.

This module provides a reusable component for selecting pseudopotentials
from configured libraries, similar to the Streamlit version.
"""

import os

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QTableWidget, QTableWidgetItem, QHeaderView
)
from PySide6.QtCore import Qt

try:
    from xespresso.pseudopotentials import (
        PseudopotentialsManager,
        DEFAULT_PSEUDOPOTENTIALS_DIR
    )
    PSEUDOPOTENTIALS_AVAILABLE = True
except ImportError:
    PSEUDOPOTENTIALS_AVAILABLE = False
    DEFAULT_PSEUDOPOTENTIALS_DIR = os.path.expanduser("~/.xespresso/pseudopotentials")


class PseudopotentialsSelectorWidget(QWidget):
    """
    Widget for selecting pseudopotentials from configured libraries.
    
    This provides two modes:
    1. Library mode: Select from configured pseudopotential configurations
    2. Manual mode: Enter pseudopotential filenames manually
    """
    
    def __init__(self, session_state, parent=None):
        super().__init__(parent)
        self.session_state = session_state
        self.elements = set()
        self.pseudo_edits = {}
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the user interface."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Info label
        self.info_label = QLabel("Load a structure to configure pseudopotentials")
        layout.addWidget(self.info_label)
        
        if PSEUDOPOTENTIALS_AVAILABLE:
            # Configuration selector
            config_layout = QHBoxLayout()
            config_layout.addWidget(QLabel("Configuration:"))
            
            self.config_combo = QComboBox()
            self.config_combo.currentTextChanged.connect(self._on_config_changed)
            config_layout.addWidget(self.config_combo, 1)
            
            refresh_btn = QPushButton("🔄")
            refresh_btn.setToolTip("Refresh configurations list")
            refresh_btn.setMaximumWidth(40)
            refresh_btn.clicked.connect(self._load_configurations)
            config_layout.addWidget(refresh_btn)
            
            layout.addLayout(config_layout)
            
            # Configuration details
            self.config_details = QLabel("")
            self.config_details.setWordWrap(True)
            self.config_details.setStyleSheet("color: gray; font-size: 11px;")
            layout.addWidget(self.config_details)
            
            # Status label
            self.status_label = QLabel("")
            self.status_label.setWordWrap(True)
            layout.addWidget(self.status_label)
        
        # Manual input section (always available as fallback)
        self.manual_group = QGroupBox("Pseudopotentials per Element")
        self.manual_layout = QFormLayout(self.manual_group)
        layout.addWidget(self.manual_group)
        
        # Load configurations
        if PSEUDOPOTENTIALS_AVAILABLE:
            self._load_configurations()
    
    def _load_configurations(self):
        """Load available pseudopotential configurations."""
        if not PSEUDOPOTENTIALS_AVAILABLE:
            return
        
        self.config_combo.blockSignals(True)
        self.config_combo.clear()
        self.config_combo.addItem("Manual Entry", None)
        
        try:
            configs = PseudopotentialsManager.list_configs(DEFAULT_PSEUDOPOTENTIALS_DIR)
            
            # Check for default configuration
            has_default = PseudopotentialsManager.has_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
            
            if has_default:
                self.config_combo.addItem("⭐ Default Configuration", "default")
            
            for config_name in configs:
                if config_name != "default":
                    self.config_combo.addItem(config_name, config_name)
            
            # Select default if available
            if has_default:
                self.config_combo.setCurrentIndex(1)
            
        except Exception as e:
            self.status_label.setText(f"⚠️ Could not load configurations: {e}")
            self.status_label.setStyleSheet("color: orange;")
        
        self.config_combo.blockSignals(False)
        self._on_config_changed(self.config_combo.currentText())
    
    def _on_config_changed(self, config_name):
        """Handle configuration selection change."""
        config_data = self.config_combo.currentData()
        
        if config_data is None:
            # Manual mode
            self.config_details.setText("Enter pseudopotential filenames manually for each element.")
            self.manual_group.setTitle("Pseudopotentials per Element (Manual)")
            self._update_manual_inputs()
            return
        
        if not PSEUDOPOTENTIALS_AVAILABLE:
            return
        
        try:
            if config_data == "default":
                config = PseudopotentialsManager.get_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
            else:
                config = PseudopotentialsManager.load_config(config_data, DEFAULT_PSEUDOPOTENTIALS_DIR)
            
            if config:
                # Show configuration details
                details = []
                if config.library:
                    details.append(f"Library: {config.library}")
                if config.version:
                    details.append(f"v{config.version}")
                if config.functional:
                    details.append(f"Functional: {config.functional}")
                details.append(f"{len(config.pseudopotentials)} elements")
                
                self.config_details.setText(" | ".join(details))
                self.manual_group.setTitle(f"Pseudopotentials ({config_name})")
                
                # Auto-fill from configuration
                self._update_inputs_from_config(config)
            else:
                self.status_label.setText(f"⚠️ Could not load configuration: {config_data}")
                self.status_label.setStyleSheet("color: orange;")
                
        except Exception as e:
            self.status_label.setText(f"❌ Error: {e}")
            self.status_label.setStyleSheet("color: red;")
    
    def set_elements(self, elements):
        """
        Set the elements that need pseudopotentials.
        
        Args:
            elements: Set or list of element symbols
        """
        self.elements = set(elements)
        self.info_label.setText(f"Configure pseudopotentials for: {', '.join(sorted(self.elements))}")
        
        # Update inputs based on current mode
        config_data = self.config_combo.currentData() if PSEUDOPOTENTIALS_AVAILABLE else None
        
        if config_data is None:
            self._update_manual_inputs()
        else:
            self._on_config_changed(self.config_combo.currentText())
    
    def _update_manual_inputs(self):
        """Update manual input fields for current elements."""
        # Clear existing inputs
        while self.manual_layout.count():
            item = self.manual_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.pseudo_edits = {}
        
        for element in sorted(self.elements):
            edit = QLineEdit()
            edit.setPlaceholderText(f"e.g., {element}.UPF or {element}.pbe-n-kjpaw_psl.1.0.0.UPF")
            self.pseudo_edits[element] = edit
            self.manual_layout.addRow(f"{element}:", edit)
        
        self.status_label.setText("")
    
    def _update_inputs_from_config(self, config):
        """Update inputs from a pseudopotential configuration."""
        # Clear existing inputs
        while self.manual_layout.count():
            item = self.manual_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        self.pseudo_edits = {}
        pseudo_dict = config.get_pseudopotentials_dict()
        
        missing = []
        for element in sorted(self.elements):
            edit = QLineEdit()
            
            if element in pseudo_dict:
                edit.setText(pseudo_dict[element])
                edit.setStyleSheet("")
            else:
                edit.setPlaceholderText(f"Not found - enter manually")
                edit.setStyleSheet("background-color: #fff3e0;")
                missing.append(element)
            
            self.pseudo_edits[element] = edit
            self.manual_layout.addRow(f"{element}:", edit)
        
        if missing:
            self.status_label.setText(f"⚠️ Missing pseudopotentials for: {', '.join(missing)}")
            self.status_label.setStyleSheet("color: orange;")
        else:
            self.status_label.setText("✅ All elements found in configuration")
            self.status_label.setStyleSheet("color: green;")
            
            # Set ESPRESSO_PSEUDO environment variable
            # Note: This affects the entire process and allows xespresso to find
            # pseudopotential files. This is the standard way to configure
            # pseudopotential paths in Quantum ESPRESSO / ASE.
            os.environ["ESPRESSO_PSEUDO"] = config.base_path
    
    def get_pseudopotentials(self):
        """
        Get the current pseudopotential configuration.
        
        Returns:
            dict: Mapping of element symbols to pseudopotential filenames
        """
        result = {}
        for element, edit in self.pseudo_edits.items():
            value = edit.text().strip()
            if value:
                result[element] = value
        return result
    
    def set_pseudopotentials(self, pseudopotentials):
        """
        Set pseudopotential values from a saved configuration.
        
        This method restores pseudopotential values from a loaded session.
        It should be called after set_elements() to ensure the input fields exist.
        
        Args:
            pseudopotentials (dict): Mapping of element symbols to pseudopotential filenames
        """
        if not pseudopotentials:
            return
        
        for element, pseudo_file in pseudopotentials.items():
            if element in self.pseudo_edits:
                self.pseudo_edits[element].setText(pseudo_file)
    
    def is_complete(self):
        """
        Check if all elements have pseudopotentials assigned.
        
        Returns:
            bool: True if all elements have pseudopotentials
        """
        for element in self.elements:
            if element not in self.pseudo_edits:
                return False
            if not self.pseudo_edits[element].text().strip():
                return False
        return True

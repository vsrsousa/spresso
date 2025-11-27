"""
Pseudopotentials Configuration Page for xespresso PyQt GUI.

This module handles the pseudopotential configuration interface.
"""

import os

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTableWidget, QTableWidgetItem,
    QHeaderView, QFileDialog
)
from PyQt5.QtCore import Qt

try:
    from xespresso.pseudopotentials import (
        get_pseudo_families, get_pseudopotentials_for_element,
        DEFAULT_PSEUDO_DIR
    )
    PSEUDO_AVAILABLE = True
except ImportError:
    PSEUDO_AVAILABLE = False


class PseudopotentialsConfigPage(QWidget):
    """Pseudopotentials configuration page widget."""
    
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
        header_label = QLabel("<h2>🧪 Pseudopotentials Configuration</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p><b>Configure pseudopotential files</b> for your calculations.</p>
<p>Pseudopotentials are selected per element in the Calculation Setup page.</p>
<p>Here you can browse available pseudopotential families and set default directories.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        # Pseudopotential Directory
        dir_group = QGroupBox("Pseudopotential Directory")
        dir_layout = QFormLayout(dir_group)
        
        dir_h_layout = QHBoxLayout()
        self.pseudo_dir_edit = QLineEdit()
        if PSEUDO_AVAILABLE:
            self.pseudo_dir_edit.setText(DEFAULT_PSEUDO_DIR)
        else:
            self.pseudo_dir_edit.setText(os.path.expanduser("~/.xespresso/pseudopotentials"))
        dir_h_layout.addWidget(self.pseudo_dir_edit)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_pseudo_dir)
        dir_h_layout.addWidget(browse_btn)
        
        dir_layout.addRow("Directory:", dir_h_layout)
        
        scan_btn = QPushButton("🔍 Scan for Pseudopotentials")
        scan_btn.clicked.connect(self._scan_pseudopotentials)
        dir_layout.addRow(scan_btn)
        
        scroll_layout.addWidget(dir_group)
        
        # Available Families
        families_group = QGroupBox("Available Pseudopotential Families")
        families_layout = QVBoxLayout(families_group)
        
        self.families_status = QLabel("")
        families_layout.addWidget(self.families_status)
        
        self.families_table = QTableWidget()
        self.families_table.setColumnCount(3)
        self.families_table.setHorizontalHeaderLabels(["Family", "Type", "Elements"])
        self.families_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        families_layout.addWidget(self.families_table)
        
        scroll_layout.addWidget(families_group)
        
        # Element Lookup
        lookup_group = QGroupBox("Element Pseudopotential Lookup")
        lookup_layout = QFormLayout(lookup_group)
        
        lookup_h_layout = QHBoxLayout()
        self.element_edit = QLineEdit()
        self.element_edit.setPlaceholderText("e.g., Fe, O, C")
        lookup_h_layout.addWidget(self.element_edit)
        
        lookup_btn = QPushButton("Find Pseudopotentials")
        lookup_btn.clicked.connect(self._lookup_element)
        lookup_h_layout.addWidget(lookup_btn)
        
        lookup_layout.addRow("Element:", lookup_h_layout)
        
        self.element_pseudos_table = QTableWidget()
        self.element_pseudos_table.setColumnCount(2)
        self.element_pseudos_table.setHorizontalHeaderLabels(["Family", "File"])
        self.element_pseudos_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        lookup_layout.addRow("Available:", self.element_pseudos_table)
        
        scroll_layout.addWidget(lookup_group)
        
        # Info
        info_label = QLabel("""
<p>💡 <b>Tips:</b></p>
<ul>
<li>Standard pseudopotential libraries are available from <a href="https://www.quantum-espresso.org/pseudopotentials">quantum-espresso.org</a></li>
<li>Popular families include: SSSP (Standard Solid-State Pseudopotentials), PSlibrary</li>
<li>Pseudopotentials are automatically assigned in Calculation Setup based on your structure</li>
</ul>
""")
        info_label.setTextFormat(Qt.RichText)
        info_label.setWordWrap(True)
        info_label.setOpenExternalLinks(True)
        scroll_layout.addWidget(info_label)
        
        # Results area
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        scroll_layout.addWidget(self.results_label)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
        
        # Initial scan
        self._scan_pseudopotentials()
    
    def _browse_pseudo_dir(self):
        """Browse for pseudopotential directory."""
        current_dir = self.pseudo_dir_edit.text() or os.path.expanduser("~")
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Pseudopotentials Directory",
            current_dir,
            QFileDialog.ShowDirsOnly
        )
        if directory:
            self.pseudo_dir_edit.setText(directory)
    
    def _scan_pseudopotentials(self):
        """Scan for available pseudopotential families."""
        if not PSEUDO_AVAILABLE:
            self.families_status.setText("⚠️ Pseudopotentials module not available")
            self.families_status.setStyleSheet("color: orange;")
            return
        
        try:
            families = get_pseudo_families()
            
            if families:
                self.families_status.setText(f"✅ Found {len(families)} pseudopotential families")
                self.families_status.setStyleSheet("color: green;")
                
                self.families_table.setRowCount(len(families))
                for i, family in enumerate(families):
                    self.families_table.setItem(i, 0, QTableWidgetItem(family.get('name', 'Unknown')))
                    self.families_table.setItem(i, 1, QTableWidgetItem(family.get('type', 'Unknown')))
                    elements = family.get('elements', [])
                    self.families_table.setItem(i, 2, QTableWidgetItem(f"{len(elements)} elements"))
            else:
                self.families_status.setText("ℹ️ No pseudopotential families found")
                self.families_status.setStyleSheet("color: blue;")
                self.families_table.setRowCount(0)
                
        except Exception as e:
            self.families_status.setText(f"⚠️ Error scanning: {e}")
            self.families_status.setStyleSheet("color: orange;")
    
    def _lookup_element(self):
        """Look up pseudopotentials for an element."""
        element = self.element_edit.text().strip()
        
        if not element:
            QMessageBox.warning(self, "Warning", "Please enter an element symbol")
            return
        
        if not PSEUDO_AVAILABLE:
            self.results_label.setText("⚠️ Pseudopotentials module not available")
            self.results_label.setStyleSheet("color: orange;")
            return
        
        try:
            pseudos = get_pseudopotentials_for_element(element)
            
            if pseudos:
                self.results_label.setText(f"✅ Found {len(pseudos)} pseudopotential(s) for {element}")
                self.results_label.setStyleSheet("color: green;")
                
                self.element_pseudos_table.setRowCount(len(pseudos))
                for i, pseudo in enumerate(pseudos):
                    self.element_pseudos_table.setItem(i, 0, QTableWidgetItem(pseudo.get('family', 'Unknown')))
                    self.element_pseudos_table.setItem(i, 1, QTableWidgetItem(pseudo.get('file', 'Unknown')))
            else:
                self.results_label.setText(f"ℹ️ No pseudopotentials found for {element}")
                self.results_label.setStyleSheet("color: blue;")
                self.element_pseudos_table.setRowCount(0)
                
        except Exception as e:
            self.results_label.setText(f"⚠️ Error: {e}")
            self.results_label.setStyleSheet("color: orange;")
    
    def refresh(self):
        """Refresh the page."""
        self._scan_pseudopotentials()

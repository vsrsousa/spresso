"""
Simplified PySide6 Main Application for xespresso GUI.

This is a streamlined version that uses standard Qt patterns:
- Simple QMainWindow with QTabWidget (not complex page navigation)
- Direct state management without complex listeners
- Lazy tab content loading
- Minimal signal connections to avoid recursion
"""

import sys
import os

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QPushButton, QFileDialog, QMessageBox,
    QStatusBar, QMenuBar, QMenu, QToolBar, QGroupBox, QFormLayout,
    QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QTextEdit,
    QScrollArea, QFrame, QSplitter, QTreeWidget, QTreeWidgetItem,
    QTableWidget, QTableWidgetItem, QHeaderView, QCheckBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction


class SimpleSessionState:
    """
    Simple session state - just a dictionary with no complex listeners.
    This avoids the recursion issues of the previous implementation.
    """
    def __init__(self):
        self._state = {
            'current_structure': None,
            'working_directory': os.path.expanduser("~"),
            'workflow_config': {},
        }
    
    def get(self, key, default=None):
        return self._state.get(key, default)
    
    def set(self, key, value):
        self._state[key] = value
    
    def __getitem__(self, key):
        return self._state.get(key)
    
    def __setitem__(self, key, value):
        self._state[key] = value


class StructureTab(QWidget):
    """Structure viewer tab - simplified."""
    
    def __init__(self, state):
        super().__init__()
        self.state = state
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Header
        header = QLabel("<h2>🔬 Structure Viewer</h2>")
        header.setTextFormat(Qt.RichText)
        layout.addWidget(header)
        
        # Status
        self.status_label = QLabel("No structure loaded")
        layout.addWidget(self.status_label)
        
        # Load structure button
        load_group = QGroupBox("Load Structure")
        load_layout = QVBoxLayout(load_group)
        
        load_btn = QPushButton("📂 Open Structure File")
        load_btn.clicked.connect(self._load_structure)
        load_layout.addWidget(load_btn)
        
        self.file_label = QLabel("")
        load_layout.addWidget(self.file_label)
        
        layout.addWidget(load_group)
        
        # Structure info
        info_group = QGroupBox("Structure Information")
        info_layout = QVBoxLayout(info_group)
        
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(200)
        info_layout.addWidget(self.info_text)
        
        layout.addWidget(info_group)
        layout.addStretch()
    
    def _load_structure(self):
        """Load a structure file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Structure File",
            os.path.expanduser("~"),
            "Structure Files (*.cif *.xyz *.pdb *.vasp);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            from ase import io as ase_io
            atoms = ase_io.read(file_path)
            self.state['current_structure'] = atoms
            
            formula = atoms.get_chemical_formula()
            natoms = len(atoms)
            
            self.status_label.setText(f"✅ Loaded: {formula} ({natoms} atoms)")
            self.status_label.setStyleSheet("color: green;")
            self.file_label.setText(file_path)
            
            # Update info
            info_lines = [
                f"Formula: {formula}",
                f"Atoms: {natoms}",
                f"Elements: {', '.join(sorted(set(atoms.get_chemical_symbols())))}",
            ]
            if atoms.cell is not None and atoms.pbc.any():
                info_lines.append(f"Volume: {atoms.get_volume():.2f} Å³")
            
            self.info_text.setText("\n".join(info_lines))
            
        except ImportError:
            QMessageBox.critical(self, "Error", "ASE library not available")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error loading file:\n{e}")


class CalculationTab(QWidget):
    """Calculation setup tab - simplified."""
    
    def __init__(self, state):
        super().__init__()
        self.state = state
        self._setup_ui()
    
    def _setup_ui(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        
        content = QWidget()
        layout = QVBoxLayout(content)
        
        # Header
        header = QLabel("<h2>📊 Calculation Setup</h2>")
        header.setTextFormat(Qt.RichText)
        layout.addWidget(header)
        
        # Basic parameters
        params_group = QGroupBox("Basic Parameters")
        params_layout = QFormLayout(params_group)
        
        self.calc_type = QComboBox()
        self.calc_type.addItems(["scf", "relax", "vc-relax"])
        params_layout.addRow("Calculation:", self.calc_type)
        
        self.ecutwfc = QDoubleSpinBox()
        self.ecutwfc.setRange(10, 200)
        self.ecutwfc.setValue(50)
        self.ecutwfc.setSuffix(" Ry")
        params_layout.addRow("Energy Cutoff:", self.ecutwfc)
        
        self.kpoints = QLineEdit("4 4 4")
        params_layout.addRow("K-points:", self.kpoints)
        
        layout.addWidget(params_group)
        
        # Pseudopotentials (simple)
        pseudo_group = QGroupBox("Pseudopotentials")
        pseudo_layout = QVBoxLayout(pseudo_group)
        
        pseudo_info = QLabel("Specify pseudopotential files for each element")
        pseudo_layout.addWidget(pseudo_info)
        
        self.pseudo_text = QTextEdit()
        self.pseudo_text.setPlaceholderText("Fe = Fe.UPF\nO = O.UPF")
        self.pseudo_text.setMaximumHeight(100)
        pseudo_layout.addWidget(self.pseudo_text)
        
        layout.addWidget(pseudo_group)
        
        # Save config button
        save_btn = QPushButton("💾 Save Configuration")
        save_btn.clicked.connect(self._save_config)
        layout.addWidget(save_btn)
        
        self.result_label = QLabel("")
        layout.addWidget(self.result_label)
        
        layout.addStretch()
        scroll.setWidget(content)
        
        main_layout = QVBoxLayout(self)
        main_layout.addWidget(scroll)
    
    def _save_config(self):
        """Save the calculation configuration."""
        config = {
            'calc_type': self.calc_type.currentText(),
            'ecutwfc': self.ecutwfc.value(),
            'kpoints': self.kpoints.text(),
            'pseudopotentials': self._parse_pseudos(),
        }
        self.state['workflow_config'] = config
        
        self.result_label.setText("✅ Configuration saved")
        self.result_label.setStyleSheet("color: green;")
    
    def _parse_pseudos(self):
        """Parse pseudopotential text."""
        pseudos = {}
        for line in self.pseudo_text.toPlainText().split('\n'):
            if '=' in line:
                parts = line.split('=')
                if len(parts) == 2:
                    elem = parts[0].strip()
                    path = parts[1].strip()
                    if elem and path:
                        pseudos[elem] = path
        return pseudos


class FileBrowserTab(QWidget):
    """Simple file browser tab."""
    
    def __init__(self, state):
        super().__init__()
        self.state = state
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Header
        header = QLabel("<h2>📂 File Browser</h2>")
        header.setTextFormat(Qt.RichText)
        layout.addWidget(header)
        
        # Working directory
        dir_layout = QHBoxLayout()
        dir_layout.addWidget(QLabel("Directory:"))
        self.dir_label = QLabel(self.state.get('working_directory', '~'))
        dir_layout.addWidget(self.dir_label, 1)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_dir)
        dir_layout.addWidget(browse_btn)
        
        refresh_btn = QPushButton("🔄 Refresh")
        refresh_btn.clicked.connect(self._refresh_files)
        dir_layout.addWidget(refresh_btn)
        
        layout.addLayout(dir_layout)
        
        # Splitter for tree and content
        splitter = QSplitter(Qt.Horizontal)
        
        # File tree
        self.file_tree = QTreeWidget()
        self.file_tree.setHeaderLabels(["Files"])
        self.file_tree.itemClicked.connect(self._on_file_clicked)
        splitter.addWidget(self.file_tree)
        
        # File content
        self.file_content = QTextEdit()
        self.file_content.setReadOnly(True)
        splitter.addWidget(self.file_content)
        
        splitter.setSizes([300, 700])
        layout.addWidget(splitter)
    
    def _browse_dir(self):
        """Browse for directory."""
        directory = QFileDialog.getExistingDirectory(
            self, "Select Directory",
            self.state.get('working_directory', os.path.expanduser("~"))
        )
        if directory:
            self.state['working_directory'] = directory
            self.dir_label.setText(directory)
            self._refresh_files()
    
    def _refresh_files(self):
        """Refresh file list - simple version with limits."""
        self.file_tree.clear()
        workdir = self.state.get('working_directory', os.path.expanduser("~"))
        
        if not os.path.isdir(workdir):
            return
        
        try:
            # Simple listing - just top level, no recursive walk
            items = []
            for entry in os.scandir(workdir):
                if entry.is_file() and not entry.name.startswith('.'):
                    items.append(entry.name)
                elif entry.is_dir() and not entry.name.startswith('.'):
                    items.append(entry.name + "/")
            
            items.sort()
            for name in items[:100]:  # Limit to 100 items
                item = QTreeWidgetItem([name])
                item.setData(0, Qt.UserRole, os.path.join(workdir, name.rstrip('/')))
                self.file_tree.addTopLevelItem(item)
                
        except (OSError, PermissionError) as e:
            self.file_content.setText(f"Error: {e}")
    
    def _on_file_clicked(self, item, column):
        """Handle file click."""
        path = item.data(0, Qt.UserRole)
        if path and os.path.isfile(path):
            try:
                with open(path, 'r', errors='replace') as f:
                    content = f.read(50000)  # Limit to 50KB
                self.file_content.setText(content)
            except Exception as e:
                self.file_content.setText(f"Error reading file: {e}")


class SimpleMainWindow(QMainWindow):
    """
    Simplified main window using standard Qt tabs.
    No complex page navigation or state listeners.
    """
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("⚛️ xespresso - Quantum ESPRESSO GUI")
        self.resize(1200, 800)
        self.setMinimumSize(800, 600)
        
        # Simple state - no listeners
        self.state = SimpleSessionState()
        
        self._setup_ui()
        self._setup_menu()
        self._setup_statusbar()
    
    def _setup_ui(self):
        """Setup the UI with a simple tab widget."""
        # Central widget with tabs
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        # Create tabs
        self.tabs.addTab(StructureTab(self.state), "🔬 Structure")
        self.tabs.addTab(CalculationTab(self.state), "📊 Calculation")
        self.tabs.addTab(FileBrowserTab(self.state), "📂 Files")
        
        # Add a simple info tab
        info_tab = QWidget()
        info_layout = QVBoxLayout(info_tab)
        info_text = QLabel("""
<h2>⚛️ xespresso GUI</h2>
<p>A simplified interface for Quantum ESPRESSO calculations.</p>

<h3>Quick Start:</h3>
<ol>
<li><b>Structure:</b> Load your atomic structure (CIF, XYZ, etc.)</li>
<li><b>Calculation:</b> Set up calculation parameters</li>
<li><b>Files:</b> Browse and view calculation files</li>
</ol>

<h3>For full features, use:</h3>
<ul>
<li>Command line: <code>python -m xespresso</code></li>
<li>Full Qt GUI: <code>python -m qtgui</code></li>
</ul>
""")
        info_text.setTextFormat(Qt.RichText)
        info_text.setWordWrap(True)
        info_layout.addWidget(info_text)
        info_layout.addStretch()
        
        self.tabs.addTab(info_tab, "ℹ️ Info")
    
    def _setup_menu(self):
        """Setup simple menu bar."""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Help menu
        help_menu = menubar.addMenu("&Help")
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)
    
    def _setup_statusbar(self):
        """Setup status bar."""
        self.statusBar().showMessage("Ready")
    
    def _show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About xespresso GUI",
            "<h2>xespresso GUI</h2>"
            "<p>Simplified PySide6 interface for Quantum ESPRESSO.</p>"
            "<p>Version: 2.0.0 (Simplified)</p>"
        )


def main():
    """Main entry point."""
    app = QApplication(sys.argv)
    app.setApplicationName("xespresso GUI")
    app.setStyle("Fusion")
    
    window = SimpleMainWindow()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

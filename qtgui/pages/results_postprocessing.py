"""
Results and Post-processing Page for xespresso PyQt GUI.

This page displays calculation results and provides post-processing tools.
"""

import os

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTextEdit,
    QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView,
    QFileDialog
)
from PyQt5.QtCore import Qt


class ResultsPostprocessingPage(QWidget):
    """Results and post-processing page widget."""
    
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
        header_label = QLabel("<h2>📈 Results & Post-Processing</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p>View calculation results and perform post-processing analysis.</p>
<p>Load output files from completed calculations to analyze:</p>
<ul>
<li>Total energies and convergence</li>
<li>Forces and stresses</li>
<li>Electronic structure (DOS, bands)</li>
<li>Relaxed structures</li>
</ul>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        # Load Results Section
        load_group = QGroupBox("📂 Load Results")
        load_layout = QVBoxLayout(load_group)
        
        dir_layout = QHBoxLayout()
        dir_layout.addWidget(QLabel("Output Directory:"))
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select calculation output directory")
        dir_layout.addWidget(self.output_dir_edit, 1)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_output_dir)
        dir_layout.addWidget(browse_btn)
        load_layout.addLayout(dir_layout)
        
        load_btn = QPushButton("📥 Load Results")
        load_btn.clicked.connect(self._load_results)
        load_layout.addWidget(load_btn)
        
        scroll_layout.addWidget(load_group)
        
        # Results Display
        results_group = QGroupBox("📊 Calculation Results")
        results_layout = QVBoxLayout(results_group)
        
        # Summary
        summary_layout = QHBoxLayout()
        
        self.energy_label = QLabel("Total Energy: --")
        self.energy_label.setStyleSheet("font-size: 14pt; font-weight: bold;")
        summary_layout.addWidget(self.energy_label)
        
        self.status_label = QLabel("Status: --")
        summary_layout.addWidget(self.status_label)
        
        results_layout.addLayout(summary_layout)
        
        # Detailed results
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setPlaceholderText("Load a calculation to view results...")
        results_layout.addWidget(self.results_text)
        
        scroll_layout.addWidget(results_group)
        
        # Convergence History
        convergence_group = QGroupBox("📉 Convergence History")
        convergence_layout = QVBoxLayout(convergence_group)
        
        self.convergence_table = QTableWidget()
        self.convergence_table.setColumnCount(4)
        self.convergence_table.setHorizontalHeaderLabels(["Iteration", "Energy (Ry)", "Delta E", "Convergence"])
        self.convergence_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        convergence_layout.addWidget(self.convergence_table)
        
        scroll_layout.addWidget(convergence_group)
        
        # Post-Processing Tools
        tools_group = QGroupBox("🔧 Post-Processing Tools")
        tools_layout = QVBoxLayout(tools_group)
        
        info_label = QLabel("""
<p>Post-processing tools available after loading results:</p>
<ul>
<li><b>DOS</b>: Density of States analysis (requires dos.x output)</li>
<li><b>Bands</b>: Band structure analysis (requires bands.x output)</li>
<li><b>Structure</b>: Extract relaxed structure from output</li>
</ul>
""")
        info_label.setTextFormat(Qt.RichText)
        info_label.setWordWrap(True)
        tools_layout.addWidget(info_label)
        
        buttons_layout = QHBoxLayout()
        
        dos_btn = QPushButton("📊 DOS Analysis")
        dos_btn.clicked.connect(self._dos_analysis)
        buttons_layout.addWidget(dos_btn)
        
        bands_btn = QPushButton("📈 Band Structure")
        bands_btn.clicked.connect(self._bands_analysis)
        buttons_layout.addWidget(bands_btn)
        
        structure_btn = QPushButton("🔬 Extract Structure")
        structure_btn.clicked.connect(self._extract_structure)
        buttons_layout.addWidget(structure_btn)
        
        tools_layout.addLayout(buttons_layout)
        
        scroll_layout.addWidget(tools_group)
        
        # Output File Viewer
        viewer_group = QGroupBox("📄 Output File Viewer")
        viewer_layout = QVBoxLayout(viewer_group)
        
        file_select_layout = QHBoxLayout()
        file_select_layout.addWidget(QLabel("Output File:"))
        self.file_combo = QComboBox()
        file_select_layout.addWidget(self.file_combo, 1)
        
        view_btn = QPushButton("View")
        view_btn.clicked.connect(self._view_file)
        file_select_layout.addWidget(view_btn)
        viewer_layout.addLayout(file_select_layout)
        
        self.file_viewer = QTextEdit()
        self.file_viewer.setReadOnly(True)
        self.file_viewer.setMaximumHeight(300)
        viewer_layout.addWidget(self.file_viewer)
        
        scroll_layout.addWidget(viewer_group)
        
        # Status
        self.status_text = QLabel("")
        self.status_text.setWordWrap(True)
        scroll_layout.addWidget(self.status_text)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
    
    def _browse_output_dir(self):
        """Browse for output directory."""
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Calculation Output Directory",
            workdir,
            QFileDialog.ShowDirsOnly
        )
        if directory:
            self.output_dir_edit.setText(directory)
    
    def _load_results(self):
        """Load calculation results."""
        output_dir = self.output_dir_edit.text()
        
        if not output_dir or not os.path.isdir(output_dir):
            QMessageBox.warning(self, "Warning", "Please select a valid output directory")
            return
        
        # Find output files
        output_files = []
        for f in os.listdir(output_dir):
            if f.endswith(('.pwo', '.out', '.log')):
                output_files.append(f)
        
        if not output_files:
            QMessageBox.warning(self, "Warning", "No output files found in directory")
            return
        
        # Update file combo
        self.file_combo.clear()
        for f in output_files:
            self.file_combo.addItem(f)
        
        # Try to parse the first output file
        output_path = os.path.join(output_dir, output_files[0])
        
        try:
            with open(output_path, 'r') as f:
                content = f.read()
            
            # Parse results
            results = self._parse_output(content)
            
            # Update display
            if results.get('energy'):
                self.energy_label.setText(f"Total Energy: {results['energy']:.6f} Ry")
            
            if results.get('converged'):
                self.status_label.setText("Status: ✅ Converged")
                self.status_label.setStyleSheet("color: green;")
            else:
                self.status_label.setText("Status: ⚠️ Not converged")
                self.status_label.setStyleSheet("color: orange;")
            
            # Results text
            results_text = []
            results_text.append(f"Output File: {output_files[0]}")
            results_text.append(f"Total Energy: {results.get('energy', 'N/A')} Ry")
            results_text.append(f"Converged: {results.get('converged', 'Unknown')}")
            results_text.append(f"SCF Iterations: {results.get('iterations', 'N/A')}")
            
            if results.get('total_force'):
                results_text.append(f"Total Force: {results['total_force']} Ry/au")
            
            self.results_text.setText("\n".join(results_text))
            
            # Update convergence table
            if results.get('scf_history'):
                history = results['scf_history']
                self.convergence_table.setRowCount(len(history))
                for i, (energy, delta) in enumerate(history):
                    self.convergence_table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
                    self.convergence_table.setItem(i, 1, QTableWidgetItem(f"{energy:.8f}"))
                    self.convergence_table.setItem(i, 2, QTableWidgetItem(f"{delta:.2e}"))
                    self.convergence_table.setItem(i, 3, QTableWidgetItem("✓" if abs(delta) < 1e-6 else ""))
            
            self.status_text.setText(f"✅ Results loaded from: {output_path}")
            self.status_text.setStyleSheet("color: green;")
            
        except Exception as e:
            self.status_text.setText(f"❌ Error loading results: {e}")
            self.status_text.setStyleSheet("color: red;")
    
    def _parse_output(self, content):
        """Parse QE output file content."""
        results = {
            'energy': None,
            'converged': False,
            'iterations': 0,
            'scf_history': [],
            'total_force': None
        }
        
        lines = content.split('\n')
        prev_energy = None
        
        for line in lines:
            # Total energy
            if '!' in line and 'total energy' in line.lower():
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        energy_str = parts[1].replace('Ry', '').strip()
                        results['energy'] = float(energy_str)
                except:
                    pass
            
            # SCF iteration
            if 'total energy' in line.lower() and '!' not in line:
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        energy_str = parts[1].replace('Ry', '').strip()
                        energy = float(energy_str)
                        delta = energy - prev_energy if prev_energy else 0
                        results['scf_history'].append((energy, delta))
                        prev_energy = energy
                        results['iterations'] += 1
                except:
                    pass
            
            # Convergence
            if 'convergence achieved' in line.lower():
                results['converged'] = True
            
            # Total force
            if 'Total force' in line:
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        force_str = parts[1].strip()
                        results['total_force'] = float(force_str)
                except:
                    pass
        
        return results
    
    def _view_file(self):
        """View selected output file."""
        output_dir = self.output_dir_edit.text()
        filename = self.file_combo.currentText()
        
        if not output_dir or not filename:
            return
        
        filepath = os.path.join(output_dir, filename)
        
        try:
            with open(filepath, 'r') as f:
                content = f.read()
            self.file_viewer.setText(content)
        except Exception as e:
            self.file_viewer.setText(f"Error reading file: {e}")
    
    def _dos_analysis(self):
        """Perform DOS analysis."""
        QMessageBox.information(
            self,
            "DOS Analysis",
            "DOS analysis requires output from dos.x calculation.\n\n"
            "Make sure you have:\n"
            "1. Completed an nscf calculation\n"
            "2. Run dos.x with appropriate input\n\n"
            "For full DOS analysis, use the xespresso.dos module."
        )
    
    def _bands_analysis(self):
        """Perform band structure analysis."""
        QMessageBox.information(
            self,
            "Band Structure Analysis",
            "Band structure analysis requires output from bands.x calculation.\n\n"
            "Make sure you have:\n"
            "1. Completed an nscf calculation along k-path\n"
            "2. Run bands.x with appropriate input\n\n"
            "For full band analysis, use the xespresso workflow modules."
        )
    
    def _extract_structure(self):
        """Extract relaxed structure from output."""
        output_dir = self.output_dir_edit.text()
        
        if not output_dir:
            QMessageBox.warning(self, "Warning", "Please load results first")
            return
        
        QMessageBox.information(
            self,
            "Extract Structure",
            "Structure extraction from relaxation output.\n\n"
            "The relaxed structure can be read using ASE:\n\n"
            "from ase.io import read\n"
            "atoms = read('espresso.pwo', format='espresso-out')\n\n"
            "This will return the final structure from the calculation."
        )
    
    def refresh(self):
        """Refresh the page."""
        pass

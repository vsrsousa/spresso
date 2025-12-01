"""
Structure Viewer Page for xespresso PySide6 GUI.

This page handles loading and visualizing atomic structures.
"""

import os
import tempfile

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTableWidget, QTableWidgetItem,
    QHeaderView, QFileDialog, QTabWidget, QTextEdit, QCheckBox,
    QDoubleSpinBox, QSpinBox, QSlider
)
from PySide6.QtCore import Qt

try:
    from ase import io as ase_io
    from ase import Atoms
    from ase.build import bulk, molecule
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

try:
    import matplotlib
    # Only set backend if not already set to avoid conflicts
    if matplotlib.get_backend() != 'QtAgg':
        matplotlib.use('QtAgg')
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    import numpy as np
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


class StructureViewerPage(QWidget):
    """Structure viewer page widget."""
    
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
        header_label = QLabel("<h2>🔬 Structure Viewer</h2>")
        header_label.setTextFormat(Qt.RichText)
        scroll_layout.addWidget(header_label)
        
        description = QLabel("""
<p>Load and visualize atomic structures.</p>
<p><b>Modular Design:</b> This page uses the structures module to handle loading and exporting.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        scroll_layout.addWidget(description)
        
        if not ASE_AVAILABLE:
            error_label = QLabel("❌ ASE not available. Structure viewing is disabled.")
            error_label.setStyleSheet("color: red; font-weight: bold;")
            scroll_layout.addWidget(error_label)
            scroll_area.setWidget(scroll_widget)
            main_layout.addWidget(scroll_area)
            return
        
        # Current Structure Info
        self.current_group = QGroupBox("📍 Currently Selected Structure")
        current_layout = QVBoxLayout(self.current_group)
        
        self.structure_info_label = QLabel("No structure loaded")
        current_layout.addWidget(self.structure_info_label)
        
        info_row = QHBoxLayout()
        self.formula_label = QLabel("")
        info_row.addWidget(self.formula_label)
        self.natoms_label = QLabel("")
        info_row.addWidget(self.natoms_label)
        self.source_label = QLabel("")
        info_row.addWidget(self.source_label)
        current_layout.addLayout(info_row)
        
        self.current_group.setVisible(False)
        scroll_layout.addWidget(self.current_group)
        
        # Tabs for different input methods
        tabs = QTabWidget()
        
        # Upload File Tab
        upload_tab = QWidget()
        upload_layout = QVBoxLayout(upload_tab)
        
        upload_btn_layout = QHBoxLayout()
        upload_btn = QPushButton("📂 Open Structure File")
        upload_btn.clicked.connect(self._open_structure_file)
        upload_btn_layout.addWidget(upload_btn)
        
        self.file_path_label = QLabel("")
        upload_btn_layout.addWidget(self.file_path_label, 1)
        upload_layout.addLayout(upload_btn_layout)
        
        upload_info = QLabel("Supported formats: CIF, XYZ, PDB, VASP, POSCAR, TRAJ, JSON")
        upload_layout.addWidget(upload_info)
        
        upload_layout.addStretch()
        tabs.addTab(upload_tab, "Upload File")
        
        # Build Structure Tab
        build_tab = QWidget()
        build_layout = QVBoxLayout(build_tab)
        
        build_type_layout = QHBoxLayout()
        build_type_layout.addWidget(QLabel("Structure Type:"))
        self.build_type_combo = QComboBox()
        self.build_type_combo.addItems(["Bulk Crystal", "Molecule"])
        self.build_type_combo.currentTextChanged.connect(self._on_build_type_changed)
        build_type_layout.addWidget(self.build_type_combo, 1)
        build_layout.addLayout(build_type_layout)
        
        # Bulk crystal options
        self.bulk_group = QGroupBox("🔷 Build Bulk Crystal")
        bulk_layout = QFormLayout(self.bulk_group)
        
        self.element_edit = QLineEdit("Fe")
        bulk_layout.addRow("Element:", self.element_edit)
        
        self.crystal_combo = QComboBox()
        self.crystal_combo.addItems(["fcc", "bcc", "hcp", "diamond", "sc"])
        bulk_layout.addRow("Crystal Structure:", self.crystal_combo)
        
        self.lattice_spin = QDoubleSpinBox()
        self.lattice_spin.setRange(0.1, 100.0)
        self.lattice_spin.setValue(3.6)
        self.lattice_spin.setSingleStep(0.1)
        self.lattice_spin.setSuffix(" Å")
        bulk_layout.addRow("Lattice Parameter:", self.lattice_spin)
        
        self.cubic_check = QCheckBox("Cubic Cell")
        self.cubic_check.setChecked(True)
        bulk_layout.addRow(self.cubic_check)
        
        build_bulk_btn = QPushButton("🔨 Build Crystal")
        build_bulk_btn.clicked.connect(self._build_crystal)
        bulk_layout.addRow(build_bulk_btn)
        
        build_layout.addWidget(self.bulk_group)
        
        # Molecule options
        self.molecule_group = QGroupBox("🧪 Build Molecule")
        molecule_layout = QFormLayout(self.molecule_group)
        
        self.molecule_edit = QLineEdit("H2O")
        self.molecule_edit.setPlaceholderText("H2O, CO2, CH4, NH3, C6H6, etc.")
        molecule_layout.addRow("Molecule Name:", self.molecule_edit)
        
        molecule_info = QLabel("💡 Tip: Try H2O, CO2, CH4, NH3, C6H6, or other common molecules")
        molecule_layout.addRow(molecule_info)
        
        build_molecule_btn = QPushButton("🔨 Build Molecule")
        build_molecule_btn.clicked.connect(self._build_molecule)
        molecule_layout.addRow(build_molecule_btn)
        
        self.molecule_group.setVisible(False)
        build_layout.addWidget(self.molecule_group)
        
        build_layout.addStretch()
        tabs.addTab(build_tab, "Build Structure")
        
        scroll_layout.addWidget(tabs)
        
        # Visualization Section
        viz_group = QGroupBox("🔬 Structure Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        # Visualization canvas
        if MATPLOTLIB_AVAILABLE:
            self.figure = Figure(figsize=(8, 6))
            self.canvas = FigureCanvas(self.figure)
            viz_layout.addWidget(self.canvas)
            
            # Visualization controls
            viz_controls = QHBoxLayout()
            
            refresh_viz_btn = QPushButton("🔄 Refresh View")
            refresh_viz_btn.clicked.connect(self._refresh_visualization)
            viz_controls.addWidget(refresh_viz_btn)
            
            viz_layout.addLayout(viz_controls)
        else:
            viz_label = QLabel("Matplotlib not available. Visualization disabled.")
            viz_layout.addWidget(viz_label)
        
        scroll_layout.addWidget(viz_group)
        
        # Structure Information
        info_group = QGroupBox("Structure Information")
        info_layout = QVBoxLayout(info_group)
        
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(200)
        info_layout.addWidget(self.info_text)
        
        scroll_layout.addWidget(info_group)
        
        # Export Section
        export_group = QGroupBox("💾 Export Options")
        export_layout = QHBoxLayout(export_group)
        
        export_layout.addWidget(QLabel("Export Format:"))
        self.export_combo = QComboBox()
        self.export_combo.addItems(["cif", "xyz", "pdb", "vasp", "json"])
        export_layout.addWidget(self.export_combo)
        
        export_btn = QPushButton("⬇️ Export Structure")
        export_btn.clicked.connect(self._export_structure)
        export_layout.addWidget(export_btn)
        
        scroll_layout.addWidget(export_group)
        
        # Results area
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        scroll_layout.addWidget(self.results_label)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        main_layout.addWidget(scroll_area)
    
    def _on_build_type_changed(self, build_type):
        """Handle build type change."""
        self.bulk_group.setVisible(build_type == "Bulk Crystal")
        self.molecule_group.setVisible(build_type == "Molecule")
    
    def _open_structure_file(self):
        """Open a structure file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Structure File",
            os.path.expanduser("~"),
            "Structure Files (*.cif *.xyz *.pdb *.vasp *.poscar *.traj *.json);;All Files (*)"
        )
        
        if file_path:
            try:
                atoms = ase_io.read(file_path)
                self._set_structure(atoms, f"File: {os.path.basename(file_path)}")
                self.file_path_label.setText(file_path)
                self.results_label.setText(f"✅ Loaded: {os.path.basename(file_path)}")
                self.results_label.setStyleSheet("color: green;")
            except Exception as e:
                self.results_label.setText(f"❌ Error loading file: {e}")
                self.results_label.setStyleSheet("color: red;")
                QMessageBox.critical(self, "Error", f"Error loading file:\n{e}")
    
    def _build_crystal(self):
        """Build a bulk crystal structure."""
        try:
            element = self.element_edit.text().strip()
            crystal = self.crystal_combo.currentText()
            a = self.lattice_spin.value()
            cubic = self.cubic_check.isChecked()
            
            atoms = bulk(element, crystal, a=a, cubic=cubic)
            self._set_structure(atoms, f"Built: {element} {crystal}")
            
            self.results_label.setText(f"✅ Built {element} {crystal} structure")
            self.results_label.setStyleSheet("color: green;")
        except Exception as e:
            self.results_label.setText(f"❌ Error building structure: {e}")
            self.results_label.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error building structure:\n{e}")
    
    def _build_molecule(self):
        """Build a molecule."""
        try:
            mol_name = self.molecule_edit.text().strip()
            
            atoms = molecule(mol_name)
            atoms.center(vacuum=5.0)
            
            self._set_structure(atoms, f"Built: {mol_name} molecule")
            
            self.results_label.setText(f"✅ Built {mol_name} molecule")
            self.results_label.setStyleSheet("color: green;")
        except Exception as e:
            self.results_label.setText(f"❌ Error building molecule: {e}")
            self.results_label.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error building molecule:\n{e}")
    
    def _set_structure(self, atoms, source):
        """Set the current structure."""
        # Prevent recursive updates
        if self._loading:
            return
        
        self._loading = True
        try:
            self.session_state['current_structure'] = atoms
            self.session_state['structure_source'] = source
            
            # Update UI
            self.current_group.setVisible(True)
            formula = atoms.get_chemical_formula()
            natoms = len(atoms)
            
            self.structure_info_label.setText(f"✅ Structure loaded: {formula}")
            self.formula_label.setText(f"<b>Formula:</b> {formula}")
            self.natoms_label.setText(f"<b>Atoms:</b> {natoms}")
            self.source_label.setText(f"<b>Source:</b> {source}")
            
            # Update structure info
            self._update_structure_info(atoms)
            
            # Update visualization
            self._refresh_visualization()
        finally:
            self._loading = False
    
    def _update_structure_info(self, atoms):
        """Update structure information display."""
        info_lines = []
        info_lines.append(f"Chemical Formula: {atoms.get_chemical_formula()}")
        info_lines.append(f"Number of Atoms: {len(atoms)}")
        
        symbols = atoms.get_chemical_symbols()
        unique_elements = list(set(symbols))
        info_lines.append(f"Elements: {', '.join(sorted(unique_elements))}")
        
        if atoms.cell is not None and atoms.pbc.any():
            info_lines.append(f"\nCell Volume: {atoms.get_volume():.2f} Å³")
            pbc_str = "".join(["T" if p else "F" for p in atoms.pbc])
            info_lines.append(f"PBC: {pbc_str}")
            
            cell_params = atoms.cell.cellpar()
            info_lines.append(f"\nCell Parameters:")
            info_lines.append(f"  a = {cell_params[0]:.3f} Å")
            info_lines.append(f"  b = {cell_params[1]:.3f} Å")
            info_lines.append(f"  c = {cell_params[2]:.3f} Å")
            info_lines.append(f"  α = {cell_params[3]:.2f}°")
            info_lines.append(f"  β = {cell_params[4]:.2f}°")
            info_lines.append(f"  γ = {cell_params[5]:.2f}°")
        
        self.info_text.setText("\n".join(info_lines))
    
    def _refresh_visualization(self):
        """Refresh the structure visualization."""
        if not MATPLOTLIB_AVAILABLE:
            return
        
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            return
        
        self.figure.clear()
        ax = self.figure.add_subplot(111, projection='3d')
        
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        # Color map for elements
        color_map = {
            'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red',
            'F': 'green', 'P': 'orange', 'S': 'yellow',
            'Cl': 'green', 'Fe': 'brown', 'Cu': 'brown',
            'Al': 'silver', 'Si': 'pink', 'Pt': 'silver'
        }
        
        colors = [color_map.get(s, 'purple') for s in symbols]
        
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                   c=colors, s=100, edgecolors='black')
        
        # Add labels
        for i, (pos, sym) in enumerate(zip(positions, symbols)):
            ax.text(pos[0], pos[1], pos[2], sym, fontsize=8)
        
        # Draw cell if present
        if atoms.cell is not None and atoms.pbc.any():
            cell = atoms.cell.array
            edges = [
                ([0, 0, 0], [1, 0, 0]), ([0, 0, 0], [0, 1, 0]), ([0, 0, 0], [0, 0, 1]),
                ([1, 0, 0], [1, 1, 0]), ([1, 0, 0], [1, 0, 1]),
                ([0, 1, 0], [1, 1, 0]), ([0, 1, 0], [0, 1, 1]),
                ([0, 0, 1], [1, 0, 1]), ([0, 0, 1], [0, 1, 1]),
                ([1, 1, 0], [1, 1, 1]), ([1, 0, 1], [1, 1, 1]), ([0, 1, 1], [1, 1, 1])
            ]
            
            for start, end in edges:
                start_pt = np.dot(start, cell)
                end_pt = np.dot(end, cell)
                ax.plot([start_pt[0], end_pt[0]], [start_pt[1], end_pt[1]], 
                        [start_pt[2], end_pt[2]], 'k-', linewidth=0.5)
        
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        
        self.figure.tight_layout()
        self.canvas.draw()
    
    def _export_structure(self):
        """Export the current structure."""
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure to export")
            return
        
        fmt = self.export_combo.currentText()
        formula = atoms.get_chemical_formula()
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Structure",
            os.path.expanduser(f"~/{formula}.{fmt}"),
            f"Structure Files (*.{fmt});;All Files (*)"
        )
        
        if file_path:
            try:
                ase_io.write(file_path, atoms, format=fmt)
                self.results_label.setText(f"✅ Exported to: {file_path}")
                self.results_label.setStyleSheet("color: green;")
            except Exception as e:
                self.results_label.setText(f"❌ Error exporting: {e}")
                self.results_label.setStyleSheet("color: red;")
    
    def refresh(self):
        """Refresh the page."""
        # Use loading guard to prevent infinite loops
        if self._loading:
            return
        self._loading = True
        try:
            atoms = self.session_state.get('current_structure')
            if atoms is not None:
                # Don't call _set_structure as it would update session state again
                # Just refresh the UI with existing structure
                self.current_group.setVisible(True)
                formula = atoms.get_chemical_formula()
                natoms = len(atoms)
                source = self.session_state.get('structure_source', 'Unknown')
                
                self.structure_info_label.setText(f"✅ Structure loaded: {formula}")
                self.formula_label.setText(f"<b>Formula:</b> {formula}")
                self.natoms_label.setText(f"<b>Atoms:</b> {natoms}")
                self.source_label.setText(f"<b>Source:</b> {source}")
                
                self._update_structure_info(atoms)
                self._refresh_visualization()
        finally:
            self._loading = False

"""
Structure Page for xespresso PySide6 GUI.

This page handles loading and visualizing atomic structures,
including ASE database operations for saving and loading structures.

The page is organized into tabs:
- Upload File: Load structures from files
- Build Structure: Build bulk crystals or molecules
- ASE Database: Load/save structures to/from database
- View Structure: Visualize the current structure
"""

import os

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTableWidget, QTableWidgetItem,
    QHeaderView, QFileDialog, QTabWidget, QTextEdit, QCheckBox,
    QDoubleSpinBox, QSpinBox, QSlider, QRadioButton, QButtonGroup
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
    from ase.db import connect as ase_db_connect
    ASE_DB_AVAILABLE = True
except ImportError:
    ASE_DB_AVAILABLE = False

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

# Default database path
DEFAULT_DB_PATH = os.path.expanduser("~/.xespresso/structures.db")


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
        
        # Header
        header_label = QLabel("<h2>🔬 Structure</h2>")
        header_label.setTextFormat(Qt.RichText)
        main_layout.addWidget(header_label)
        
        description = QLabel("""
<p>Load, build, and visualize atomic structures. Use the tabs below to load or create structures, 
then view them in the "View Structure" tab.</p>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        main_layout.addWidget(description)
        
        if not ASE_AVAILABLE:
            error_label = QLabel("❌ ASE not available. Structure viewing is disabled.")
            error_label.setStyleSheet("color: red; font-weight: bold;")
            main_layout.addWidget(error_label)
            return
        
        # Current Structure Info (always visible at the top)
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
        main_layout.addWidget(self.current_group)
        
        # Main tabs
        self.tabs = QTabWidget()
        
        # Upload File Tab
        self._create_upload_tab()
        
        # Build Structure Tab
        self._create_build_tab()
        
        # ASE Database Tab
        self._create_database_tab()
        
        # View Structure Tab (visualization)
        self._create_view_tab()
        
        main_layout.addWidget(self.tabs)
        
        # Results area (at the bottom, outside tabs)
        self.results_label = QLabel("")
        self.results_label.setWordWrap(True)
        main_layout.addWidget(self.results_label)
    
    def _create_upload_tab(self):
        """Create the upload file tab."""
        upload_tab = QWidget()
        upload_layout = QVBoxLayout(upload_tab)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # File upload section
        upload_group = QGroupBox("📂 Upload Structure File")
        upload_group_layout = QVBoxLayout(upload_group)
        
        upload_btn_layout = QHBoxLayout()
        upload_btn = QPushButton("📂 Open Structure File")
        upload_btn.clicked.connect(self._open_structure_file)
        upload_btn_layout.addWidget(upload_btn)
        
        self.file_path_label = QLabel("")
        upload_btn_layout.addWidget(self.file_path_label, 1)
        upload_group_layout.addLayout(upload_btn_layout)
        
        upload_info = QLabel("Supported formats: CIF, XYZ, PDB, VASP, POSCAR, TRAJ, JSON")
        upload_info.setStyleSheet("color: gray;")
        upload_group_layout.addWidget(upload_info)
        
        scroll_layout.addWidget(upload_group)
        
        # Save to database section
        if ASE_DB_AVAILABLE:
            self._create_save_to_db_section(scroll_layout, "upload")
        
        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_widget)
        upload_layout.addWidget(scroll_area)
        
        self.tabs.addTab(upload_tab, "📂 Upload File")
    
    def _create_build_tab(self):
        """Create the build structure tab."""
        build_tab = QWidget()
        build_layout = QVBoxLayout(build_tab)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Structure type selector
        type_group = QGroupBox("🔧 Structure Type")
        type_layout = QHBoxLayout(type_group)
        type_layout.addWidget(QLabel("Structure Type:"))
        self.build_type_combo = QComboBox()
        self.build_type_combo.addItems(["Bulk Crystal", "Molecule"])
        self.build_type_combo.currentTextChanged.connect(self._on_build_type_changed)
        type_layout.addWidget(self.build_type_combo, 1)
        scroll_layout.addWidget(type_group)
        
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
        
        scroll_layout.addWidget(self.bulk_group)
        
        # Molecule options
        self.molecule_group = QGroupBox("🧪 Build Molecule")
        molecule_layout = QFormLayout(self.molecule_group)
        
        self.molecule_edit = QLineEdit("H2O")
        self.molecule_edit.setPlaceholderText("H2O, CO2, CH4, NH3, C6H6, etc.")
        molecule_layout.addRow("Molecule Name:", self.molecule_edit)
        
        molecule_info = QLabel("💡 Tip: Try H2O, CO2, CH4, NH3, C6H6, or other common molecules")
        molecule_info.setStyleSheet("color: gray;")
        molecule_layout.addRow(molecule_info)
        
        build_molecule_btn = QPushButton("🔨 Build Molecule")
        build_molecule_btn.clicked.connect(self._build_molecule)
        molecule_layout.addRow(build_molecule_btn)
        
        self.molecule_group.setVisible(False)
        scroll_layout.addWidget(self.molecule_group)
        
        # Save to database section
        if ASE_DB_AVAILABLE:
            self._create_save_to_db_section(scroll_layout, "build")
        
        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_widget)
        build_layout.addWidget(scroll_area)
        
        self.tabs.addTab(build_tab, "🔨 Build Structure")
    
    def _create_save_to_db_section(self, parent_layout, prefix):
        """Create the save to database section (reusable for upload and build tabs)."""
        save_group = QGroupBox("💾 Save to ASE Database")
        save_layout = QFormLayout(save_group)
        
        save_info = QLabel("After loading/building a structure, you can save it to the database for quick access.")
        save_info.setWordWrap(True)
        save_info.setStyleSheet("color: gray;")
        save_layout.addRow(save_info)
        
        # Database path
        db_path_layout = QHBoxLayout()
        db_path_edit = QLineEdit()
        db_path_edit.setText(DEFAULT_DB_PATH)
        db_path_edit.setToolTip("Path to ASE database file")
        db_path_layout.addWidget(db_path_edit)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(lambda: self._browse_db_path_for(db_path_edit))
        db_path_layout.addWidget(browse_btn)
        save_layout.addRow("Database Path:", db_path_layout)
        
        # Structure name
        name_edit = QLineEdit()
        name_edit.setPlaceholderText("e.g., my_structure, bulk_Fe")
        save_layout.addRow("Structure Name:", name_edit)
        
        # Tags
        tags_edit = QLineEdit()
        tags_edit.setPlaceholderText("e.g., bulk, metal, test (comma-separated)")
        save_layout.addRow("Tags:", tags_edit)
        
        # Save button
        save_btn = QPushButton("💾 Save Current Structure to Database")
        save_btn.clicked.connect(lambda: self._save_structure_to_db(db_path_edit, name_edit, tags_edit))
        save_layout.addRow(save_btn)
        
        # Status label
        status_label = QLabel("")
        status_label.setWordWrap(True)
        save_layout.addRow(status_label)
        
        # Store references for later use
        setattr(self, f'{prefix}_db_path_edit', db_path_edit)
        setattr(self, f'{prefix}_name_edit', name_edit)
        setattr(self, f'{prefix}_tags_edit', tags_edit)
        setattr(self, f'{prefix}_save_status', status_label)
        
        parent_layout.addWidget(save_group)
    
    def _create_database_tab(self):
        """Create the ASE database tab."""
        db_tab = QWidget()
        db_layout = QVBoxLayout(db_tab)
        
        if not ASE_DB_AVAILABLE:
            db_unavailable = QLabel("❌ ASE database module not available.")
            db_unavailable.setStyleSheet("color: red;")
            db_layout.addWidget(db_unavailable)
            db_layout.addStretch()
            self.tabs.addTab(db_tab, "📚 ASE Database")
            return
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        db_desc = QLabel("""
<p><b>ASE Database:</b> Load structures from your database library.</p>
<p>💡 <b>Tip:</b> Use the Upload and Build tabs to add new structures to your database.</p>
""")
        db_desc.setTextFormat(Qt.RichText)
        db_desc.setWordWrap(True)
        scroll_layout.addWidget(db_desc)
        
        # Database path configuration
        db_path_group = QGroupBox("📁 Database Location")
        db_path_layout = QHBoxLayout(db_path_group)
        db_path_layout.addWidget(QLabel("Database Path:"))
        self.db_path_edit = QLineEdit()
        self.db_path_edit.setText(DEFAULT_DB_PATH)
        self.db_path_edit.setToolTip("Path to ASE database file")
        db_path_layout.addWidget(self.db_path_edit, 1)
        
        db_browse_btn = QPushButton("Browse...")
        db_browse_btn.clicked.connect(self._browse_db_path)
        db_path_layout.addWidget(db_browse_btn)
        scroll_layout.addWidget(db_path_group)
        
        # Load from database section
        load_group = QGroupBox("📥 Load from Database")
        load_layout = QVBoxLayout(load_group)
        
        refresh_db_btn = QPushButton("🔄 Refresh Database List")
        refresh_db_btn.clicked.connect(self._refresh_db_list)
        load_layout.addWidget(refresh_db_btn)
        
        self.db_status_label = QLabel("")
        load_layout.addWidget(self.db_status_label)
        
        self.db_structures_table = QTableWidget()
        self.db_structures_table.setColumnCount(5)
        self.db_structures_table.setHorizontalHeaderLabels(["ID", "Name", "Formula", "Atoms", "Tags"])
        self.db_structures_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.db_structures_table.setMaximumHeight(250)
        self.db_structures_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.db_structures_table.setSelectionMode(QTableWidget.SingleSelection)
        load_layout.addWidget(self.db_structures_table)
        
        load_selected_btn = QPushButton("📥 Load Selected Structure")
        load_selected_btn.clicked.connect(self._load_from_database)
        load_layout.addWidget(load_selected_btn)
        
        self.db_load_result = QLabel("")
        self.db_load_result.setWordWrap(True)
        load_layout.addWidget(self.db_load_result)
        
        scroll_layout.addWidget(load_group)
        
        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_widget)
        db_layout.addWidget(scroll_area)
        
        self.tabs.addTab(db_tab, "📚 ASE Database")
    
    def _create_view_tab(self):
        """Create the view structure tab (visualization)."""
        view_tab = QWidget()
        view_layout = QVBoxLayout(view_tab)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        view_desc = QLabel("""
<p><b>Structure Visualization:</b> View the currently loaded structure.</p>
<p>Load or build a structure using the other tabs, then come here to visualize it.</p>
""")
        view_desc.setTextFormat(Qt.RichText)
        view_desc.setWordWrap(True)
        scroll_layout.addWidget(view_desc)
        
        # Visualization Section
        viz_group = QGroupBox("🔬 3D Structure Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        if MATPLOTLIB_AVAILABLE:
            self.figure = Figure(figsize=(8, 6))
            self.canvas = FigureCanvas(self.figure)
            viz_layout.addWidget(self.canvas)
            
            # Visualization controls
            viz_controls = QHBoxLayout()
            
            refresh_viz_btn = QPushButton("🔄 Refresh View")
            refresh_viz_btn.clicked.connect(self._refresh_visualization)
            viz_controls.addWidget(refresh_viz_btn)
            
            viz_controls.addStretch()
            viz_layout.addLayout(viz_controls)
        else:
            viz_label = QLabel("Matplotlib not available. Visualization disabled.")
            viz_label.setStyleSheet("color: red;")
            viz_layout.addWidget(viz_label)
        
        scroll_layout.addWidget(viz_group)
        
        # Structure Information
        info_group = QGroupBox("📋 Structure Information")
        info_layout = QVBoxLayout(info_group)
        
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(200)
        info_layout.addWidget(self.info_text)
        
        scroll_layout.addWidget(info_group)
        
        # Export Section
        export_group = QGroupBox("💾 Export Structure")
        export_layout = QHBoxLayout(export_group)
        
        export_layout.addWidget(QLabel("Export Format:"))
        self.export_combo = QComboBox()
        self.export_combo.addItems(["cif", "xyz", "pdb", "vasp", "json"])
        export_layout.addWidget(self.export_combo)
        
        export_btn = QPushButton("⬇️ Export Structure")
        export_btn.clicked.connect(self._export_structure)
        export_layout.addWidget(export_btn)
        
        export_layout.addStretch()
        scroll_layout.addWidget(export_group)
        
        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_widget)
        view_layout.addWidget(scroll_area)
        
        # Store reference to view tab for later use
        self.view_tab = view_tab
        self.tabs.addTab(view_tab, "🔬 View Structure")
    
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
                # Save the full file path for session restoration
                self.session_state['structure_file_path'] = file_path
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
            
            # Update structure info in view tab
            self._update_structure_info(atoms)
            
            # Update visualization in view tab
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
    
    def _browse_db_path(self):
        """Browse for database file path."""
        current_path = self.db_path_edit.text() or DEFAULT_DB_PATH
        current_dir = os.path.dirname(current_path)
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Database File",
            current_dir,
            "ASE Database Files (*.db);;All Files (*)"
        )
        
        if file_path:
            self.db_path_edit.setText(file_path)
    
    def _browse_db_path_for(self, line_edit):
        """Browse for database file path for a specific line edit."""
        current_path = line_edit.text() or DEFAULT_DB_PATH
        current_dir = os.path.dirname(current_path)
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Database File",
            current_dir,
            "ASE Database Files (*.db);;All Files (*)"
        )
        
        if file_path:
            line_edit.setText(file_path)
    
    def _refresh_db_list(self):
        """Refresh the database structures list."""
        if not ASE_DB_AVAILABLE:
            return
        
        db_path = self.db_path_edit.text().strip()
        if not db_path:
            self.db_status_label.setText("⚠️ Please enter a database path")
            self.db_status_label.setStyleSheet("color: orange;")
            return
        
        # Normalize the path
        db_path = os.path.abspath(os.path.expanduser(db_path))
        
        if not os.path.exists(db_path):
            self.db_status_label.setText("ℹ️ Database does not exist yet. It will be created when you save your first structure.")
            self.db_status_label.setStyleSheet("color: blue;")
            self.db_structures_table.setRowCount(0)
            return
        
        try:
            db = ase_db_connect(db_path)
            rows = list(db.select())
            
            if rows:
                self.db_status_label.setText(f"✅ Found {len(rows)} structure(s) in database")
                self.db_status_label.setStyleSheet("color: green;")
                
                self.db_structures_table.setRowCount(len(rows))
                for i, row in enumerate(rows):
                    self.db_structures_table.setItem(i, 0, QTableWidgetItem(str(row.id)))
                    
                    # Get name from key_value_pairs if it exists
                    name = row.key_value_pairs.get('name', '') if row.key_value_pairs else ''
                    self.db_structures_table.setItem(i, 1, QTableWidgetItem(name))
                    
                    self.db_structures_table.setItem(i, 2, QTableWidgetItem(row.formula))
                    self.db_structures_table.setItem(i, 3, QTableWidgetItem(str(row.natoms)))
                    
                    # Get tags from key_value_pairs (excluding 'name' as it's now in its own column)
                    tags_list = [k for k in row.key_value_pairs.keys() if k != 'name'] if row.key_value_pairs else []
                    tags = ", ".join(tags_list)
                    self.db_structures_table.setItem(i, 4, QTableWidgetItem(tags))
            else:
                self.db_status_label.setText("ℹ️ Database is empty. Use Upload or Build tabs to add structures.")
                self.db_status_label.setStyleSheet("color: blue;")
                self.db_structures_table.setRowCount(0)
                
        except Exception as e:
            self.db_status_label.setText(f"❌ Error reading database: {e}")
            self.db_status_label.setStyleSheet("color: red;")
            self.db_structures_table.setRowCount(0)
    
    def _load_from_database(self):
        """Load a structure from the database."""
        if not ASE_DB_AVAILABLE:
            return
        
        selected_rows = self.db_structures_table.selectedItems()
        if not selected_rows:
            QMessageBox.warning(self, "Warning", "Please select a structure to load")
            return
        
        # Get the ID from the first column of the selected row
        selected_row = self.db_structures_table.currentRow()
        id_item = self.db_structures_table.item(selected_row, 0)
        if not id_item:
            return
        
        try:
            structure_id = int(id_item.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid structure ID")
            return
        
        db_path = os.path.abspath(os.path.expanduser(self.db_path_edit.text().strip()))
        
        try:
            db = ase_db_connect(db_path)
            row = db.get(id=structure_id)
            atoms = row.toatoms()
            
            self._set_structure(atoms, f"Database: ID {structure_id}")
            # Save the database path for session restoration
            self.session_state['structure_db_path'] = db_path
            
            self.db_load_result.setText(f"✅ Loaded structure ID {structure_id}: {atoms.get_chemical_formula()}")
            self.db_load_result.setStyleSheet("color: green;")
            self.results_label.setText(f"✅ Loaded from database: {atoms.get_chemical_formula()}")
            self.results_label.setStyleSheet("color: green;")
            
            # Switch to the view tab after loading
            if hasattr(self, 'view_tab'):
                self.tabs.setCurrentWidget(self.view_tab)
            
        except Exception as e:
            self.db_load_result.setText(f"❌ Error loading structure: {e}")
            self.db_load_result.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error loading structure:\n{e}")
    
    def _save_structure_to_db(self, db_path_edit, name_edit, tags_edit):
        """Save the current structure to the database."""
        if not ASE_DB_AVAILABLE:
            return
        
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded. Load or build a structure first.")
            return
        
        db_path = db_path_edit.text().strip()
        if not db_path:
            QMessageBox.warning(self, "Warning", "Please enter a database path")
            return
        
        # Normalize the path
        db_path = os.path.abspath(os.path.expanduser(db_path))
        
        try:
            # Create database directory if it doesn't exist
            db_dir = os.path.dirname(db_path)
            if db_dir:
                os.makedirs(db_dir, exist_ok=True)
            
            db = ase_db_connect(db_path)
            
            # Prepare key-value pairs
            key_value_pairs = {}
            
            # Add name
            save_name = name_edit.text().strip()
            if save_name:
                key_value_pairs['name'] = save_name
            
            # Add source info
            source = self.session_state.get('structure_source', '')
            if source:
                key_value_pairs['source'] = source
            
            # Parse and add tags
            tags = tags_edit.text().strip()
            if tags:
                for tag in tags.split(','):
                    tag = tag.strip()
                    if tag:
                        key_value_pairs[tag] = True
            
            # Save to database
            db.write(atoms, **key_value_pairs)
            
            self.results_label.setText(f"✅ Structure saved to database: {db_path}")
            self.results_label.setStyleSheet("color: green;")
            
            QMessageBox.information(
                self, "Success",
                f"Structure saved to database!\n\n"
                f"Formula: {atoms.get_chemical_formula()}\n"
                f"Database: {db_path}\n\n"
                "💡 Go to the 'ASE Database' tab to see the updated list."
            )
            
            # Clear the name and tags
            name_edit.clear()
            tags_edit.clear()
            
        except Exception as e:
            self.results_label.setText(f"❌ Error saving to database: {e}")
            self.results_label.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error saving to database:\n{e}")

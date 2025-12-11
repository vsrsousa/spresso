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
    QDoubleSpinBox, QSpinBox, QSlider, QRadioButton, QButtonGroup,
    QDialog, QDialogButtonBox
)
from PySide6.QtCore import Qt

# Required imports - these are always available as they're in requirements.txt
from ase import io as ase_io
from ase import Atoms
from ase.build import bulk, molecule
from ase.visualize import view as ase_view

import matplotlib
# Only set backend if not already set to avoid conflicts
if matplotlib.get_backend() != 'QtAgg':
    matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np

# Optional import - ase.db may not be available in all installations
try:
    from ase.db import connect as ase_db_connect
    ASE_DB_AVAILABLE = True
except ImportError:
    ASE_DB_AVAILABLE = False

# Default database path
DEFAULT_DB_PATH = os.path.expanduser("~/.xespresso/structures.db")


class StructureViewerPage(QWidget):
    """Structure viewer page widget."""
    
    # Viewer type constants
    VIEWER_INTERACTIVE = "Interactive 3D"
    VIEWER_SIMPLE = "Simple 3D"
    
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
        
        db_buttons_layout = QHBoxLayout()
        
        load_selected_btn = QPushButton("📥 Load Selected Structure")
        load_selected_btn.clicked.connect(self._load_from_database)
        db_buttons_layout.addWidget(load_selected_btn)
        
        edit_selected_btn = QPushButton("✏️ Edit Properties")
        edit_selected_btn.clicked.connect(self._edit_database_entry)
        db_buttons_layout.addWidget(edit_selected_btn)
        
        delete_selected_btn = QPushButton("🗑️ Delete Entry")
        delete_selected_btn.clicked.connect(self._delete_database_entry)
        db_buttons_layout.addWidget(delete_selected_btn)
        
        load_layout.addLayout(db_buttons_layout)
        
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
        
        # Viewer selection and external viewer button
        viewer_select_layout = QHBoxLayout()
        viewer_select_layout.addWidget(QLabel("Embedded Viewer:"))
        self.viewer_combo = QComboBox()
        self.viewer_combo.addItems([self.VIEWER_INTERACTIVE, self.VIEWER_SIMPLE])
        self.viewer_combo.currentTextChanged.connect(self._on_viewer_changed)
        viewer_select_layout.addWidget(self.viewer_combo)
        
        # Add external ASE GUI button
        self.ase_gui_btn = QPushButton("🚀 Open in ASE GUI")
        self.ase_gui_btn.setToolTip("Open structure in external ASE GUI window")
        self.ase_gui_btn.clicked.connect(self._open_ase_gui)
        viewer_select_layout.addWidget(self.ase_gui_btn)
        
        viewer_select_layout.addStretch()
        scroll_layout.addLayout(viewer_select_layout)
        
        # Visualization Section
        viz_group = QGroupBox("🔬 3D Structure Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        # Container for different viewers
        self.viewer_container = QWidget()
        self.viewer_container_layout = QVBoxLayout(self.viewer_container)
        self.viewer_container_layout.setContentsMargins(0, 0, 0, 0)
        
        # Create matplotlib viewers (always available)
        # Interactive matplotlib viewer with toolbar
        self.interactive_matplotlib_widget = QWidget()
        interactive_layout = QVBoxLayout(self.interactive_matplotlib_widget)
        interactive_layout.setContentsMargins(0, 0, 0, 0)
        
        self.interactive_figure = Figure(figsize=(8, 6))
        self.interactive_canvas = FigureCanvas(self.interactive_figure)
        
        # Add navigation toolbar for interactivity
        self.toolbar = NavigationToolbar(self.interactive_canvas, self.interactive_matplotlib_widget)
        interactive_layout.addWidget(self.toolbar)
        interactive_layout.addWidget(self.interactive_canvas)
        
        # Visualization controls
        viz_controls = QHBoxLayout()
        refresh_viz_btn = QPushButton("🔄 Refresh View")
        refresh_viz_btn.clicked.connect(self._refresh_visualization)
        viz_controls.addWidget(refresh_viz_btn)
        viz_controls.addStretch()
        interactive_layout.addLayout(viz_controls)
        
        # Simple matplotlib viewer without toolbar
        self.simple_matplotlib_widget = QWidget()
        simple_layout = QVBoxLayout(self.simple_matplotlib_widget)
        simple_layout.setContentsMargins(0, 0, 0, 0)
        
        self.simple_figure = Figure(figsize=(8, 6))
        self.simple_canvas = FigureCanvas(self.simple_figure)
        simple_layout.addWidget(self.simple_canvas)
        
        # Simple controls
        simple_controls = QHBoxLayout()
        simple_refresh_btn = QPushButton("🔄 Refresh View")
        simple_refresh_btn.clicked.connect(self._refresh_visualization)
        simple_controls.addWidget(simple_refresh_btn)
        simple_controls.addStretch()
        simple_layout.addLayout(simple_controls)
        
        # Add container to viz_layout
        viz_layout.addWidget(self.viewer_container)
        
        # Initialize with default viewer
        self._on_viewer_changed(self.viewer_combo.currentText())
        
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
                from ..utils import read_structure
                atoms = read_structure(file_path)
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
    
    def _on_viewer_changed(self, viewer_name):
        """Handle viewer selection change."""
        atoms = self.session_state.get('current_structure')
        
        # Clear current viewer
        while self.viewer_container_layout.count():
            item = self.viewer_container_layout.takeAt(0)
            if item.widget():
                item.widget().setParent(None)
        
        # Show selected viewer
        if viewer_name == self.VIEWER_INTERACTIVE:
            self._show_interactive_matplotlib_viewer(atoms)
        elif viewer_name == self.VIEWER_SIMPLE:
            self._show_simple_matplotlib_viewer(atoms)
    
    def _show_interactive_matplotlib_viewer(self, atoms):
        """Show interactive matplotlib viewer with navigation toolbar."""
        if self.interactive_matplotlib_widget:
            self.viewer_container_layout.addWidget(self.interactive_matplotlib_widget)
            if atoms is not None:
                self._refresh_interactive_matplotlib_visualization()
    
    def _show_simple_matplotlib_viewer(self, atoms):
        """Show simple matplotlib viewer without toolbar."""
        if self.simple_matplotlib_widget:
            self.viewer_container_layout.addWidget(self.simple_matplotlib_widget)
            if atoms is not None:
                self._refresh_simple_matplotlib_visualization()
    
    def _open_ase_gui(self):
        """Open structure in external ASE GUI window."""
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded. Load a structure first.")
            return
        
        try:
            # Open ASE GUI in external window
            ase_view(atoms)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error opening ASE GUI:\n{e}")
    
    def _refresh_visualization(self):
        """Refresh the current visualization."""
        if not hasattr(self, 'viewer_combo') or self.viewer_combo is None:
            return
        
        viewer_name = self.viewer_combo.currentText()
        
        if viewer_name == self.VIEWER_INTERACTIVE:
            self._refresh_interactive_matplotlib_visualization()
        elif viewer_name == self.VIEWER_SIMPLE:
            self._refresh_simple_matplotlib_visualization()
        else:
            # Default to first available viewer
            self._on_viewer_changed(viewer_name)
    
    def _refresh_interactive_matplotlib_visualization(self):
        """Refresh the interactive matplotlib structure visualization."""
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            return
        
        self._draw_structure_on_axes(self.interactive_figure, self.interactive_canvas, atoms)
    
    def _refresh_simple_matplotlib_visualization(self):
        """Refresh the simple matplotlib structure visualization."""
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            return
        
        self._draw_structure_on_axes(self.simple_figure, self.simple_canvas, atoms)
    
    def _draw_structure_on_axes(self, figure, canvas, atoms):
        """Draw structure on matplotlib axes."""
        figure.clear()
        ax = figure.add_subplot(111, projection='3d')
        
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
        
        figure.tight_layout()
        canvas.draw()
    
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
    
    def _edit_database_entry(self):
        """Edit properties (name, tags) of a database entry."""
        if not ASE_DB_AVAILABLE:
            return
        
        selected_rows = self.db_structures_table.selectedItems()
        if not selected_rows:
            QMessageBox.warning(self, "Warning", "Please select a structure to edit")
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
        
        if not os.path.exists(db_path):
            QMessageBox.warning(self, "Warning", "Database file not found")
            return
        
        try:
            db = ase_db_connect(db_path)
            row = db.get(id=structure_id)
            
            # Get current values
            current_name = row.key_value_pairs.get('name', '') if row.key_value_pairs else ''
            current_tags = [k for k in row.key_value_pairs.keys() if k != 'name' and k != 'source'] if row.key_value_pairs else []
            current_tags_str = ", ".join(current_tags)
            
            # Create edit dialog
            dialog = QDialog(self)
            dialog.setWindowTitle(f"Edit Database Entry (ID: {structure_id})")
            dialog.setMinimumWidth(400)
            
            layout = QVBoxLayout(dialog)
            
            info_label = QLabel(f"<b>Formula:</b> {row.formula}<br><b>Atoms:</b> {row.natoms}")
            info_label.setTextFormat(Qt.RichText)
            layout.addWidget(info_label)
            
            form_layout = QFormLayout()
            
            name_edit = QLineEdit(current_name)
            name_edit.setPlaceholderText("Structure name")
            form_layout.addRow("Name:", name_edit)
            
            tags_edit = QLineEdit(current_tags_str)
            tags_edit.setPlaceholderText("comma-separated tags")
            form_layout.addRow("Tags:", tags_edit)
            
            layout.addLayout(form_layout)
            
            button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            button_box.accepted.connect(dialog.accept)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)
            
            if dialog.exec() == QDialog.Accepted:
                # Update the entry
                new_name = name_edit.text().strip()
                new_tags = tags_edit.text().strip()
                
                # Start with existing key-value pairs to preserve unmodified data
                key_value_pairs = dict(row.key_value_pairs) if row.key_value_pairs else {}
                
                # Remove old tags (but keep name and source)
                keys_to_remove = [k for k in key_value_pairs.keys() if k not in ['name', 'source']]
                for k in keys_to_remove:
                    del key_value_pairs[k]
                
                # Update name
                if new_name:
                    key_value_pairs['name'] = new_name
                elif 'name' in key_value_pairs:
                    del key_value_pairs['name']
                
                # Add new tags
                if new_tags:
                    for tag in new_tags.split(','):
                        tag = tag.strip()
                        if tag:
                            key_value_pairs[tag] = True
                
                # Update the entry by updating its key_value_pairs
                db.update(id=structure_id, **key_value_pairs)
                
                self.db_load_result.setText(f"✅ Updated entry ID {structure_id}")
                self.db_load_result.setStyleSheet("color: green;")
                
                # Refresh the list
                self._refresh_db_list()
                
        except Exception as e:
            self.db_load_result.setText(f"❌ Error editing entry: {e}")
            self.db_load_result.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error editing entry:\n{e}")
    
    def _delete_database_entry(self):
        """Delete a structure from the database."""
        if not ASE_DB_AVAILABLE:
            return
        
        selected_rows = self.db_structures_table.selectedItems()
        if not selected_rows:
            QMessageBox.warning(self, "Warning", "Please select a structure to delete")
            return
        
        # Get the ID from the first column of the selected row
        selected_row = self.db_structures_table.currentRow()
        id_item = self.db_structures_table.item(selected_row, 0)
        name_item = self.db_structures_table.item(selected_row, 1)
        formula_item = self.db_structures_table.item(selected_row, 2)
        
        if not id_item:
            return
        
        try:
            structure_id = int(id_item.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid structure ID")
            return
        
        # Confirm deletion
        name = name_item.text() if name_item else ""
        formula = formula_item.text() if formula_item else ""
        display_name = name if name else formula
        
        reply = QMessageBox.question(
            self,
            "Confirm Deletion",
            f"Are you sure you want to delete this structure?\n\n"
            f"ID: {structure_id}\n"
            f"Name: {display_name}\n"
            f"Formula: {formula}\n\n"
            f"This action cannot be undone.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        
        if reply != QMessageBox.Yes:
            return
        
        db_path = os.path.abspath(os.path.expanduser(self.db_path_edit.text().strip()))
        
        if not os.path.exists(db_path):
            QMessageBox.warning(self, "Warning", "Database file not found")
            return
        
        try:
            db = ase_db_connect(db_path)
            db.delete([structure_id])
            
            self.db_load_result.setText(f"✅ Deleted entry ID {structure_id}")
            self.db_load_result.setStyleSheet("color: green;")
            
            # Refresh the list
            self._refresh_db_list()
            
        except Exception as e:
            self.db_load_result.setText(f"❌ Error deleting entry: {e}")
            self.db_load_result.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error deleting entry:\n{e}")

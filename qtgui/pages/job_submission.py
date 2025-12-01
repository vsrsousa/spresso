"""
Job Submission Page for xespresso PySide6 GUI.

This page handles file browsing, dry run, and job submission.
"""

import os

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTextEdit,
    QTabWidget, QTreeWidget, QTreeWidgetItem, QFileDialog,
    QHeaderView, QSplitter
)
from PySide6.QtCore import Qt

from qtgui.utils import validate_path_under_base, safe_makedirs

try:
    from ase import Atoms
    from ase import io as ase_io
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


class JobSubmissionPage(QWidget):
    """Job submission page widget."""
    
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self._loading = False  # Guard to prevent recursive updates
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the user interface."""
        main_layout = QVBoxLayout(self)
        
        # Header
        header_label = QLabel("<h2>🚀 Job Submission & File Management</h2>")
        header_label.setTextFormat(Qt.RichText)
        main_layout.addWidget(header_label)
        
        description = QLabel("""
<p><b>Generate calculation files, browse directories, and run calculations using xespresso.</b></p>
<ul>
<li><b>Generate Files (Dry Run)</b>: Creates Espresso calculator and generates input files</li>
<li><b>File Browser</b>: Browse, view, and edit existing calculation files</li>
<li><b>Run Calculation</b>: Creates Espresso calculator and runs with get_potential_energy()</li>
</ul>
""")
        description.setTextFormat(Qt.RichText)
        description.setWordWrap(True)
        main_layout.addWidget(description)
        
        # Tabs
        tabs = QTabWidget()
        
        # File Browser Tab
        browser_tab = QWidget()
        browser_layout = QVBoxLayout(browser_tab)
        self._setup_browser_tab(browser_layout)
        tabs.addTab(browser_tab, "📂 File Browser")
        
        # Dry Run Tab
        dry_run_tab = QWidget()
        dry_run_layout = QVBoxLayout(dry_run_tab)
        self._setup_dry_run_tab(dry_run_layout)
        tabs.addTab(dry_run_tab, "🧪 Generate Files (Dry Run)")
        
        # Run Calculation Tab
        run_tab = QWidget()
        run_layout = QVBoxLayout(run_tab)
        self._setup_run_tab(run_layout)
        tabs.addTab(run_tab, "🚀 Run Calculation")
        
        main_layout.addWidget(tabs)
    
    def _setup_browser_tab(self, layout):
        """Setup the file browser tab."""
        # Working directory
        workdir_layout = QHBoxLayout()
        workdir_layout.addWidget(QLabel("📁 Working Directory:"))
        self.workdir_label = QLabel("")
        workdir_layout.addWidget(self.workdir_label, 1)
        
        refresh_btn = QPushButton("🔄 Refresh")
        refresh_btn.clicked.connect(self._refresh_browser)
        workdir_layout.addWidget(refresh_btn)
        layout.addLayout(workdir_layout)
        
        # Splitter for tree and content
        splitter = QSplitter(Qt.Horizontal)
        
        # File tree
        tree_widget = QWidget()
        tree_layout = QVBoxLayout(tree_widget)
        tree_layout.setContentsMargins(0, 0, 0, 0)
        
        tree_layout.addWidget(QLabel("<b>Calculation Folders</b>"))
        self.file_tree = QTreeWidget()
        self.file_tree.setHeaderLabels(["Name"])
        self.file_tree.itemClicked.connect(self._on_file_selected)
        tree_layout.addWidget(self.file_tree)
        
        splitter.addWidget(tree_widget)
        
        # File content
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        content_layout.setContentsMargins(0, 0, 0, 0)
        
        content_layout.addWidget(QLabel("<b>File Content</b>"))
        self.file_content = QTextEdit()
        self.file_content.setReadOnly(True)
        content_layout.addWidget(self.file_content)
        
        # File info
        self.file_info_label = QLabel("")
        content_layout.addWidget(self.file_info_label)
        
        splitter.addWidget(content_widget)
        splitter.setSizes([300, 700])
        
        layout.addWidget(splitter)
        
        # Update browser
        self._refresh_browser()
    
    def _setup_dry_run_tab(self, layout):
        """Setup the dry run tab."""
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Description
        desc = QLabel("""
<p><b>Generate Calculation Files (Dry Run)</b></p>
<p>Generate input files and job scripts for your calculation <b>without submitting the job</b>.</p>
<p>This is useful for testing your configuration before submission.</p>
""")
        desc.setTextFormat(Qt.RichText)
        desc.setWordWrap(True)
        scroll_layout.addWidget(desc)
        
        # Status
        self.dry_run_status = QLabel("")
        self.dry_run_status.setWordWrap(True)
        scroll_layout.addWidget(self.dry_run_status)
        
        # Configuration Summary
        config_group = QGroupBox("📋 Current Configuration")
        config_layout = QVBoxLayout(config_group)
        
        self.config_summary = QTextEdit()
        self.config_summary.setReadOnly(True)
        self.config_summary.setMaximumHeight(150)
        config_layout.addWidget(self.config_summary)
        
        scroll_layout.addWidget(config_group)
        
        # Output Directory
        output_group = QGroupBox("📁 Output Directory")
        output_layout = QFormLayout(output_group)
        
        self.output_dir_label = QLabel("")
        output_layout.addRow("Base Directory:", self.output_dir_label)
        
        self.label_edit = QLineEdit()
        self.label_edit.setPlaceholderText("e.g., scf/Al")
        output_layout.addRow("Calculation Label:", self.label_edit)
        
        self.full_path_label = QLabel("")
        output_layout.addRow("Full Path:", self.full_path_label)
        
        scroll_layout.addWidget(output_group)
        
        # Generate Button
        generate_btn = QPushButton("🧪 Generate Files (Dry Run)")
        generate_btn.clicked.connect(self._generate_files)
        scroll_layout.addWidget(generate_btn)
        
        # Generated Files Preview
        preview_group = QGroupBox("📄 Generated Files")
        preview_layout = QVBoxLayout(preview_group)
        
        self.generated_files_list = QTextEdit()
        self.generated_files_list.setReadOnly(True)
        self.generated_files_list.setMaximumHeight(100)
        preview_layout.addWidget(self.generated_files_list)
        
        self.input_preview = QTextEdit()
        self.input_preview.setReadOnly(True)
        preview_layout.addWidget(self.input_preview)
        
        scroll_layout.addWidget(preview_group)
        
        # Results
        self.dry_run_results = QLabel("")
        self.dry_run_results.setWordWrap(True)
        scroll_layout.addWidget(self.dry_run_results)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)
        
        # Update configuration display
        self._update_dry_run_config()
    
    def _setup_run_tab(self, layout):
        """Setup the run calculation tab."""
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.NoFrame)
        
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # Description
        desc = QLabel("""
<p><b>Run Calculation</b></p>
<p>Run a calculation using xespresso by calling <code>calc.get_potential_energy()</code>.</p>
<p>This will:</p>
<ul>
<li>Create an Espresso calculator from your configuration</li>
<li>Check for previous calculation results and convergence</li>
<li>Automatically generate input files if they don't exist</li>
<li>Run the calculation and return the energy</li>
</ul>
""")
        desc.setTextFormat(Qt.RichText)
        desc.setWordWrap(True)
        scroll_layout.addWidget(desc)
        
        # Status
        self.run_status = QLabel("")
        self.run_status.setWordWrap(True)
        scroll_layout.addWidget(self.run_status)
        
        # Configuration Summary
        config_group = QGroupBox("📋 Current Configuration")
        config_layout = QVBoxLayout(config_group)
        
        self.run_config_summary = QTextEdit()
        self.run_config_summary.setReadOnly(True)
        self.run_config_summary.setMaximumHeight(150)
        config_layout.addWidget(self.run_config_summary)
        
        scroll_layout.addWidget(config_group)
        
        # Output Location
        output_group = QGroupBox("📁 Output Location")
        output_layout = QFormLayout(output_group)
        
        self.run_output_dir_label = QLabel("")
        output_layout.addRow("Base Directory:", self.run_output_dir_label)
        
        self.run_label_edit = QLineEdit()
        self.run_label_edit.setPlaceholderText("e.g., scf/Al")
        output_layout.addRow("Calculation Label:", self.run_label_edit)
        
        self.run_full_path_label = QLabel("")
        output_layout.addRow("Full Path:", self.run_full_path_label)
        
        scroll_layout.addWidget(output_group)
        
        # Run Button
        run_btn = QPushButton("🚀 Run Calculation")
        run_btn.clicked.connect(self._run_calculation)
        scroll_layout.addWidget(run_btn)
        
        # Results
        results_group = QGroupBox("📊 Results")
        results_layout = QVBoxLayout(results_group)
        
        self.energy_label = QLabel("")
        self.energy_label.setStyleSheet("font-size: 16pt; font-weight: bold;")
        results_layout.addWidget(self.energy_label)
        
        self.run_results = QTextEdit()
        self.run_results.setReadOnly(True)
        results_layout.addWidget(self.run_results)
        
        scroll_layout.addWidget(results_group)
        
        scroll_layout.addStretch()
        
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)
        
        # Update configuration display
        self._update_run_config()
    
    def _refresh_browser(self):
        """Refresh the file browser.
        
        This method scans the working directory for calculation folders.
        To prevent blocking the main thread during large directory scans,
        we limit the depth and number of directories scanned.
        """
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        self.workdir_label.setText(workdir)
        
        self.file_tree.clear()
        
        if not os.path.exists(workdir):
            return
        
        # Find calculation folders with limits to prevent blocking
        input_extensions = (".in", ".pwi", ".phi", ".ppi", ".bandi")
        output_extensions = (".sh", ".slurm", ".out", ".pwo")
        max_dirs_scanned = 100  # Limit to prevent blocking on large directories
        dirs_scanned = 0
        
        try:
            for root, dirs, files in os.walk(workdir, topdown=True):
                # Count directories to prevent scanning too many
                dirs_scanned += 1
                if dirs_scanned > max_dirs_scanned:
                    self.file_info_label.setText(
                        f"⚠️ Directory scan limited to {max_dirs_scanned} folders. "
                        "Navigate to a specific calculation folder for better results."
                    )
                    break
                
                depth = root[len(workdir):].count(os.sep)
                if depth < 4:
                    # Prune hidden directories and common large directories to speed up walk
                    dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ('node_modules', '__pycache__', '.git', 'venv', 'env')]
                    
                    has_input = any(
                        f.endswith(input_extensions) or 
                        f == "job_file" or 
                        f.endswith((".sh", ".slurm"))
                        for f in files
                    )
                    
                    if has_input:
                        rel_path = os.path.relpath(root, workdir)
                        item = QTreeWidgetItem([rel_path])
                        item.setData(0, Qt.UserRole, root)
                        
                        # Add files as children
                        for f in files:
                            if (f.endswith(input_extensions) or 
                                f == "job_file" or 
                                f.endswith(output_extensions)):
                                child = QTreeWidgetItem([f])
                                child.setData(0, Qt.UserRole, os.path.join(root, f))
                                item.addChild(child)
                        
                        self.file_tree.addTopLevelItem(item)
                else:
                    # Don't descend into directories beyond depth 4
                    dirs[:] = []
        except (OSError, IOError) as e:
            self.file_info_label.setText(f"Error scanning: {e}")
    
    def _on_file_selected(self, item, column):
        """Handle file selection."""
        path = item.data(0, Qt.UserRole)
        
        if path and os.path.isfile(path):
            try:
                with open(path, 'r') as f:
                    content = f.read()
                
                self.file_content.setText(content)
                
                stat = os.stat(path)
                info = f"Size: {stat.st_size} bytes | Lines: {len(content.splitlines())}"
                self.file_info_label.setText(info)
            except Exception as e:
                self.file_content.setText(f"Error reading file: {e}")
    
    def _update_dry_run_config(self):
        """Update dry run configuration display."""
        self._update_status_and_config()
        
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        self.output_dir_label.setText(workdir)
        
        config = self.session_state.get('workflow_config', {})
        atoms = self.session_state.get('current_structure')
        
        if config.get('label'):
            self.label_edit.setText(config['label'])
        elif atoms is not None and ASE_AVAILABLE and isinstance(atoms, Atoms):
            calc_type = config.get('calc_type', 'scf')
            formula = atoms.get_chemical_formula()
            self.label_edit.setText(f"{calc_type}/{formula}")
        
        # Update full path
        label = self.label_edit.text() or "scf/structure"
        self.full_path_label.setText(os.path.join(workdir, label))
    
    def _update_run_config(self):
        """Update run configuration display."""
        self._update_status_and_config(run_tab=True)
        
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        self.run_output_dir_label.setText(workdir)
        
        config = self.session_state.get('workflow_config', {})
        atoms = self.session_state.get('current_structure')
        
        if config.get('label'):
            self.run_label_edit.setText(config['label'])
        elif atoms is not None and ASE_AVAILABLE and isinstance(atoms, Atoms):
            calc_type = config.get('calc_type', 'scf')
            formula = atoms.get_chemical_formula()
            self.run_label_edit.setText(f"{calc_type}/{formula}")
        
        # Update full path
        label = self.run_label_edit.text() or "scf/structure"
        self.run_full_path_label.setText(os.path.join(workdir, label))
    
    def _update_status_and_config(self, run_tab=False):
        """Update status and configuration summary."""
        atoms = self.session_state.get('current_structure')
        config = self.session_state.get('workflow_config', {})
        
        status_label = self.run_status if run_tab else self.dry_run_status
        summary_widget = self.run_config_summary if run_tab else self.config_summary
        
        # Check structure
        if atoms is None:
            status_label.setText("⚠️ No structure loaded. Please load a structure first.")
            status_label.setStyleSheet("color: orange;")
            return
        
        if not ASE_AVAILABLE or not isinstance(atoms, Atoms):
            status_label.setText("⚠️ Invalid structure type.")
            status_label.setStyleSheet("color: orange;")
            return
        
        # Check configuration
        if not config.get('pseudopotentials'):
            status_label.setText("⚠️ No calculation configured. Please configure in Calculation Setup.")
            status_label.setStyleSheet("color: orange;")
            return
        
        status_label.setText(f"✅ Structure loaded: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
        status_label.setStyleSheet("color: green;")
        
        # Update summary
        import json
        summary_widget.setText(json.dumps(config, indent=2, default=str))
    
    def _generate_files(self):
        """Generate calculation files (dry run)."""
        atoms = self.session_state.get('current_structure')
        config = self.session_state.get('workflow_config', {})
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded.")
            return
        
        if not config.get('pseudopotentials'):
            QMessageBox.warning(self, "Warning", "No calculation configured.")
            return
        
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        label = self.label_edit.text()
        
        if not label:
            QMessageBox.warning(self, "Warning", "Please enter a calculation label.")
            return
        
        full_path = os.path.join(workdir, label)
        
        # Validate path using centralized validation utility
        is_valid, full_path, error_msg = validate_path_under_base(full_path, workdir)
        if not is_valid:
            QMessageBox.critical(self, "Error", f"Invalid calculation label: {error_msg}")
            return
        
        try:
            # Create output directory using safe utility
            safe_makedirs(full_path)
            
            # Save structure file
            if ASE_AVAILABLE:
                structure_filename = f"{atoms.get_chemical_formula()}.cif"
                structure_path = os.path.join(full_path, structure_filename)
                ase_io.write(structure_path, atoms)
            
            # Generate a simple input file preview
            input_content = self._generate_input_preview(atoms, config, label)
            
            # Save input file
            input_path = os.path.join(full_path, "espresso.pwi")
            with open(input_path, 'w') as f:
                f.write(input_content)
            
            # Update displays
            generated = ["espresso.pwi"]
            if ASE_AVAILABLE:
                generated.append(structure_filename)
            
            self.generated_files_list.setText("\n".join([f"✅ {f}" for f in generated]))
            self.input_preview.setText(input_content)
            
            self.dry_run_results.setText(f"""
✅ <b>Files generated successfully!</b>

Files created in: <code>{full_path}</code>

<b>Next Steps:</b>
<ul>
<li>Review the files using the File Browser tab</li>
<li>Edit the files if needed</li>
<li>Run the calculation using the Run Calculation tab</li>
</ul>
""")
            self.dry_run_results.setTextFormat(Qt.RichText)
            self.dry_run_results.setStyleSheet("color: green;")
            
            # Refresh browser
            self._refresh_browser()
            
        except Exception as e:
            self.dry_run_results.setText(f"❌ Error generating files: {e}")
            self.dry_run_results.setStyleSheet("color: red;")
    
    def _generate_input_preview(self, atoms, config, label):
        """Generate a preview of the input file."""
        lines = []
        lines.append("&CONTROL")
        lines.append(f"  calculation = '{config.get('calc_type', 'scf')}'")
        lines.append(f"  prefix = '{label.split('/')[-1] if '/' in label else label}'")
        lines.append("  pseudo_dir = './pseudo'")
        lines.append("  outdir = './tmp'")
        lines.append("/")
        lines.append("")
        lines.append("&SYSTEM")
        lines.append(f"  ecutwfc = {config.get('ecutwfc', 50.0)}")
        lines.append(f"  ecutrho = {config.get('ecutrho', 400.0)}")
        lines.append(f"  occupations = '{config.get('occupations', 'smearing')}'")
        if config.get('occupations') == 'smearing':
            lines.append(f"  smearing = '{config.get('smearing', 'gaussian')}'")
            lines.append(f"  degauss = {config.get('degauss', 0.02)}")
        lines.append("/")
        lines.append("")
        lines.append("&ELECTRONS")
        lines.append(f"  conv_thr = {config.get('conv_thr', 1.0e-8)}")
        lines.append("/")
        lines.append("")
        lines.append("ATOMIC_SPECIES")
        for elem, pseudo in config.get('pseudopotentials', {}).items():
            # Use a generic mass
            lines.append(f"  {elem}  1.0  {pseudo}")
        lines.append("")
        lines.append("! Note: This is a preview. Full input requires atomic positions and cell.")
        lines.append(f"! Structure: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
        
        return "\n".join(lines)
    
    def _run_calculation(self):
        """Run the calculation."""
        atoms = self.session_state.get('current_structure')
        config = self.session_state.get('workflow_config', {})
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded.")
            return
        
        if not config.get('pseudopotentials'):
            QMessageBox.warning(self, "Warning", "No calculation configured.")
            return
        
        # Note: In a real implementation, this would use xespresso to run the calculation
        # For now, we show a message that the calculation framework is not available
        
        QMessageBox.information(
            self,
            "Run Calculation",
            "To run calculations, please use the command line:\n\n"
            "python -m xespresso\n\n"
            "Or use the Streamlit GUI which has full calculation support.\n\n"
            "This PyQt GUI currently supports configuration and file generation."
        )
        
        self.run_results.setText("""
ℹ️ <b>Calculation Execution</b>

Full calculation execution requires xespresso to be properly installed
with all dependencies (Quantum ESPRESSO, MPI, etc.).

<b>To run a calculation:</b>
<ol>
<li>Use the "Generate Files" tab to create input files</li>
<li>Review and edit files as needed</li>
<li>Run manually: <code>pw.x -in espresso.pwi > espresso.pwo</code></li>
<li>Or use the Streamlit GUI for integrated execution</li>
</ol>
""")
        self.run_results.setTextFormat(Qt.RichText)
    
    def refresh(self):
        """Refresh the page."""
        self._refresh_browser()
        self._update_dry_run_config()
        self._update_run_config()

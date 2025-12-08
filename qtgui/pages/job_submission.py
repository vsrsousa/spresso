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
    QHeaderView, QSplitter, QApplication
)
from PySide6.QtCore import Qt

from qtgui.utils import validate_path_under_base, safe_makedirs, ASE_ESPRESSO_COMMAND_TEMPLATE

try:
    from ase import Atoms
    from ase import io as ase_io
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

try:
    from xespresso.xio import write_espresso_in
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False


# Default orbitals for common elements in DFT+U calculations
DEFAULT_HUBBARD_ORBITALS = {
    # 3d transition metals
    'Ti': '3d', 'V': '3d', 'Cr': '3d', 'Mn': '3d', 'Fe': '3d',
    'Co': '3d', 'Ni': '3d', 'Cu': '3d', 'Zn': '3d',
    # 4d transition metals
    'Zr': '4d', 'Nb': '4d', 'Mo': '4d', 'Tc': '4d', 'Ru': '4d',
    'Rh': '4d', 'Pd': '4d', 'Ag': '4d', 'Cd': '4d',
    # 4f rare earths
    'Ce': '4f', 'Pr': '4f', 'Nd': '4f', 'Pm': '4f', 'Sm': '4f',
    'Eu': '4f', 'Gd': '4f', 'Tb': '4f', 'Dy': '4f', 'Ho': '4f',
    'Er': '4f', 'Tm': '4f', 'Yb': '4f', 'Lu': '4f',
    # p-block elements (often used with oxygen)
    'O': '2p', 'S': '3p', 'Se': '4p',
    'N': '2p', 'P': '3p', 'As': '4p',
}


def _get_default_orbital(element):
    """Get default orbital for element in Hubbard calculations.
    
    Args:
        element: Chemical symbol (e.g., 'Fe', 'Mn')
    
    Returns:
        str: Default orbital (e.g., '3d', '4f', '2p')
    """
    return DEFAULT_HUBBARD_ORBITALS.get(element, '3d')


def _get_input_extension(code):
    """Get input file extension for the given QE code.
    
    Following ASE_ESPRESSO_COMMAND pattern: PREFIX.PACKAGEi
    
    Args:
        code: QE code name (e.g., 'pw', 'dos', 'projwfc', 'bands', 'pp', 'ph')
    
    Returns:
        str: File extension including dot (e.g., '.pwi', '.dosi', '.projwfci')
    """
    # Following xespresso's naming convention: {prefix}.{package}i
    # This matches the .PACKAGEi pattern in ASE_ESPRESSO_COMMAND
    if not code:
        return '.pwi'  # Default to pw.x
    
    # For most codes, it's just .{code}i
    return f'.{code}i'


def _get_output_extension(code):
    """Get output file extension for the given QE code.
    
    Following ASE_ESPRESSO_COMMAND pattern: PREFIX.PACKAGEo
    
    Args:
        code: QE code name (e.g., 'pw', 'dos', 'projwfc', 'bands', 'pp', 'ph')
    
    Returns:
        str: File extension including dot (e.g., '.pwo', '.doso', '.projwfco')
    """
    # Following xespresso's naming convention: {prefix}.{package}o
    # This matches the .PACKAGEo pattern in ASE_ESPRESSO_COMMAND
    if not code:
        return '.pwo'  # Default to pw.x
    
    # For most codes, it's just .{code}o
    return f'.{code}o'


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
        self.label_edit.setPlaceholderText("e.g., Al/scf")
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
<p>Run a calculation using xespresso by calling <code>atoms.get_potential_energy()</code>.</p>
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
        self.run_label_edit.setPlaceholderText("e.g., Al/scf")
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
        
        This method scans the working directory and displays its contents.
        Strategy:
        1. Show all first-level directories and files
        2. For directories, recursively scan for calculation files
        3. Display calculation files as children of their parent directories
        """
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        self.workdir_label.setText(workdir)
        
        self.file_tree.clear()
        # Clear the file content viewer when refreshing
        self.file_content.clear()
        self.file_info_label.setText("")
        
        if not os.path.exists(workdir):
            return
        
        # QE input/output file extensions
        # Pattern: PREFIX.PACKAGEi for input, PREFIX.PACKAGEo for output
        # Common extensions: .pwi/.pwo, .phi/.pho, .dosi/.doso, .bandi/.bando, etc.
        input_extensions = (".in", ".pwi", ".phi", ".ppi", ".bandi", ".dosi", ".projwfci", 
                           ".nebi", ".hpi", ".tddfpti", ".cpi")
        output_extensions = (".out", ".pwo", ".pho", ".ppo", ".bando", ".doso", ".projwfco",
                            ".nebo", ".hpo", ".tddfpto", ".cpo", ".sh", ".slurm")
        max_items_per_dir = 100  # Limit items per directory to prevent UI slowdown
        
        try:
            # Get all items in the working directory
            try:
                items = sorted(os.listdir(workdir))
            except (OSError, IOError) as e:
                self.file_info_label.setText(f"Error reading directory: {e}")
                return
            
            # Process each top-level item
            for item_name in items:
                item_path = os.path.join(workdir, item_name)
                
                # Skip hidden items
                if item_name.startswith('.'):
                    continue
                
                # Skip common large directories
                if item_name in ('node_modules', '__pycache__', '.git', 'venv', 'env'):
                    continue
                
                if os.path.isdir(item_path):
                    # Add directory to tree
                    dir_item = QTreeWidgetItem([item_name + "/"])
                    dir_item.setData(0, Qt.UserRole, item_path)
                    
                    # Scan directory recursively for calculation files
                    calc_folders = []
                    try:
                        for root, dirs, files in os.walk(item_path, topdown=True):
                            # Limit depth to 3 levels below each top-level directory
                            depth = root[len(item_path):].count(os.sep)
                            if depth >= 3:
                                dirs[:] = []
                                continue
                            
                            # Prune hidden and large directories
                            dirs[:] = [d for d in dirs if not d.startswith('.') and 
                                      d not in ('node_modules', '__pycache__', '.git', 'venv', 'env')]
                            
                            # Check if this directory has calculation files
                            calc_files = [f for f in files if (
                                f.endswith(input_extensions) or 
                                f == "job_file" or 
                                f.endswith(output_extensions)
                            )]
                            
                            if calc_files:
                                rel_path = os.path.relpath(root, item_path)
                                calc_folders.append((rel_path, root, calc_files))
                                
                                # Limit items to prevent slowdown
                                if len(calc_folders) >= max_items_per_dir:
                                    break
                    except (OSError, IOError):
                        pass  # Skip directories we can't read
                    
                    # Add calculation folders as children
                    if calc_folders:
                        for rel_path, abs_path, calc_files in calc_folders:
                            if rel_path == ".":
                                # Files are directly in this directory
                                subfolder_item = dir_item
                            else:
                                # Files are in a subdirectory
                                subfolder_item = QTreeWidgetItem([rel_path])
                                subfolder_item.setData(0, Qt.UserRole, abs_path)
                                dir_item.addChild(subfolder_item)
                            
                            # Add files as children
                            for f in calc_files[:max_items_per_dir]:
                                file_item = QTreeWidgetItem([f])
                                file_item.setData(0, Qt.UserRole, os.path.join(abs_path, f))
                                subfolder_item.addChild(file_item)
                    
                    self.file_tree.addTopLevelItem(dir_item)
                    
                elif os.path.isfile(item_path):
                    # Add files at the root level if they're relevant
                    if (item_name.endswith(input_extensions) or 
                        item_name == "job_file" or 
                        item_name.endswith(output_extensions)):
                        file_item = QTreeWidgetItem([item_name])
                        file_item.setData(0, Qt.UserRole, item_path)
                        self.file_tree.addTopLevelItem(file_item)
                    
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
            self.label_edit.setText(f"{formula}/{calc_type}")
        
        # Update full path
        label = self.label_edit.text() or "structure/scf"
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
            self.run_label_edit.setText(f"{formula}/{calc_type}")
        
        # Update full path
        label = self.run_label_edit.text() or "structure/scf"
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
        
        # Check configuration - be more lenient to avoid false negatives after session save
        # Only show "not configured" if config is completely empty or missing key fields
        if not config or (not config.get('pseudopotentials') and not config.get('calc_type')):
            status_label.setText("⚠️ No calculation configured. Please configure in Calculation Setup.")
            status_label.setStyleSheet("color: orange;")
            return
        
        # Show warning if pseudopotentials are missing but other config exists
        if not config.get('pseudopotentials'):
            status_label.setText("⚠️ Pseudopotentials not configured. Please configure in Calculation Setup.")
            status_label.setStyleSheet("color: orange;")
            return
        
        status_label.setText(f"✅ Structure loaded: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
        status_label.setStyleSheet("color: green;")
        
        # Update summary
        import json
        summary_widget.setText(json.dumps(config, indent=2, default=str))
    
    def _get_prefix_from_label(self, label):
        """Extract prefix from label.
        
        Args:
            label: Calculation label (may include path separators)
            
        Returns:
            str: Prefix for calculation (basename of label)
        """
        return label.split('/')[-1] if '/' in label else label
    
    def _setup_xespresso_environment(self):
        """Set up environment variables required by xespresso.
        
        This sets ASE_ESPRESSO_COMMAND which tells xespresso how to execute
        Quantum ESPRESSO commands following the pattern:
        LAUNCHER PACKAGE.x PARALLEL -in PREFIX.PACKAGEi > PREFIX.PACKAGEo
        """
        os.environ['ASE_ESPRESSO_COMMAND'] = ASE_ESPRESSO_COMMAND_TEMPLATE
    
    def _generate_files(self):
        """Generate calculation files (dry run).
        
        Uses the pre-created Espresso calculator and prepared_atoms from Calculation Setup page.
        Simply updates the label and calls write_input().
        """
        atoms = self.session_state.get('current_structure')
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded.")
            return
        
        # Check if calculator was prepared in Calculation Setup
        calc = self.session_state.get('espresso_calculator')
        prepared_atoms = self.session_state.get('prepared_atoms', atoms)
        
        if calc is None:
            QMessageBox.warning(
                self, 
                "Calculator Not Prepared",
                "Please go to Calculation Setup page and click 'Prepare Calculation' first.\n\n"
                "This will create the Espresso calculator object that's needed for file generation."
            )
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
            # Set up environment for xespresso
            self._setup_xespresso_environment()
            
            # Create output directory using safe utility
            safe_makedirs(full_path)
            
            # Save structure file using the PREPARED atoms (may have magnetic/Hubbard config)
            structure_filename = None
            if ASE_AVAILABLE:
                structure_filename = f"{prepared_atoms.get_chemical_formula()}.cif"
                structure_path = os.path.join(full_path, structure_filename)
                ase_io.write(structure_path, prepared_atoms)
            
            # Get configuration to determine the QE code being used
            config = self.session_state.get('workflow_config', {})
            selected_code = config.get('selected_code', 'pw')
            
            # Prepare input parameters for xespresso
            prefix = self._get_prefix_from_label(label)
            
            # Get correct input file extension based on the QE code
            input_extension = _get_input_extension(selected_code)
            input_filename = f"{prefix}{input_extension}"
            
            # Set the directory and prefix for the calculation
            # This tells xespresso where to write the files
            calc.directory = full_path
            calc.prefix = prefix
            
            # Call write_input to generate input file AND job_file via scheduler
            # xespresso's write_input method:
            # 1. Writes the input file (e.g., .pwi, .dosi) to calc.directory
            # 2. Calls set_queue() to set up the scheduler
            # 3. Scheduler writes job_file to calc.directory
            calc.write_input(prepared_atoms)
            
            input_path = os.path.join(full_path, input_filename)
            
            # Read the generated input file for preview
            with open(input_path, 'r') as f:
                input_content = f.read()
            
            # Check if job_file was created
            job_file_path = os.path.join(full_path, "job_file")
            job_file_exists = os.path.exists(job_file_path)
            
            # Update displays
            generated = [input_filename]
            if job_file_exists:
                generated.append("job_file")
            if structure_filename:
                generated.append(structure_filename)
            
            self.generated_files_list.setText("\n".join([f"✅ {f}" for f in generated]))
            self.input_preview.setText(input_content)
            
            self.dry_run_results.setText(f"""
✅ <b>Files generated successfully!</b>

Files created in: <code>{full_path}</code>

<b>Generated Files:</b>
<ul>
<li><b>{input_filename}</b> - QE input file with atomic positions and cell</li>
{f'<li><b>job_file</b> - Execution script</li>' if job_file_exists else ''}
{f'<li><b>{structure_filename}</b> - Structure file</li>' if structure_filename else ''}
</ul>

<b>Next Steps:</b>
<ul>
<li>Review the files using the File Browser tab</li>
<li>Edit the files if needed</li>
<li>Run the calculation using: <code>bash job_file</code></li>
</ul>
""")
            self.dry_run_results.setTextFormat(Qt.RichText)
            self.dry_run_results.setStyleSheet("color: green;")
            
            # Process pending events to ensure filesystem changes are reflected
            QApplication.processEvents()
            
            # Refresh browser
            self._refresh_browser()
            
        except Exception as e:
            import traceback
            traceback.print_exc()
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
    
    def _build_input_data(self, config, prefix):
        """Build input_data dictionary for xespresso from GUI config.
        
        Args:
            config: Configuration dictionary from the GUI
            prefix: Calculation prefix (used for outdir)
            
        Returns:
            dict: Input data dictionary formatted for xespresso's write_espresso_in
        """
        def ensure_input_ntyp(input_data):
            """Ensure input_ntyp key exists in input_data."""
            if 'input_ntyp' not in input_data:
                input_data['input_ntyp'] = {}
        
        input_data = {
            'CONTROL': {
                'calculation': config.get('calc_type', 'scf'),
                'prefix': prefix,
                'pseudo_dir': './pseudo',
                'outdir': './tmp',
                'verbosity': 'high',
            },
            'SYSTEM': {
                'ecutwfc': config.get('ecutwfc', 50.0),
                'ecutrho': config.get('ecutrho', 400.0),
                'occupations': config.get('occupations', 'smearing'),
            },
            'ELECTRONS': {
                'conv_thr': config.get('conv_thr', 1.0e-8),
            },
        }
        
        # Add smearing parameters if using smearing occupations
        if config.get('occupations') == 'smearing':
            input_data['SYSTEM']['smearing'] = config.get('smearing', 'gaussian')
            input_data['SYSTEM']['degauss'] = config.get('degauss', 0.02)
        
        # Add magnetic configuration if enabled
        if config.get('enable_magnetism') and config.get('magnetic_config'):
            input_data['SYSTEM']['nspin'] = 2
            # Use lowercase 'input_ntyp' as required by xespresso
            ensure_input_ntyp(input_data)
            input_data['input_ntyp']['starting_magnetization'] = {}
            for element, mag_values in config.get('magnetic_config', {}).items():
                if isinstance(mag_values, list):
                    input_data['input_ntyp']['starting_magnetization'][element] = mag_values[0]
                else:
                    input_data['input_ntyp']['starting_magnetization'][element] = mag_values
        
        # Add Hubbard U configuration if enabled
        if config.get('enable_hubbard') and config.get('hubbard_u'):
            # Determine Hubbard format
            hubbard_format = config.get('hubbard_format', 'old')
            qe_version = config.get('qe_version', '')
            
            # Auto-detect format from QE version if not explicitly set
            use_new_format = False
            if hubbard_format == 'new':
                use_new_format = True
            elif hubbard_format == 'old':
                use_new_format = False
            elif qe_version:
                # Parse version and determine format
                try:
                    major = int(qe_version.split('.')[0])
                    use_new_format = (major >= 7)
                except (ValueError, IndexError):
                    pass
            
            if use_new_format:
                # NEW FORMAT (QE >= 7.0): Use 'hubbard' dictionary with HUBBARD card
                # NOTE: lda_plus_u should NOT be set when using HUBBARD card (new format)
                hubbard_dict = {
                    'projector': config.get('hubbard_projector', 'atomic'),
                    'u': {},
                    'v': []
                }
                
                # Build U parameters with orbital specifications
                for element, u_value in config.get('hubbard_u', {}).items():
                    if u_value > 0:
                        # Get orbital from config, default to common orbitals
                        orbital = config.get('hubbard_orbitals', {}).get(element, _get_default_orbital(element))
                        hubbard_dict['u'][f"{element}-{orbital}"] = u_value
                
                # Add V parameters if present
                if config.get('hubbard_v'):
                    hubbard_dict['v'] = config['hubbard_v']
                
                input_data['hubbard'] = hubbard_dict
                if qe_version:
                    input_data['qe_version'] = qe_version
            else:
                # OLD FORMAT (QE < 7.0): Use 'input_ntyp' with Hubbard_U
                # Only set lda_plus_u for old format
                input_data['SYSTEM']['lda_plus_u'] = True
                ensure_input_ntyp(input_data)
                input_data['input_ntyp']['Hubbard_U'] = {}
                for element, u_value in config.get('hubbard_u', {}).items():
                    if u_value > 0:
                        input_data['input_ntyp']['Hubbard_U'][element] = u_value
        
        # Add relaxation parameters if doing relaxation
        if config.get('calc_type') in ('relax', 'vc-relax'):
            input_data['IONS'] = {}
            if config.get('forc_conv_thr'):
                input_data['CONTROL']['forc_conv_thr'] = config.get('forc_conv_thr')
        
        return input_data
    
    def _create_job_file(self, job_file_path, prefix, config):
        """Create a job_file script for executing the calculation.
        
        Args:
            job_file_path: Full path to the job_file
            prefix: Calculation prefix
            config: Configuration dictionary from the GUI
        """
        lines = ["#!/bin/bash", ""]
        
        # Add comment header
        lines.append("# Job script generated by xespresso GUI")
        lines.append(f"# Calculation: {config.get('calc_type', 'scf')}")
        lines.append("")
        
        # Module loading from codes configuration
        modules = config.get('modules')
        if modules:
            lines.append("# Load required modules")
            if isinstance(modules, list):
                for module in modules:
                    lines.append(f"module load {module}")
            else:
                lines.append(f"module load {modules}")
            lines.append("")
        
        # Additional environment setup (user can customize)
        lines.append("# Additional environment setup (uncomment and modify as needed)")
        lines.append("# source /path/to/espresso/env.sh")
        if not modules:
            lines.append("# module load <your-quantum-espresso-module>")
        lines.append("")
        
        # Get launcher from machine configuration
        launcher = ""
        nprocs = config.get('nprocs', 1)
        machine = self.session_state.get('calc_machine')
        
        if machine and hasattr(machine, 'launcher') and machine.launcher:
            # Use launcher from machine config
            launcher = machine.launcher
            # Replace {nprocs} placeholder if present
            if '{nprocs}' in launcher:
                launcher = launcher.replace('{nprocs}', str(nprocs))
            launcher = launcher + " "
        elif nprocs > 1:
            # Fallback to default mpirun launcher
            launcher = f"mpirun -np {nprocs} "
        
        # Execution command
        lines.append(f"{launcher}pw.x -in {prefix}.pwi > {prefix}.pwo")
        lines.append("")
        
        # Write the job file
        with open(job_file_path, 'w') as f:
            f.write("\n".join(lines))
        
        # Make it executable
        os.chmod(job_file_path, 0o755)
    
    def _run_calculation(self):
        """Run the calculation.
        
        Uses the pre-created Espresso calculator and prepared_atoms from Calculation Setup page.
        Simply updates the label and calls atoms.get_potential_energy().
        """
        atoms = self.session_state.get('current_structure')
        
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded.")
            return
        
        # Check if calculator was prepared in Calculation Setup
        calc = self.session_state.get('espresso_calculator')
        prepared_atoms = self.session_state.get('prepared_atoms', atoms)
        
        if calc is None:
            QMessageBox.warning(
                self, 
                "Calculator Not Prepared",
                "Please go to Calculation Setup page and click 'Prepare Calculation' first.\n\n"
                "This will create the Espresso calculator object that's needed for execution."
            )
            return
        
        if not XESPRESSO_AVAILABLE:
            QMessageBox.warning(
                self,
                "xespresso Not Available",
                "xespresso module is not installed or not available.\n\n"
                "Install it with: pip install xespresso"
            )
            return
        
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        label = self.run_label_edit.text()
        
        if not label:
            QMessageBox.warning(self, "Warning", "Please enter a calculation label.")
            return
        
        full_path = os.path.join(workdir, label)
        
        # Validate path using centralized validation utility
        from qtgui.utils import validate_path_under_base, safe_makedirs
        is_valid, full_path, error_msg = validate_path_under_base(full_path, workdir)
        if not is_valid:
            QMessageBox.critical(self, "Error", f"Invalid calculation label: {error_msg}")
            return
        
        try:
            # Set up environment for xespresso
            self._setup_xespresso_environment()
            
            # Create output directory
            safe_makedirs(full_path)
            
            # Set the directory and prefix for the calculation
            # This tells xespresso where to write and read files
            prefix = self._get_prefix_from_label(label)
            calc.directory = full_path
            calc.prefix = prefix
            
            # Update status
            self.run_status.setText("⏳ Running calculation... (this may take a while)")
            self.run_status.setStyleSheet("color: blue;")
            self.run_results.setHtml(f"""
<b>Calculation started</b>

Working directory: <code>{full_path}</code>

The calculation is now running. This may take several minutes to hours
depending on system size and computational resources.

<b>Note:</b> The GUI may appear frozen during calculation execution.
This is normal for local calculations.
""")
            QApplication.processEvents()  # Update UI
            
            # Ensure calculator is attached to atoms (following xespresso pattern)
            # This is required for atoms.get_potential_energy() to work
            # Defensive check for robustness (e.g., restored from old session state)
            if prepared_atoms.calc is None or prepared_atoms.calc != calc:
                prepared_atoms.calc = calc
            
            # Run the calculation using ASE/xespresso pattern
            # xespresso's get_potential_energy() will:
            # 1. Check for previous results
            # 2. Generate input files if needed using write_input()
            # 3. Submit the job via the configured scheduler
            # 4. Wait for completion and parse output
            energy = prepared_atoms.get_potential_energy()
            
            # Success! Display results
            self.run_status.setText("✅ Calculation completed successfully!")
            self.run_status.setStyleSheet("color: green;")
            
            self.energy_label.setText(f"Energy: {energy:.6f} eV")
            
            # Get prefix for output files
            prefix = self._get_prefix_from_label(label)
            
            # Get additional results if available
            results_text = f"""
<b>✅ Calculation Completed Successfully!</b>

<b>Energy:</b> {energy:.6f} eV

<b>Output Directory:</b> <code>{full_path}</code>

<b>Output Files:</b>
<ul>
<li>Input file: <code>{prefix}.pwi</code></li>
<li>Output file: <code>{prefix}.pwo</code></li>
<li>Job script: <code>job_file</code> (if using scheduler)</li>
</ul>

<b>Next Steps:</b>
<ul>
<li>View output files in the File Browser tab</li>
<li>Analyze results in Results & Post-Processing</li>
<li>Continue workflow (e.g., bands, DOS calculation)</li>
</ul>
"""
            
            # Add forces if available
            # Some calculation types may not have forces available
            if hasattr(prepared_atoms, 'get_forces'):
                try:
                    forces = prepared_atoms.get_forces()
                    max_force = float((forces**2).sum(axis=1).max()**0.5)
                    results_text += f"\n<b>Max Force:</b> {max_force:.6f} eV/Å\n"
                except (RuntimeError, KeyError) as e:
                    # Forces not available for this calculation type (e.g., SCF)
                    # This is expected and not an error
                    import logging
                    logging.debug(f"Forces not available: {e}")
                except Exception as e:
                    # Unexpected error getting forces
                    import logging
                    logging.warning(f"Unexpected error getting forces: {e}")
            
            self.run_results.setHtml(results_text)
            
            # Process pending events to ensure filesystem changes are reflected
            QApplication.processEvents()
            
            # Refresh browser to show new files
            self._refresh_browser()
            
            QMessageBox.information(
                self,
                "Calculation Complete",
                f"Calculation completed successfully!\n\nEnergy: {energy:.6f} eV"
            )
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            
            self.run_status.setText("❌ Calculation failed!")
            self.run_status.setStyleSheet("color: red;")
            
            self.energy_label.setText("")
            
            self.run_results.setHtml(f"""
<b>❌ Calculation Failed</b>

<b>Error:</b> {str(e)}

<b>Possible Causes:</b>
<ul>
<li>Quantum ESPRESSO not installed or not in PATH</li>
<li>Missing pseudopotential files</li>
<li>Invalid calculation parameters</li>
<li>Insufficient computational resources</li>
<li>MPI not configured correctly (for parallel runs)</li>
</ul>

<b>Troubleshooting:</b>
<ol>
<li>Check that QE is installed: <code>which pw.x</code></li>
<li>Verify pseudopotential paths are correct</li>
<li>Review input file in: <code>{full_path}</code></li>
<li>Try running manually: <code>bash job_file</code></li>
<li>Check QE output for detailed error messages</li>
</ol>

<b>Detailed Error:</b>
<pre>{error_traceback}</pre>
""")
            
            QMessageBox.critical(
                self,
                "Calculation Error",
                f"Calculation failed with error:\n\n{str(e)}\n\n"
                f"See details in the Results section below."
            )
    
    def refresh(self):
        """Refresh the page."""
        self._refresh_browser()
        self._update_dry_run_config()
        self._update_run_config()

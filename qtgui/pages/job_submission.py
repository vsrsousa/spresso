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

try:
    from xespresso.xio import write_espresso_in
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False


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
        
        This method scans the working directory for calculation folders.
        To prevent blocking the main thread during large directory scans,
        we limit the depth and number of folders found.
        """
        workdir = self.session_state.get('working_directory', os.path.expanduser("~"))
        self.workdir_label.setText(workdir)
        
        self.file_tree.clear()
        
        if not os.path.exists(workdir):
            return
        
        # Find calculation folders with limits to prevent blocking
        input_extensions = (".in", ".pwi", ".phi", ".ppi", ".bandi")
        output_extensions = (".sh", ".slurm", ".out", ".pwo")
        max_folders_found = 100  # Limit number of calculation folders to display
        folders_found = 0
        
        try:
            for root, dirs, files in os.walk(workdir, topdown=True):
                depth = root[len(workdir):].count(os.sep)
                
                # Limit search depth
                if depth >= 4:
                    dirs[:] = []  # Don't descend further
                    continue
                
                # Prune hidden directories and common large directories to speed up walk
                dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ('node_modules', '__pycache__', '.git', 'venv', 'env')]
                
                has_input = any(
                    f.endswith(input_extensions) or 
                    f == "job_file" or 
                    f.endswith((".sh", ".slurm"))
                    for f in files
                )
                
                if has_input:
                    folders_found += 1
                    if folders_found > max_folders_found:
                        self.file_info_label.setText(
                            f"⚠️ Found {max_folders_found}+ calculation folders. "
                            "Navigate to a specific folder for better results."
                        )
                        break
                    
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
    
    def _generate_files(self):
        """Generate calculation files (dry run).
        
        Uses xespresso's Espresso calculator to generate a proper QE input file
        with atomic positions and cell parameters from the loaded structure,
        and creates the job_file script via xespresso's scheduler.
        """
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
            structure_filename = None
            if ASE_AVAILABLE:
                structure_filename = f"{atoms.get_chemical_formula()}.cif"
                structure_path = os.path.join(full_path, structure_filename)
                ase_io.write(structure_path, atoms)
            
            # Prepare input parameters for xespresso
            prefix = label.split('/')[-1] if '/' in label else label
            input_filename = f"{prefix}.pwi"
            
            # Build input_data dictionary for xespresso
            input_data = self._build_input_data(config, prefix)
            
            # Get k-points configuration
            kspacing = config.get('kspacing')
            kpts = config.get('kpts')
            
            # Use xespresso's Espresso calculator for proper dry run
            if XESPRESSO_AVAILABLE:
                try:
                    from xespresso import Espresso
                    
                    # Build queue configuration for job_file generation
                    # First check if machine is configured and use it
                    machine = self.session_state.get('calc_machine')
                    if machine:
                        # Machine object is configured, convert to queue dict using to_queue()
                        try:
                            queue = machine.to_queue()
                        except (AttributeError, TypeError) as e:
                            # Fallback if machine doesn't have to_queue() or it fails
                            print(f"Warning: Could not convert machine to queue: {e}")
                            queue = {
                                'execution': 'local',
                                'scheduler': 'direct',
                            }
                    else:
                        # Build queue configuration manually
                        queue = {
                            'execution': 'local',
                            'scheduler': 'direct',
                        }
                        
                        # Add resources if configured
                        if config.get('adjust_resources'):
                            resources = config.get('resources', {})
                            queue.update(resources)
                    
                    # Create Espresso calculator
                    calc_kwargs = {
                        'label': full_path,
                        'pseudopotentials': config.get('pseudopotentials', {}),
                        'input_data': input_data,
                        'queue': queue,
                    }
                    
                    if kspacing:
                        calc_kwargs['kspacing'] = kspacing
                    elif kpts:
                        calc_kwargs['kpts'] = kpts
                    else:
                        calc_kwargs['kpts'] = (4, 4, 4)  # Default k-points
                    
                    # Set parallel command if multiple processors
                    nprocs = config.get('nprocs', 1)
                    if nprocs > 1:
                        calc_kwargs['parallel'] = f'-np {nprocs}'
                    
                    calc = Espresso(**calc_kwargs)
                    
                    # Make a copy of atoms to avoid modifying the original
                    atoms_copy = atoms.copy()
                    atoms_copy.calc = calc
                    
                    # Call write_input to generate input file AND job_file via scheduler
                    calc.write_input(atoms_copy)
                    
                    input_path = f"{calc.label}.pwi"
                    
                except Exception as e:
                    # If xespresso calculator fails, fall back to direct write
                    import traceback
                    traceback.print_exc()
                    input_path = os.path.join(full_path, input_filename)
                    # Pass either kspacing or kpts, not both
                    write_kwargs = {
                        'input_data': input_data,
                        'pseudopotentials': config.get('pseudopotentials', {}),
                    }
                    if kspacing:
                        write_kwargs['kspacing'] = kspacing
                    elif kpts:
                        write_kwargs['kpts'] = kpts
                    write_espresso_in(input_path, atoms, **write_kwargs)
                    # Create job_file manually as fallback
                    job_file_path = os.path.join(full_path, "job_file")
                    self._create_job_file(job_file_path, prefix, config)
            else:
                # Fallback to simple preview if xespresso not available
                input_path = os.path.join(full_path, input_filename)
                input_content = self._generate_input_preview(atoms, config, label)
                with open(input_path, 'w') as f:
                    f.write(input_content)
                # Create job_file manually
                job_file_path = os.path.join(full_path, "job_file")
                self._create_job_file(job_file_path, prefix, config)
            
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
            if 'input_ntyp' not in input_data:
                input_data['input_ntyp'] = {}
            input_data['input_ntyp']['starting_magnetization'] = {}
            for element, mag_values in config.get('magnetic_config', {}).items():
                if isinstance(mag_values, list):
                    input_data['input_ntyp']['starting_magnetization'][element] = mag_values[0]
                else:
                    input_data['input_ntyp']['starting_magnetization'][element] = mag_values
        
        # Add Hubbard U configuration if enabled
        if config.get('enable_hubbard') and config.get('hubbard_u'):
            input_data['SYSTEM']['lda_plus_u'] = True
            # Use lowercase 'input_ntyp' as required by xespresso
            if 'input_ntyp' not in input_data:
                input_data['input_ntyp'] = {}
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
        
        # Environment setup (user can customize)
        lines.append("# Environment setup (uncomment and modify as needed)")
        lines.append("# source /path/to/espresso/env.sh")
        lines.append("# module load quantum-espresso")
        lines.append("")
        
        # Execution command
        nprocs = config.get('nprocs', 1)
        if nprocs > 1:
            lines.append(f"mpirun -np {nprocs} pw.x -in {prefix}.pwi > {prefix}.pwo")
        else:
            lines.append(f"pw.x -in {prefix}.pwi > {prefix}.pwo")
        lines.append("")
        
        # Write the job file
        with open(job_file_path, 'w') as f:
            f.write("\n".join(lines))
        
        # Make it executable
        os.chmod(job_file_path, 0o755)
    
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

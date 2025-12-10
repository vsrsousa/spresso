"""
Results and Post-processing Page for xespresso PySide6 GUI.

This page displays calculation results and provides post-processing tools.
"""

import os
import re

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QComboBox, QPushButton, QGroupBox, QFormLayout,
    QMessageBox, QScrollArea, QFrame, QTextEdit,
    QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView,
    QFileDialog
)
from PySide6.QtCore import Qt

# Unit conversion constants
# Energy: 1 Rydberg = 13.605693122994 eV (CODATA 2018)
RY_TO_EV = 13.605693122994
EV_TO_RY = 1.0 / RY_TO_EV

# Force: 1 Ry/Bohr = 25.71104309541616 eV/Angstrom
RYBOHR_TO_EVANG = 25.71104309541616
EVANG_TO_RYBOHR = 1.0 / RYBOHR_TO_EVANG

# Stress/Pressure: 1 eV/Angstrom^3 = 160.21766208 GPa = 1602.1766208 kbar
EVANG3_TO_KBAR = 1602.1766208


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
        self.convergence_table.setColumnCount(3)
        self.convergence_table.setHorizontalHeaderLabels(["Iteration", "Energy (Ry)", "Delta E"])
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
            # First, try to use ASE to read the output (preferred method)
            results = self._parse_with_ase(output_path)
            
            # If ASE parsing fails or doesn't get all info, supplement with manual parsing
            if not results or results.get('energy') is None:
                with open(output_path, 'r') as f:
                    content = f.read()
                results = self._parse_output(content)
            else:
                # Supplement ASE results with manual parsing for additional fields
                with open(output_path, 'r') as f:
                    content = f.read()
                manual_results = self._parse_output(content)
                # Merge results, preferring ASE values when available
                for key in ['scf_history', 'is_magnetic', 'pressure', 'stress_tensor']:
                    if key in manual_results and manual_results[key]:
                        results[key] = manual_results[key]
            
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
            
            if results.get('fermi_energy') is not None:
                results_text.append(f"Fermi Energy: {results['fermi_energy']:.4f} eV")
            
            # Only show magnetic properties if calculation is magnetic
            if results.get('is_magnetic'):
                if results.get('total_magnetization') is not None:
                    results_text.append(f"Total Magnetization: {results['total_magnetization']:.4f} Bohr mag/cell")
                
                if results.get('absolute_magnetization') is not None:
                    results_text.append(f"Absolute Magnetization: {results['absolute_magnetization']:.4f} Bohr mag/cell")
            
            if results.get('total_force'):
                results_text.append(f"Total Force: {results['total_force']} Ry/au")
            
            if results.get('pressure') is not None:
                results_text.append(f"Pressure: {results['pressure']:.2f} kbar")
            
            # Only show magnetic moments if calculation is magnetic
            if results.get('is_magnetic') and results.get('magnetic_moments'):
                results_text.append("\nMagnetic Moments per Atom:")
                for atom_info in results['magnetic_moments']:
                    # Use atom symbol if available, otherwise just the number
                    atom_label = f"{atom_info.get('symbol', 'X')}{atom_info['atom']}" if 'symbol' in atom_info else f"Atom {atom_info['atom']}"
                    results_text.append(f"  {atom_label}: charge={atom_info['charge']:.4f}, magn={atom_info['magn']:.4f}")
            
            # Show forces per atom if available
            if results.get('forces'):
                results_text.append("\nForces per Atom (Ry/au):")
                for force_info in results['forces']:
                    results_text.append(f"  Atom {force_info['atom']}: fx={force_info['fx']:10.6f}, fy={force_info['fy']:10.6f}, fz={force_info['fz']:10.6f}")
            
            # Show stress tensor if available
            if results.get('stress_tensor'):
                results_text.append("\nStress Tensor (kbar):")
                for i, row in enumerate(results['stress_tensor']):
                    results_text.append(f"  [{row[0]:10.2f}  {row[1]:10.2f}  {row[2]:10.2f}]")
            
            self.results_text.setText("\n".join(results_text))
            
            # Update convergence table
            if results.get('scf_history'):
                history = results['scf_history']
                self.convergence_table.setRowCount(len(history))
                for i, (energy, delta) in enumerate(history):
                    self.convergence_table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
                    self.convergence_table.setItem(i, 1, QTableWidgetItem(f"{energy:.8f}"))
                    self.convergence_table.setItem(i, 2, QTableWidgetItem(f"{delta:.2e}"))
            
            self.status_text.setText(f"✅ Results loaded from: {output_path}")
            self.status_text.setStyleSheet("color: green;")
            
        except Exception as e:
            self.status_text.setText(f"❌ Error loading results: {e}")
            self.status_text.setStyleSheet("color: red;")
    
    def _parse_with_ase(self, output_path):
        """Parse QE output using ASE (primary method).
        
        ASE is always available since xespresso depends on it.
        Returns a results dictionary compatible with _parse_output.
        """
        try:
            from ase import io
            
            # Read the output file with ASE
            atoms = io.read(output_path, format='espresso-out')
            
            if atoms.calc is None:
                return None
            
            calc_results = atoms.calc.results
            
            results = {
                'energy': None,
                'converged': True,  # If ASE read it, it likely converged
                'iterations': 0,
                'scf_history': [],
                'total_force': None,
                'fermi_energy': None,
                'total_magnetization': None,
                'absolute_magnetization': None,
                'magnetic_moments': [],
                'forces': [],
                'pressure': None,
                'stress_tensor': None,
                'is_magnetic': False
            }
            
            # Extract energy (in eV, convert to Ry)
            if 'energy' in calc_results:
                results['energy'] = calc_results['energy'] * EV_TO_RY
            
            # Extract forces
            if 'forces' in calc_results and calc_results['forces'] is not None:
                forces_array = calc_results['forces']
                # Convert from eV/Angstrom to Ry/Bohr
                for i, force in enumerate(forces_array):
                    results['forces'].append({
                        'atom': i + 1,
                        'fx': force[0] * EVANG_TO_RYBOHR,
                        'fy': force[1] * EVANG_TO_RYBOHR,
                        'fz': force[2] * EVANG_TO_RYBOHR
                    })
                # Calculate total force magnitude
                total_f = sum([f[0]**2 + f[1]**2 + f[2]**2 for f in forces_array])**0.5
                results['total_force'] = total_f * EVANG_TO_RYBOHR
            
            # Extract stress
            if 'stress' in calc_results and calc_results['stress'] is not None:
                stress_array = calc_results['stress']
                # Convert from eV/Angstrom^3 to kbar
                # Stress is in Voigt notation: [xx, yy, zz, yz, xz, xy]
                results['stress_tensor'] = [
                    [stress_array[0] * EVANG3_TO_KBAR, stress_array[5] * EVANG3_TO_KBAR, stress_array[4] * EVANG3_TO_KBAR],
                    [stress_array[5] * EVANG3_TO_KBAR, stress_array[1] * EVANG3_TO_KBAR, stress_array[3] * EVANG3_TO_KBAR],
                    [stress_array[4] * EVANG3_TO_KBAR, stress_array[3] * EVANG3_TO_KBAR, stress_array[2] * EVANG3_TO_KBAR]
                ]
                # Calculate pressure (negative trace / 3)
                results['pressure'] = -(stress_array[0] + stress_array[1] + stress_array[2]) * EVANG3_TO_KBAR / 3.0
            
            # Extract magnetic moments
            if 'magmoms' in calc_results and calc_results['magmoms'] is not None:
                magmoms_array = calc_results['magmoms']
                if any(abs(m) > 1e-6 for m in magmoms_array):
                    results['is_magnetic'] = True
                    # Get atom symbols for labeling (atoms are in same order as input)
                    symbols = atoms.get_chemical_symbols()
                    # Show all atoms' magnetic moments, regardless of magnitude
                    for i, magmom in enumerate(magmoms_array):
                        results['magnetic_moments'].append({
                            'atom': i + 1,
                            'symbol': symbols[i] if i < len(symbols) else 'X',
                            'charge': 0.0,  # Not available from ASE
                            'magn': magmom
                        })
                    # Calculate total magnetization
                    results['total_magnetization'] = sum(magmoms_array)
                    results['absolute_magnetization'] = sum(abs(m) for m in magmoms_array)
            
            # Extract Fermi energy
            try:
                if hasattr(atoms.calc, 'get_fermi_level'):
                    results['fermi_energy'] = atoms.calc.get_fermi_level()
            except Exception:
                pass
            
            return results
            
        except Exception as e:
            # ASE parsing failed, will fall back to manual parsing
            return None
    
    def _parse_output(self, content):
        """Parse QE output file content.
        
        The exclamation point in '!    total energy' is the definitive indicator
        that an SCF calculation has converged in Quantum ESPRESSO.
        """
        results = {
            'energy': None,
            'converged': False,
            'iterations': 0,
            'scf_history': [],
            'total_force': None,
            'fermi_energy': None,
            'total_magnetization': None,
            'absolute_magnetization': None,
            'magnetic_moments': [],
            'forces': [],
            'pressure': None,
            'stress_tensor': None,
            'is_magnetic': False
        }
        
        lines = content.split('\n')
        prev_energy = None
        in_forces_section = False
        in_stress_section = False
        
        for i, line in enumerate(lines):
            # Total energy with '!' indicates CONVERGED calculation
            # This is the definitive convergence indicator in QE output
            if '!' in line and 'total energy' in line.lower():
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        energy_str = parts[1].replace('Ry', '').strip()
                        results['energy'] = float(energy_str)
                        # The presence of '!' in the total energy line means converged
                        results['converged'] = True
                except (ValueError, IndexError):
                    pass
            
            # SCF iteration energies (without '!')
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
                except (ValueError, IndexError):
                    pass
            
            # Additional convergence indicators (backup checks)
            # Also extract iteration count from "convergence has been achieved in X iterations"
            if 'convergence achieved' in line.lower() or \
               'convergence has been achieved' in line.lower():
                results['converged'] = True
                # Extract iteration count if present
                # Example: "convergence has been achieved in  10 iterations"
                if 'in' in line.lower() and 'iteration' in line.lower():
                    try:
                        # Find the number between "in" and "iteration(s)"
                        match = re.search(r'in\s+(\d+)\s+iterations?', line.lower())
                        if match:
                            iterations = int(match.group(1))
                            # Only update if we haven't counted SCF iterations yet
                            if results['iterations'] == 0:
                                results['iterations'] = iterations
                    except (ValueError, AttributeError):
                        pass
            
            # Total force
            if 'Total force' in line:
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        force_str = parts[1].strip()
                        results['total_force'] = float(force_str)
                except (ValueError, IndexError):
                    pass
            
            # Fermi energy
            if 'the Fermi energy is' in line.lower():
                try:
                    parts = line.split('is')
                    if len(parts) > 1:
                        fermi_str = parts[1].replace('ev', '').replace('eV', '').strip()
                        results['fermi_energy'] = float(fermi_str)
                except (ValueError, IndexError):
                    pass
            
            # Total magnetization - marks calculation as magnetic
            if 'total magnetization' in line.lower() and '=' in line:
                results['is_magnetic'] = True
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        mag_str = parts[1].split('Bohr')[0].strip()
                        results['total_magnetization'] = float(mag_str)
                except (ValueError, IndexError):
                    pass
            
            # Absolute magnetization
            if 'absolute magnetization' in line.lower() and '=' in line:
                results['is_magnetic'] = True
                try:
                    parts = line.split('=')
                    if len(parts) > 1:
                        mag_str = parts[1].split('Bohr')[0].strip()
                        results['absolute_magnetization'] = float(mag_str)
                except (ValueError, IndexError):
                    pass
            
            # Magnetic moment per site
            # Show ALL magnetic moments regardless of value
            # Note: Atom symbols will be added by ASE parsing; manual parsing doesn't extract symbols
            if 'atom:' in line.lower() and 'charge:' in line.lower() and 'magn:' in line.lower():
                results['is_magnetic'] = True
                try:
                    # Parse line like: "     atom:    1    charge:   14.5678    magn:    1.9876    constr:    0.0000"
                    parts = line.split()
                    atom_idx = None
                    charge = None
                    magn = None
                    
                    for j, part in enumerate(parts):
                        if part.lower() == 'atom:' and j + 1 < len(parts):
                            atom_idx = int(parts[j + 1])
                        elif part.lower() == 'charge:' and j + 1 < len(parts):
                            charge = float(parts[j + 1])
                        elif part.lower() == 'magn:' and j + 1 < len(parts):
                            magn = float(parts[j + 1])
                    
                    # Include ALL atoms with magnetic moments, no filtering by magnitude
                    if atom_idx is not None and charge is not None and magn is not None:
                        results['magnetic_moments'].append({
                            'atom': atom_idx,
                            'charge': charge,
                            'magn': magn
                            # 'symbol' will be added when merging with ASE results
                        })
                except (ValueError, IndexError):
                    pass
            
            # Forces section
            if 'Forces acting on atoms' in line:
                in_forces_section = True
                continue
            
            if in_forces_section:
                # Parse force line: "     atom    1 type  1   force =     0.00012345    0.00023456   -0.00034567"
                if 'atom' in line.lower() and 'force' in line.lower() and '=' in line:
                    try:
                        parts = line.split('=')
                        if len(parts) > 1:
                            force_parts = parts[1].split()
                            if len(force_parts) >= 3:
                                atom_line = parts[0].split()
                                atom_idx = None
                                for j, part in enumerate(atom_line):
                                    if part.lower() == 'atom' and j + 1 < len(atom_line):
                                        atom_idx = int(atom_line[j + 1])
                                        break
                                
                                if atom_idx is not None:
                                    results['forces'].append({
                                        'atom': atom_idx,
                                        'fx': float(force_parts[0]),
                                        'fy': float(force_parts[1]),
                                        'fz': float(force_parts[2])
                                    })
                    except (ValueError, IndexError):
                        pass
                
                # End of forces section
                if 'Total force' in line:
                    in_forces_section = False
            
            # Stress/Pressure section
            if 'entering subroutine stress' in line.lower():
                in_stress_section = True
                continue
            
            if in_stress_section:
                # Parse pressure line: "          total   stress  (Ry/bohr**3)                   (kbar)     P=   12.34"
                if 'P=' in line or 'P =' in line:
                    try:
                        # Extract pressure value
                        if 'P=' in line:
                            pressure_str = line.split('P=')[1].strip().split()[0]
                        else:
                            pressure_str = line.split('P =')[1].strip().split()[0]
                        results['pressure'] = float(pressure_str)
                    except (ValueError, IndexError):
                        pass
                
                # Parse stress tensor (3 lines of 3 values each in kbar)
                # We look for lines with numbers in the stress section
                if line.strip() and not any(keyword in line.lower() for keyword in ['total', 'stress', 'kbar', 'ry/bohr']):
                    try:
                        parts = line.split()
                        # Check if we have numeric values (at least 3 for the first half, or 3 for the kbar half)
                        if len(parts) >= 6:
                            # Extract kbar values (usually the last 3 columns)
                            kbar_values = [float(parts[-3]), float(parts[-2]), float(parts[-1])]
                            if results['stress_tensor'] is None:
                                results['stress_tensor'] = []
                            results['stress_tensor'].append(kbar_values)
                            
                            # After 3 rows, we're done
                            if len(results['stress_tensor']) >= 3:
                                in_stress_section = False
                    except (ValueError, IndexError):
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
        """Refresh the page and auto-load output directory from session."""
        # Auto-load output directory from working directory and calculation label
        working_dir = self.session_state.get('working_directory', os.path.expanduser("~"))
        config = self.session_state.get('workflow_config', {})
        label = config.get('label', '')
        
        if working_dir and label:
            # Construct the expected output directory
            output_dir = os.path.join(working_dir, label)
            if os.path.exists(output_dir):
                self.output_dir_edit.setText(output_dir)
            else:
                # If full path doesn't exist, just set working directory
                self.output_dir_edit.setText(working_dir)

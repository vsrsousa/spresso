"""
Simplified workflow for common Quantum ESPRESSO calculations.

This module provides an easy-to-use interface for running common calculations
like SCF and structure relaxation from CIF files or ASE Atoms objects.
"""

import os
import logging
import numpy as np
from typing import Dict, Optional, Union, Tuple, List
from pathlib import Path
from ase import Atoms
from ase.io import read
from ase.io.espresso import kspacing_to_grid
from xespresso import Espresso
from xespresso.tools import setup_magnetic_config
from xespresso.machines import load_machine
from xespresso.pseudopotentials import load_pseudopotentials_config
from xespresso.codes import load_codes_config
from xespresso.schedulers import RemoteJobMonitor


logger = logging.getLogger(__name__)


# Preset configurations for different calculation protocols
PRESETS = {
    'fast': {
        'ecutwfc': 30.0,
        'ecutrho': 240.0,
        'conv_thr': 1.0e-6,
        'kspacing': 0.5,  # Angstrom^-1
        'mixing_beta': 0.7,
        'electron_maxstep': 100,
        'occupations': 'smearing',  # Required for smearing
        'smearing': 'cold',  # Better than gaussian for metals
        'degauss': 0.02,  # Ry
    },
    'moderate': {
        'ecutwfc': 50.0,
        'ecutrho': 400.0,
        'conv_thr': 1.0e-8,
        'kspacing': 0.3,  # Angstrom^-1
        'mixing_beta': 0.5,
        'electron_maxstep': 200,
        'occupations': 'smearing',  # Required for smearing
        'smearing': 'cold',  # Better than gaussian for metals
        'degauss': 0.015,  # Ry
    },
    'accurate': {
        'ecutwfc': 80.0,
        'ecutrho': 640.0,
        'conv_thr': 1.0e-10,
        'kspacing': 0.15,  # Angstrom^-1
        'mixing_beta': 0.3,
        'electron_maxstep': 300,
        'occupations': 'smearing',  # Required for smearing
        'smearing': 'cold',  # Better than gaussian for metals
        'degauss': 0.01,  # Ry (tighter for accurate)
    }
}


class CalculationWorkflow:
    """
    Simplified workflow for running Quantum ESPRESSO calculations.
    
    This class provides an easy interface for:
    - Reading structures from CIF files
    - Setting up calculations with protocol presets
    - Running SCF or relaxation calculations
    - Using k-spacing instead of explicit k-points
    
    Examples:
        >>> # SCF calculation from CIF file
        >>> workflow = CalculationWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     protocol='moderate'
        ... )
        >>> workflow.run_scf(label='scf/silicon')
        
        >>> # Relax calculation with custom k-spacing
        >>> workflow = CalculationWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     protocol='fast',
        ...     kspacing=0.4
        ... )
        >>> workflow.run_relax(label='relax/silicon')
    """
    
    def __init__(
        self,
        atoms: Atoms,
        protocol: str = 'moderate',
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        kspacing: Optional[float] = None,
        input_data: Optional[Dict] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        expand_cell: bool = False,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize a calculation workflow.
        
        Args:
            atoms: ASE Atoms object representing the structure
            protocol: Protocol preset: 'fast', 'moderate', or 'accurate' (default: 'moderate')
            pseudopotentials: Optional dictionary mapping element symbols to pseudopotential filenames.
                            Either this or pseudopotentials_config must be provided.
                            If pseudopotentials_config is also provided, it takes precedence.
                            Example: {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe-n.UPF'}
            pseudopotentials_config: Optional name of pseudopotentials configuration to load.
                                   The workflow automatically extracts only the pseudopotentials
                                   needed for elements present in the structure.
                                   Configuration must exist in ~/.xespresso/pseudopotentials/
                                   Example: 'SSSP_efficiency' loads from SSSP_efficiency.json
            kspacing: K-point spacing in Angstrom^-1 (physical units). If None, uses preset value.
                     The workflow automatically handles the 2π normalization when converting to k-points.
                     Example: kspacing=0.20 will give the same k-points as
                     ase.io.espresso.kspacing_to_grid(atoms, 0.20/(2*np.pi))
            input_data: Additional input parameters (merged with preset)
            magnetic_config: Magnetic configuration. Can be:
                           - 'ferro' or 'ferromagnetic': All atoms ferromagnetic
                           - 'antiferro' or 'antiferromagnetic': Alternating spin
                           - Dict: Element-based config, e.g. {'Fe': [1, -1], 'O': [0]}
                           Also supports Hubbard parameters in the dict format
            expand_cell: If True, expand cell to accommodate magnetic configuration
            queue: Queue configuration dictionary for job submission (local or remote).
                   This is directly passed to the Espresso calculator.
            machine: Name of a machine configuration to load from ~/.xespresso/machines/.
                    If provided, the machine configuration is loaded and converted to a queue dict.
                    Cannot be used together with 'queue' parameter.
            code_version: Optional Quantum ESPRESSO version to use (e.g., '7.2', '6.8').
                         If provided with machine parameter, automatically loads code configuration
                         for that version from ~/.xespresso/codes/ and extracts modules to add to queue.
                         This enables using different QE versions on the same machine.
                         Example: code_version='7.2' loads modules like 'quantum-espresso/7.2'
            **kwargs: Additional parameters passed to Espresso calculator
        """
        self.atoms = atoms.copy()  # Work with a copy to avoid modifying original
        self.protocol = protocol
        self.extra_kwargs = kwargs
        self.expand_cell = expand_cell
        self.pseudopotentials_base_path = None  # Will be set if loading from config
        self._pseudo_config = None  # Will store config object if loaded from config
        
        # Handle pseudopotentials: either config name or explicit dict (config takes precedence)
        if pseudopotentials_config is not None:
            # Load from config file and extract only needed elements
            pseudopotentials = self._load_pseudopotentials_from_config(pseudopotentials_config)
        elif pseudopotentials is None:
            raise ValueError(
                "Must provide either 'pseudopotentials' dictionary or "
                "'pseudopotentials_config' name to load from ~/.xespresso/pseudopotentials/"
            )
        
        self.original_pseudopotentials = pseudopotentials
        
        # Handle machine and queue configuration
        if queue is not None and machine is not None:
            raise ValueError(
                "Cannot specify both 'queue' and 'machine' parameters. "
                "Use 'queue' for direct configuration or 'machine' to load from config."
            )
        
        if machine is not None:
            # Load machine configuration
            self.queue = load_machine(machine_name=machine)
            
            # Load code configuration and extract modules for specified version
            if code_version is not None:
                self._merge_code_modules_into_queue(machine, code_version)
        else:
            self.queue = queue
        
        # Get preset configuration
        if protocol not in PRESETS:
            raise ValueError(
                f"Protocol must be one of {list(PRESETS.keys())}, got '{protocol}'"
            )
        
        self.preset = PRESETS[protocol].copy()
        
        # Override k-spacing if provided
        if kspacing is not None:
            self.preset['kspacing'] = kspacing
        
        # Initialize input_data early so it can be used in magnetic config
        self.input_data = self.preset.copy()
        if input_data:
            self.input_data.update(input_data)
        
        # Handle magnetic configuration if provided
        if magnetic_config is not None:
            self._apply_magnetic_config(magnetic_config)
        else:
            self.pseudopotentials = pseudopotentials
        
        # Remove kspacing from input_data as it will be converted to kpts
        self.kspacing = self.input_data.pop('kspacing', None)
        
        # Ensure ecutwfc meets pseudopotential recommendations
        self._adjust_ecutwfc_for_pseudos()
    
    def _load_pseudopotentials_from_config(self, config_name: str) -> Dict[str, str]:
        """
        Load pseudopotentials from a configuration file and extract only the
        elements present in the structure.
        
        Also stores the base_path for later use in finding pseudopotential files.
        
        Args:
            config_name: Name of the pseudopotentials configuration to load
                        (e.g., 'SSSP_efficiency' loads from ~/.xespresso/pseudopotentials/SSSP_efficiency.json)
        
        Returns:
            Dictionary mapping element symbols to pseudopotential filenames,
            filtered to only include elements in self.atoms
        
        Raises:
            ValueError: If configuration not found or required elements not in config
        """
        config = load_pseudopotentials_config(config_name)
        
        if config is None:
            raise ValueError(
                f"Pseudopotentials configuration '{config_name}' not found. "
                f"Please save it first using create_pseudopotentials_config(). "
                f"Expected location: ~/.xespresso/pseudopotentials/{config_name}.json"
            )
        
        # Store base_path for remote execution and pseudopotential lookup
        if hasattr(config, 'base_path'):
            self.pseudopotentials_base_path = config.base_path
            logger.info(f"DEBUG: Loaded pseudopotentials_base_path = {self.pseudopotentials_base_path}")
        else:
            self.pseudopotentials_base_path = None
            logger.warning(f"DEBUG: config has no base_path attribute!")
            logger.warning(f"DEBUG: config object = {config}")
            logger.warning(f"DEBUG: config type = {type(config)}")
        
        # Get unique elements in the atomic structure
        elements_in_atoms = set(self.atoms.get_chemical_symbols())
        
        # Extract pseudopotentials for required elements
        pseudopotentials = {}
        missing_elements = []
        
        for element in elements_in_atoms:
            pseudo_obj = config.get_pseudopotential(element)
            if pseudo_obj is not None:
                # Store ONLY filename for transfer logic in _transfer_pseudopotentials
                # The full path will be handled separately via ESPRESSO_PSEUDO env var
                filename = pseudo_obj.filename
                pseudopotentials[element] = filename
                logger.info(f"DEBUG: Extracted {element} -> {filename}")
            else:
                missing_elements.append(element)
        
        if missing_elements:
            available = config.list_elements()
            raise ValueError(
                f"Pseudopotentials configuration '{config_name}' is missing "
                f"the following elements: {missing_elements}. "
                f"Available elements: {available}"
            )
        
        # Store config object for later access to suggested_ecutwfc
        self._pseudo_config = config
        
        return pseudopotentials
    
    def _get_min_ecut_from_pseudos(self) -> Optional[float]:
        """
        Get the minimum suggested ecutwfc from all pseudopotentials in this workflow.
        
        Returns the MAXIMUM suggested_ecutwfc from all elements in the structure,
        ensuring compatibility with all pseudopotentials.
        
        Returns:
            Maximum suggested ecutwfc in Ry, or None if not available from any pseudo
        """
        max_suggested_ecut = None
        
        # Try to get from config object first (if loaded from config)
        if self._pseudo_config is not None:
            elements_in_atoms = set(self.atoms.get_chemical_symbols())
            for element in elements_in_atoms:
                pseudo_obj = self._pseudo_config.get_pseudopotential(element)
                if pseudo_obj and hasattr(pseudo_obj, 'suggested_ecutwfc') and pseudo_obj.suggested_ecutwfc:
                    ecut = pseudo_obj.suggested_ecutwfc
                    if max_suggested_ecut is None or ecut > max_suggested_ecut:
                        max_suggested_ecut = ecut
                        logger.debug(f"Element {element}: suggested_ecutwfc = {ecut} Ry")
        
        if max_suggested_ecut is not None:
            logger.info(f"Minimum required ecutwfc from pseudopotentials: {max_suggested_ecut} Ry")
        else:
            logger.debug("No suggested_ecutwfc information available from pseudopotentials")
        
        return max_suggested_ecut
    
    def _get_ecutrho_ratio_for_pseudos(self) -> float:
        """
        Determine the appropriate ecutrho/ecutwfc ratio based on pseudopotential types.
        
        - Norm-Conserving (NC): ratio = 4.0
        - Ultrasoft (US): ratio >= 8.0
        - PAW: ratio >= 8.0
        
        Returns the MAXIMUM ratio needed to ensure compatibility with all pseudopotentials.
        
        Returns:
            Ratio (ecutrho/ecutwfc) - default 4.0 if no pseudo_config
        """
        if self._pseudo_config is None:
            return 4.0  # Default ratio for NC pseudos
        
        max_ratio = 4.0  # Start with NC ratio
        
        try:
            elements_in_atoms = set(self.atoms.get_chemical_symbols())
            for element in elements_in_atoms:
                pseudo_obj = self._pseudo_config.get_pseudopotential(element)
                if pseudo_obj and hasattr(pseudo_obj, 'type') and pseudo_obj.type:
                    pseudo_type = pseudo_obj.type.upper()
                    
                    # Determine ratio based on pseudopotential type
                    if 'PAW' in pseudo_type or 'PROJECTOR' in pseudo_type:
                        ratio = 8.0
                        type_label = 'PAW'
                    elif 'ULTRASOFT' in pseudo_type or 'US' in pseudo_type:
                        ratio = 8.0
                        type_label = 'Ultrasoft'
                    elif 'NORM-CONSERVING' in pseudo_type or 'NC' in pseudo_type or 'ONCV' in pseudo_type:
                        ratio = 4.0
                        type_label = 'Norm-Conserving'
                    else:
                        # Unknown type, use conservative ratio (8.0 is safer)
                        ratio = 8.0
                        type_label = f"Unknown ({pseudo_obj.type})"
                    
                    max_ratio = max(max_ratio, ratio)
                    logger.debug(f"Element {element}: type={type_label}, ratio={ratio}")
                else:
                    logger.debug(f"Element {element}: no type info available")
        except Exception as e:
            logger.debug(f"Could not determine ecutrho ratio from pseudos: {e}")
        
        if max_ratio > 4.0:
            logger.info(f"Using ecutrho/ecutwfc ratio: {max_ratio} (Ultrasoft/PAW pseudopotentials detected)")
        
        return max_ratio
    
    def _adjust_ecutwfc_for_pseudos(self):
        """
        Adjust ecutwfc in input_data to ensure it meets pseudopotential recommendations.
        
        If the current ecutwfc is below the maximum suggested value from pseudopotentials,
        increase it to meet the requirement and log a warning.
        Also adjusts ecutrho based on the pseudopotential type.
        """
        min_ecut = self._get_min_ecut_from_pseudos()
        
        if min_ecut is None:
            # No suggested values available, use preset as-is
            return
        
        current_ecut = self.input_data.get('ecutwfc', self.preset.get('ecutwfc'))
        
        if current_ecut < min_ecut:
            logger.warning(
                f"Protocol '{self.protocol}' specifies ecutwfc={current_ecut} Ry, "
                f"but pseudopotentials require at least {min_ecut} Ry. "
                f"Adjusting ecutwfc to {min_ecut} Ry to ensure physical correctness."
            )
            self.input_data['ecutwfc'] = min_ecut
            
            # Adjust ecutrho based on pseudopotential type
            ratio = self._get_ecutrho_ratio_for_pseudos()
            adjusted_ecutrho = min_ecut * ratio
            logger.info(f"Adjusting ecutrho to {adjusted_ecutrho} Ry (ratio={ratio})")
            self.input_data['ecutrho'] = adjusted_ecutrho
        else:
            logger.debug(
                f"ecutwfc={current_ecut} Ry meets pseudopotential requirement (min: {min_ecut} Ry)"
            )
    
    def _merge_code_modules_into_queue(self, machine_name: str, code_version: str):
        """
        Load code configuration for a specific machine and version,
        extract modules, and merge them into the queue configuration.
        
        This enables using different QE versions on the same machine by
        automatically loading version-specific modules.
        
        Args:
            machine_name: Name of the machine configuration
            code_version: QE version string (e.g., '7.2', '6.8')
        
        Raises:
            ValueError: If code configuration not found or version not available
        """
        try:
            codes_config = load_codes_config(machine_name, version=code_version, verbose=False)
        except Exception as e:
            raise ValueError(
                f"Could not load codes configuration for machine '{machine_name}' "
                f"with version '{code_version}': {e}"
            )
        
        if codes_config is None:
            raise ValueError(
                f"Codes configuration not found for machine '{machine_name}'. "
                f"Please create it using create_codes_config() or detect_qe_codes()."
            )
        
        # Validate that the requested version is available
        available_versions = codes_config.list_versions() if hasattr(codes_config, 'list_versions') else []
        if available_versions and code_version not in available_versions:
            raise ValueError(
                f"QE version '{code_version}' not available for machine '{machine_name}'. "
                f"Available versions: {', '.join(available_versions)}"
            )
        
        # Extract modules for this version
        modules = None
        
        # Try to get modules from version-specific configuration
        if hasattr(codes_config, 'versions') and codes_config.versions:
            if code_version in codes_config.versions:
                version_config = codes_config.versions[code_version]
                modules = version_config.get('modules')
        
        # Fallback to top-level modules if no version-specific ones found
        if modules is None and hasattr(codes_config, 'modules'):
            modules = codes_config.modules
        
        # Merge modules into queue configuration if found
        if modules:
            if self.queue is None:
                self.queue = {}
            
            # Ensure queue is a dict
            if not isinstance(self.queue, dict):
                # Convert to dict if it's an object with to_queue() method
                if hasattr(self.queue, 'to_queue'):
                    self.queue = self.queue.to_queue()
                else:
                    self.queue = {}
            
            # Add modules to queue
            self.queue['use_modules'] = True
            self.queue['modules'] = modules
    
    def _apply_magnetic_config(self, magnetic_config: Union[str, Dict]):
        """Apply magnetic configuration using setup_magnetic_config."""
        from xespresso.tools import set_ferromagnetic, set_antiferromagnetic
        
        if isinstance(magnetic_config, str):
            magnetic_config = magnetic_config.lower()
            if magnetic_config in ['ferro', 'ferromagnetic']:
                # Simple ferromagnetic configuration
                config = set_ferromagnetic(
                    self.atoms, 
                    magnetic_moment=1.0, 
                    pseudopotentials=self.original_pseudopotentials
                )
                # atoms modified in-place, just update config
                self.pseudopotentials = config.get('pseudopotentials', self.original_pseudopotentials.copy())
                if 'input_ntyp' in config:
                    if 'input_ntyp' not in self.input_data:
                        self.input_data['input_ntyp'] = {}
                    self.input_data['input_ntyp'].update(config['input_ntyp'])
                # Set nspin=2 for magnetic calculation
                self.input_data['nspin'] = 2
            elif magnetic_config in ['antiferro', 'antiferromagnetic']:
                # Simple antiferromagnetic configuration
                # For antiferromagnetic, we need to determine sublattices
                # Check if we have enough atoms for alternating configuration
                n_atoms = len(self.atoms)
                symbols = self.atoms.get_chemical_symbols()
                unique_symbols = set(symbols)
                
                # Count atoms of each element
                symbol_counts = {s: symbols.count(s) for s in unique_symbols}
                
                # Check if antiferro is feasible with current structure
                # If we have only 1 atom of a kind, we need to expand to proper BCC structure
                needs_expansion = any(count < 2 for count in symbol_counts.values()) and len(unique_symbols) == 1
                
                if needs_expansion:
                    # For BCC structure with antiferro, we need 2 atoms at proper positions
                    # Átomo 1: (0, 0, 0) - corner
                    # Átomo 2: (a/2, a/2, a/2) - body center
                    logger.info(f"Antiferromagnetic BCC requires proper 2-atom structure. Reconstructing...")
                    
                    # Get cell parameters from current (primitive) cell
                    cell = self.atoms.get_cell()
                    cell_volume = self.atoms.get_volume()
                    
                    # For BCC, the conventional cell is 2x the primitive cell
                    # Create proper BCC with 2 atoms
                    symbol = symbols[0]  # Get element symbol
                    
                    # Create new atoms with proper BCC structure
                    bcc_atoms = Atoms(
                        symbols=[symbol, symbol],
                        positions=[[0, 0, 0], [cell[0, 0]/2, cell[1, 1]/2, cell[2, 2]/2]],
                        cell=cell,
                        pbc=True
                    )
                    self.atoms = bcc_atoms
                    logger.info(f"  BCC structure reconstructed: {symbol} at (0,0,0) and ({cell[0, 0]/2:.4f}, {cell[1, 1]/2:.4f}, {cell[2, 2]/2:.4f})")
                    n_atoms = len(self.atoms)
                
                # Now create sublattices from alternating atoms
                sublattice1 = list(range(0, n_atoms, 2))
                sublattice2 = list(range(1, n_atoms, 2))
                
                config = set_antiferromagnetic(
                    self.atoms,
                    sublattice_indices=[sublattice1, sublattice2],
                    magnetic_moment=1.0,
                    pseudopotentials=self.original_pseudopotentials
                )
                # atoms modified in-place, just update config
                self.pseudopotentials = config.get('pseudopotentials', self.original_pseudopotentials.copy())
                if 'input_ntyp' in config:
                    if 'input_ntyp' not in self.input_data:
                        self.input_data['input_ntyp'] = {}
                    self.input_data['input_ntyp'].update(config['input_ntyp'])
                # Set nspin=2 for magnetic calculation
                self.input_data['nspin'] = 2
            else:
                raise ValueError(
                    f"Unknown magnetic configuration: '{magnetic_config}'. "
                    "Use 'ferro', 'antiferro', or a dict with element-based config."
                )
        elif isinstance(magnetic_config, dict):
            # Element-based configuration with possible Hubbard parameters
            config = setup_magnetic_config(
                self.atoms,
                magnetic_config,
                pseudopotentials=self.original_pseudopotentials,
                expand_cell=self.expand_cell
            )
            self.atoms = config['atoms']
            self.pseudopotentials = config.get('pseudopotentials', self.original_pseudopotentials)
            
            # Merge special input_data from magnetic config
            if 'input_ntyp' in config:
                if 'input_ntyp' not in self.input_data:
                    self.input_data['input_ntyp'] = {}
                self.input_data['input_ntyp'].update(config['input_ntyp'])
            
            # Handle Hubbard parameters in new format
            if 'hubbard' in config:
                self.input_data['hubbard'] = config['hubbard']
            if 'hubbard_v' in config:
                self.input_data['hubbard_v'] = config['hubbard_v']
            if 'qe_version' in config:
                self.input_data['qe_version'] = config.get('qe_version')
            if 'lda_plus_u' in config:
                self.input_data['lda_plus_u'] = config['lda_plus_u']
            
            # Set nspin=2 for polarized magnetic calculation
            self.input_data['nspin'] = 2
        else:
            raise TypeError(
                f"magnetic_config must be str or dict, got {type(magnetic_config)}"
            )
    
    @classmethod
    def from_cif(
        cls,
        cif_file: Union[str, Path],
        protocol: str = 'moderate',
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        kspacing: Optional[float] = None,
        input_data: Optional[Dict] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        expand_cell: bool = False,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        **kwargs
    ) -> 'CalculationWorkflow':
        """
        Create a workflow from a CIF file.
        
        Args:
            cif_file: Path to CIF file
            protocol: Protocol preset: 'fast', 'moderate', or 'accurate' (default: 'moderate')
            pseudopotentials: Optional dictionary mapping element symbols to pseudopotential filenames.
                            Either this or pseudopotentials_config must be provided.
            pseudopotentials_config: Optional name of pseudopotentials configuration to load.
                                   The workflow automatically extracts elements needed from the CIF file.
            kspacing: K-point spacing in Angstrom^-1 (physical units)
            input_data: Additional input parameters
            magnetic_config: Magnetic configuration ('ferro', 'antiferro', or element dict)
            expand_cell: If True, expand cell to accommodate magnetic configuration
            queue: Queue configuration dictionary for job submission
            machine: Name of a machine configuration to load
            code_version: Optional Quantum ESPRESSO version to use with this machine
            **kwargs: Additional parameters passed to Espresso calculator
            
        Returns:
            CalculationWorkflow: Initialized workflow object
        """
        atoms = read(str(cif_file))
        return cls(
            atoms,
            protocol=protocol,
            pseudopotentials=pseudopotentials,
            pseudopotentials_config=pseudopotentials_config,
            kspacing=kspacing,
            input_data=input_data,
            magnetic_config=magnetic_config,
            expand_cell=expand_cell,
            queue=queue,
            machine=machine,
            code_version=code_version,
            **kwargs
        )
    
    def _get_kpts(self) -> Union[Tuple[int, int, int], str]:
        """
        Calculate k-points from k-spacing using ase.io.espresso.kspacing_to_grid.
        
        Returns:
            Tuple of k-points or 'gamma'
        """
        if self.kspacing is not None:
            # Convert kspacing to k-point grid
            # Note: kspacing_to_grid expects spacing in units of 2*pi/Angstrom
            # So we need to convert from Angstrom^-1
            kpts = kspacing_to_grid(self.atoms, self.kspacing / (2 * np.pi))
            return tuple(kpts)
        else:
            # Default to gamma point if no k-spacing specified
            return (1, 1, 1)
    
    def _check_convergence(self, calc: Espresso, calculation_type: str = 'scf') -> dict:
        """
        Check if calculation converged and inform user.
        
        Args:
            calc: Espresso calculator object
            calculation_type: 'scf', 'relax', 'vc-relax', etc
            
        Returns:
            dict: Convergence status information with keys:
                - job_done: bool, JOB DONE found in output
                - scf_converged: bool, convergence achieved
                - scf_iterations: int or None, number of SCF cycles
                - final_energy: float or None, in Ry
                - message: str, formatted status message
        """
        import re
        
        convergence_info = {
            'job_done': False,
            'scf_converged': False,
            'scf_iterations': None,
            'final_energy': None,
            'message': None,
        }
        
        try:
            output = calc.results.get('output', '')
            if not output:
                convergence_info['message'] = "⚠ No output available to check convergence"
                print("\n" + "="*70)
                print(f"CALCULATION STATUS: {calculation_type.upper()}")
                print("="*70)
                print(convergence_info['message'])
                print("="*70)
                return convergence_info
            
            # Check basic markers
            convergence_info['job_done'] = 'JOB DONE' in output
            convergence_info['scf_converged'] = 'convergence has been achieved' in output
            
            # Extract SCF iterations - try multiple patterns
            # Pattern 1: "convergence has been achieved in X iterations"
            scf_match = re.search(r'convergence has been achieved in\s+(\d+)\s+iterations?', output)
            if scf_match:
                convergence_info['scf_iterations'] = int(scf_match.group(1))
            else:
                # Pattern 2: "number of scf cycles = X" (fallback)
                scf_match = re.search(r'number of scf cycles\s*=\s*(\d+)', output)
                if scf_match:
                    convergence_info['scf_iterations'] = int(scf_match.group(1))
            
            # Extract final energy - try multiple patterns
            # Pattern 1: "!    total energy              =     X Ry"
            energy_match = re.search(r'!\s+total energy\s*=\s*(-?\d+\.\d+)\s*Ry', output)
            if energy_match:
                convergence_info['final_energy'] = float(energy_match.group(1))
            else:
                # Pattern 2: "Final energy = X Ry" (fallback)
                energy_match = re.search(r'Final energy\s*=\s*(-?\d+\.\d+)\s*Ry', output)
                if energy_match:
                    convergence_info['final_energy'] = float(energy_match.group(1))
            
            # Generate status message
            if convergence_info['job_done'] and convergence_info['scf_converged']:
                convergence_info['message'] = (
                    f"✓ {calculation_type.upper()} calculation CONVERGED successfully"
                )
            elif convergence_info['job_done'] and not convergence_info['scf_converged']:
                convergence_info['message'] = (
                    f"⚠ WARNING: {calculation_type.upper()} completed but SCF did NOT converge!\n"
                    f"  Iterations: {convergence_info['scf_iterations']} | "
                    f"Try increasing electron_maxstep or decreasing conv_thr"
                )
            else:
                convergence_info['message'] = (
                    f"✗ {calculation_type.upper()} calculation FAILED or incomplete"
                )
            
            # Print status to user
            print("\n" + "="*70)
            print(f"CALCULATION STATUS: {calculation_type.upper()}")
            print("="*70)
            print(f"Job completed (JOB DONE):    {convergence_info['job_done']}")
            print(f"SCF converged:               {convergence_info['scf_converged']}")
            print(f"SCF iterations:              {convergence_info['scf_iterations']}")
            print(f"Final energy (Ry):           {convergence_info['final_energy']}")
            print("-"*70)
            print(convergence_info['message'])
            print("="*70)
            
        except Exception as e:
            logger.warning(f"Could not parse convergence info: {e}")
            convergence_info['message'] = f"⚠ Could not verify convergence: {e}"
            print(f"\n⚠ Warning: {convergence_info['message']}")
        
        return convergence_info

    def _monitor_remote_job(self, calc: Espresso, job_id: str, timeout: int = 3600, poll_interval: int = 30) -> dict:
        """
        Monitor SLURM job status and detect stuck/failed jobs.
        
        Args:
            calc: Espresso calculator object
            job_id: SLURM job ID
            timeout: Maximum time to wait in seconds (default 3600s = 1 hour)
            poll_interval: Time between status checks in seconds (default 30s)
            
        Returns:
            dict: Job status info with keys:
                - job_id: str, SLURM job ID
                - state: str, final job state (RUNNING, COMPLETED, FAILED, etc)
                - elapsed_time: int, seconds elapsed
                - pending_time: int, seconds spent in PENDING state
                - reason: str, status reason or error message
                - success: bool, whether job completed successfully
                - message: str, formatted status message
        """
        import time
        import subprocess
        
        job_status = {
            'job_id': job_id,
            'state': None,
            'elapsed_time': 0,
            'pending_time': 0,
            'reason': None,
            'success': False,
            'message': None,
        }
        
        start_time = time.time()
        pending_start = None
        
        print(f"\n{'='*70}")
        print(f"REMOTE JOB MONITORING: {job_id}")
        print(f"{'='*70}")
        print(f"Timeout: {timeout}s | Poll interval: {poll_interval}s")
        print(f"{'-'*70}")
        
        try:
            while True:
                elapsed = int(time.time() - start_time)
                
                # Check if timeout exceeded
                if elapsed > timeout:
                    job_status['elapsed_time'] = elapsed
                    job_status['message'] = (
                        f"✗ JOB TIMEOUT: {job_id} exceeded {timeout}s limit\n"
                        f"  Job may still be running on remote. Check manually with: squeue -j {job_id}"
                    )
                    print(f"\n{job_status['message']}")
                    return job_status
                
                # Query SLURM status on remote system
                try:
                    # Get remote connection from calc
                    remote_conn = getattr(calc, 'remote', None)
                    if remote_conn is None:
                        raise ValueError("No remote connection available for job monitoring")
                    
                    # Execute squeue on remote system
                    stdout, stderr = remote_conn.run_command(
                        f"squeue -j {job_id} -h -o '%T,%r,%M'"
                    )
                    
                    # Check if job is still in queue
                    if not stdout.strip():
                        # Job not found in queue (probably completed) - check final status with sacct
                        try:
                            stdout_sacct, stderr_sacct = remote_conn.run_command(
                                f"sacct -j {job_id} --format=State -n -P"
                            )
                            if stdout_sacct.strip():
                                # Parse sacct output - get the last line (most recent state)
                                lines = stdout_sacct.strip().split('\n')
                                last_line = lines[-1] if lines else ""
                                state = last_line.split('|')[0].strip() if '|' in last_line else last_line.strip()
                                
                                job_status['state'] = state
                                job_status['elapsed_time'] = elapsed
                                
                                if state == 'COMPLETED':
                                    job_status['success'] = True
                                    job_status['message'] = f"✓ JOB {job_id} completed successfully"
                                elif state in ['FAILED', 'TIMEOUT', 'CANCELLED', 'OUT_OF_MEMORY']:
                                    job_status['success'] = False
                                    job_status['message'] = f"✗ JOB {job_id} {state}"
                                else:
                                    # Other states (RUNNING, PENDING shouldn't happen here)
                                    job_status['success'] = False
                                    job_status['message'] = f"? JOB {job_id} finished with state: {state}"
                            else:
                                # Could not get sacct status, assume completed for backward compatibility
                                job_status['state'] = 'COMPLETED'
                                job_status['success'] = True
                                job_status['message'] = f"✓ JOB {job_id} completed (status unknown)"
                        except Exception as sacct_e:
                            logger.warning(f"Could not check job status with sacct: {sacct_e}")
                            # Assume completed for backward compatibility
                            job_status['state'] = 'COMPLETED'
                            job_status['success'] = True
                            job_status['message'] = f"✓ JOB {job_id} completed (sacct failed)"
                        
                        print(f"\n{job_status['message']}")
                        return job_status
                    
                    # Parse squeue output: STATE,REASON,ELAPSED
                    output = stdout.strip()
                    if output:
                        parts = output.split(',')
                        state = parts[0].strip() if len(parts) > 0 else 'UNKNOWN'
                        reason = parts[1].strip() if len(parts) > 1 else ''
                        elapsed_str = parts[2].strip() if len(parts) > 2 else ''
                        
                        job_status['state'] = state
                        job_status['reason'] = reason
                        job_status['elapsed_time'] = elapsed
                        
                        # Print status update
                        status_line = f"[{elapsed:5d}s] State: {state:10s} | Reason: {reason:20s}"
                        print(f"\r{status_line}", end='', flush=True)
                        
                        # Check for problematic states
                        if state == 'FAILED':
                            job_status['message'] = (
                                f"✗ JOB FAILED: {job_id}\n"
                                f"  Reason: {reason}\n"
                                f"  Check output for details: {reason}"
                            )
                            print(f"\n\n{job_status['message']}")
                            return job_status
                        
                        elif state == 'CANCELLED':
                            job_status['message'] = (
                                f"✗ JOB CANCELLED: {job_id}\n"
                                f"  Reason: {reason}"
                            )
                            print(f"\n\n{job_status['message']}")
                            return job_status
                        
                        elif state == 'PENDING':
                            if pending_start is None:
                                pending_start = time.time()
                            pending_time = int(time.time() - pending_start)
                            job_status['pending_time'] = pending_time
                            
                            # Alert if pending for too long
                            if pending_time > 120:  # 2 minutes
                                print(f"\n\n⚠ WARNING: Job {job_id} stuck in PENDING for {pending_time}s")
                                print(f"  Reason: {reason}")
                                print(f"  Possible causes:")
                                print(f"    - Compute node is down")
                                print(f"    - Resource limit reached")
                                print(f"    - Queue misconfiguration")
                                print(f"  Manual check: squeue -j {job_id}")
                                print(f"  To cancel: scancel {job_id}")
                                
                                # Continue monitoring but alert user
                                user_input = input(f"\nContinue waiting? (y/n): ")
                                if user_input.lower() != 'y':
                                    # Try to cancel job remotely
                                    try:
                                        remote_conn.run_command(f"scancel {job_id}")
                                    except Exception as cancel_e:
                                        logger.warning(f"Failed to cancel job remotely: {cancel_e}")
                                    job_status['message'] = f"✗ JOB CANCELLED BY USER: {job_id}"
                                    return job_status
                        
                        elif state == 'RUNNING':
                            pending_start = None  # Reset pending timer
                        
                        elif state == 'COMPLETED':
                            job_status['success'] = True
                            job_status['message'] = f"✓ JOB {job_id} completed successfully"
                            print(f"\n\n{job_status['message']}")
                            return job_status
                
                except subprocess.TimeoutExpired:
                    print(f"\n⚠ Warning: squeue command timed out")
                except Exception as e:
                    logger.warning(f"Error querying job status: {e}")
                
                # Wait before next poll
                time.sleep(poll_interval)
        
        except KeyboardInterrupt:
            print(f"\n\n⚠ Monitoring interrupted by user")
            job_status['message'] = f"User interrupted monitoring for job {job_id}"
            return job_status
        except Exception as e:
            logger.error(f"Error during remote job monitoring: {e}")
            job_status['message'] = f"Error monitoring job: {e}"
            return job_status
    
    def submit_scf_batch(
        self,
        label: str = 'scf',
        wait_for_completion: bool = False,
        timeout: int = 3600,
        poll_interval: int = 30,
        **calc_kwargs
    ) -> Dict:
        """
        Submit a single SCF calculation without blocking (non-blocking mode).
        
        This method is designed for batch/parallel execution on HPC systems.
        It submits the job and returns immediately with job metadata.
        
        Args:
            label: Directory/label for the calculation
            wait_for_completion: If True, block until job completes. If False, return immediately after submission.
            timeout: Maximum time to wait for job completion (seconds), only if wait_for_completion=True
            poll_interval: Time between status checks (seconds), only if wait_for_completion=True
            **calc_kwargs: Additional parameters for the Espresso calculator
            
        Returns:
            Dict with keys:
                - 'calc': Espresso calculator object
                - 'job_id': str, SLURM job ID (if remote), None otherwise
                - 'label': str, calculation label
                - 'submitted': bool, True if submitted successfully
                - 'completed': bool, True if job completed (only if wait_for_completion=True)
        """
        if not self.queue or self.queue.get('execution') != 'remote':
            raise ValueError(
                "submit_scf_batch() is only supported for remote queue systems (SLURM). "
                "Use run_scf() for local calculations."
            )
        
        # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': 'scf',
            'input_data': self.input_data.copy(),
            'kpts': self._get_kpts(),
        }
        
        params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = self.input_data.get('ecutrho', 400.0)
        
        if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
            params['input_data']['pseudo_dir'] = './pseudo'
        
        if self.queue is not None:
            params['queue'] = self.queue
        
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create calculator
        calc = Espresso(**params)
        self.atoms.calc = calc
        self.last_calc = calc
        
        # Step 1: Write input
        calc.write_input(self.atoms)
        calc.atoms = self.atoms
        
        # Step 2: Execute (submits remotely, returns immediately)
        logger.info(f"Submitting SCF calculation: {label}")
        calc.execute()
        
        job_id = getattr(calc, 'last_job_id', None)
        logger.info(f"SCF calculation submitted with job ID: {job_id}")
        
        # Store remote connection for later monitoring
        if hasattr(calc, 'scheduler') and hasattr(calc.scheduler, 'remote'):
            calc.remote = calc.scheduler.remote
        
        result = {
            'calc': calc,
            'job_id': job_id,
            'label': label,
            'submitted': True,
            'completed': False,
        }
        
        # If requested, wait for job to complete
        if wait_for_completion:
            if job_id:
                monitor_result = self._monitor_remote_job(calc, job_id, timeout=timeout, poll_interval=poll_interval)
                result['completed'] = monitor_result['success']
                if monitor_result['success']:
                    # Retrieve output and read results
                    monitor = RemoteJobMonitor(calc)
                    if monitor.wait(timeout=60, poll_interval=5):
                        monitor.retrieve_output()
                        calc.read_results()
                        logger.info(f"Job {job_id} completed and output retrieved")
            else:
                logger.warning(f"Could not track job completion: no job ID available")
        
        return result
    
    def submit_scf_batch_multiple(
        self,
        parameter_sets: List[Dict],
        verbose: bool = True
    ) -> List[Dict]:
        """
        Submit multiple SCF calculations in batch mode (non-blocking).
        
        All jobs are submitted quickly without waiting for completion.
        Use wait_for_batch_jobs() to monitor and collect results.
        
        Args:
            parameter_sets: List of dicts, each containing:
                - 'label': str, unique label for this calculation
                - 'ecutwfc': float (optional), override ecutwfc
                - 'kspacing': float (optional), override kspacing
                - 'other_params': other parameters to pass to Espresso
            verbose: Print submission status
            
        Returns:
            List of dicts with submission info for each job:
                - 'calc': Espresso calculator
                - 'job_id': SLURM job ID
                - 'label': calculation label
                - 'submitted': bool
                - 'completed': bool
        """
        results = []
        
        if verbose:
            print(f"\nSubmitting {len(parameter_sets)} SCF calculations in batch mode...")
            print(f"{'-'*70}")
        
        for i, params in enumerate(parameter_sets):
            label = params.get('label', f'scf_{i}')
            
            # Create a modified workflow with the specific parameters
            input_data_override = self.input_data.copy()
            if 'ecutwfc' in params:
                input_data_override['ecutwfc'] = params['ecutwfc']
            if 'ecutrho' in params:
                input_data_override['ecutrho'] = params['ecutrho']
            
            # Add any other input data params from the parameter set
            for key, value in params.items():
                if key not in ['label', 'ecutwfc', 'ecutrho', 'kspacing']:
                    input_data_override[key] = value
            
            # Create new workflow with overridden parameters
            try:
                temp_workflow = CalculationWorkflow(
                    self.atoms,
                    protocol=self.protocol,
                    pseudopotentials=self.pseudopotentials,
                    kspacing=params.get('kspacing', self.preset.get('kspacing')),
                    input_data=input_data_override,
                    queue=self.queue,
                    **self.extra_kwargs
                )
                
                # Copy the pseudopotentials_base_path if it exists (for remote transfer)
                if hasattr(self, 'pseudopotentials_base_path'):
                    temp_workflow.pseudopotentials_base_path = self.pseudopotentials_base_path
                
                # Submit this calculation
                result = temp_workflow.submit_scf_batch(label=label, wait_for_completion=False)
                results.append(result)
                
                if verbose:
                    job_id_str = f" (job_id: {result['job_id']})" if result['job_id'] else ""
                    print(f"  [{i+1}/{len(parameter_sets)}] {label}{job_id_str}")
                    
            except Exception as e:
                logger.error(f"Error submitting calculation {label}: {e}")
                results.append({
                    'calc': None,
                    'job_id': None,
                    'label': label,
                    'submitted': False,
                    'completed': False,
                    'error': str(e),
                })
        
        if verbose:
            print(f"{'-'*70}")
            print(f"Submitted {sum(1 for r in results if r['submitted'])}/{len(parameter_sets)} calculations\n")
        
        return results
    
    def wait_for_batch_jobs(
        self,
        batch_results: List[Dict],
        timeout: int = 3600,
        poll_interval: int = 30,
        verbose: bool = True
    ) -> List[Dict]:
        """
        Monitor and wait for multiple submitted jobs to complete.
        
        Args:
            batch_results: List of dicts returned from submit_scf_batch_multiple()
            timeout: Maximum time to wait (seconds)
            poll_interval: Time between status checks (seconds)
            verbose: Print monitoring progress
            
        Returns:
            List of dicts with completion status for each job:
                - 'label': str
                - 'job_id': str
                - 'completed': bool
                - 'success': bool
                - 'energy': float (if successful and readable)
                - 'error': str (if failed)
        """
        import time
        
        results = []
        active_jobs = {r['job_id']: r for r in batch_results if r['job_id']}
        
        if verbose:
            print(f"\nMonitoring {len(active_jobs)} remote jobs (timeout: {timeout}s, poll: {poll_interval}s)...")
            print(f"{'-'*70}")
        
        start_time = time.time()
        completed_count = 0
        
        while active_jobs and (time.time() - start_time) < timeout:
            # Query status of all jobs at once
            if verbose:
                elapsed = int(time.time() - start_time)
                print(f"[{elapsed}s] Checking {len(active_jobs)} jobs...")
            
            jobs_to_remove = []
            
            for job_id, job_info in active_jobs.items():
                calc = job_info['calc']
                label = job_info['label']
                
                try:
                    # Get job status
                    remote_conn = getattr(calc, 'remote', None)
                    if not remote_conn:
                        logger.warning(f"No remote connection for {label}")
                        jobs_to_remove.append(job_id)
                        continue
                    
                    # Check SLURM status
                    stdout, _ = remote_conn.run_command(f"squeue -j {job_id} -h")
                    
                    if not stdout.strip():
                        # Job completed, check final status
                        stdout_sacct, _ = remote_conn.run_command(
                            f"sacct -j {job_id} --format=State -n -P"
                        )
                        
                        lines = stdout_sacct.strip().split('\n') if stdout_sacct.strip() else []
                        state = lines[-1].split('|')[0].strip() if lines else 'UNKNOWN'
                        
                        success = state == 'COMPLETED'
                        
                        # Try to retrieve energy
                        energy = None
                        try:
                            monitor = RemoteJobMonitor(calc)
                            if monitor.wait(timeout=30, poll_interval=5):
                                monitor.retrieve_output()
                                calc.read_results()
                                energy = calc.results.get('energy')
                        except Exception as e:
                            logger.debug(f"Could not retrieve results for {label}: {e}")
                        
                        results.append({
                            'label': label,
                            'job_id': job_id,
                            'completed': True,
                            'success': success,
                            'state': state,
                            'energy': energy,
                        })
                        
                        jobs_to_remove.append(job_id)
                        completed_count += 1
                        
                        if verbose:
                            status = "✓" if success else "✗"
                            print(f"  {status} {label}: {state}")
                
                except Exception as e:
                    logger.warning(f"Error checking status for {label}: {e}")
            
            # Remove completed jobs
            for job_id in jobs_to_remove:
                del active_jobs[job_id]
            
            # Wait before next check
            if active_jobs:
                time.sleep(poll_interval)
        
        # Timeout or all completed
        if active_jobs:
            if verbose:
                print(f"\n⚠ Timeout reached with {len(active_jobs)} jobs still running")
            
            for job_id, job_info in active_jobs.items():
                results.append({
                    'label': job_info['label'],
                    'job_id': job_id,
                    'completed': False,
                    'success': False,
                    'error': 'Timeout',
                })
        
        if verbose:
            print(f"{'-'*70}")
            print(f"Completed {completed_count + len([r for r in results if r['completed']])}/{len(batch_results)} jobs\n")
        
        return results

    def run_scf(
        self,
        label: str = 'scf',
        **calc_kwargs
    ) -> Espresso:
        """
        Run a self-consistent field (SCF) calculation.
        
        Args:
            label: Directory/label for the calculation
            **calc_kwargs: Additional parameters for the Espresso calculator
            
        Returns:
            Espresso: Calculator object with results
        """
        # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
        # This allows remote_mixin to find pseudopotentials
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': 'scf',
            'input_data': self.input_data.copy(),
            'kpts': self._get_kpts(),
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = self.input_data.get('ecutrho', 400.0)
        
        # Set pseudo_dir when using pseudopotentials_config
        if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
            params['input_data']['pseudo_dir'] = './pseudo'
        
        # Add queue configuration if provided
        if self.queue is not None:
            params['queue'] = self.queue
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create calculator
        calc = Espresso(**params)
        self.atoms.calc = calc
        self.last_calc = calc  # Track last calculator for monitoring
        
        # If remote non-blocking: control execution steps to avoid retry loop
        if self.queue and self.queue.get('execution') == 'remote' and not self.queue.get('wait_for_completion', False):
            logger.info("Remote non-blocking: executing with automatic job monitoring...")
            
            # Step 1: Write input (with atoms, so _transfer_pseudopotentials won't need to call it again)
            calc.write_input(self.atoms)
            
            # IMPORTANT: Set calc.atoms so _transfer_pseudopotentials() can use it if needed
            calc.atoms = self.atoms
            
            # Step 2: Execute (submits job remotely)
            calc.execute()
            
            # IMPORTANT: Store remote connection on calc for RemoteJobMonitor to access
            if hasattr(calc, 'scheduler') and hasattr(calc.scheduler, 'remote'):
                calc.remote = calc.scheduler.remote
            
            # Step 3: Get job ID and monitor SLURM job status
            job_id = getattr(calc, 'last_job_id', None)
            
            if job_id is None:
                raise RuntimeError(
                    "Remote job submission failed: No job ID returned from scheduler. "
                    "Please check scheduler configuration and job submission logs."
                )
            
            logger.info(f"Remote job {job_id} submitted. Monitoring SLURM status...")
            timeout = self.queue.get('job_timeout', 3600)
            job_monitor_result = self._monitor_remote_job(calc, job_id, timeout=timeout, poll_interval=30)
            
            if not job_monitor_result['success']:
                raise RuntimeError(f"Remote job {job_id} failed: {job_monitor_result['message']}")
            
            # Step 4: Job completed in queue, now fetch output using RemoteJobMonitor
            monitor = RemoteJobMonitor(calc)
            if monitor.wait(timeout=60, poll_interval=5):  # Short timeout since job already completed
                monitor.retrieve_output()
                logger.info("Remote job output retrieved.")
                # Step 5: Read results
                calc.read_results()
            else:
                raise RuntimeError(f"Failed to retrieve output for job {job_id}")
        else:
            # Local or remote blocking: use normal run() with retry logic
            calc.run(atoms=self.atoms)
        
        # Check convergence and inform user
        self._check_convergence(calc, calculation_type='scf')
        
        return calc
    
    def run_nscf(
        self,
        label: str = 'nscf',
        kpts: tuple = None,
        nbnd: int = None,
        wf_collect: bool = True,
        npools: int = None,
        **calc_kwargs
    ) -> Espresso:
        """
        Run a non-self-consistent field (NSCF) calculation for band structure.
        
        NSCF reads the charge density from a previous SCF calculation and computes
        electronic structure on a denser k-point mesh without updating electron density.
        Essential for Wannier interpolation and band structure analysis.
        
        Args:
            label: Directory/label for the calculation
            kpts: K-point mesh tuple (e.g., (12, 12, 12) for dense mesh)
                  If None, uses self.kpts from initialization
            nbnd: Number of bands to compute (must be > n_electrons/2)
                  If None, uses preset value or estimates from pseudopotentials
            wf_collect: If True, collect wavefunctions on each k-point (required for Wannier)
            npools: Number of k-point pools for parallelization (e.g., -npools 4)
                    Allows distributing k-points across processes
            **calc_kwargs: Additional Espresso calculator parameters
            
        Returns:
            Espresso: Calculator object with NSCF results
            
        Notes:
            - Must run SCF first to generate charge density
            - Generates prefix.save/wavefunction.* files (can be large!)
            - Set wf_collect=True for Wannier calculations
            - High nbnd increases memory but necessary for accurate interpolation
        """
        if kpts is None:
            kpts = self._get_kpts()
        
        if nbnd is None:
            # Estimate nbnd from pseudopotentials if not specified
            nbnd = self._estimate_nbnd()
        
        # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Prepare input_data (copy from preset)
        input_data = self.input_data.copy()
        
        # NSCF-specific parameters
        input_data['nbnd'] = nbnd  # Override with larger value for band structure
        
        if wf_collect:
            input_data['wf_collect'] = True  # Collect wavefunctions for post-processing
        
        if npools is not None:
            # Note: -npools is a command-line argument, not in &control
            # Will be passed as extra kwargs to Espresso
            calc_kwargs['npools'] = npools
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': 'nscf',
            'input_data': input_data,
            'kpts': kpts,
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = input_data.get('ecutrho', 400.0)
        
        # Set pseudo_dir when using pseudopotentials_config
        if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
            params['input_data']['pseudo_dir'] = './pseudo'
        
        # Add queue configuration if provided
        if self.queue is not None:
            params['queue'] = self.queue
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create calculator
        calc = Espresso(**params)
        self.atoms.calc = calc
        self.last_calc = calc  # Track last calculator for monitoring
        
        # If remote non-blocking: control execution steps to avoid retry loop
        if self.queue and self.queue.get('execution') == 'remote' and not self.queue.get('wait_for_completion', False):
            logger.info("Remote non-blocking: executing NSCF with automatic job monitoring...")
            
            # Step 1: Write input
            calc.write_input(self.atoms)
            calc.atoms = self.atoms
            
            # Step 2: Execute (submits job remotely)
            calc.execute()
            
            # IMPORTANT: Store remote connection on calc for RemoteJobMonitor to access
            if hasattr(calc, 'scheduler') and hasattr(calc.scheduler, 'remote'):
                calc.remote = calc.scheduler.remote
            
            # Step 3: Monitor SLURM job status
            logger.info(f"Remote job {calc.last_job_id} submitted. Monitoring SLURM status...")
            job_id = calc.last_job_id
            timeout = self.queue.get('job_timeout', 7200)  # NSCF may take longer
            job_monitor_result = self._monitor_remote_job(calc, job_id, timeout=timeout, poll_interval=30)
            
            if not job_monitor_result['success']:
                raise RuntimeError(f"Remote job {job_id} failed: {job_monitor_result['message']}")
            
            # Step 4: Fetch output
            monitor = RemoteJobMonitor(calc)
            if monitor.wait(timeout=60, poll_interval=5):
                monitor.retrieve_output()
                logger.info("Remote NSCF output retrieved.")
                calc.read_results()
            else:
                raise RuntimeError(f"Failed to retrieve NSCF output for job {job_id}")
        else:
            # Local or remote blocking: use normal run()
            calc.run(atoms=self.atoms)
        
        # Check convergence and inform user
        self._check_convergence(calc, calculation_type='nscf')
        
        logger.info(f"NSCF calculation completed. Wavefunctions saved in {label}/")
        if wf_collect:
            logger.info(f"  wf_collect=True: Wavefunctions available for post-processing (Wannier, bands, etc)")
        
        return calc
    
    def run_dos(
        self,
        nscf_label: str = 'nscf',
        dos_label: str = 'dos',
        Emin: Optional[float] = None,
        Emax: Optional[float] = None,
        DeltaE: float = 0.01,
        degauss: Optional[float] = None,
        ngauss: int = 0,
        pdos: bool = False,
    ):
        """
        Run Density of States (DOS) calculation with spin polarization support.
        
        This is a post-processing step that requires a prior NSCF calculation.
        The DOS is calculated from the electron density converged in NSCF.
        
        For magnetic systems (nspin=2 or nspin=4), automatically produces 
        spin-polarized DOS showing up and down electron contributions separately,
        allowing analysis of magnetic ordering and site-projected properties.
        
        Note: The package changes from 'pw' (SCF/NSCF) to 'dos' for this calculation.
        
        Args:
            nscf_label: Directory of the NSCF calculation (default 'nscf')
            dos_label: Save results in this directory (default 'dos')
            Emin: Minimum energy for DOS (eV relative to Fermi). If None, uses -30 eV
            Emax: Maximum energy for DOS (eV relative to Fermi). If None, uses +10 eV
            DeltaE: Energy grid spacing (eV, default 0.01)
            degauss: Gaussian broadening (eV). If None, uses preset value
            ngauss: Gaussian broadening type (0=cold, 1=Fermi-Dirac). Default 0
            pdos: If True, compute local/projected DOS by atomic site (not orbital)
            
        Returns:
            EspressoDos: Post-processing calculator with DOS results
            
        Notes:
            - For magnetic systems: DOS includes spin-polarized contributions
            - PDOS useful for analyzing magnetic ordering in transition metals
            - Total DOS = DOS(up) + DOS(down) for magnetic systems
            
        Example:
            >>> # SCF + NSCF for magnetic system
            >>> scf_calc = workflow.run_scf(label='scf')  # nspin=2
            >>> workflow.atoms = scf_calc.atoms
            >>> nscf_calc = workflow.run_nscf(label='nscf', kpts=(12, 12, 12))
            
            >>> # DOS post-processing (separates spin-up and spin-down)
            >>> dos_result = workflow.run_dos(
            ...     nscf_label='nscf',
            ...     Emin=-30,
            ...     Emax=10,
            ...     pdos=False  # Set True for projected DOS
            ... )
            
            >>> # Plot spin-polarized DOS
            >>> from xespresso.dos import DOS
            >>> dos = DOS(label='nscf', prefix='nscf')
            >>> dos.read_dos()
            >>> dos.plot_dos(Emin=-30, Emax=10, smearing=[0.01])
        """
        from xespresso.post.dos import EspressoDos
        import re
        
        logger.info(f"Starting DOS post-processing from {nscf_label}...")
        
        # Check if NSCF calculation directory exists
        nscf_path = Path(nscf_label)
        if not nscf_path.exists():
            raise FileNotFoundError(
                f"NSCF calculation directory '{nscf_label}' not found. "
                f"Please run NSCF first: workflow.run_nscf(label='{nscf_label}')"
            )
        
        # Extract prefix from NSCF directory path
        nscf_prefix = nscf_label.split('/')[-1] if '/' in nscf_label else nscf_label
        
        # Detect if system is magnetic by checking input_data
        nspin = self.input_data.get('nspin', 1)
        is_magnetic = nspin > 1
        
        # Set default energy windows if not specified
        if Emin is None:
            Emin = -30.0  # 30 eV below Fermi level
        if Emax is None:
            Emax = 10.0   # 10 eV above Fermi level
        
        # Use degauss from preset if not specified
        if degauss is None:
            degauss = self.input_data.get('degauss', 0.01)
        
        # Prepare DOS parameters
        dos_params = {
            'Emin': Emin,
            'Emax': Emax,
            'DeltaE': DeltaE,
            'degauss': degauss,
            'ngauss': ngauss,
        }
        
        # Log system analysis
        print("\n" + "="*70)
        print("DOS CALCULATION - MAGNETIC SYSTEM ANALYSIS")
        print("="*70)
        print(f"NSCF source: {nscf_label}")
        print(f"Magnetic system (nspin={nspin}): {is_magnetic}")
        if is_magnetic:
            if nspin == 2:
                print(f"  → Collinear magnetism: spin-up and spin-down electrons")
            elif nspin == 4:
                print(f"  → Non-collinear magnetism: full spinor calculation")
        print("\nDOS Parameters:")
        print(f"  Energy window: [{Emin}, {Emax}] eV (relative to Fermi)")
        print(f"  Grid spacing (DeltaE): {DeltaE} eV")
        print(f"  Broadening (degauss): {degauss} eV")
        print(f"  Broadening type (ngauss): {ngauss} (0=Methfessel-Paxton)")
        print(f"  PDOS (projected): {pdos}")
        
        if is_magnetic:
            print("\n" + "-"*70)
            print("Spin-Polarized Analysis:")
            print("  DOS file will contain separate contributions for:")
            print("    • Spin-UP electrons")
            print("    • Spin-DOWN electrons")
            print("    • Total DOS = DOS(↑) + DOS(↓)")
            if pdos:
                print("  Site-projected contributions available for each atom")
        print("="*70 + "\n")
        
        logger.info(
            f"DOS parameters: Emin={Emin} eV, Emax={Emax} eV, DeltaE={DeltaE} eV, "
            f"degauss={degauss} eV, ngauss={ngauss}, "
            f"magnetic={is_magnetic}, pdos={pdos}"
        )
        
        # Create DOS post-processor (note: package changes from 'pw' to 'dos')
        dos_calc = EspressoDos(
            parent_directory=nscf_label,
            prefix=nscf_prefix,
            queue=self.queue,
            parallel=self.queue.get('parallel', '') if self.queue else '',
            **dos_params
        )
        
        # Execute DOS calculation
        dos_calc.run()
        
        logger.info(f"DOS calculation completed. Results in {dos_calc.directory}/")
        
        # Provide next steps information
        if is_magnetic:
            logger.info(
                f"Spin-polarized DOS calculated. Use xespresso.dos.DOS class to analyze:\n"
                f"  - Compare DOS(up) vs DOS(down) to validate magnetic ordering\n"
                f"  - Site projections show local magnetization\n"
                f"  - Orbital decomposition available with PDOS"
            )
        else:
            logger.info(f"DOS output can be plotted using xespresso.dos.DOS class")
        
        return dos_calc
    
    def run_bands(
        self,
        label: str = 'bands',
        bandpath_type: str = 'auto',
        mode: str = 'explicit',
        **calc_kwargs
    ):
        """
        Run band structure calculation along high-symmetry k-path.
        
        This calculates the electronic band structure using the charge density 
        converged from SCF/NSCF. Automatically generates or uses a high-symmetry 
        k-point path based on crystal symmetry.
        
        Args:
            label: Directory/label for the calculation (default 'bands')
            bandpath_type: How to generate k-point path:
                - 'auto': Automatic from cell.bandpath() (default, recommended)
                - 'custom': User provides custom k-path (not yet implemented)
            mode: Calculation mode:
                - 'explicit': Calculate at explicit k-points along path (default)
                - 'interpolated': Interpolate from NSCF (not yet implemented)
            **calc_kwargs: Additional parameters for the calculator
            
        Returns:
            Espresso: Calculator with band structure results
            
        Notes:
            - Reuses charge density from SCF/NSCF (same label without /bands)
            - For magnetic systems: respects nspin, produces magnetized bands
            - Automatic detection of high-symmetry points (Γ, X, W, L, K, U, etc)
            - Standard path: Al example gives GXWKGLUWLK,UX (50 k-points)
            
        Example:
            >>> # SCF calculation
            >>> scf = workflow.run_scf(label='scf')
            >>> 
            >>> # Band structure along auto-generated path
            >>> bands = workflow.run_bands(label='bands')
            >>> 
            >>> # For magnetic system (automatically respects nspin=2)
            >>> bands_mag = workflow.run_bands(label='bands_mag')
            >>> # Shows band splitting from magnetic moment
            
            >>> # Extract and plot
            >>> try:
            ...     bs = bands.band_structure()
            ...     bs.reference = bands.get_fermi_level()
            ...     bs.plot()
            ... except:
            ...     print("Use xespresso plotting utilities")
        """
        logger.info(f"Starting band structure calculation...")
        
        # Generate band path from crystal symmetry
        if bandpath_type == 'auto':
            bandpath = self.atoms.cell.bandpath()
            logger.info(f"Auto-generated band path from crystal symmetry:")
            logger.info(f"  Path: {bandpath.path}")
            logger.info(f"  High-symmetry points: {list(bandpath.special_points.keys())}")
            logger.info(f"  Total k-points: {len(bandpath.kpts)}")
            kpts = bandpath
        else:
            raise NotImplementedError(
                f"bandpath_type='{bandpath_type}' not yet implemented. "
                f"Use 'auto' for automatic generation from cell symmetry."
            )
        
        # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Prepare input_data (copy from preset)
        input_data = self.input_data.copy()
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': 'bands',  # High-symmetry k-path calculation
            'input_data': input_data,
            'kpts': kpts,
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = input_data.get('ecutrho', 400.0)
        
        # Set pseudo_dir when using pseudopotentials_config
        if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
            params['input_data']['pseudo_dir'] = './pseudo'
        
        # Add queue configuration if provided
        if self.queue is not None:
            params['queue'] = self.queue
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create calculator for band structure
        calc = Espresso(**params)
        self.atoms.calc = calc
        self.last_calc = calc  # Track last calculator for monitoring
        
        # Run calculation (local or remote)
        if self.queue and self.queue.get('execution') == 'remote' and not self.queue.get('wait_for_completion', False):
            logger.info("Remote non-blocking: executing band structure with automatic job monitoring...")
            
            # Step 1: Write input
            calc.write_input(self.atoms)
            calc.atoms = self.atoms
            
            # Step 2: Execute (submits job remotely)
            calc.execute()
            
            # Store remote connection on calc for RemoteJobMonitor to access
            if hasattr(calc, 'scheduler') and hasattr(calc.scheduler, 'remote'):
                calc.remote = calc.scheduler.remote
            
            # Step 3: Monitor SLURM job status
            logger.info(f"Remote job {calc.last_job_id} submitted. Monitoring SLURM status...")
            job_id = calc.last_job_id
            timeout = self.queue.get('job_timeout', 3600)  # 1 hour for bands
            job_monitor_result = self._monitor_remote_job(calc, job_id, timeout=timeout, poll_interval=30)
            
            if not job_monitor_result['success']:
                raise RuntimeError(f"Remote job {job_id} failed: {job_monitor_result['message']}")
            
            # Step 4: Fetch output
            monitor = RemoteJobMonitor(calc)
            if monitor.wait(timeout=60, poll_interval=5):
                monitor.retrieve_output()
                logger.info("Remote band structure output retrieved.")
                calc.read_results()
            else:
                raise RuntimeError(f"Failed to retrieve band structure output for job {job_id}")
        else:
            # Local or remote blocking: use normal run()
            calc.run(atoms=self.atoms)
        
        # Check convergence
        self._check_convergence(calc, calculation_type='bands')
        
        # Report magnetic band structure info
        nspin = self.input_data.get('nspin', 1)
        is_magnetic = nspin > 1
        
        print("\n" + "="*70)
        print("BAND STRUCTURE CALCULATION COMPLETED")
        print("="*70)
        print(f"Calculation: {label}/")
        print(f"High-symmetry path: {bandpath.path}")
        print(f"K-points calculated: {len(bandpath.kpts)}")
        if is_magnetic:
            if nspin == 2:
                print(f"Spin-polarized bands: YES (nspin=2, collinear)")
                print(f"  → Separate band structures for spin-up and spin-down")
            elif nspin == 4:
                print(f"Non-collinear bands: YES (nspin=4)")
                print(f"  → Full spinor band structure")
        else:
            print(f"Magnetic bands: No (nspin=1)")
        print("="*70)
        
        logger.info(f"Band structure calculation completed in {label}/")
        logger.info(f"Use Espresso.band_structure() to extract BandPath object")
        
        return calc
    
    def run_projwfc(
        self,
        nscf_label: str = 'nscf',
        projwfc_label: str = 'projwfc',
        Emin: Optional[float] = None,
        Emax: Optional[float] = None,
        DeltaE: float = 0.01,
        degauss: Optional[float] = None,
        ngauss: int = 0,
        lsym: int = 1,
        pawproj: int = 0,
        filpdos: Optional[str] = None,
        lowdin: bool = False,
    ):
        """
        Run Projections on Atomic Wavefunctions (PROJWFC) post-processing.
        
        This computes the local and orbital-projected density of states (PDOS)
        from a prior NSCF calculation. Useful for understanding which atoms
        and orbitals contribute to the electronic structure, especially
        important for magnetic systems and strongly correlated materials.
        
        Note: The package changes from 'pw' (SCF/NSCF) to 'projwfc' for this calculation.
        
        Args:
            nscf_label: Directory of the NSCF calculation (default 'nscf')
            projwfc_label: Output directory label (default 'projwfc')
            Emin: Minimum energy for PDOS (eV relative to Fermi). If None, uses -30 eV
            Emax: Maximum energy for PDOS (eV relative to Fermi). If None, uses +10 eV
            DeltaE: Energy grid spacing (eV, default 0.01)
            degauss: Gaussian broadening (eV). If None, uses preset value
            ngauss: Gaussian broadening type (0=Methfessel-Paxton, 1=Fermi-Dirac)
            lsym: Symmetrize projections (0=no, 1=yes). Default 1
            pawproj: PAW projector type (0=Rydberg, 1=m_j dependent). Default 0
            filpdos: Prefix for output PDOS files. If None, uses nscf prefix
            lowdin: If True, calculate and output Lowdin charges to file (default False)
            
        Returns:
            EspressoProjwfc: Post-processing calculator with PDOS results
            
        Notes:
            - For magnetic systems: Projections include spin-polarized contributions
            - PDOS shows which atoms/orbitals contribute at each energy
            - Essential for validating magnetic orderings and orbital occupations
            - Output files: {prefix}.pdos_* (one per orbital symmetry type)
            - If lowdin=True: generates lowdin_charges.dat with atomic population analysis
            
        Example:
            >>> # SCF + NSCF for transition metal oxide
            >>> scf_calc = workflow.run_scf(label='scf')
            >>> workflow.atoms = scf_calc.atoms
            >>> nscf_calc = workflow.run_nscf(label='nscf', kpts=(12, 12, 12))
            
            >>> # PROJWFC for orbital-resolved analysis
            >>> projwfc_result = workflow.run_projwfc(
            ...     nscf_label='nscf',
            ...     Emin=-30,
            ...     Emax=10,
            ...     DeltaE=0.01
            ... )
            
            >>> # PROJWFC with Lowdin charges
            >>> projwfc_result = workflow.run_projwfc(
            ...     nscf_label='nscf',
            ...     Emin=-30,
            ...     Emax=10,
            ...     lowdin=True  # Generate lowdin_charges.dat
            ... )
            
            >>> # Analyze PDOS
            >>> from xespresso.dos import DOS
            >>> dos = DOS(label='nscf', prefix='nscf')
            >>> dos.read_pdos()
            >>> dos.plot_pdos(Emin=-30, Emax=10, smearing=[0.01])
        """
        from xespresso.post.projwfc import EspressoProjwfc
        
        logger.info(f"Starting PROJWFC post-processing from {nscf_label}...")
        
        # Check if NSCF calculation directory exists
        nscf_path = Path(nscf_label)
        if not nscf_path.exists():
            raise FileNotFoundError(
                f"NSCF calculation directory '{nscf_label}' not found. "
                f"Please run NSCF first: workflow.run_nscf(label='{nscf_label}')"
            )
        
        # Extract prefix from NSCF directory path
        nscf_prefix = nscf_label.split('/')[-1] if '/' in nscf_label else nscf_label
        
        # Detect if system is magnetic by checking input_data
        nspin = self.input_data.get('nspin', 1)
        is_magnetic = nspin > 1
        
        # Set default energy windows if not specified
        if Emin is None:
            Emin = -30.0  # 30 eV below Fermi level
        if Emax is None:
            Emax = 10.0   # 10 eV above Fermi level
        
        # Use degauss from preset if not specified
        if degauss is None:
            degauss = self.input_data.get('degauss', 0.01)
        
        # Set output filename if not specified
        if filpdos is None:
            filpdos = nscf_prefix
        
        # Set Lowdin charges output filename if requested
        filowdin = None
        if lowdin:
            filowdin = f"{nscf_prefix}_lowdin_charges.dat"
        
        # Prepare PROJWFC parameters
        projwfc_params = {
            'Emin': Emin,
            'Emax': Emax,
            'DeltaE': DeltaE,
            'degauss': degauss,
            'ngauss': ngauss,
            'lsym': lsym,
            'pawproj': pawproj,
            'filpdos': filpdos,
            'filowdin': filowdin,
        }
        
        # Log system analysis
        print("\n" + "="*70)
        print("PROJWFC (PDOS) CALCULATION - ORBITAL ANALYSIS")
        print("="*70)
        print(f"NSCF source: {nscf_label}")
        print(f"Magnetic system (nspin={nspin}): {is_magnetic}")
        if is_magnetic:
            if nspin == 2:
                print(f"  → Collinear magnetism: separate spin-up/down projections")
            elif nspin == 4:
                print(f"  → Non-collinear magnetism: spinor projections")
        print("\nPROJWFC Parameters:")
        print(f"  Energy window: [{Emin}, {Emax}] eV (relative to Fermi)")
        print(f"  Grid spacing (DeltaE): {DeltaE} eV")
        print(f"  Broadening (degauss): {degauss} eV")
        print(f"  Broadening type (ngauss): {ngauss}")
        print(f"  Symmetry projection (lsym): {lsym}")
        print(f"  PAW projector (pawproj): {pawproj}")
        print(f"  Output prefix (filpdos): {filpdos}")
        print(f"  Lowdin charges: {lowdin}")
        if lowdin:
            print(f"    → Output file: {filowdin}")
        
        if is_magnetic:
            print("\n" + "-"*70)
            print("Orbital-Projected Analysis:")
            print("  PDOS will show contributions from:")
            print("    • Individual atoms (site-projected)")
            print("    • Different orbitals (s, p, d, f, etc.)")
            print("    • Spin-up and spin-down for magnetic systems")
            if lowdin:
                print("    • Lowdin charges: atomic population analysis")
            print("  Useful for validating magnetic orderings in transition metals")
        else:
            print("\n" + "-"*70)
            print("Orbital-Projected Analysis:")
            print("  PDOS will decompose electronic structure by:")
            print("    • Atomic sites")
            print("    • Orbital angular momentum (s, p, d, f)")
            if lowdin:
                print("    • Lowdin charges: atomic population analysis")
        print("="*70 + "\n")
        
        logger.info(
            f"PROJWFC parameters: Emin={Emin} eV, Emax={Emax} eV, DeltaE={DeltaE} eV, "
            f"degauss={degauss} eV, ngauss={ngauss}, lsym={lsym}, "
            f"pawproj={pawproj}, filpdos={filpdos}, magnetic={is_magnetic}"
        )
        
        # Create PROJWFC post-processor (note: package changes from 'pw' to 'projwfc')
        projwfc_calc = EspressoProjwfc(
            parent_directory=nscf_label,
            prefix=nscf_prefix,
            queue=self.queue,
            parallel=self.queue.get('parallel', '') if self.queue else '',
            **projwfc_params
        )
        
        # Execute PROJWFC calculation
        projwfc_calc.run()
        
        logger.info(f"PROJWFC calculation completed. Results in {projwfc_calc.directory}/")
        
        # Provide next steps information
        analysis_msg = f"Projected DOS (PDOS) calculated. Use xespresso.dos.DOS class to analyze:\n"
        analysis_msg += f"  - read_pdos() to load projection data\n"
        analysis_msg += f"  - plot_pdos() to visualize orbital contributions\n"
        analysis_msg += f"  - Compare site/orbital contributions to validate electronic structure\n"
        analysis_msg += f"  - For magnetic systems: analyze spin-up vs spin-down projections"
        
        if lowdin:
            analysis_msg += f"\nLowdin charges calculated in {filowdin}:\n"
            analysis_msg += f"  - Contains atomic population analysis using Lowdin transformation\n"
            analysis_msg += f"  - Shows charge distribution by atom and orbital\n"
            analysis_msg += f"  - Includes spilling parameter and magnetic moments"
        
        logger.info(analysis_msg)
        
        return projwfc_calc
    
    def _estimate_nbnd(self) -> int:
        """
        Estimate number of bands needed for band structure calculations.
        
        Uses pseudopotential valence electrons + buffer for unoccupied states.
        
        Returns:
            int: Recommended nbnd value
        """
        from xespresso.workflow.wannier_workflow import suggest_nbnd_from_pseudos
        
        # Use the existing helper from wannier_workflow
        nbnd = suggest_nbnd_from_pseudos(self.pseudopotentials, buffer=20)
        logger.info(f"Estimated nbnd: {nbnd}")
        
        return nbnd
    
    def run_relax(
        self,
        label: str = 'relax',
        relax_type: str = 'relax',
        **calc_kwargs
    ) -> Espresso:
        """
        Run a structure relaxation calculation.
        
        Args:
            label: Directory/label for the calculation
            relax_type: Type of relaxation: 'relax' (ions only) or 'vc-relax' (ions + cell)
            **calc_kwargs: Additional parameters for the Espresso calculator
            
        Returns:
            Espresso: Calculator object with results
        """
        if relax_type not in ['relax', 'vc-relax']:
            raise ValueError(
                f"relax_type must be 'relax' or 'vc-relax', got '{relax_type}'"
            )
        
        # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
        # This allows remote_mixin to find pseudopotentials
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': relax_type,
            'input_data': self.input_data.copy(),
            'kpts': self._get_kpts(),
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = self.input_data.get('ecutrho', 400.0)
        
        # Set pseudo_dir when using pseudopotentials_config
        if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
            params['input_data']['pseudo_dir'] = './pseudo'
        
        # Add queue configuration if provided
        if self.queue is not None:
            params['queue'] = self.queue
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create calculator
        calc = Espresso(**params)
        self.atoms.calc = calc
        self.last_calc = calc  # Track last calculator for monitoring
        
        # If remote non-blocking: control execution steps to avoid retry loop
        if self.queue and self.queue.get('execution') == 'remote' and not self.queue.get('wait_for_completion', False):
            logger.info("Remote non-blocking: executing with automatic job monitoring...")
            
            # Step 1: Write input (with atoms, so _transfer_pseudopotentials won't need to call it again)
            calc.write_input(self.atoms)
            
            # IMPORTANT: Set calc.atoms so _transfer_pseudopotentials() can use it if needed
            calc.atoms = self.atoms
            
            # Step 2: Execute (submits job remotely)
            calc.execute()
            
            # IMPORTANT: Store remote connection on calc for RemoteJobMonitor to access
            if hasattr(calc, 'scheduler') and hasattr(calc.scheduler, 'remote'):
                calc.remote = calc.scheduler.remote
            
            # Step 3: Monitor SLURM job status (detect stuck jobs)
            job_id = getattr(calc, 'last_job_id', None)
            
            if job_id is None:
                raise RuntimeError(
                    "Remote job submission failed: No job ID returned from scheduler. "
                    "Check scheduler configuration and job submission logs."
                )
            
            logger.info(f"Remote job {job_id} submitted. Monitoring SLURM status...")
            timeout = self.queue.get('job_timeout', 3600)
            job_monitor_result = self._monitor_remote_job(calc, job_id, timeout=timeout, poll_interval=30)
            
            if not job_monitor_result['success']:
                raise RuntimeError(f"Remote job {job_id} failed: {job_monitor_result['message']}")
            
            # Step 4: Job completed in queue, now fetch output using RemoteJobMonitor
            monitor = RemoteJobMonitor(calc)
            if monitor.wait(timeout=60, poll_interval=5):  # Short timeout since job already completed
                monitor.retrieve_output()
                logger.info("Remote job output retrieved.")
                # Step 5: Read results
                calc.read_results()
            else:
                raise RuntimeError(f"Failed to retrieve output for job {job_id}")
        else:
            # Local or remote blocking: use normal run() with retry logic
            calc.run(atoms=self.atoms)
        
        # Check convergence and inform user
        self._check_convergence(calc, calculation_type=relax_type)
        
        return calc
    
    def write_input(
        self,
        atoms: Optional[Atoms] = None,
        label: str = 'scf',
        calculation: str = 'scf',
        **calc_kwargs
    ) -> Espresso:
        """
        Generate input files without running the calculation (dry run).
        
        Args:
            atoms: Optional ASE Atoms object. If None, uses workflow atoms.
            label: Directory/label for the calculation
            calculation: Type of calculation: 'scf', 'relax', 'vc-relax'
            **calc_kwargs: Additional parameters for the Espresso calculator
            
        Returns:
            Espresso: Calculator object with input files written
            
        Example:
            >>> workflow.write_input(label='scf/si-test')
            >>> # Files generated to scf/si-test/.pwi and scf/si-test/job_file
        """
        if atoms is None:
            atoms = self.atoms
        
        # Set ESPRESSO_PSEUDO if we have pseudopotentials_config
        # This allows remote_mixin to find pseudopotentials
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
        
        # Prepare parameters
        params = {
            'pseudopotentials': self.pseudopotentials,
            'label': label,
            'calculation': calculation,
            'input_data': self.input_data.copy(),
            'kpts': self._get_kpts(),
        }
        
        # Add ecutwfc and ecutrho at top level
        params['ecutwfc'] = self.input_data.get('ecutwfc', 50.0)
        params['ecutrho'] = self.input_data.get('ecutrho', 400.0)
        
        # Set pseudo_dir when using pseudopotentials_config
        if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
            params['input_data']['pseudo_dir'] = './pseudo'
        
        # Note: write_input() is a dry run, never use queue (always local)
        
        # Merge with extra kwargs
        params.update(self.extra_kwargs)
        params.update(calc_kwargs)
        
        # Create calculator and write input files (no execution)
        calc = Espresso(**params)
        calc.write_input(atoms)
        self.last_calc = calc  # Track last calculator for monitoring
        
        return calc
    
    def get_monitor(self) -> RemoteJobMonitor:
        """
        Get a monitor for the last remote job execution.
        
        Use this to track non-blocking remote jobs started with run_scf() or run_relax().
        
        Returns:
            RemoteJobMonitor: Monitor object for the last executed calculator
            
        Raises:
            ValueError: If no remote job was executed or job lacks required tracking info
            
        Examples:
            >>> calc = workflow.run_scf(label='scf/si-test')
            >>> monitor = workflow.get_monitor()
            >>> print(monitor.status())  # 'running', 'completed', 'failed'
            >>> if monitor.wait(timeout=3600):
            ...     local_path, output = monitor.retrieve_output()
        """
        if not hasattr(self, 'last_calc') or self.last_calc is None:
            raise ValueError("No calculation has been executed yet")
        
        return RemoteJobMonitor(self.last_calc)
    
    def get_atoms(self) -> Atoms:
        """Get the current atoms object."""
        return self.atoms
        """Get information about the current protocol preset."""
        return {
            'protocol': self.protocol,
            'preset': self.preset,
            'kpts': self._get_kpts(),
            'kspacing': self.kspacing,
        }


def quick_scf(
    structure: Union[str, Path, Atoms],
    pseudopotentials: Dict[str, str],
    label: str = 'scf',
    protocol: str = 'moderate',
    kspacing: Optional[float] = None,
    magnetic_config: Optional[Union[str, Dict]] = None,
    expand_cell: bool = False,
    queue: Optional[Dict] = None,
    machine: Optional[str] = None,
    **kwargs
) -> Espresso:
    """
    Quick SCF calculation helper function.
    
    Args:
        structure: CIF file path or ASE Atoms object
        pseudopotentials: Dictionary mapping element symbols to pseudopotential files
        label: Directory/label for the calculation
        protocol: Protocol preset: 'fast', 'moderate', or 'accurate'
        kspacing: K-point spacing in Angstrom^-1 (physical units)
        magnetic_config: Magnetic configuration ('ferro', 'antiferro', or element dict)
        expand_cell: If True, expand cell to accommodate magnetic configuration
        queue: Queue configuration dictionary for job submission (local or remote)
        machine: Name of a machine configuration to load from ~/.xespresso/machines/
        **kwargs: Additional parameters for the calculator
        
    Returns:
        Espresso: Calculator object with results
        
    Example:
        >>> calc = quick_scf(
        ...     'structure.cif',
        ...     {'Si': 'Si.pbe.UPF'},
        ...     protocol='fast'
        ... )
        >>> # With magnetic configuration
        >>> calc = quick_scf(
        ...     atoms,
        ...     {'Fe': 'Fe.pbe-spn.UPF'},
        ...     magnetic_config='antiferro',
        ...     protocol='moderate'
        ... )
        >>> # With remote execution
        >>> calc = quick_scf(
        ...     'structure.cif',
        ...     {'Fe': 'Fe.pbe-spn.UPF'},
        ...     protocol='moderate',
        ...     machine='cluster1'  # Load from ~/.xespresso/machines/cluster1.json
        ... )
    """
    if isinstance(structure, (str, Path)):
        workflow = CalculationWorkflow.from_cif(
            structure, pseudopotentials, protocol, kspacing, 
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    else:
        workflow = CalculationWorkflow(
            structure, pseudopotentials, protocol, kspacing,
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    
    return workflow.run_scf(label=label)


def quick_relax(
    structure: Union[str, Path, Atoms],
    pseudopotentials: Dict[str, str],
    label: str = 'relax',
    protocol: str = 'moderate',
    kspacing: Optional[float] = None,
    relax_type: str = 'relax',
    magnetic_config: Optional[Union[str, Dict]] = None,
    expand_cell: bool = False,
    queue: Optional[Dict] = None,
    machine: Optional[str] = None,
    **kwargs
) -> Espresso:
    """
    Quick structure relaxation helper function.
    
    Args:
        structure: CIF file path or ASE Atoms object
        pseudopotentials: Dictionary mapping element symbols to pseudopotential files
        label: Directory/label for the calculation
        protocol: Protocol preset: 'fast', 'moderate', or 'accurate'
        kspacing: K-point spacing in Angstrom^-1 (physical units)
        relax_type: Type of relaxation: 'relax' or 'vc-relax'
        magnetic_config: Magnetic configuration ('ferro', 'antiferro', or element dict)
        expand_cell: If True, expand cell to accommodate magnetic configuration
        queue: Queue configuration dictionary for job submission (local or remote)
        machine: Name of a machine configuration to load from ~/.xespresso/machines/
        **kwargs: Additional parameters for the calculator
        
    Returns:
        Espresso: Calculator object with results
        
    Example:
        >>> calc = quick_relax(
        ...     'structure.cif',
        ...     {'Si': 'Si.pbe.UPF'},
        ...     protocol='moderate',
        ...     relax_type='vc-relax'
        ... )
        >>> # With Hubbard parameters
        >>> calc = quick_relax(
        ...     atoms,
        ...     {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
        ...     magnetic_config={'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}},
        ...     protocol='accurate'
        ... )
        >>> # With remote execution on SLURM cluster
        >>> calc = quick_relax(
        ...     'structure.cif',
        ...     {'Fe': 'Fe.pbe-spn.UPF'},
        ...     protocol='moderate',
        ...     machine='slurm_cluster'  # Load from config
        ... )
    """
    if isinstance(structure, (str, Path)):
        workflow = CalculationWorkflow.from_cif(
            structure, pseudopotentials, protocol, kspacing,
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    else:
        workflow = CalculationWorkflow(
            structure, pseudopotentials, protocol, kspacing,
            magnetic_config=magnetic_config, expand_cell=expand_cell,
            queue=queue, machine=machine, **kwargs
        )
    
    return workflow.run_relax(label=label, relax_type=relax_type)

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
    },
    'moderate': {
        'ecutwfc': 50.0,
        'ecutrho': 400.0,
        'conv_thr': 1.0e-8,
        'kspacing': 0.3,  # Angstrom^-1
        'mixing_beta': 0.5,
        'electron_maxstep': 200,
    },
    'accurate': {
        'ecutwfc': 80.0,
        'ecutrho': 640.0,
        'conv_thr': 1.0e-10,
        'kspacing': 0.15,  # Angstrom^-1
        'mixing_beta': 0.3,
        'electron_maxstep': 300,
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
                pseudopotentials[element] = pseudo_obj.filename
                logger.info(f"DEBUG: Extracted {element} -> {pseudo_obj.filename}")
            else:
                missing_elements.append(element)
        
        if missing_elements:
            available = config.list_elements()
            raise ValueError(
                f"Pseudopotentials configuration '{config_name}' is missing "
                f"the following elements: {missing_elements}. "
                f"Available elements: {available}"
            )
        
        return pseudopotentials
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
            elif magnetic_config in ['antiferro', 'antiferromagnetic']:
                # Simple antiferromagnetic configuration
                # For antiferromagnetic, we need to determine sublattices
                # Simple approach: alternate atoms
                n_atoms = len(self.atoms)
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
            
            # Step 3: Wait for remote job completion
            logger.info(f"Remote job {calc.last_job_id} submitted. Waiting for completion...")
            monitor = RemoteJobMonitor(calc)
            timeout = self.queue.get('job_timeout', 3600)
            if monitor.wait(timeout=timeout, poll_interval=10):
                monitor.retrieve_output()
                logger.info("Remote job completed and output retrieved.")
                # Step 4: Read results
                calc.read_results()
            else:
                raise RuntimeError(f"Remote job {calc.last_job_id} timed out after {timeout}s")
        else:
            # Local or remote blocking: use normal run() with retry logic
            calc.run(atoms=self.atoms)
        
        return calc
    
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
            
            # Step 3: Wait for remote job completion
            logger.info(f"Remote job {calc.last_job_id} submitted. Waiting for completion...")
            monitor = RemoteJobMonitor(calc)
            timeout = self.queue.get('job_timeout', 3600)
            if monitor.wait(timeout=timeout, poll_interval=10):
                monitor.retrieve_output()
                logger.info("Remote job completed and output retrieved.")
                # Step 4: Read results
                calc.read_results()
            else:
                raise RuntimeError(f"Remote job {calc.last_job_id} timed out after {timeout}s")
        else:
            # Local or remote blocking: use normal run() with retry logic
            calc.run(atoms=self.atoms)
        
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

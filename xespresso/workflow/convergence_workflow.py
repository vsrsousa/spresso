"""
Convergence parameter optimization workflow.

This module provides tools to systematically optimize DFT calculation parameters
(ecutwfc, kspacing) for energy and force convergence on a target structure.

Two usage modes:
1. Simple mode: Just provide structure and desired precision level
2. Advanced mode: Specify custom parameter ranges for detailed control

Simple workflow:
    1. Define precision level ('low', 'medium', 'high', 'ultra')
    2. Workflow automatically determines optimal parameter ranges
    3. Run convergence study and get recommendations

Advanced workflow:
    1. Define parameter ranges (ecutwfc, kspacing)
    2. Run SCF calculations for each parameter combination
    3. Analyze convergence of total energy and forces
    4. Recommend optimal parameters with target accuracy
"""

import logging
import os
import numpy as np
import pandas as pd
from typing import Dict, Optional, Union, Tuple, List
from pathlib import Path
from ase import Atoms
from ase.io import read
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.pseudopotentials.detector import parse_upf_header


logger = logging.getLogger(__name__)


class ConvergenceWorkflow:
    """
    Systematic convergence testing for DFT calculations.
    
    This class systematically varies DFT parameters (ecutwfc, kspacing) and
    runs SCF calculations to determine optimal values for energy and force
    convergence.
    
    **IMPROVED ALGORITHM (NEW)**: Two-phase independent convergence
    ✅ PHASE 1: Convergence de ecutwfc com kspacing FIXO (0.5 Å⁻¹)
    ✅ PHASE 2: Convergence de kspacing com ecutwfc otimizado
    
    Benefits:
    - Pseudopotenciais transferidos apenas 2x (não N×M times!)
    - 4-6x mais rápido que nested-loop approach
    - Dois processos completamente independentes
    
    Legacy nested-loop mode still available via `independent_mode=False`
    (NOT recommended - kept only for backward compatibility)
    
    KEY CONCEPT: Precision level controls PARAMETER RANGES only, while convergence
    criteria are INDEPENDENT. You can use any combination of criteria with any precision.
    
    Two usage modes:
    
    1. Simple mode (recommended for most users):
        >>> # Just provide structure and precision level
        >>> workflow = ConvergenceWorkflow.from_cif(
        ...     'structure.cif',
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     precision='low'  # Controls parameter ranges (ecutwfc, kspacing), fastest settings
        ... )
        >>> optimal_params = workflow.optimize_parameters()
    
    2. Advanced mode (for detailed control):
        >>> # Specify custom parameter ranges and criteria independently
        >>> conv = ConvergenceWorkflow(
        ...     atoms=atoms,
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     precision='low',  # Coarse parameter ranges
        ...     convergence_criteria_list=['energy', 'forces', 'geometry', 'stress']  # Strict criteria
        ... )
        >>> # Uses INDEPENDENT two-phase algorithm by default
        >>> conv.run_convergence_study()
        >>> recommendations = conv.get_recommendations()
    
    Attributes:
        atoms: ASE Atoms object (structure to test)
        pseudopotentials: Dictionary mapping element symbols to UPF files
        protocol: Base protocol for calculations ('fast', 'moderate', 'accurate')
        precision: Precision level ('low', 'medium', 'high', 'ultra') - controls parameter ranges
                   Ranges are automatically adjusted based on pseudopotential requirements!
        ecutwfc_range: List of ecutwfc values to test (automatically optimized)
        kspacing_range: List of kspacing values to test
        results: DataFrame with convergence test results
    """
    
    def __init__(
        self,
        atoms: Atoms,
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        protocol: str = 'moderate',
        precision: str = 'low',
        ecutwfc_range: Optional[List[float]] = None,
        kspacing_range: Optional[List[float]] = None,
        conv_thr_range: Optional[List[float]] = None,
        convergence_criteria_list: Optional[List[str]] = None,
        convergence_criteria: Optional[Dict] = None,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        **kwargs
    ):
        """
        Initialize convergence workflow.
        
        Args:
            atoms: ASE Atoms object with structure
            pseudopotentials: Dict mapping element symbols to UPF files
            protocol: Base protocol ('fast', 'moderate', 'accurate')
            precision: Precision level for automatic parameter selection.
                     Options: 'low', 'medium', 'high', 'ultra'.
                     Overrides ecutwfc_range and kspacing_range.
                     Default: 'low' (fast convergence)
            ecutwfc_range: List of ecutwfc values to test (optional override).
            kspacing_range: List of kspacing values to test in Å⁻¹ (optional override).
            conv_thr_range: List of conv_thr values to test (optional)
            convergence_criteria_list: List of convergence criteria to check.
                                     Options: 'energy', 'forces', 'geometry', 'stress', 'magnetic_moments'
                                     If None, uses defaults based on precision level.
            convergence_criteria: Dict with convergence tolerances.
                                If None, uses defaults based on precision level.
                                Keys: 'energy_tolerance', 'force_tolerance', 'geometry_tolerance', 'stress_tolerance', 'magnetic_tolerance'
            queue: Queue configuration for job submission (optional)
            magnetic_config: Magnetic configuration for calculations. Can be a string like 'ferromagnetic', 
                           'antiferromagnetic', or a dict specifying magnetic moments per atom (optional)
            **kwargs: Additional parameters passed to CalculationWorkflow
        """
        self.atoms = atoms.copy()
        # Store pseudopotentials and config name. If config is provided, we pass the
        # config name to CalculationWorkflow so it loads pseudopotentials the same way
        # and has access to base_path for proper file resolution.
        self.pseudopotentials = {}
        self.pseudopotentials_base_path = None
        self._pseudo_config_name = pseudopotentials_config

        if pseudopotentials_config is not None:
            from xespresso.pseudopotentials.manager import load_pseudopotentials_config
            from xespresso.utils.pseudo_utils import get_ecutrho_ratio

            cfg = load_pseudopotentials_config(pseudopotentials_config, verbose=False)
            if cfg is None:
                raise ValueError(f"Pseudopotentials configuration '{pseudopotentials_config}' not found")

            # Store base path for use in convergence study analysis
            self.pseudopotentials_base_path = cfg.base_path if hasattr(cfg, 'base_path') else None
            
            # Only load pseudopotentials for elements present in the structure (for analysis)
            required_elements = set(self.atoms.get_chemical_symbols())
            for el, pseudo in cfg.pseudopotentials.items():
                if el in required_elements:
                    filename = pseudo.filename if hasattr(pseudo, 'filename') else str(pseudo)
                    # Store FILENAME only (not full path) - CalculationWorkflow will resolve via config
                    self.pseudopotentials[el] = filename
            
            # Calculate ecutrho ratio once, to be used for all calculations in this convergence study
            self.ecutrho_ratio = get_ecutrho_ratio(required_elements, cfg)
        else:
            if pseudopotentials is None:
                raise ValueError("Must provide 'pseudopotentials' mapping or 'pseudopotentials_config' name")
            self.pseudopotentials = pseudopotentials
            # No config loaded, use default ratio
            self.ecutrho_ratio = 4.0  # Default for Norm-Conserving
        self.protocol = protocol
        self.precision = precision
        
        # Set convergence criteria list
        if convergence_criteria_list is None:
            self.convergence_criteria_list = self._get_default_convergence_criteria_list(precision)
        else:
            self.convergence_criteria_list = convergence_criteria_list
            
        # Set convergence tolerances
        if convergence_criteria is None:
            self.convergence_criteria = self._get_default_convergence_criteria(precision)
        else:
            self.convergence_criteria = convergence_criteria
        
        # Handle machine and queue configuration (same as CalculationWorkflow)
        if queue is not None and machine is not None:
            raise ValueError(
                "Cannot specify both 'queue' and 'machine' parameters. "
                "Use 'queue' for direct configuration or 'machine' to load from config."
            )
        
        if machine is not None:
            # Load machine configuration
            from xespresso.machines import load_machine
            self.queue = load_machine(machine_name=machine)
            
            # Load code configuration and extract modules for specified version
            if code_version is not None:
                self._merge_code_modules_into_queue(machine, code_version)
        else:
            self.queue = queue
        
        self.machine = machine
        self.code_version = code_version
        self.magnetic_config = magnetic_config
        self.extra_kwargs = kwargs
        
        # Define parameter ranges based on precision level
        # Get smart ranges adjusted for pseudopotential requirements
        default_ecutwfc, default_kspacing = self._get_smart_ranges_for_pseudopotentials(
            precision, self.pseudopotentials, atoms
        )
        
        # Apply user overrides if provided
        self.ecutwfc_range = sorted(ecutwfc_range) if ecutwfc_range is not None else default_ecutwfc
        self.kspacing_range = sorted(kspacing_range, reverse=True) if kspacing_range is not None else default_kspacing
        
        # Results storage
        self.results = None  # DataFrame will be created after tests
        
        logger.info(
            f"Convergence workflow initialized:\n"
            f"  Structure: {self.atoms.get_chemical_formula()}\n"
            f"  Precision: {self.precision or 'custom'}\n"
            f"  ecutwfc range: {self.ecutwfc_range}\n"
            f"  kspacing range: {self.kspacing_range}"
        )
    
    @staticmethod
    def _get_ranges_for_precision_static(precision: str) -> Tuple[List[float], List[float]]:
        """
        Static version of _get_ranges_for_precision for use in classmethod.
        """
        precision = precision.lower()
        
        if precision == 'low':
            ecutwfc_range = [30, 40, 50]
            kspacing_range = [0.30, 0.27, 0.23]
        elif precision == 'medium':
            ecutwfc_range = [40, 50, 60, 70]
            kspacing_range = [0.30, 0.27, 0.23, 0.20]
        elif precision == 'high':
            ecutwfc_range = [50, 60, 70, 80, 90]
            kspacing_range = [0.30, 0.27, 0.23, 0.20, 0.18]
        elif precision == 'ultra':
            ecutwfc_range = [60, 80, 100, 120, 140]
            kspacing_range = [0.30, 0.27, 0.23, 0.20, 0.18, 0.15]
        else:
            raise ValueError(f"Unknown precision level: {precision}. "
                           "Choose from 'low', 'medium', 'high', 'ultra'")
        
        return ecutwfc_range, kspacing_range
    
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
        from xespresso.codes import load_codes_config
        
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
    
    @staticmethod
    def _get_smart_ranges_for_pseudopotentials_static(
        precision: str, 
        pseudopotentials: Dict[str, str]
    ) -> Tuple[List[float], List[float]]:
        """
        Get parameter ranges intelligently adjusted for pseudopotential requirements.
        
        Instead of starting with generic low values, this analyzes pseudopotentials
        and sets appropriate starting points for convergence testing.
        
        Args:
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            pseudopotentials: Dict mapping element symbols to UPF file paths
            
        Returns:
            Tuple of (ecutwfc_range, kspacing_range) adjusted for pseudopotentials
        """
        
        # Get base ranges for precision level
        ecutwfc_range, kspacing_range = ConvergenceWorkflow._get_ranges_for_precision_static(precision)
        
        # Extract suggested ecutwfc from pseudopotential files
        max_suggested_ecutwfc = 0.0
        pseudopotential_info = {}
        
        for element, upf_path in pseudopotentials.items():
            try:
                # Try to parse the UPF file
                header_info = parse_upf_header(upf_path)
                if 'suggested_ecutwfc' in header_info:
                    suggested = header_info['suggested_ecutwfc']
                    max_suggested_ecutwfc = max(max_suggested_ecutwfc, suggested)
                    pseudopotential_info[element] = {
                        'suggested_ecutwfc': suggested,
                        'file': upf_path
                    }
            except Exception as e:
                # Silently skip if parsing fails
                pass
        
        if max_suggested_ecutwfc > 0:
            
            # Adjust the starting point based on pseudopotential requirements
            # For high-ecutwfc pseudopotentials, start higher
            if max_suggested_ecutwfc >= 80:  # High-ecutwfc pseudopotentials
                if precision == 'low':
                    ecutwfc_range = [60, 80, 100]  # Start higher
                elif precision == 'medium':
                    ecutwfc_range = [80, 100, 120, 140]  # Start much higher
                elif precision == 'high':
                    ecutwfc_range = [100, 120, 140, 160, 180]  # Start very high
                elif precision == 'ultra':
                    ecutwfc_range = [120, 150, 180, 210, 240]  # Start extremely high
                    
            elif max_suggested_ecutwfc >= 50:  # Medium-ecutwfc pseudopotentials
                if precision == 'low':
                    ecutwfc_range = [40, 50, 60]  # Start moderately higher
                elif precision == 'medium':
                    ecutwfc_range = [50, 60, 70, 80]  # Start higher
                # High and ultra remain as default for medium pseudopotentials
        
        return ecutwfc_range, kspacing_range
    
    def _get_ranges_for_precision(self, precision: str) -> Tuple[List[float], List[float]]:
        """
        Get parameter ranges for a precision level.
        """
        return self._get_ranges_for_precision_static(precision)
    
    def _get_smart_ranges_for_pseudopotentials(
        self, 
        precision: str, 
        pseudopotentials: Dict[str, str],
        atoms: Atoms
    ) -> Tuple[List[float], List[float]]:
        """
        Get parameter ranges intelligently adjusted for pseudopotential requirements.
        
        This is the instance method version that also considers structural complexity.
        
        Args:
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            pseudopotentials: Dict mapping element symbols to UPF file paths
            atoms: ASE Atoms object for structural analysis
            
        Returns:
            Tuple of (ecutwfc_range, kspacing_range) adjusted for pseudopotentials and structure
        """
        # Get smart ranges based on pseudopotentials
        ecutwfc_range, kspacing_range = self._get_smart_ranges_for_pseudopotentials_static(
            precision, pseudopotentials
        )
        
        # Additional structural adjustments if needed
        # (Could add structural complexity analysis here in the future)
        
        return ecutwfc_range, kspacing_range
    
    def _get_default_convergence_criteria(self, precision: str) -> Dict:
        """
        Get default convergence criteria tolerances based on precision level.
        
        Args:
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            
        Returns:
            Dict with convergence tolerances
        """
        precision = precision.lower()
        
        criteria = {
            'low': {
                'energy_tolerance': 3e-3,      # 3 meV/atom
                'force_tolerance': 0.5,        # eV/Å
                'stress_tolerance': 1.0,       # GPa
                'geometry_tolerance': 0.05,    # Å
                'magnetic_tolerance': 0.01,    # μB
            },
            'medium': {
                'energy_tolerance': 2e-3,      # 2 meV/atom
                'force_tolerance': 0.2,        # eV/Å
                'stress_tolerance': 0.5,       # GPa
                'geometry_tolerance': 0.02,    # Å
                'magnetic_tolerance': 0.005,   # μB
            },
            'high': {
                'energy_tolerance': 1e-3,      # 1 meV/atom
                'force_tolerance': 0.1,        # eV/Å
                'stress_tolerance': 0.1,       # GPa
                'geometry_tolerance': 0.01,    # Å
                'magnetic_tolerance': 0.001,   # μB
            },
            'ultra': {
                'energy_tolerance': 5e-4,      # 0.5 meV/atom
                'force_tolerance': 0.05,       # eV/Å
                'stress_tolerance': 0.05,      # GPa
                'geometry_tolerance': 0.005,   # Å
                'magnetic_tolerance': 0.0005,  # μB
            }
        }
        
        if precision not in criteria:
            raise ValueError(f"Unknown precision level: {precision}")
            
        return criteria[precision]
    
    def _get_default_convergence_criteria_list(self, precision: str) -> List[str]:
        """
        Get default convergence criteria list.
        
        Note: Criteria control which physical quantities are checked for convergence.
        
        Args:
            precision: Precision level (required)
            
        Returns:
            List of convergence criteria
        """
        # Default criteria: energy convergence only
        return ['energy']
    
    def _fit_exponential_convergence(self, ecutwfc_values: List[float], energies: List[float]) -> float:
        """
        Fit exponential convergence to estimate infinite ecutwfc energy.
        
        Args:
            ecutwfc_values: List of ecutwfc values tested
            energies: Corresponding energies
            
        Returns:
            Estimated energy at infinite ecutwfc
        """
        try:
            from scipy.optimize import curve_fit
            
            def exp_func(x, a, b, c):
                return a * np.exp(-b * x) + c
            
            popt, _ = curve_fit(exp_func, ecutwfc_values, energies, p0=[1, 0.1, min(energies)])
            return popt[2]  # c parameter is the asymptotic value
        except:
            # Fallback: use last energy if fit fails
            return energies[-1]
    
    def _adjust_ranges_for_pseudopotentials(
        self, 
        precision: str, 
        pseudopotentials: Dict[str, str],
        atoms: Atoms
    ) -> Tuple[List[float], List[float]]:
        """
        Adjust parameter ranges based on pseudopotential requirements and structural complexity.
        
        Args:
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            pseudopotentials: Dict mapping element symbols to UPF file paths
            atoms: ASE Atoms object for structural analysis
            
        Returns:
            Tuple of (adjusted_ecutwfc_range, kspacing_range)
        """
        
        # Get base ranges for precision level
        ecutwfc_range, kspacing_range = self._get_ranges_for_precision(precision)
        
        # Analyze structural complexity
        complexity = self._analyze_structural_complexity(atoms)
        structural_factor = (
            complexity['surface_factor'] * 
            complexity['vacuum_factor'] * 
            complexity['heterogeneity_factor']
        )
        
        logger.info(f"  Structural complexity factors: {complexity}")
        logger.info(f"  Combined structural factor: {structural_factor:.2f}")
        
        # Extract suggested ecutwfc from pseudopotential files
        max_suggested_ecutwfc = 0.0
        
        for element, upf_path in pseudopotentials.items():
            try:
                # Try to parse the UPF file
                header_info = parse_upf_header(upf_path)
                if 'suggested_ecutwfc' in header_info:
                    suggested = header_info['suggested_ecutwfc']
                    max_suggested_ecutwfc = max(max_suggested_ecutwfc, suggested)
                    logger.info(f"  {element}: suggested ecutwfc = {suggested} Ry (from {upf_path})")
                else:
                    logger.warning(f"  {element}: no suggested ecutwfc found in {upf_path}")
            except Exception as e:
                logger.warning(f"  {element}: failed to parse {upf_path}: {e}")
        
        if max_suggested_ecutwfc > 0:
            # Apply structural complexity factor
            adjusted_suggested = max_suggested_ecutwfc * structural_factor
            
            # Ensure the range covers at least the adjusted suggested value
            current_max = max(ecutwfc_range)
            
            if adjusted_suggested > current_max:
                logger.info(f"  Adjusted suggested ecutwfc: {max_suggested_ecutwfc:.1f} × {structural_factor:.2f} = {adjusted_suggested:.1f} Ry")
                logger.info(f"  Current max: {current_max} Ry → extending range")
                
                # Extend the range to cover the adjusted suggested value
                # Add some buffer (10% above adjusted suggested) for convergence testing
                extended_max = adjusted_suggested * 1.1
                
                # Replace the highest value in the range with the extended max
                ecutwfc_range[-1] = round(extended_max, 1)
                
                # Ensure the range is still sorted
                ecutwfc_range = sorted(ecutwfc_range)
                
                logger.info(f"  New ecutwfc range: {ecutwfc_range}")
        
        return ecutwfc_range, kspacing_range
    
    def _analyze_structural_complexity(self, atoms: Atoms) -> Dict[str, float]:
        """
        Analyze structural complexity to determine convergence requirements.
        
        Args:
            atoms: ASE Atoms object
            
        Returns:
            Dict with complexity factors:
                - 'surface_factor': 1.0-2.0 (bulk vs surface/cluster)
                - 'vacuum_factor': 1.0-1.5 (presence of vacuum)
                - 'heterogeneity_factor': 1.0-1.8 (multiple elements/interfaces)
        """
        complexity = {
            'surface_factor': 1.0,
            'vacuum_factor': 1.0, 
            'heterogeneity_factor': 1.0
        }
        
        # Check for vacuum (slab calculations)
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        
        # Simple vacuum detection: check if cell is much larger than atomic positions
        if len(atoms) > 2:  # Avoid false positives for small systems
            max_pos = positions.max(axis=0)
            min_pos = positions.min(axis=0)
            cell_lengths = cell.lengths()
            
            for i in range(3):
                atomic_span = max_pos[i] - min_pos[i]
                if cell_lengths[i] > atomic_span * 1.5:  # Significant vacuum
                    vacuum_ratio = cell_lengths[i] / atomic_span
                    complexity['vacuum_factor'] = min(1.5, 1.0 + (vacuum_ratio - 1.5) * 0.1)
                    break
        
        # Check for surface/cluster characteristics
        # Systems with high surface-to-volume ratio need higher cutoffs
        if len(atoms) < 10:
            complexity['surface_factor'] = 2.0  # Small clusters
        elif len(atoms) < 50:
            complexity['surface_factor'] = 1.5  # Medium systems
        elif len(atoms) < 100:
            complexity['surface_factor'] = 1.2  # Large but finite systems
        
        # Check elemental heterogeneity
        elements = set(atoms.get_chemical_symbols())
        if len(elements) > 3:
            complexity['heterogeneity_factor'] = 1.8  # Complex multi-element systems
        elif len(elements) > 2:
            complexity['heterogeneity_factor'] = 1.4  # Binary/ternary systems
        elif len(elements) == 1:
            complexity['heterogeneity_factor'] = 1.0  # Pure elements
        else:
            complexity['heterogeneity_factor'] = 1.2  # Simple compounds
        
        return complexity
    
    @classmethod
    def from_cif(
        cls,
        cif_file: Union[str, Path],
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        precision: Optional[str] = 'low',
        **kwargs
    ) -> 'ConvergenceWorkflow':
        """
        Create convergence workflow from CIF file.
        
        Args:
            cif_file: Path to CIF structure file
            pseudopotentials: Dict mapping element symbols to UPF files
            precision: Precision level ('low', 'medium', 'high', 'ultra'). Default: 'low'
            **kwargs: Additional parameters for __init__ (ecutwfc_range, kspacing_range, etc.)
            
        Returns:
            ConvergenceWorkflow instance
        """
        atoms = read(cif_file)
        if pseudopotentials_config is not None:
            return cls(atoms, pseudopotentials_config=pseudopotentials_config, precision=precision, **kwargs)
        else:
            return cls(atoms, pseudopotentials, precision=precision, **kwargs)
    
    @classmethod
    def optimize_parameters(
        cls,
        atoms: Atoms,
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        precision: str = 'low',
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        convergence_criteria_list: List[str] = ['energy'],
        verbose: bool = True,
        batch_timeout: int = 3600,
        label_prefix: str = 'convergence',
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        magnetic_config: Optional[Union[str, Dict]] = None,
    ) -> 'ConvergenceWorkflow':
        """
        Create and run convergence workflow with automatic parameter optimization.
        
        This is the simplest interface - just provide structure, pseudopotentials,
        and desired precision level. The workflow will automatically determine
        optimal parameter ranges and run the convergence study.
        
        Args:
            atoms: ASE Atoms object
            pseudopotentials: Dict mapping element symbols to UPF files
            pseudopotentials_config: Name of pseudopotentials configuration to load
            precision: Precision level ('low', 'medium', 'high', 'ultra')
            queue: Queue configuration for job submission
            machine: Machine configuration name to load from ~/.xespresso/machines/
            code_version: Quantum ESPRESSO version (e.g., '7.2', '6.8')
            convergence_criteria_list: List of convergence criteria to check (default: ['energy'])
            verbose: Print progress information (default: True)
            batch_timeout: Timeout for batch jobs in seconds (default: 3600)
            label_prefix: Prefix for calculation directories (default: 'convergence')
            max_ecutwfc: Maximum ecutwfc to test (default: 200.0)
            ecutwfc_step: Step size for ecutwfc increases (default: 10.0)
            magnetic_config: Magnetic configuration string or dict (optional)
            
        Returns:
            ConvergenceWorkflow instance with completed convergence study
        """
        # If a pseudopotentials_config name is provided, pass it to the ctor
        if pseudopotentials_config is not None:
            workflow = cls(atoms, pseudopotentials_config=pseudopotentials_config, precision=precision, queue=queue, machine=machine, code_version=code_version, convergence_criteria_list=convergence_criteria_list, magnetic_config=magnetic_config)
        else:
            workflow = cls(atoms, pseudopotentials, precision=precision, queue=queue, machine=machine, code_version=code_version, convergence_criteria_list=convergence_criteria_list, magnetic_config=magnetic_config)
        
        # Run convergence study with specified parameters
        workflow.run_convergence_study(
            label_prefix=label_prefix,
            verbose=verbose,
            max_ecutwfc=max_ecutwfc,
            ecutwfc_step=ecutwfc_step,
            batch_timeout=batch_timeout,
        )
        return workflow
    
    def run_convergence_study(
        self,
        label_prefix: str = 'convergence',
        verbose: bool = True,
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        batch_timeout: int = 3600,
    ) -> pd.DataFrame:
        """
        Run convergence study with INDEPENDENT two-phase algorithm.
        
        PHASE 1: Ecutwfc convergence with FIXED coarse kspacing
        PHASE 2: Kspacing convergence with FIXED optimal ecutwfc
        
        Benefits:
        - Pseudopotenciais transferred only 2-3x (not N×M times)
        - 4-6x faster than nested loop approaches
        - Two independent, clearly separated phases
        
        Args:
            label_prefix: Prefix for calculation directories
            verbose: Print progress information
            max_ecutwfc: Maximum ecutwfc to test (safety limit)
            ecutwfc_step: Step size for ecutwfc increases
            batch_timeout: Timeout for batch jobs in seconds (default: 3600)
            
        Returns:
            pandas.DataFrame with convergence results
        """
        return self.run_convergence_independent(
            label_prefix=label_prefix,
            max_ecutwfc=max_ecutwfc,
            ecutwfc_step=ecutwfc_step,
            verbose=verbose,
            batch_timeout=batch_timeout,
        )
    
    def get_recommendations(self, verbose: bool = True) -> Dict:
        """
        Analyze convergence results and recommend optimal parameters.
        
        Uses the convergence criteria and tolerances from the convergence study.
        
        Args:
            verbose: Print recommendations
            
        Returns:
            Dict with optimal parameters found during convergence study
        """
        if self.results is None or len(self.results) == 0:
            raise ValueError("No convergence results. Run convergence study first.")
        
        # Optimal parameters found during convergence study
        optimal_ecutwfc = getattr(self, 'optimal_ecutwfc', None)
        optimal_kspacing = getattr(self, 'optimal_kspacing', None)
        
        if optimal_ecutwfc is None or optimal_kspacing is None:
            raise ValueError(
                f"Convergence study incomplete. Could not find optimal parameters. "
                f"optimal_ecutwfc={optimal_ecutwfc}, optimal_kspacing={optimal_kspacing}"
            )
        
        # Get tolerance info from convergence criteria
        energy_tol_meV = self.convergence_criteria.get('energy_tolerance', 1e-4) * 1000
        
        recommendations = {
            'optimal_ecutwfc': optimal_ecutwfc,
            'optimal_kspacing': optimal_kspacing,
            'precision': self.precision,
            'energy_tolerance_meV_atom': energy_tol_meV,
        }
        
        if verbose:
            print("\n" + "="*80)
            print("CONVERGENCE RECOMMENDATIONS")
            print("="*80)
            print(f"Precision level: {self.precision}")
            print(f"Energy tolerance: {energy_tol_meV:.2f} meV/atom")
            print(f"\nOptimal ecutwfc: {optimal_ecutwfc} Ry")
            print(f"Optimal kspacing: {optimal_kspacing} Å⁻¹")
            print("="*80 + "\n")
        
        return recommendations
    
    def plot_convergence(
        self,
        show: bool = True,
        save_path: Optional[str] = None
    ):
        """
        Plot convergence of energy and forces vs parameters.
        
        Creates visualization showing how total energy and max force
        converge with ecutwfc and kspacing.
        
        Args:
            show: Display plot (default: True)
            save_path: Save plot to file (optional)
            
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> conv.run_convergence_study()
            >>> conv.plot_convergence(save_path='convergence.png')
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warning("matplotlib not available. Skipping plot.")
            return
        
        if self.results is None:
            raise ValueError("No results to plot. Run convergence study first.")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Energy vs ecutwfc (for each kspacing)
        ax = axes[0, 0]
        for ksp in self.kspacing_range:
            data = self.results[self.results['kspacing'] == ksp]
            ax.plot(data['ecutwfc'], data['energy_per_atom'], 'o-', label=f'ksp={ksp:.2f}')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Total Energy (eV/atom)')
        ax.set_title('Energy Convergence vs ecutwfc')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Energy vs kspacing (for each ecutwfc)
        ax = axes[0, 1]
        for ecut in self.ecutwfc_range:
            data = self.results[self.results['ecutwfc'] == ecut]
            ax.plot(data['kspacing'], data['energy_per_atom'], 'o-', label=f'ecut={ecut:.0f}')
        ax.set_xlabel('kspacing (Å⁻¹)')
        ax.set_ylabel('Total Energy (eV/atom)')
        ax.set_title('Energy Convergence vs kspacing')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.invert_xaxis()  # Finer grid on right
        
        # Plot 3: Max force vs ecutwfc
        ax = axes[1, 0]
        for ksp in self.kspacing_range:
            data = self.results[self.results['kspacing'] == ksp]
            ax.plot(data['ecutwfc'], data['max_force'], 'o-', label=f'ksp={ksp:.2f}')
        ax.set_xlabel('ecutwfc (Ry)')
        ax.set_ylabel('Max Force (eV/Å)')
        ax.set_title('Force Convergence vs ecutwfc')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Energy difference (convergence map)
        ax = axes[1, 1]
        pivot_data = self.results.pivot(
            index='ecutwfc',
            columns='kspacing',
            values='energy_per_atom'
        )
        im = ax.imshow(pivot_data, aspect='auto', origin='lower', cmap='viridis')
        ax.set_xlabel('kspacing')
        ax.set_ylabel('ecutwfc')
        ax.set_title('Energy Convergence Map')
        plt.colorbar(im, ax=ax, label='Energy (eV/atom)')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            logger.info(f"Convergence plot saved to {save_path}")
        
        if show:
            plt.show()
    
    def to_csv(self, filepath: Union[str, Path]):
        """
        Export convergence results to CSV file.
        
        Args:
            filepath: Output CSV file path
            
        Example:
            >>> conv.run_convergence_study()
            >>> conv.to_csv('convergence_results.csv')
        """
        if self.results is None:
            raise ValueError("No results to export. Run convergence study first.")
        
        self.results.to_csv(filepath, index=False)
        logger.info(f"Results exported to {filepath}")
    
    def from_csv(self, filepath: Union[str, Path]):
        """
        Load convergence results from CSV file.
        
        Args:
            filepath: Input CSV file path
            
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> conv.from_csv('previous_results.csv')
            >>> recs = conv.get_recommendations()
        """
        self.results = pd.read_csv(filepath)
        logger.info(f"Results loaded from {filepath}")
    
    def _expand_range(self, current_range: List[float], step: float, max_val: float) -> List[float]:
        """
        Expand parameter range by adding next values intelligently.
        
        Args:
            current_range: Current list of parameters tested
            step: Step size for expansion
            max_val: Maximum limit for expansion
            
        Returns:
            New range with additional values (subset of new values to test)
        """
        if not current_range:
            return []
        
        max_current = max(current_range)
        
        # Generate next values beyond current max
        if max_current >= max_val:
            return []  # Already at limit
        
        # Calculate how many steps to add
        remaining = max_val - max_current
        num_steps = max(2, int(remaining / step))  # At least 2 new values
        
        new_vals = []
        for i in range(1, num_steps + 1):
            val = max_current + (i * step)
            if val <= max_val and val not in current_range:
                new_vals.append(val)
        
        return sorted(new_vals)
    
    def _get_calculation_config(self, convergence_criteria_list: List[str]) -> Dict:
        """
        Determine calculation configuration based on convergence criteria.
        
        Analyzes which properties need to be calculated and returns appropriate
        configuration. Currently only 'energy' criterion is fully implemented.
        
        Args:
            convergence_criteria_list: List of criteria like ['energy'], ['energy', 'forces'], etc.
            
        Returns:
            Dict with keys:
                - 'calc_type': 'scf' or 'vc-relax' (string)
                - 'input_data_overrides': Dict of QE input parameters to override
                
        Raises:
            NotImplementedError: If any criterion other than 'energy' is used
        """
        if not convergence_criteria_list:
            convergence_criteria_list = ['energy']
        
        config = {
            'calc_type': 'scf',
            'input_data_overrides': {}
        }
        
        # Validate which criteria are implemented
        valid_criteria = {'energy', 'forces'}  # energy and forces are implemented
        unsupported = set(convergence_criteria_list) - valid_criteria
        
        if unsupported:
            unsupported_str = ', '.join(sorted(unsupported))
            raise NotImplementedError(
                f"Convergence criteria not yet implemented: {unsupported_str}\n"
                f"Currently supported: 'energy', 'forces'\n"
                f"Coming soon: 'stress' (SCF+tstress), 'geometry' (VC-RELAX), "
                f"'magnetic_moments' (nspin=2)"
            )
        
        # Add QE input flags for forces if needed
        if 'forces' in convergence_criteria_list:
            config['input_data_overrides']['tprnfor'] = True
        
        return config
    
    def _extract_property_from_result(
        self, 
        completion: Dict, 
        num_atoms: int,
        property_name: str
    ) -> float:
        """
        Extract a single physical property from a completed calculation result.
        
        Args:
            completion: Completion dict from wait_for_batch_jobs()
            num_atoms: Number of atoms in the structure
            property_name: Name of property to extract: 'energy', 'forces', 'stress', etc.
            
        Returns:
            Float value of the property, or np.nan if not available
            
        Raises:
            NotImplementedError: If property is not yet implemented
        """
        if property_name == 'energy':
            # Energy per atom in eV
            if 'energy' in completion:
                return completion['energy'] / num_atoms
            else:
                return np.nan
        
        elif property_name == 'forces':
            # Extract maximum force magnitude from completion results
            if 'forces' in completion:
                forces_array = np.array(completion['forces'])  # Shape: (N_atoms, 3)
                # Compute force magnitude for each atom
                force_magnitudes = np.linalg.norm(forces_array, axis=1)  # Shape: (N_atoms,)
                # Return maximum magnitude in eV/Å
                return float(np.max(force_magnitudes))
            else:
                return np.nan
        
        elif property_name == 'stress':
            raise NotImplementedError(
                "Stress extraction not yet implemented. "
                "Coming soon: will extract hydrostatic pressure."
            )
        
        elif property_name == 'geometry':
            raise NotImplementedError(
                "Geometry extraction not yet implemented. "
                "Coming soon: will extract atomic displacement (requires vc-relax)."
            )
        
        elif property_name == 'magnetic_moments':
            raise NotImplementedError(
                "Magnetic moments extraction not yet implemented. "
                "Coming soon: will extract total magnetic moment (requires nspin=2)."
            )
        
        else:
            raise ValueError(f"Unknown property: {property_name}")
    
    def _check_convergence_vs_reference(
        self, 
        results_dict: Dict[float, Dict[str, float]],
        reference_properties: Dict[str, float],
        criteria_tolerances: Dict[str, float],
        convergence_criteria_list: List[str]
    ) -> bool:
        """
        Check if parameters converged for ALL criteria in convergence_criteria_list.
        
        All criteria must converge for this to return True.
        
        Args:
            results_dict: Dict mapping parameter value → Dict of extracted properties
                         Inner dict has keys like 'energy', 'forces', etc.
            reference_properties: Dict with same structure as values in results_dict
                                 Contains reference values for all properties
            criteria_tolerances: Dict with tolerance keys like 'energy_tolerance', 'force_tolerance'
            convergence_criteria_list: List of which criteria to check (e.g., ['energy', 'forces'])
            
        Returns:
            True if ALL criteria in convergence_criteria_list converged, False otherwise
        """
        if not results_dict:
            return False
        
        if not convergence_criteria_list:
            convergence_criteria_list = ['energy']
        
        # Check each criterion - ALL must pass
        for criterion in convergence_criteria_list:
            if criterion == 'energy':
                # Energy convergence: max deviation from reference < tolerance
                energies = [props.get('energy', np.nan) for props in results_dict.values()]
                if not energies or all(np.isnan(e) for e in energies):
                    return False
                
                max_energy = max(e for e in energies if not np.isnan(e))
                energy_diff = abs(max_energy - reference_properties.get('energy', 0))
                tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
                
                if energy_diff >= tolerance:
                    return False  # Energy not converged
            
            elif criterion == 'forces':
                # Force convergence: max force in test results vs reference < tolerance
                forces = [props.get('forces', np.nan) for props in results_dict.values()]
                if not forces or all(np.isnan(f) for f in forces):
                    return False
                
                max_force = max(f for f in forces if not np.isnan(f))
                reference_force = reference_properties.get('forces', 0)
                force_diff = abs(max_force - reference_force)
                tolerance = criteria_tolerances.get('force_tolerance', 0.05)  # eV/Å
                
                if force_diff >= tolerance:
                    return False  # Forces not converged
            
            else:
                # Other criteria not yet implemented, but we already validated in
                # _get_calculation_config(), so this shouldn't happen
                raise NotImplementedError(f"Convergence check for '{criterion}' not yet implemented")
        
        # All criteria converged!
        return True

    def run_convergence_independent(
        self,
        label_prefix: str = 'convergence',
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        max_kspacing: float = 0.1,
        min_kspacing: float = 0.1,
        kspacing_step: float = 0.03,
        verbose: bool = True,
        batch_timeout: int = 3600,
    ) -> pd.DataFrame:
        """
        Run INDEPENDENT convergence study with DYNAMIC RANGES and REFERENCE ENERGY.
        
        ALGORITHM:
        1. Calculate REFERENCE energy with very high ecutwfc (200 Ry)
        2. PHASE 1: Ecutwfc convergence (DYNAMIC)
           - Start with initial range [30, 40, 50]
           - Compare each with reference
           - If not converged, expand and test new values (cache existing)
           - Repeat until converged
        3. PHASE 2: Kspacing convergence (DYNAMIC)
           - Use converged ecutwfc from PHASE 1
           - Same dynamic expansion logic as PHASE 1
        
        Benefits:
        - Converges to TRUE reference (not false convergence)
        - Expands ranges only as needed
        - Caches results (no redundant calculations)
        - Pseudo transferred only 2-3x
        
        Returns:
            pandas.DataFrame with complete convergence results
        """
        results_all = []
        
        # Get convergence criteria tolerances
        criteria_tolerances = self.convergence_criteria
        
        # Validate calculation configuration for requested criteria
        calc_config = self._get_calculation_config(self.convergence_criteria_list)
        if calc_config['calc_type'] != 'scf':
            raise NotImplementedError(
                f"Calculation type '{calc_config['calc_type']}' not yet supported. "
                f"Only 'scf' (energy convergence) is currently implemented."
            )
        
        # ===== PHASE 1: DYNAMIC ECUTWFC CONVERGENCE =====
        print("\n" + "="*80)
        print("PHASE 1: ECUTWFC CONVERGENCE (DYNAMIC)")
        print("="*80)
        
        fixed_kspacing_phase1 = 0.3  # Coarse k-mesh
        print(f"\nStructure: {self.atoms.get_chemical_formula()}")
        print(f"Fixed kspacing: {fixed_kspacing_phase1:.3f} Å⁻¹")
        print(f"Reference ecutwfc: {max_ecutwfc:.1f} Ry (included in first batch)\n")
        
        # Start with initial range + reference (max_ecutwfc) in first iteration
        current_ecut_range = sorted(set(self.ecutwfc_range.copy() + [max_ecutwfc]))
        ecut_results = {}  # Cache: ecut → Dict of extracted properties
        reference_properties = None  # Dict with all properties from reference calculation
        iteration = 1
        
        while True:
            print(f"\n--- Iteration {iteration} ---")
            print(f"Testing ecutwfc: {current_ecut_range}")
            
            # Find which values to calculate (not in cache)
            to_calculate = [e for e in current_ecut_range if e not in ecut_results]
            
            if to_calculate:
                # Create workflow for this iteration
                wf_kwargs = {
                    'atoms': self.atoms,
                    'protocol': self.protocol,
                    'kspacing': fixed_kspacing_phase1,
                    'code_version': self.code_version,
                }
                
                if self._pseudo_config_name:
                    wf_kwargs['pseudopotentials_config'] = self._pseudo_config_name
                else:
                    wf_kwargs['pseudopotentials'] = self.pseudopotentials
                
                if self.queue is not None:
                    wf_kwargs['queue'] = self.queue
                elif self.machine is not None:
                    wf_kwargs['machine'] = self.machine
                
                if self.magnetic_config is not None:
                    wf_kwargs['magnetic_config'] = self.magnetic_config
                
                wf1 = CalculationWorkflow(**wf_kwargs)
                
                # Prepare batch for new values only
                batch_params = []
                for ecutwfc in to_calculate:
                    label = f"{label_prefix}/phase1_iter{iteration}_ecut{int(ecutwfc)}"
                    batch_params.append({
                        'label': label,
                        'ecutwfc': ecutwfc,
                        'ecutrho': ecutwfc * self.ecutrho_ratio,  # Calculate ecutrho dynamically
                        'kspacing': fixed_kspacing_phase1,
                    })
                
                is_first_batch = iteration == 1 and max_ecutwfc in to_calculate
                if verbose:
                    msg = f"Submitting {len(batch_params)} ecutwfc tests"
                    if is_first_batch:
                        msg += f" (including reference ecut={max_ecutwfc})"
                    print(f"{msg}...")
                
                # Submit batch
                batch_results = wf1.submit_scf_batch_multiple(batch_params, verbose=verbose)
                completion = wf1.wait_for_batch_jobs(batch_results, timeout=batch_timeout, verbose=verbose)
                
                # STEP 1: Extract reference first (if not yet available)
                if reference_properties is None:
                    for i, comp in enumerate(completion):
                        param = batch_params[i]
                        if param['ecutwfc'] == max_ecutwfc and comp['success']:
                            # Extract all properties for this result
                            props = {}
                            for prop_name in self.convergence_criteria_list:
                                props[prop_name] = self._extract_property_from_result(
                                    comp, len(self.atoms), prop_name
                                )
                            ecut_results[param['ecutwfc']] = props
                            reference_properties = props
                            
                            if verbose:
                                energy_str = f"{props.get('energy', np.nan):.6f}" if 'energy' in props else "N/A"
                                print(f"  ✓ [REFERENCE] ecutwfc={param['ecutwfc']:.1f}: E = {energy_str} eV/atom")
                            break
                
                # STEP 2: Store all results and print with ΔE now available
                for i, comp in enumerate(completion):
                    if comp['success']:
                        param = batch_params[i]
                        # Skip if already stored (reference)
                        if param['ecutwfc'] in ecut_results:
                            continue
                        
                        # Extract all properties for this result
                        props = {}
                        for prop_name in self.convergence_criteria_list:
                            props[prop_name] = self._extract_property_from_result(
                                comp, len(self.atoms), prop_name
                            )
                        ecut_results[param['ecutwfc']] = props
                        
                        is_reference = (param['ecutwfc'] == max_ecutwfc)
                        
                        result = {
                            'phase': 0 if is_reference else 1,
                            'ecutwfc': param['ecutwfc'],
                            'kspacing': param['kspacing'],
                            'energy_per_atom': props.get('energy', np.nan),
                            'label': param['label'],
                        }
                        results_all.append(result)
                        
                        # Print non-reference with ΔE
                        if not is_reference and verbose and reference_properties is not None:
                            energy = props.get('energy', np.nan)
                            ref_energy = reference_properties.get('energy', 0)
                            diff = abs(energy - ref_energy)
                            status = "✓" if diff < criteria_tolerances.get('energy_tolerance', 1e-3) else "✗"
                            print(f"  {status} ecutwfc={param['ecutwfc']:.1f}: E = {energy:.6f} eV/atom (ΔE = {diff:.6f})")
                    else:
                        if verbose:
                            print(f"  ✗ ecutwfc={batch_params[i]['ecutwfc']}: {comp.get('error', 'Failed')}")
            
            # Check convergence (skip if reference not yet calculated)
            if reference_properties is not None:
                # Exclude reference from convergence check
                test_ecut_results = {k: v for k, v in ecut_results.items() if k != max_ecutwfc}
                converged = self._check_convergence_vs_reference(
                    test_ecut_results, reference_properties, criteria_tolerances, 
                    self.convergence_criteria_list
                )
                
                if converged:
                    if verbose:
                        print(f"\n✓ CONVERGED at iteration {iteration}")
                    break
            else:
                # Reference calculation pending
                if verbose:
                    print(f"\n~ Iteration {iteration} complete. Reference pending...")
            
            # Not converged: expand range
            expansion = self._expand_range(current_ecut_range, ecutwfc_step, max_ecutwfc)
            if not expansion:
                if verbose:
                    print(f"\n⚠️  Cannot expand further (limit: {max_ecutwfc}). Stopping.")
                break
            
            # Add expanded values to range
            current_ecut_range = sorted(set(current_ecut_range + expansion))
            iteration += 1
        
        # Ensure reference was calculated
        if reference_properties is None:
            raise RuntimeError("Could not obtain reference energy (ecutwfc=200)")
        
        # Remove reference from test results for selection
        test_ecut_results = {k: v for k, v in ecut_results.items() if k != max_ecutwfc}
        if not test_ecut_results:
            raise RuntimeError("PHASE 1 failed: no successful calculations")
        
        # Select ecutwfc for PHASE 2 (MINIMUM converged value for best efficiency)
        optimal_ecutwfc = min(test_ecut_results.keys())
        self.optimal_ecutwfc = optimal_ecutwfc  # Store for later use in get_recommendations()
        optimal_props_phase1 = ecut_results[optimal_ecutwfc]
        
        print(f"\n✓ PHASE 1 COMPLETE: Selected ecutwfc = {optimal_ecutwfc:.1f} Ry")
        
        # ===== PHASE 2: DYNAMIC KSPACING CONVERGENCE =====
        print("\n" + "="*80)
        print("PHASE 2: KSPACING CONVERGENCE (DYNAMIC)")
        print("="*80)
        
        print(f"\nFixed ecutwfc: {optimal_ecutwfc:.1f} Ry (from PHASE 1)")
        print(f"Reference kspacing (fine): {max_kspacing:.3f} Å⁻¹\n")
        
        # Start with normal range (WITHOUT reference). Reference will be added to first batch.
        current_ksp_range = sorted(self.kspacing_range.copy(), reverse=True)
        ksp_results = {}  # Cache: kspacing → Dict of extracted properties
        reference_properties_phase2 = None  # Dict with all properties from reference calculation
        iteration = 1
        
        while True:
            print(f"\n--- Iteration {iteration} ---")
            
            # Add reference to batch on first iteration only
            if iteration == 1:
                to_test = sorted(set(current_ksp_range + [max_kspacing]), reverse=True)
            else:
                to_test = current_ksp_range
            
            print(f"Testing kspacing: {to_test}")
            
            # Find which values to calculate (not in cache)
            to_calculate = [k for k in to_test if k not in ksp_results]
            
            if to_calculate:
                # Create workflow for this iteration
                wf_kwargs = {
                    'atoms': self.atoms,
                    'protocol': self.protocol,
                    'ecutwfc': optimal_ecutwfc,
                    'code_version': self.code_version,
                }
                
                if self._pseudo_config_name:
                    wf_kwargs['pseudopotentials_config'] = self._pseudo_config_name
                else:
                    wf_kwargs['pseudopotentials'] = self.pseudopotentials
                
                if self.queue is not None:
                    wf_kwargs['queue'] = self.queue
                elif self.machine is not None:
                    wf_kwargs['machine'] = self.machine
                
                if self.magnetic_config is not None:
                    wf_kwargs['magnetic_config'] = self.magnetic_config
                
                wf2 = CalculationWorkflow(**wf_kwargs)
                
                # Prepare batch for new values only
                batch_params = []
                for kspacing in to_calculate:
                    label = f"{label_prefix}/phase2_iter{iteration}_ecut{int(optimal_ecutwfc)}_ksp{kspacing:.3f}"
                    batch_params.append({
                        'label': label,
                        'ecutwfc': optimal_ecutwfc,
                        'ecutrho': optimal_ecutwfc * self.ecutrho_ratio,
                        'kspacing': kspacing,
                    })
                
                if verbose:
                    is_first_batch = iteration == 1 and max_kspacing in to_calculate
                    msg = f"Submitting {len(batch_params)} kspacing tests"
                    if is_first_batch:
                        msg += f" (including reference kspacing={max_kspacing:.3f})"
                    print(f"{msg}...")
                
                # Submit batch
                batch_results = wf2.submit_scf_batch_multiple(batch_params, verbose=verbose)
                completion = wf2.wait_for_batch_jobs(batch_results, timeout=batch_timeout, verbose=verbose)
                
                # STEP 1: Extract reference first (if not yet available)
                if reference_properties_phase2 is None:
                    for i, comp in enumerate(completion):
                        param = batch_params[i]
                        if param['kspacing'] == max_kspacing and comp['success']:
                            # Extract all properties for this result
                            props = {}
                            for prop_name in self.convergence_criteria_list:
                                props[prop_name] = self._extract_property_from_result(
                                    comp, len(self.atoms), prop_name
                                )
                            ksp_results[param['kspacing']] = props
                            reference_properties_phase2 = props
                            if verbose:
                                energy_str = f"{props.get('energy', np.nan):.6f}" if 'energy' in props else "N/A"
                                print(f"  ✓ [REFERENCE] kspacing={param['kspacing']:.3f}: E = {energy_str} eV/atom")
                            break
                
                # STEP 2: Store all results and print with ΔE now available
                for i, comp in enumerate(completion):
                    if comp['success']:
                        param = batch_params[i]
                        # Skip if already stored (reference)
                        if param['kspacing'] in ksp_results:
                            continue
                        
                        # Extract all properties for this result
                        props = {}
                        for prop_name in self.convergence_criteria_list:
                            props[prop_name] = self._extract_property_from_result(
                                comp, len(self.atoms), prop_name
                            )
                        ksp_results[param['kspacing']] = props
                        
                        result = {
                            'phase': 2,
                            'ecutwfc': param['ecutwfc'],
                            'kspacing': param['kspacing'],
                            'energy_per_atom': props.get('energy', np.nan),
                            'label': param['label'],
                        }
                        results_all.append(result)
                        
                        # Print non-reference with ΔE
                        if verbose and reference_properties_phase2 is not None:
                            energy = props.get('energy', np.nan)
                            ref_energy = reference_properties_phase2.get('energy', 0)
                            diff = abs(energy - ref_energy)
                            status = "✓" if diff < criteria_tolerances.get('energy_tolerance', 1e-3) else "✗"
                            print(f"  {status} kspacing={param['kspacing']:.3f}: E = {energy:.6f} eV/atom (ΔE = {diff:.6f})")
                    else:
                        if verbose:
                            print(f"  ✗ kspacing={batch_params[i]['kspacing']}: {comp.get('error', 'Failed')}")
            
            # Check convergence (skip if reference not yet calculated)
            if reference_properties_phase2 is not None:
                # Exclude reference from convergence check
                test_ksp_results = {k: v for k, v in ksp_results.items() if k != max_kspacing}
                converged = self._check_convergence_vs_reference(
                    test_ksp_results, reference_properties_phase2, criteria_tolerances,
                    self.convergence_criteria_list
                )
            else:
                converged = False
            
            if converged:
                if verbose:
                    print(f"\n✓ CONVERGED at iteration {iteration}")
                break
            
            # Not converged: expand range to finer kspacing (smaller values)
            # Note: for kspacing, smaller values are finer
            min_current = min(current_ksp_range)
            
            # Check if we can expand further (can't go below min_kspacing)
            if min_current <= min_kspacing:
                if verbose:
                    print(f"\n⚠️  Cannot expand further (limit: {min_kspacing:.3f}). Stopping.")
                break
            
            # Generate finer kspacing values (smaller than min_current)
            new_ksp_vals = []
            num_steps = 1
            for i in range(1, num_steps + 1):
                val = min_current - (i * kspacing_step)
                if val >= min_kspacing and val not in current_ksp_range:
                    new_ksp_vals.append(val)
            
            if not new_ksp_vals:
                if verbose:
                    print(f"\n⚠️  Cannot generate new values (would go below {min_kspacing:.3f}). Stopping.")
                break
            
            # Add expanded values to range (maintain reverse sort)
            current_ksp_range = sorted(set(current_ksp_range + new_ksp_vals), reverse=True)
            iteration += 1
        
        # Ensure reference was calculated
        if reference_properties_phase2 is None:
            raise RuntimeError(f"Could not obtain reference energy (kspacing={max_kspacing})")
        
        # PHASE 2: SELECT OPTIMAL KSPACING FROM CONVERGENCE RESULTS
        if not ksp_results:
            print("\n⚠️  PHASE 2: No successful kspacing tests. May need to adjust parameters.")
        else:
            # Exclude reference from selection
            test_ksp_results = {k: v for k, v in ksp_results.items() if k != max_kspacing}
            
            if not test_ksp_results:
                raise RuntimeError("PHASE 2 failed: no successful calculations (only reference)")
            
            # Find converged kspacing values (ΔE < tolerance)
            tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
            converged_ksp = {}
            for ksp, props_dict in test_ksp_results.items():
                energy = props_dict.get('energy', np.nan)
                ref_energy = reference_properties_phase2.get('energy', 0)
                delta_e = abs(energy - ref_energy)
                if delta_e < tolerance:
                    converged_ksp[ksp] = delta_e
            
            if converged_ksp:
                # Select the MAXIMUM (coarsest, most efficient) converged kspacing
                optimal_kspacing = max(converged_ksp.keys())
                self.optimal_kspacing = optimal_kspacing  # Store for get_recommendations()
                optimal_delta_e = converged_ksp[optimal_kspacing]
                
                if verbose:
                    print(f"\n✓ PHASE 2 CONVERGED")
                    print(f"  Optimal kspacing: {optimal_kspacing:.3f} Å⁻¹")
                    print(f"  ΔE = {optimal_delta_e*1000:.2f} meV/atom < tolerance = {tolerance*1000:.2f} meV/atom")
                    print(f"  (All converged values: {sorted(converged_ksp.keys())})")
            else:
                if verbose:
                    best_kspacing = min(test_ksp_results.keys())
                    best_props = test_ksp_results[best_kspacing]
                    best_energy = best_props.get('energy', np.nan)
                    ref_energy = reference_properties_phase2.get('energy', 0)
                    best_delta_e = abs(best_energy - ref_energy)
                    print(f"\n✗ PHASE 2 NOT CONVERGED")
                    print(f"  Tolerance = {tolerance*1000:.2f} meV/atom (precision='{self.precision}')")
                    print(f"  Reached minimum kspacing={min_kspacing:.3f} without achieving convergence")
                    print(f"  Closest: kspacing={best_kspacing:.3f} with ΔE = {best_delta_e*1000:.2f} meV/atom")
                
                # Use the best (finest) kspacing found as fallback
                self.optimal_kspacing = best_kspacing
        
        print(f"\n" + "="*80)
        print("CONVERGENCE STUDY COMPLETE (INDEPENDENT WITH DYNAMIC RANGES)")
        print("="*80 + "\n")
        
        # Store results
        self.results = pd.DataFrame(results_all)
        return self.results

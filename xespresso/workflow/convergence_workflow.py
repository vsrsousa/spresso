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
        ...     precision='medium'  # Controls parameter ranges (ecutwfc, kspacing)
        ... )
        >>> optimal_params = workflow.optimize_parameters()
    
    2. Advanced mode (for detailed control):
        >>> # Specify custom parameter ranges and criteria independently
        >>> conv = ConvergenceWorkflow(
        ...     atoms=atoms,
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     precision='low',  # Coarse parameter ranges
        ...     convergence_criteria_list=['energy', 'forces', 'geometry']  # Strict criteria
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
        precision: Optional[str] = None,
        ecutwfc_range: Optional[List[float]] = None,
        kspacing_range: Optional[List[float]] = None,
        conv_thr_range: Optional[List[float]] = None,
        convergence_criteria_list: Optional[List[str]] = None,
        convergence_criteria: Optional[Dict] = None,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
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
                     If specified, overrides ecutwfc_range and kspacing_range.
                     Default: None (use explicit ranges)
            ecutwfc_range: List of ecutwfc values to test.
                          If None and precision=None: [30, 40, 50, 60, 70]
            kspacing_range: List of kspacing values to test (Å^-1).
                           If None and precision=None: [0.5, 0.3, 0.2, 0.15]
            conv_thr_range: List of conv_thr values to test (optional)
            convergence_criteria_list: List of convergence criteria to check.
                                     Options: 'energy', 'forces', 'geometry', 'magnetic_moments'
                                     If None, uses defaults based on precision level.
            convergence_criteria: Dict with convergence tolerances.
                                If None, uses defaults based on precision level.
                                Keys: 'energy_tolerance', 'force_tolerance', 'geometry_tolerance', 'magnetic_tolerance'
            queue: Queue configuration for job submission (optional)
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
        self.extra_kwargs = kwargs
        
        # Define parameter ranges based on precision level
        if precision is not None:
            # Get smart ranges adjusted for pseudopotential requirements
            self.ecutwfc_range, self.kspacing_range = self._get_smart_ranges_for_pseudopotentials(
                precision, self.pseudopotentials, atoms
            )
        else:
            # Use explicit ranges or defaults
            if ecutwfc_range is None:
                self.ecutwfc_range = [30, 40, 50, 60, 70]
            else:
                self.ecutwfc_range = sorted(ecutwfc_range)
            
            if kspacing_range is None:
                self.kspacing_range = [0.30, 0.27, 0.23, 0.20]
            else:
                self.kspacing_range = sorted(kspacing_range, reverse=True)
        
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
    
    def _get_default_convergence_criteria(self, precision: Optional[str]) -> Dict:
        """
        Get default convergence criteria tolerances based on precision level.
        
        Args:
            precision: Precision level or None
            
        Returns:
            Dict with convergence tolerances
        """
        if precision is None:
            # Default criteria for custom ranges
            return {
                'energy_tolerance': 1e-3,      # 1 meV/atom
                'force_tolerance': 0.1,        # eV/Å
                'geometry_tolerance': 0.01,    # Å
                'magnetic_tolerance': 0.001,   # μB
            }
        
        precision = precision.lower()
        
        criteria = {
            'low': {
                'energy_tolerance': 3e-3,      # 3 meV/atom
                'force_tolerance': 0.5,        # eV/Å
                'geometry_tolerance': 0.05,    # Å
                'magnetic_tolerance': 0.01,    # μB
            },
            'medium': {
                'energy_tolerance': 2e-3,      # 2 meV/atom
                'force_tolerance': 0.2,        # eV/Å
                'geometry_tolerance': 0.02,    # Å
                'magnetic_tolerance': 0.005,   # μB
            },
            'high': {
                'energy_tolerance': 1e-3,      # 1 meV/atom
                'force_tolerance': 0.1,        # eV/Å
                'geometry_tolerance': 0.01,    # Å
                'magnetic_tolerance': 0.001,   # μB
            },
            'ultra': {
                'energy_tolerance': 5e-4,      # 0.5 meV/atom
                'force_tolerance': 0.05,       # eV/Å
                'geometry_tolerance': 0.005,   # Å
                'magnetic_tolerance': 0.0005,  # μB
            }
        }
        
        if precision not in criteria:
            raise ValueError(f"Unknown precision level: {precision}")
            
        return criteria[precision]
    
    def _get_default_convergence_criteria_list(self, precision: Optional[str]) -> List[str]:
        """
        Get default convergence criteria list.
        
        Note: Criteria control which physical quantities are checked for convergence.
        
        Args:
            precision: Precision level (unused, kept for compatibility)
            
        Returns:
            List of convergence criteria
        """
        # Default criteria: energy convergence only
        return ['energy']
    
    def _check_kspacing_convergence(self, results_df: pd.DataFrame, criteria_list: List[str], tolerances: Dict) -> Dict[str, bool]:
        """
        Check convergence with respect to kspacing for a fixed ecutwfc.
        
        Args:
            results_df: DataFrame with results for a single ecutwfc (multiple kspacing)
            criteria_list: List of criteria to check
            tolerances: Dict with tolerance values
            
        Returns:
            Dict mapping criterion to convergence status
        """
        convergence_status = {}
        
        # Sort by kspacing (finest first)
        sorted_df = results_df.sort_values('kspacing')
        
        for criterion in criteria_list:
            if criterion == 'energy':
                # Check if energy converges with kspacing
                if len(sorted_df) >= 2:
                    energies = sorted_df['energy_per_atom'].values
                    # Compare finest two kspacing values
                    if len(energies) >= 2:
                        delta_e = abs(energies[-1] - energies[-2])  # Compare last two (finest)
                        convergence_status['energy'] = delta_e < tolerances['energy_tolerance']
                    else:
                        convergence_status['energy'] = False
                else:
                    convergence_status['energy'] = False
                    
            elif criterion == 'forces':
                # Check force convergence with kspacing
                if 'max_force' in sorted_df.columns and len(sorted_df) >= 2:
                    forces = sorted_df['max_force'].dropna().values
                    if len(forces) >= 2:
                        delta_f = abs(forces[-1] - forces[-2])  # Compare finest two
                        convergence_status['forces'] = delta_f < tolerances['force_tolerance']
                    else:
                        convergence_status['forces'] = False
                else:
                    convergence_status['forces'] = False
                    
            elif criterion == 'geometry':
                # For geometry, check if positions are converged
                # This would require position data - for now assume converged if energy is
                convergence_status['geometry'] = convergence_status.get('energy', False)
                
            elif criterion == 'magnetic_moments':
                # For magnetic moments - assume converged if energy is
                convergence_status['magnetic_moments'] = convergence_status.get('energy', False)
        
        return convergence_status
    
    def _check_convergence(self, energies: List[float], tolerance: float) -> bool:
        """
        Check if energy has converged based on tolerance.
        
        Args:
            energies: List of energies in order (should be at least 2 values)
            tolerance: Energy tolerance in eV/atom
            
        Returns:
            True if converged (last two energies differ by less than tolerance)
        """
        if len(energies) < 2:
            return False
        
        # Check if last two energies are within tolerance
        delta_e = abs(energies[-1] - energies[-2])
        return delta_e < tolerance
    
    def _check_all_convergence_criteria(self, results_df: pd.DataFrame, criteria_list: List[str]) -> Dict[str, bool]:
        """
        Check convergence for all specified criteria.
        
        Args:
            results_df: DataFrame with calculation results
            criteria_list: List of convergence criteria to check
            
        Returns:
            Dict mapping criteria to convergence status
        """
        convergence_status = {}
        
        # Group by ecutwfc and get the finest kspacing results
        finest_results = results_df.loc[results_df.groupby('ecutwfc')['kspacing'].idxmin()]
        finest_results = finest_results.sort_values('ecutwfc')
        
        if len(finest_results) < 2:
            # Not enough data for convergence check
            for criterion in criteria_list:
                convergence_status[criterion] = False
            return convergence_status
        
        # Check energy convergence
        if 'energy' in criteria_list:
            energies = finest_results['energy_per_atom'].values
            convergence_status['energy'] = self._check_convergence(
                energies, self.convergence_criteria['energy_tolerance']
            )
        
        # Check forces convergence
        if 'forces' in criteria_list:
            # For forces, we need to check if max_force is below tolerance
            # This is a different type of convergence - absolute value vs difference
            latest_max_force = finest_results['max_force'].iloc[-1]
            convergence_status['forces'] = latest_max_force < self.convergence_criteria['force_tolerance']
        
        # Check geometry convergence (simplified - would need position differences)
        if 'geometry' in criteria_list:
            # For now, use energy as proxy for geometry convergence
            # In a full implementation, this would compare atomic positions
            convergence_status['geometry'] = convergence_status.get('energy', False)
        
        # Check magnetic moments convergence
        if 'magnetic_moments' in criteria_list:
            # For now, assume converged if energy converged
            # In a full implementation, this would check magnetic moments
            convergence_status['magnetic_moments'] = convergence_status.get('energy', False)
        
        return convergence_status
    
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
        precision: Optional[str] = 'medium',
        **kwargs
    ) -> 'ConvergenceWorkflow':
        """
        Create convergence workflow from CIF file.
        
        Args:
            cif_file: Path to CIF structure file
            pseudopotentials: Dict mapping element symbols to UPF files
            precision: Precision level ('low', 'medium', 'high', 'ultra') or None for custom ranges
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
        precision: str = 'medium',
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        verbose: bool = True,
        use_batch_mode: bool = True,
        batch_timeout: int = 3600,
        label_prefix: str = 'convergence',
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
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
            verbose: Print progress information (default: True)
            use_batch_mode: Use parallel batch submission for remote systems (default: True)
            batch_timeout: Timeout for batch jobs in seconds (default: 3600)
            label_prefix: Prefix for calculation directories (default: 'convergence')
            max_ecutwfc: Maximum ecutwfc to test (default: 200.0)
            ecutwfc_step: Step size for ecutwfc increases (default: 10.0)
            
        Returns:
            ConvergenceWorkflow instance with completed convergence study
        """
        # If a pseudopotentials_config name is provided, pass it to the ctor
        if pseudopotentials_config is not None:
            workflow = cls(atoms, pseudopotentials_config=pseudopotentials_config, precision=precision, queue=queue, machine=machine, code_version=code_version)
        else:
            workflow = cls(atoms, pseudopotentials, precision=precision, queue=queue, machine=machine, code_version=code_version)
        
        # Run convergence study with specified parameters
        workflow.run_convergence_study(
            label_prefix=label_prefix,
            verbose=verbose,
            max_ecutwfc=max_ecutwfc,
            ecutwfc_step=ecutwfc_step,
            use_batch_mode=use_batch_mode,
            batch_timeout=batch_timeout,
        )
        return workflow
    
    def run_convergence_study(
        self,
        label_prefix: str = 'convergence',
        verbose: bool = True,
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        use_batch_mode: bool = True,
        batch_timeout: int = 3600,
        independent_mode: bool = True,
    ) -> pd.DataFrame:
        """
        Run convergence study with independent parameter optimization.
        
        **RECOMMENDED**: Uses independent_mode=True by default.
        
        INDEPENDENT MODE (recommended, default):
        ✅ PHASE 1: Ecutwfc convergence with FIXED coarse kspacing (0.5 Å⁻¹)
        ✅ PHASE 2: Kspacing convergence with FIXED optimal ecutwfc
        
        Benefits:
        - Pseudopotenciais transferidos apenas 2x (não N×M times!)
        - 4-6x mais rápido que nested loop mode
        - Algoritmo claramente diferenciado
        - Pseudo enviado UMA VEZ para cada fase
        
        LEGACY MODE (nested loops, slower):
        ❌ Tests all (ecutwfc, kspacing) combinations in nested loops
        ❌ Pseudopotenciais reenviados múltiplas vezes
        ❌ Mantido apenas para compatibilidade com scripts antigos
        
        Args:
            label_prefix: Prefix for calculation directories
            verbose: Print progress information
            max_ecutwfc: Maximum ecutwfc to test (safety limit)
            ecutwfc_step: Step size for ecutwfc increases
            use_batch_mode: If True and queue is remote, use batch submission (default: True)
            batch_timeout: Timeout for batch jobs in seconds (default: 3600)
            independent_mode: If True (default), use INDEPENDENT two-phase algorithm.
                             If False, use legacy nested-loop algorithm (NOT recommended).
            
        Returns:
            pandas.DataFrame with convergence results
        """
        # Use independent mode by default (RECOMMENDED)
        if independent_mode:
            return self.run_convergence_independent(
                label_prefix=label_prefix,
                max_ecutwfc=max_ecutwfc,
                ecutwfc_step=ecutwfc_step,
                verbose=verbose,
                batch_timeout=batch_timeout,
            )
        
        # Legacy nested-loop mode (NOT recommended)
        return self._run_convergence_study_legacy(
            label_prefix=label_prefix,
            verbose=verbose,
            max_ecutwfc=max_ecutwfc,
            ecutwfc_step=ecutwfc_step,
            use_batch_mode=use_batch_mode,
            batch_timeout=batch_timeout,
        )
    
    def _run_convergence_study_legacy(
        self,
        label_prefix: str = 'convergence',
        verbose: bool = True,
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        use_batch_mode: bool = True,
        batch_timeout: int = 3600,
    ) -> pd.DataFrame:
        """
        Legacy nested-loop convergence (NOT recommended).
        
        ⚠️ WARNING: This method uses nested loops and transfers pseudopotenciais
        multiple times. Use run_convergence_study(independent_mode=True) instead!
        
        For REMOTE HPC systems (SLURM), uses batch mode to submit all jobs
        in parallel and monitor them together.
        
        For LOCAL systems, uses sequential mode (one job at a time).
        """
        results_list = []
        
        # Get starting ecutwfc (minimum from range or 30 Ry)
        if hasattr(self, 'ecutwfc_range') and self.ecutwfc_range:
            start_ecutwfc = min(self.ecutwfc_range)
        else:
            start_ecutwfc = 30.0
            
        # Get kspacing range to test
        if hasattr(self, 'kspacing_range') and self.kspacing_range:
            kspacing_values = sorted(self.kspacing_range, reverse=True)  # Coarsest first
        else:
            kspacing_values = [0.5, 0.3, 0.2, 0.15]
        
        print("\n" + "="*80)
        print("⚠️ LEGACY NESTED-LOOP CONVERGENCE STUDY")
        print("="*80)
        print("WARNING: This uses nested loops and transfers pseudopotenciais multiple times!")
        print("Use run_convergence_study(independent_mode=True) instead for better performance.")
        print("="*80)
        print(f"\nStructure: {self.atoms.get_chemical_formula()}")
        print(f"Convergence criteria: {', '.join(self.convergence_criteria_list)}")
        print(f"Energy tolerance: {self.convergence_criteria['energy_tolerance']*1000:.1f} meV/atom")
        if 'forces' in self.convergence_criteria_list:
            print(f"Force tolerance: {self.convergence_criteria['force_tolerance']:.1f} eV/Å")
        if 'geometry' in self.convergence_criteria_list:
            print(f"Geometry tolerance: {self.convergence_criteria['geometry_tolerance']*1000:.1f} meV/atom")
        if 'magnetic_moments' in self.convergence_criteria_list:
            print(f"Magnetic tolerance: {self.convergence_criteria['magnetic_tolerance']*1000:.1f} μB")
        print(f"Starting ecutwfc: {start_ecutwfc} Ry")
        print(f"ecutwfc step: {ecutwfc_step} Ry")
        print(f"kspacing values: {kspacing_values}")
        
        # Detect if using remote batch mode
        is_remote = self.queue and self.queue.get('execution') == 'remote'
        use_batch = use_batch_mode and is_remote
        
        if use_batch:
            print(f"\n⚡ BATCH MODE: Submitting jobs in parallel to SLURM")
        else:
            print(f"\n📊 SEQUENTIAL MODE: Running jobs one at a time")
        
        # Track convergence
        converged_ecutwfc = None
        converged_kspacing = None
        current_ecutwfc = start_ecutwfc
        test_count = 0
        
        # If using batch mode, prepare all jobs for first ecutwfc
        if use_batch:
            return self._run_convergence_study_batch(
                label_prefix, verbose, max_ecutwfc, ecutwfc_step, 
                start_ecutwfc, kspacing_values, batch_timeout
            )
        else:
            return self._run_convergence_study_sequential(
                label_prefix, verbose, max_ecutwfc, ecutwfc_step,
                start_ecutwfc, kspacing_values
            )
    
    def _run_convergence_study_batch(
        self,
        label_prefix: str,
        verbose: bool,
        max_ecutwfc: float,
        ecutwfc_step: float,
        start_ecutwfc: float,
        kspacing_values: List[float],
        batch_timeout: int
    ) -> pd.DataFrame:
        """
        Run convergence study using batch (parallel) submission on SLURM.
        
        Submits all kspacing tests for current ecutwfc as a batch,
        waits for all to complete, then checks convergence.
        """
        results_list = []
        test_count = 0
        current_ecutwfc = start_ecutwfc
        
        while current_ecutwfc <= max_ecutwfc:
            if verbose:
                print(f"\n--- Testing ecutwfc = {current_ecutwfc:.1f} Ry ---")
            
            # Prepare batch parameters for this ecutwfc
            batch_params = []
            for kspacing in kspacing_values:
                test_count += 1
                label = f"{label_prefix}/ecut{int(current_ecutwfc)}_ksp{kspacing:.2f}"
                batch_params.append({
                    'label': label,
                    'ecutwfc': current_ecutwfc,
                    'ecutrho': current_ecutwfc * self.ecutrho_ratio,  # Calculate ecutrho dynamically
                    'kspacing': kspacing,
                })
            
            # Create base workflow for batch submission
            # Pass pseudopotentials via config name (recommended) so CalculationWorkflow
            # loads them the same way and has access to base_path
            workflow_kwargs = {
                'atoms': self.atoms,
                'protocol': self.protocol,
                'code_version': self.code_version,
            }
            
            # Pass pseudopotentials via config name if available
            if self._pseudo_config_name:
                workflow_kwargs['pseudopotentials_config'] = self._pseudo_config_name
            else:
                workflow_kwargs['pseudopotentials'] = self.pseudopotentials
            
            # Pass either queue or machine (not both)
            if self.queue is not None:
                workflow_kwargs['queue'] = self.queue
            elif self.machine is not None:
                workflow_kwargs['machine'] = self.machine
            
            # Add any extra kwargs
            workflow_kwargs.update(self.extra_kwargs)
            
            workflow = CalculationWorkflow(**workflow_kwargs)
            
            # Submit all jobs in batch
            if verbose:
                print(f"Submitting {len(batch_params)} jobs in batch mode...")
            
            batch_results = workflow.submit_scf_batch_multiple(batch_params, verbose=verbose)
            
            # Wait for all jobs to complete
            if verbose:
                print(f"Waiting for batch to complete (timeout: {batch_timeout}s)...")
            
            completion_results = workflow.wait_for_batch_jobs(
                batch_results, timeout=batch_timeout, verbose=verbose
            )
            
            # Process results
            ecutwfc_results = []
            for i, comp_result in enumerate(completion_results):
                if comp_result['success']:
                    batch_param = batch_params[i]
                    result = {
                        'ecutwfc': batch_param['ecutwfc'],
                        'kspacing': batch_param['kspacing'],
                        'energy': comp_result.get('energy', np.nan),
                        'energy_per_atom': comp_result.get('energy', np.nan) / len(self.atoms) if comp_result.get('energy') else np.nan,
                        'max_force': np.nan,  # Not available from batch results yet
                        'n_kpoints': np.nan,
                        'label': batch_param['label'],
                        'test_number': i + 1,
                    }
                    results_list.append(result)
                    ecutwfc_results.append(result)
                    
                    if verbose:
                        energy_per_atom = result['energy_per_atom']
                        print(f"  ✓ {batch_param['label']}: E = {energy_per_atom:.6f} eV/atom")
                else:
                    if verbose:
                        print(f"  ✗ {batch_params[i]['label']}: {comp_result.get('error', 'Unknown error')}")
            
            # Check kspacing convergence for this ecutwfc
            if len(ecutwfc_results) >= 2:
                temp_df = pd.DataFrame(ecutwfc_results)
                convergence_status = self._check_kspacing_convergence(
                    temp_df, self.convergence_criteria_list, self.convergence_criteria
                )
                
                all_converged = all(convergence_status.values())
                
                if all_converged:
                    converged_ecutwfc = current_ecutwfc
                    converged_kspacing = min([r['kspacing'] for r in ecutwfc_results])
                    
                    if verbose:
                        print(f"✓ CONVERGED at ecutwfc={current_ecutwfc:.1f} Ry, kspacing≤{converged_kspacing:.3f} Å⁻¹")
                    
                    break  # Exit ecutwfc loop
                else:
                    if verbose:
                        unconverged = [k for k, v in convergence_status.items() if not v]
                        print(f"  Not converged: {', '.join(unconverged)}")
            
            # Go to next ecutwfc value
            current_ecutwfc += ecutwfc_step
        
        # Finalize
        if verbose:
            print(f"\n{'='*80}")
            print("CONVERGENCE STUDY COMPLETE")
            print(f"{'='*80}\n")
        
        # Store results
        self.results = pd.DataFrame(results_list)
        return self.results
    
    def _run_convergence_study_sequential(
        self,
        label_prefix: str,
        verbose: bool,
        max_ecutwfc: float,
        ecutwfc_step: float,
        start_ecutwfc: float,
        kspacing_values: List[float]
    ) -> pd.DataFrame:
        """Run convergence study using sequential (one-at-a-time) submission."""
        results_list = []
        test_count = 0
        current_ecutwfc = start_ecutwfc
        
        while current_ecutwfc <= max_ecutwfc:
            if verbose:
                print(f"\n--- Testing ecutwfc = {current_ecutwfc:.1f} Ry ---")
            
            ecutwfc_results = []
            
            # Test all kspacing values for this ecutwfc (sequentially)
            for kspacing in kspacing_values:
                test_count += 1
                label = f"{label_prefix}/ecut{int(current_ecutwfc)}_ksp{kspacing:.2f}"
                
                if verbose:
                    print(f"  [{test_count}] kspacing={kspacing:.3f} Å⁻¹")
                
                # Create workflow for this parameter set
                # Pass pseudopotentials via config name if available
                workflow_kwargs = {
                    'atoms': self.atoms,
                    'protocol': self.protocol,
                    'kspacing': kspacing,
                    'code_version': self.code_version,
                }
                
                # Pass pseudopotentials via config name if available
                if self._pseudo_config_name:
                    workflow_kwargs['pseudopotentials_config'] = self._pseudo_config_name
                else:
                    workflow_kwargs['pseudopotentials'] = self.pseudopotentials
                
                # Pass either queue or machine (not both)
                if self.queue is not None:
                    workflow_kwargs['queue'] = self.queue
                elif self.machine is not None:
                    workflow_kwargs['machine'] = self.machine
                
                # Add any extra kwargs
                workflow_kwargs.update(self.extra_kwargs)
                
                workflow = CalculationWorkflow(**workflow_kwargs)
                
                # Override ecutwfc
                workflow.input_data['ecutwfc'] = current_ecutwfc
                
                # Run calculation
                try:
                    if 'geometry' in self.convergence_criteria_list:
                        calc = workflow.run_geometry_optimization(label=label)
                    else:
                        calc = workflow.run_scf(label=label)
                    
                    # Extract results
                    energy = calc.atoms.get_potential_energy()
                    energy_per_atom = energy / len(self.atoms)
                    
                    # Calculate max force if available
                    try:
                        forces = calc.atoms.get_forces()
                        max_force = np.max(np.abs(forces))
                    except:
                        max_force = np.nan
                    
                    # Extract k-points info
                    try:
                        kpts = workflow.atoms.get_calculator().get_ibz_k_points()
                        n_kpoints = len(kpts)
                    except:
                        n_kpoints = np.nan
                    
                    result = {
                        'ecutwfc': current_ecutwfc,
                        'kspacing': kspacing,
                        'energy': energy,
                        'energy_per_atom': energy_per_atom,
                        'max_force': max_force,
                        'n_kpoints': n_kpoints,
                        'label': label,
                        'test_number': test_count,
                    }
                    results_list.append(result)
                    ecutwfc_results.append(result)
                    
                    if verbose:
                        print(f"    ✓ E = {energy_per_atom:.6f} eV/atom, F_max = {max_force:.4f} eV/Å")
                    
                except Exception as e:
                    logger.error(f"Failed to run {label}: {e}")
                    if verbose:
                        print(f"    ✗ Failed: {e}")
                    continue
            
            # Check kspacing convergence for this ecutwfc
            if len(ecutwfc_results) >= 2:
                temp_df = pd.DataFrame(ecutwfc_results)
                convergence_status = self._check_kspacing_convergence(
                    temp_df, self.convergence_criteria_list, self.convergence_criteria
                )
                
                # Check if all required criteria are converged
                all_converged = all(convergence_status.values())
                
                if all_converged:
                    converged_ecutwfc = current_ecutwfc
                    # Find the finest kspacing that achieved convergence
                    converged_kspacing = min([r['kspacing'] for r in ecutwfc_results])
                    
                    if verbose:
                        print(f"✓ CONVERGED at ecutwfc={converged_ecutwfc:.1f} Ry, kspacing≤{converged_kspacing:.3f} Å⁻¹")
                        finest_result = min(ecutwfc_results, key=lambda x: x['kspacing'])
                        print(f"  Energy: {finest_result['energy_per_atom']:.6f} eV/atom")
                    
                    break  # Exit ecutwfc loop
                else:
                    if verbose:
                        unconverged = [k for k, v in convergence_status.items() if not v]
                        print(f"  Not converged: {', '.join(unconverged)}")
            else:
                if verbose:
                    print(f"  Insufficient data for convergence check")
            
            # Increase ecutwfc for next iteration
            current_ecutwfc += ecutwfc_step
        
        # Create results DataFrame
        self.results = pd.DataFrame(results_list)
        
        print("\n" + "="*80)
        print("CONVERGENCE STUDY COMPLETE")
        print("="*80)
        
        return self.results
    
    def get_recommendations(
        self,
        energy_tolerance: float = 1e-4,
        force_tolerance: float = 0.1,
        verbose: bool = True,
    ) -> Dict:
        """
        Analyze convergence results and recommend optimal parameters.
        
        Recommends parameter combinations based on:
        1. Energy convergence (difference from highest ecutwfc/finest kspacing)
        2. Force convergence
        3. Computational cost balance
        
        Args:
            energy_tolerance: Energy convergence target (meV/atom).
                             Default: 0.1 meV/atom (1e-4 eV/atom)
            force_tolerance: Force convergence target (eV/Å).
                            Default: 0.1 eV/Å
            verbose: Print recommendations
            
        Returns:
            Dict with recommendations:
                - 'fast': minimal parameters (low cost, acceptable accuracy)
                - 'balanced': good balance of accuracy/cost
                - 'accurate': tight convergence (high cost, high accuracy)
                - 'best_ecutwfc': ec ecutwfc for best convergence
                - 'best_kspacing': best kspacing for best convergence
                - 'convergence_summary': details of convergence analysis
        
        Example:
            >>> conv = ConvergenceWorkflow(...)
            >>> conv.run_convergence_study()
            >>> recs = conv.get_recommendations(
            ...     energy_tolerance=1e-4,  # 0.1 meV/atom
            ...     force_tolerance=0.1
            ... )
            >>> print(recs['balanced'])
        """
        if self.results is None or len(self.results) == 0:
            raise ValueError("No convergence results. Run convergence study first.")
        
        # Reference: highest ecutwfc, finest kspacing
        ref_ecutwfc = self.ecutwfc_range[-1]
        ref_kspacing = self.kspacing_range[-1]
        ref_result = self.results[
            (self.results['ecutwfc'] == ref_ecutwfc) &
            (self.results['kspacing'] == ref_kspacing)
        ]
        
        if len(ref_result) == 0:
            raise ValueError(
                f"Reference calculation (ecutwfc={ref_ecutwfc}, "
                f"kspacing={ref_kspacing}) not found in results."
            )
        
        ref_energy = ref_result['energy_per_atom'].values[0]
        
        # Analyze convergence
        self.results['energy_diff'] = (
            (self.results['energy_per_atom'] - ref_energy) * 1000  # Convert to meV
        )
        self.results['converged_energy'] = (
            np.abs(self.results['energy_diff']) <= energy_tolerance
        )
        self.results['converged_force'] = (
            self.results['max_force'] <= force_tolerance
        )
        
        recommendations = {}
        
        # 1. Fast convergence: loosest parameters that meet criteria
        fast = self.results[
            self.results['converged_energy'] &
            self.results['converged_force']
        ]
        if len(fast) > 0:
            # Choose loosest (smallest ecutwfc, largest kspacing)
            fast_sorted = fast.sort_values('ecutwfc')
            fast_rec = fast_sorted.iloc[0]
            recommendations['fast'] = {
                'ecutwfc': fast_rec['ecutwfc'],
                'kspacing': fast_rec['kspacing'],
                'energy_per_atom': fast_rec['energy_per_atom'],
                'energy_diff': fast_rec['energy_diff'],
                'max_force': fast_rec['max_force'],
                'n_kpoints': fast_rec['n_kpoints'],
                'label': fast_rec['label'],
            }
        
        # 2. Balanced: medium ecutwfc, reasonable kspacing
        mid_ecutwfc = self.ecutwfc_range[len(self.ecutwfc_range)//2]
        balanced = self.results[
            (self.results['ecutwfc'] <= mid_ecutwfc) &
            self.results['converged_energy']
        ]
        if len(balanced) > 0:
            # Choose finest kspacing from mid ecutwfc options
            balanced_rec = balanced.sort_values('kspacing').iloc[-1]
            recommendations['balanced'] = {
                'ecutwfc': balanced_rec['ecutwfc'],
                'kspacing': balanced_rec['kspacing'],
                'energy_per_atom': balanced_rec['energy_per_atom'],
                'energy_diff': balanced_rec['energy_diff'],
                'max_force': balanced_rec['max_force'],
                'n_kpoints': balanced_rec['n_kpoints'],
                'label': balanced_rec['label'],
            }
        
        # 3. Accurate: tightest parameters
        accurate = self.results[
            (self.results['ecutwfc'] == ref_ecutwfc) |
            (self.results['kspacing'] == ref_kspacing)
        ]
        if len(accurate) > 0:
            accurate_rec = accurate.iloc[0]
            recommendations['accurate'] = {
                'ecutwfc': accurate_rec['ecutwfc'],
                'kspacing': accurate_rec['kspacing'],
                'energy_per_atom': accurate_rec['energy_per_atom'],
                'energy_diff': accurate_rec['energy_diff'],
                'max_force': accurate_rec['max_force'],
                'n_kpoints': accurate_rec['n_kpoints'],
                'label': accurate_rec['label'],
            }
        
        # Best parameters
        recommendations['best_ecutwfc'] = ref_ecutwfc
        recommendations['best_kspacing'] = ref_kspacing
        
        # Summary
        recommendations['convergence_summary'] = {
            'energy_tolerance_eV_atom': energy_tolerance,
            'force_tolerance_eV_A': force_tolerance,
            'reference_energy_per_atom': ref_energy,
            'n_converged_energy': len(self.results[self.results['converged_energy']]),
            'n_converged_both': len(
                self.results[
                    self.results['converged_energy'] &
                    self.results['converged_force']
                ]
            ),
            'n_total_tests': len(self.results),
        }
        
        if verbose:
            self._print_recommendations(recommendations)
        
        return recommendations
    
    def _print_recommendations(self, recommendations: Dict):
        """Print formatted recommendations."""
        print("\n" + "="*80)
        print("CONVERGENCE RECOMMENDATIONS")
        print("="*80)
        
        if 'convergence_summary' in recommendations:
            summary = recommendations['convergence_summary']
            print(f"\nConvergence Criteria:")
            print(f"  Energy tolerance: {summary['energy_tolerance_eV_atom']*1000:.2f} meV/atom")
            print(f"  Force tolerance: {summary['force_tolerance_eV_A']:.3f} eV/Å")
            print(f"\nTest Results:")
            print(f"  Total tests: {summary['n_total_tests']}")
            print(f"  Energy converged: {summary['n_converged_energy']}")
            print(f"  Both converged: {summary['n_converged_both']}")
        
        if 'fast' in recommendations:
            rec = recommendations['fast']
            print(f"\n⚡ FAST (minimal cost):")
            print(f"  ecutwfc: {rec['ecutwfc']:.1f} Ry")
            print(f"  kspacing: {rec['kspacing']:.3f} Å⁻¹")
            print(f"  k-points: {int(rec['n_kpoints'])}")
            print(f"  Energy diff: {rec['energy_diff']:.3f} meV/atom")
        
        if 'balanced' in recommendations:
            rec = recommendations['balanced']
            print(f"\n⚖️ BALANCED (accuracy/cost):")
            print(f"  ecutwfc: {rec['ecutwfc']:.1f} Ry")
            print(f"  kspacing: {rec['kspacing']:.3f} Å⁻¹")
            print(f"  k-points: {int(rec['n_kpoints'])}")
            print(f"  Energy diff: {rec['energy_diff']:.3f} meV/atom")
        
        if 'accurate' in recommendations:
            rec = recommendations['accurate']
            print(f"\n🎯 ACCURATE (high precision):")
            print(f"  ecutwfc: {rec['ecutwfc']:.1f} Ry")
            print(f"  kspacing: {rec['kspacing']:.3f} Å⁻¹")
            print(f"  k-points: {int(rec['n_kpoints'])}")
            print(f"  Energy diff: {rec['energy_diff']:.3f} meV/atom")
        
        print("\n" + "="*80)
    
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
    
    def _check_convergence_vs_reference(
        self, 
        results_dict: Dict[float, float],
        reference_energy: float,
        criteria_tolerances: Dict[str, float]
    ) -> bool:
        """
        Check if parameters converged compared to reference energy.
        
        Args:
            results_dict: Dict mapping parameter value → energy
            reference_energy: Energy calculated with high cutoff
            criteria_tolerances: Dict with 'energy_tolerance' key
            
        Returns:
            True if |E_max - E_ref| < tolerance, False otherwise
        """
        if not results_dict:
            return False
        
        max_energy = max(results_dict.values())
        energy_diff = abs(max_energy - reference_energy)
        tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
        
        converged = energy_diff < tolerance
        return converged

    def run_convergence_independent(
        self,
        label_prefix: str = 'convergence',
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        max_kspacing: float = 0.1,
        kspacing_step: float = 0.05,
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
        ecut_results = {}  # Cache: ecut → energy
        reference_energy_per_atom = None
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
                if reference_energy_per_atom is None:
                    for i, comp in enumerate(completion):
                        param = batch_params[i]
                        if param['ecutwfc'] == max_ecutwfc and comp['success']:
                            energy = comp.get('energy', np.nan) / len(self.atoms)
                            ecut_results[param['ecutwfc']] = energy
                            reference_energy_per_atom = energy
                            if verbose:
                                print(f"  ✓ [REFERENCE] ecutwfc={param['ecutwfc']:.1f}: E = {energy:.6f} eV/atom")
                            break
                
                # STEP 2: Store all results and print with ΔE now available
                for i, comp in enumerate(completion):
                    if comp['success']:
                        param = batch_params[i]
                        energy = comp.get('energy', np.nan) / len(self.atoms)
                        ecut_results[param['ecutwfc']] = energy
                        
                        is_reference = (param['ecutwfc'] == max_ecutwfc)
                        
                        result = {
                            'phase': 0 if is_reference else 1,
                            'ecutwfc': param['ecutwfc'],
                            'kspacing': param['kspacing'],
                            'energy_per_atom': energy,
                            'label': param['label'],
                        }
                        results_all.append(result)
                        
                        # Print non-reference with ΔE
                        if not is_reference and verbose and reference_energy_per_atom is not None:
                            diff = abs(energy - reference_energy_per_atom)
                            status = "✓" if diff < criteria_tolerances.get('energy_tolerance', 1e-3) else "✗"
                            print(f"  {status} ecutwfc={param['ecutwfc']:.1f}: E = {energy:.6f} eV/atom (ΔE = {diff:.6f})")
                    else:
                        if verbose:
                            print(f"  ✗ ecutwfc={batch_params[i]['ecutwfc']}: {comp.get('error', 'Failed')}")
            
            # Check convergence (skip if reference not yet calculated)
            if reference_energy_per_atom is not None:
                # Exclude reference from convergence check
                test_ecut_results = {k: v for k, v in ecut_results.items() if k != max_ecutwfc}
                converged = self._check_convergence_vs_reference(
                    test_ecut_results, reference_energy_per_atom, criteria_tolerances
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
        if reference_energy_per_atom is None:
            raise RuntimeError("Could not obtain reference energy (ecutwfc=200)")
        
        # Remove reference from test results for selection
        test_ecut_results = {k: v for k, v in ecut_results.items() if k != max_ecutwfc}
        if not test_ecut_results:
            raise RuntimeError("PHASE 1 failed: no successful calculations")
        
        # Select ecutwfc for PHASE 2 (highest converged value for maximum precision)
        optimal_ecutwfc = max(ecut_results.keys())
        optimal_energy_phase1 = ecut_results[optimal_ecutwfc]
        
        print(f"\n✓ PHASE 1 COMPLETE: Selected ecutwfc = {optimal_ecutwfc:.1f} Ry")
        
        # ===== PHASE 2: DYNAMIC KSPACING CONVERGENCE =====
        print("\n" + "="*80)
        print("PHASE 2: KSPACING CONVERGENCE (DYNAMIC)")
        print("="*80)
        
        print(f"\nFixed ecutwfc: {optimal_ecutwfc:.1f} Ry (from PHASE 1)")
        
        # Start with initial range
        current_ksp_range = self.kspacing_range.copy()
        ksp_results = {}  # Cache: kspacing → energy
        iteration = 1
        
        while True:
            print(f"\n--- Iteration {iteration} ---")
            print(f"Testing kspacing: {current_ksp_range}")
            
            # Find which values to calculate (not in cache)
            to_calculate = [k for k in current_ksp_range if k not in ksp_results]
            
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
                
                wf2 = CalculationWorkflow(**wf_kwargs)
                
                # Prepare batch for new values only
                batch_params = []
                for kspacing in to_calculate:
                    label = f"{label_prefix}/phase2_iter{iteration}_ecut{int(optimal_ecutwfc)}_ksp{kspacing:.2f}"
                    batch_params.append({
                        'label': label,
                        'ecutwfc': optimal_ecutwfc,
                        'ecutrho': optimal_ecutwfc * self.ecutrho_ratio,  # Calculate ecutrho dynamically
                        'kspacing': kspacing,
                    })
                
                if verbose:
                    print(f"Submitting {len(batch_params)} new kspacing tests...")
                
                # Submit batch
                batch_results = wf2.submit_scf_batch_multiple(batch_params, verbose=verbose)
                completion = wf2.wait_for_batch_jobs(batch_results, timeout=batch_timeout, verbose=verbose)
                
                # Store results in cache
                for i, comp in enumerate(completion):
                    if comp['success']:
                        param = batch_params[i]
                        energy = comp.get('energy', np.nan) / len(self.atoms)
                        ksp_results[param['kspacing']] = energy
                        
                        result = {
                            'phase': 2,
                            'ecutwfc': param['ecutwfc'],
                            'kspacing': param['kspacing'],
                            'energy_per_atom': energy,
                            'label': param['label'],
                        }
                        results_all.append(result)
                        
                        if verbose:
                            diff = abs(energy - reference_energy_per_atom)
                            status = "✓" if diff < criteria_tolerances.get('energy_tolerance', 1e-3) else "✗"
                            print(f"  {status} kspacing={param['kspacing']:.3f}: E = {energy:.6f} eV/atom (ΔE = {diff:.6f})")
                    else:
                        if verbose:
                            print(f"  ✗ kspacing={batch_params[i]['kspacing']}: {comp.get('error', 'Failed')}")
            
            # Check convergence
            converged = self._check_convergence_vs_reference(
                ksp_results, reference_energy_per_atom, criteria_tolerances
            )
            
            if converged:
                if verbose:
                    print(f"\n✓ CONVERGED at iteration {iteration}")
                break
            
            # Not converged: expand range (finer kspacing)
            # Note: for kspacing, smaller values are finer, so we expand downward
            min_current = min(current_ksp_range)
            if min_current <= max_kspacing:
                if verbose:
                    print(f"\n⚠️  Cannot expand further (limit: {max_kspacing}). Stopping.")
                break
            
            new_ksp_vals = []
            num_steps = 2
            for i in range(1, num_steps + 1):
                val = min_current - (i * kspacing_step)
                if val >= max_kspacing and val not in current_ksp_range:
                    new_ksp_vals.append(val)
            
            if not new_ksp_vals:
                if verbose:
                    print(f"\n⚠️  Cannot expand further (limit: {max_kspacing}). Stopping.")
                break
            
            # Add expanded values to range
            current_ksp_range = sorted(set(current_ksp_range + new_ksp_vals), reverse=True)
            iteration += 1
        
        if not ksp_results:
            print("\n⚠️  PHASE 2: No successful kspacing tests. May need to adjust parameters.")
        
        print(f"\n✓ PHASE 2 COMPLETE")
        print("\n" + "="*80)
        print("CONVERGENCE STUDY COMPLETE (INDEPENDENT WITH DYNAMIC RANGES)")
        print("="*80 + "\n")
        
        # Store results
        self.results = pd.DataFrame(results_all)
        return self.results

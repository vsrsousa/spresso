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
from xespresso import kpts_from_spacing
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.pseudopotentials.detector import parse_upf_header, get_suggested_min_ecutwfc_from_pseudos, get_ecutrho_ratio_from_pseudos
from xespresso.utils.pseudo_utils import discover_pseudopotential_directory, get_ecutrho_ratio


# Minimum ecutwfc default: 30.0 Ry
# User can override with min_ecutwfc parameter if needed
DEFAULT_MIN_ECUTWFC = 30.0
DEFAULT_MAX_ECUTWFC = 200.0
DEFAULT_INITIAL_KSPACING = 0.3  # Coarse k-mesh for Phase 1 (configurable)
DEFAULT_MIN_KSPACING = 0.1      # Minimum k-spacing limit (convergence stops here)
logger = logging.getLogger(__name__)

# ═════════════════════════════════════════════════════════════════════════════════
# PSEUDOPOTENTIAL DATABASE FOR AUTO-DETECTION OF ecutwfc_min_for_fit
# ═════════════════════════════════════════════════════════════════════════════════
#
# This database stores the minimum ecutwfc where the basis set becomes complete
# and exponential convergence begins for different pseudopotential types.
#
# Two convergence regions:
#   1. BASIS-INCOMPLETE REGION (ecut < ecut_min_for_fit):
#      - Behavior: Non-exponential, large energy jumps (100-2000 meV)
#      - Cause: Pseudopotential basis set is too small to represent density
#      - Physics: Different regime, not statistical outliers!
#      - Action: EXCLUDE from exponential fit
#   
#   2. CONVERGENCE REGION (ecut ≥ ecut_min_for_fit):
#      - Behavior: Exponential decay E(ecut) = E_inf + A*exp(-B*ecut)
#      - Cause: Sufficient basis for smooth energy convergence
#      - Physics: Standard exponential basis convergence
#      - Action: USE for exponential fit to determine E_inf, A, B
#
# Default values based on extensive testing:
#   - Norm-Conserving (NC) with moderate augmentation: 40-50 Ry
#   - Ultra-Soft (US): 25-35 Ry
#   - PAW (hardest): 65-100 Ry
#   - Extended basis NC: 50-70 Ry
#
PSEUDOPOTENTIAL_ECUT_MIN_DATABASE = {
    # Format: 'pattern_in_filename' → ecut_min_for_fit
    # Checked in order of specificity (longest match = best)
    
    # PSL (Pseudo Dojo) pseudopotentials
    'psl.1.0.0.nc-xl': 60,       # Extended basis NC
    'psl.1.0.0.n-rrkjusxc': 60,  # Extended basis NC with XC
    'psl.1.0.0.n-rrkjus': 50,    # Standard NC (your Au pseudo!)
    'psl.1.0.0.rrkjus': 45,      # Standard NC
    'psl.1.0.0.nc': 40,          # Slim basis NC
    'psl.1.0.0.us': 30,          # Ultra-Soft
    'psl.1.0.0.paw': 75,         # PAW
    
    # Older PSL versions
    'psl.0.3.1.n-rrkjus': 50,
    'psl.0.3.1.us': 30,
    'psl.0.3.1.paw': 75,
    
    # GBRV pseudopotentials
    'gbrv_bh': 35,               # GBRV Bestow Hartree
    'gbrv_bh_us': 30,            # GBRV BH Ultra-Soft
    'gbrv_bh_paw': 60,           # GBRV BH PAW
    
    # SG15 (Simple Generic) pseudopotentials
    'sg15_oncvpsp': 50,          # Standard SG15
    'sg15_oncvpsp_fr': 55,       # SG15 Full Relativistic
    
    # SSSP (Standard Solid State Pseudopotentials)
    'sssp_pbe_fr_v1.2': 55,      # SSSP efficiency variant
    'sssp_pbe_1.1': 50,          # SSSP standard
    
    # Default fallback by pattern
    'paw': 70,                   # Generic PAW
    'us': 30,                    # Generic Ultra-Soft
    'nc': 45,                    # Generic Norm-Conserving
}

# Inverse database: given a pseudopotential filename, find the minimum ecutwfc
# Matching logic: longest matching substring wins (more specific pattern first)
def _get_ecut_min_from_pseudo_filename(pseudo_filename: str) -> Optional[float]:
    """Look up ecut_min from pseudopotential filename in database."""
    if not pseudo_filename:
        return None
    
    pseudo_lower = pseudo_filename.lower()
    best_match = None
    best_score = 0
    
    for pattern, ecut_min in PSEUDOPOTENTIAL_ECUT_MIN_DATABASE.items():
        if pattern in pseudo_lower:
            score = len(pattern)  # Longer match = more specific = better
            if score > best_score:
                best_match = ecut_min
                best_score = score
    
    return best_match


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
        precision: Precision level ('low', 'medium', 'high', 'ultra') - controls convergence tolerance
        results: DataFrame with convergence test results
    """
    
    def __init__(
        self,
        atoms: Atoms,
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        protocol: str = 'moderate',
        precision: str = 'low',
        min_ecutwfc: Optional[float] = None,
        max_ecutwfc: float = DEFAULT_MAX_ECUTWFC,
        initial_kspacing: float = DEFAULT_INITIAL_KSPACING,
        ecut_vals: Optional[List[float]] = None,
        ecut_range: Optional[Dict[str, float]] = None,
        conv_thr_range: Optional[List[float]] = None,
        convergence_criteria_list: Optional[List[str]] = None,
        convergence_criteria: Optional[Dict] = None,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        hubbard_config: Optional[Union[str, Dict]] = None,
        enhance_nbands: bool = False,
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
                     Overrides ecut_range and kspacing_range.
                     Default: 'low' (fast convergence)
            min_ecutwfc: Minimum ecutwfc cutoff for convergence study (default: 30.0 Ry).
                        Convergence will start from this value and increase until convergence is reached.
                        If None, uses DEFAULT_MIN_ECUTWFC (30.0 Ry).
            max_ecutwfc: Maximum ecutwfc cutoff for convergence study (default: 200.0 Ry).
                        Upper bound for convergence range.
            ecut_vals: Optional explicit list of ecutwfc values to test (e.g., [30, 40, 50, 60, 70]).
                          If provided, these exact values will be tested in order (HIGHEST PRIORITY).
                          Mutually exclusive with ecut_range.
                          Example: ecut_vals=[30, 40, 50, 60, 70]
            ecut_range: Optional dict for range-based specification (HIGHER PRIORITY than min/max/step).
                          Supports two modes:
                          
                          LINSPACE MODE (uniform spacing):
                            {'min': 30, 'max': 70, 'n_points': 5}
                            → Tests 5 evenly-spaced values
                            → Useful for comparison studies, papers
                          
                          ARANGE MODE (step-based):
                            {'min': 30, 'max': 70, 'step': 10}
                            → Tests values with fixed step size
                            → Useful for production runs with controlled spacing
                          
                          Mutually exclusive with ecut_vals.
            initial_kspacing: Initial k-spacing value for Phase 1 (and Phase 2 starting point) in Å⁻¹.
                            Default: 0.3 Å⁻¹ (coarse k-mesh). User can customize for specific needs.
                            Phase 1: Fixed at this value while converging ecutwfc
                            Phase 2: Starts from this value and refines downward

            conv_thr_range: List of conv_thr values to test (optional)
            convergence_criteria_list: List of convergence criteria to check.
                                     Options: 'energy', 'forces', 'geometry', 'stress', 'magnetic_moments'
                                     If None, uses defaults based on precision level.
            convergence_criteria: Dict with convergence tolerances.
                                If None, uses defaults based on precision level.
                                Keys: 'energy_tolerance', 'force_tolerance', 'geometry_tolerance', 'stress_tolerance', 'magnetic_tolerance'
            queue: Queue configuration for job submission (optional)
            magnetic_config: Magnetic configuration for calculations. Can be a string like 'ferromagnetic', 
                           'antiferromagnetic', or a dict specifying magnetic moments per atom (optional).
                           Can also include Hubbard U parameters with 'U' key:
                           Example: {'Fe': {'mag': [1, -1], 'U': 4.3}}
            hubbard_config: Hubbard parameter configuration (alternative to including in magnetic_config).
                          For QE >= 7.0, can be a dict with orbital-specific U values:
                          Example: {'Fe': {'3d': 4.3, '4s': 0.0}}
                          (optional)
            enhance_nbands: If True, automatically calculates nbnd as the exact total number of 
                           valence electrons in the structure. Default False uses traditional 
                           estimation with buffer.
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
            
            #  Discover the base directory for pseudopotentials if not absolute paths
            try:
                resolved_pseudos, self.pseudopotentials_base_path = discover_pseudopotential_directory(pseudopotentials)
                # Extract FILENAMES from resolved absolute paths (same as pseudoconfig logic)
                # This matches the pseudoconfig behavior: store filenames, not full paths
                self.pseudopotentials = {}
                for element, full_path in resolved_pseudos.items():
                    filename = os.path.basename(full_path)
                    self.pseudopotentials[element] = filename
                    logger.info(f"  Discovered {element}: {filename} from {self.pseudopotentials_base_path}")
            except FileNotFoundError as e:
                raise FileNotFoundError(str(e))
            
            # Auto-detect ecutrho_ratio from pseudopotential types
            self.ecutrho_ratio = get_ecutrho_ratio_from_pseudos(
                self.pseudopotentials,
                self.pseudopotentials_base_path
            )
            if self.ecutrho_ratio == 4.0:
                logger.info(f"Auto-detected ecutrho_ratio = {self.ecutrho_ratio:.1f} (Norm-Conserving pseudos)")
            else:
                logger.info(f"Auto-detected ecutrho_ratio = {self.ecutrho_ratio:.1f} (Ultrasoft/PAW or mixed pseudos)")
        
        # Set ESPRESSO_PSEUDO environment variable if we discovered a base path
        # This ensures CalculationWorkflow can find pseudopotentials via env var
        if self.pseudopotentials_base_path:
            os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
            logger.info(f"Set ESPRESSO_PSEUDO={self.pseudopotentials_base_path}")
        
        self.protocol = protocol
        self.precision = precision
        self.initial_kspacing = initial_kspacing
        self.magnetic_config = magnetic_config
        self.hubbard_config = hubbard_config
        self.queue = queue
        self.machine = machine
        self.code_version = code_version
        self.enhance_nbands = enhance_nbands
        self.workflow_kwargs = kwargs
        
        # Auto-detect min_ecutwfc from pseudopotentials if not provided by user
        if min_ecutwfc is None:
            detected_min = get_suggested_min_ecutwfc_from_pseudos(
                self.pseudopotentials, 
                self.pseudopotentials_base_path
            )
            if detected_min is not None:
                min_ecutwfc = detected_min
                logger.info(f"Auto-detected min_ecutwfc = {detected_min:.1f} Ry from UPF headers")
        
        # Store ecutwfc range parameters (can be overridden per run)
        self.min_ecutwfc = min_ecutwfc if min_ecutwfc is not None else DEFAULT_MIN_ECUTWFC
        self.max_ecutwfc = max_ecutwfc
        self.ecut_vals = ecut_vals  # Explicit list mode
        self.ecut_range = ecut_range    # Linspace/arange mode
        
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
        self.hubbard_config = hubbard_config
        self.extra_kwargs = kwargs
        
        # Store code_version to be passed to CalculationWorkflow instances
        # It will be auto-converted to qe_version in input_data for Hubbard format detection
        if code_version is not None:
            logger.info(f"Code version set to {code_version} - will be used for Hubbard format selection")
        
        # Results storage
        self.results = None  # DataFrame will be created after tests
        
        # Persistent cache for convergence results (reused between multiple runs)
        self.ecut_results_cache = {}  # {ecutwfc: {property_name: value}}
        self.kspacing_results_cache = {}  # {kspacing: {property_name: value}}
        self.reference_ecut_result = None  # Store reference calculation result
        
        logger.info(
            f"Convergence workflow initialized:\n"
            f"  Structure: {self.atoms.get_chemical_formula()}\n"
            f"  Precision: {self.precision or 'custom'}\n"
            f"  Convergence criteria: {self.convergence_criteria_list}\n"
            f"  Algorithm: Dynamic range expansion based on convergence"
        )
    
    def _build_ecutwfc_range(self, ecut_step: float = 10.0) -> List[float]:
        """
        Build ecutwfc range based on user specification.
        
        Supports 4 modes with clear precedence:
        1. Explicit list (self.ecut_vals) - HIGHEST priority
        2. Dict-based range (self.ecut_range) - LINSPACE or ARANGE mode
        3. Step-based range (default) - Uses min/max/step parameters (ARANGE mode)
        
        Args:
            ecut_step: Step size for default ARANGE mode (default: 10.0 Ry)
        
        Returns:
            List of ecutwfc values to test (sorted)
            
        Examples:
            # Mode 1: Explicit list
            >>> wf.ecut_vals = [30, 40, 50, 60, 70]
            >>> wf._build_ecutwfc_range()
            [30, 40, 50, 60, 70]
            
            # Mode 2: Linspace (5 uniform points)
            >>> wf.ecut_range = {'min': 30, 'max': 70, 'n_points': 5}
            >>> wf._build_ecutwfc_range()
            [30.0, 40.0, 50.0, 60.0, 70.0]
            
            # Mode 3: Arange (step-based) - EXPLICIT
            >>> wf.ecut_range = {'min': 30, 'max': 70, 'step': 10}
            >>> wf._build_ecutwfc_range()
            [30, 40, 50, 60, 70]
            
            # Mode 3: Arange (step-based) - DEFAULT (NEW!)
            >>> wf.ecut_vals = None; wf.ecut_range = None
            >>> wf._build_ecutwfc_range(ecut_step=10)
            [30, 40, 50, 60, 70, 80, ...]  # Automatic step-based range up to max_ecutwfc
        """
        
        # PRIORITY 1: Explicit list (highest priority)
        if self.ecut_vals is not None:
            if not isinstance(self.ecut_vals, (list, tuple)):
                raise TypeError(f"ecut_vals must be list or tuple, got {type(self.ecut_vals)}")
            ecutwfc_list = sorted(list(self.ecut_vals))
            if len(ecutwfc_list) < 1:
                raise ValueError("ecut_vals must contain at least 1 value")
            return ecutwfc_list
        
        # PRIORITY 2: Dict-based range specification
        if self.ecut_range is not None:
            if not isinstance(self.ecut_range, dict):
                raise TypeError(f"ecut_range must be dict, got {type(self.ecut_range)}")
            
            # Check for linspace mode (n_points specified)
            if 'n_points' in self.ecut_range:
                min_val = self.ecut_range.get('min')
                max_val = self.ecut_range.get('max')
                n_points = self.ecut_range.get('n_points')
                
                if min_val is None or max_val is None or n_points is None:
                    raise ValueError(
                        "Linspace mode requires 'min', 'max', and 'n_points' keys. "
                        f"Got: {self.ecut_range}"
                    )
                
                if n_points < 2:
                    raise ValueError(f"n_points must be >= 2, got {n_points}")
                
                ecutwfc_list = np.linspace(min_val, max_val, n_points).tolist()
                return sorted(ecutwfc_list)
            
            # Check for arange mode (step specified)
            elif 'step' in self.ecut_range:
                min_val = self.ecut_range.get('min')
                max_val = self.ecut_range.get('max')
                step = self.ecut_range.get('step')
                
                if min_val is None or max_val is None or step is None:
                    raise ValueError(
                        "Arange mode requires 'min', 'max', and 'step' keys. "
                        f"Got: {self.ecut_range}"
                    )
                
                if step <= 0:
                    raise ValueError(f"step must be positive, got {step}")
                
                # Use arange: include max value if it's a multiple of step
                ecutwfc_list = np.arange(min_val, max_val + step/2, step).tolist()
                return sorted(ecutwfc_list)
            
            else:
                raise ValueError(
                    "ecut_range dict must have either 'n_points' (linspace mode) "
                    "or 'step' (arange mode). "
                    f"Got keys: {list(self.ecut_range.keys())}"
                )
        
        # PRIORITY 3: Default ARANGE mode - automatic step-based range (NEW DEFAULT!)
        # Use arange from min to max with specified step
        ecutwfc_list = np.arange(self.min_ecutwfc, self.max_ecutwfc + ecut_step/2, ecut_step).tolist()
        return sorted(ecutwfc_list)

    
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
            **kwargs: Additional parameters for __init__ (ecut_range, kspacing_range, etc.)
            
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
        ecut_min: Optional[float] = None,
        ecut_max: float = 200.0,
        ecut_step: float = 10.0,
        ecut_vals: Optional[List[float]] = None,
        ecut_range: Optional[Dict[str, float]] = None,
        n_ecut: Optional[int] = None,
        phases: str = 'both',
        magnetic_config: Optional[Union[str, Dict]] = None,
        hubbard_config: Optional[Union[str, Dict]] = None,
        **kwargs
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
            ecut_min: Minimum ecutwfc to test (default: None, uses DEFAULT_MIN_ECUTWFC=30.0)
            ecut_max: Maximum ecutwfc to test (default: 200.0) - ignored if ecut_range provided
            ecut_step: Step size for ecutwfc increases (default: 10.0) - ignored if ecut_range provided
            ecut_vals: Optional explicit list of ecutwfc values to test (e.g., [30, 40, 50, 60, 70])
                          Highest priority - overrides other range specifications
            ecut_range: Optional dict for range-based specification.
                          Supports two modes:
                          - {'min': 30, 'max': 150, 'n_points': 5} - Linspace (uniform spacing)
                          - {'min': 30, 'max': 150, 'step': 10} - Arange (step-based)
            n_ecut: Shorthand for linspace mode - equivalent to 
                            ecut_range={'min': min_ecutwfc, 'max': max_ecutwfc, 'n_points': n_ecut}
            phases: Which phases to run ('ecut', 'kpt', or 'both', default: 'both')
            magnetic_config: Magnetic configuration string or dict (optional)
            hubbard_config: Hubbard parameter configuration (optional).
                          Can be dict with U values or string for new/old format selection.
            **kwargs: Additional parameters passed to CalculationWorkflow (e.g., nbnd, conv_thr)
            
        Returns:
            ConvergenceWorkflow instance with completed convergence study
        """
        # If a pseudopotentials_config name is provided, pass it to the ctor
        if pseudopotentials_config is not None:
            workflow = cls(
                atoms,
                pseudopotentials_config=pseudopotentials_config,
                precision=precision,
                queue=queue,
                machine=machine,
                code_version=code_version,
                convergence_criteria_list=convergence_criteria_list,
                magnetic_config=magnetic_config,
                hubbard_config=hubbard_config,
                min_ecutwfc=ecut_min,
                max_ecutwfc=ecut_max,
                ecut_vals=ecut_vals,
                ecut_range=ecut_range,
                **kwargs
            )
        else:
            workflow = cls(
                atoms,
                pseudopotentials,
                precision=precision,
                queue=queue,
                machine=machine,
                code_version=code_version,
                convergence_criteria_list=convergence_criteria_list,
                magnetic_config=magnetic_config,
                hubbard_config=hubbard_config,
                min_ecutwfc=ecut_min,
                max_ecutwfc=ecut_max,
                ecut_vals=ecut_vals,
                ecut_range=ecut_range,
                **kwargs
            )
        
        # Run convergence study with specified parameters
        workflow.run_convergence_study(
            label_prefix=label_prefix,
            verbose=verbose,
            ecut_min=ecut_min,
            ecut_max=ecut_max,
            ecut_step=ecut_step,
            ecut_vals=ecut_vals,
            ecut_range=ecut_range,
            n_ecut=n_ecut,
            batch_timeout=batch_timeout,
            phases=phases,
        )
        return workflow
    
    def run_convergence_study(
        self,
        label_prefix: str = 'convergence',
        verbose: bool = True,
        ecut_min: Optional[float] = None,
        ecut_max: float = 200.0,
        ecut_step: float = 10.0,
        ecut_vals: Optional[List[float]] = None,
        ecut_range: Optional[Dict[str, float]] = None,
        n_ecut: Optional[int] = None,
        batch_timeout: int = 3600,
        phases: str = 'both',
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
            ecut_min: Minimum ecutwfc to test (ignored if ecut_range provided, default: uses __init__ value)
            ecut_max: Maximum ecutwfc to test (ignored if ecut_range provided, default: 200.0)
            ecut_step: Step size for ecutwfc increases (ignored if ecut_range provided, default: 10.0)
            ecut_vals: Optional explicit list of ecutwfc values to test (highest priority)
            ecut_range: Optional dict for range specification ({'min': ..., 'max': ..., 'step': ...} or {'min': ..., 'max': ..., 'n_points': ...})
            n_ecut: Shorthand for linspace mode - creates ecut_range={'min': ..., 'max': ..., 'n_points': n_ecut}
            batch_timeout: Timeout for batch jobs in seconds (default: 3600)
            phases: Which phases to run ('ecut', 'kpt', or 'both', default: 'both')
            
        Returns:
            pandas.DataFrame with convergence results
        """
        return self.run_convergence(
            label_prefix=label_prefix,
            ecut_min=ecut_min,
            ecut_max=ecut_max,
            ecut_step=ecut_step,
            ecut_vals=ecut_vals,
            ecut_range=ecut_range,
            n_ecut=n_ecut,
            verbose=verbose,
            batch_timeout=batch_timeout,
            phases=phases,
        )
    
    def get_recommendations(self, verbose: bool = True) -> Dict:
        """
        Analyze convergence results and recommend optimal parameters.
        
        Uses the convergence criteria and tolerances from the convergence study,
        including information from exponential fit analysis in Phase 1.
        
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
        
        # Add exponential fit information if available
        fit_result = getattr(self, 'phase1_fit_result', None)
        if fit_result and fit_result.get('success'):
            recommendations['exponential_fit'] = {
                'E_inf': fit_result['E_inf'],  # Extrapolated asymptotic energy
                'A': fit_result['A'],           # Exponential amplitude
                'B': fit_result['B'],           # Decay constant (Ry⁻¹)
                'R_squared': fit_result['R_squared'],  # Goodness of fit
                'min_ecutwfc_for_tolerance': fit_result['min_ecutwfc_for_tolerance'],  # Extrapolated ecut for tolerance
                'method': 'exponential_decay'
            }
        
        if verbose:
            print("\n" + "="*80)
            print("CONVERGENCE RECOMMENDATIONS")
            print("="*80)
            print(f"Precision level: {self.precision}")
            print(f"Energy tolerance: {energy_tol_meV:.2f} meV/atom")
            print(f"\nOptimal ecutwfc: {optimal_ecutwfc} Ry (TESTED)")
            print(f"Optimal kspacing: {optimal_kspacing} Å⁻¹ (TESTED)")
            
            # Show exponential fit extrapolation if available
            if 'exponential_fit' in recommendations:
                fit = recommendations['exponential_fit']
                print(f"\n📊 Exponential Fit Analysis (Phase 1):")
                print(f"  Asymptotic energy E_inf = {fit['E_inf']:.8f} eV")
                print(f"  R² = {fit['R_squared']:.6f}")
                print(f"  Estimated ecutwfc for ΔE < {energy_tol_meV:.2f} meV: {fit['min_ecutwfc_for_tolerance']:.1f} Ry")
                
                # Comparison
                if fit['min_ecutwfc_for_tolerance'] > optimal_ecutwfc:
                    margin = ((fit['min_ecutwfc_for_tolerance'] - optimal_ecutwfc) / optimal_ecutwfc * 100)
                    print(f"  \n  ✓ Tested ecutwfc {optimal_ecutwfc} Ry is {margin:.1f}% BELOW estimated value")
                    print(f"    → Provides safety margin for numerical stability")
                else:
                    margin = ((optimal_ecutwfc - fit['min_ecutwfc_for_tolerance']) / fit['min_ecutwfc_for_tolerance'] * 100)
                    print(f"  \n  ⚠ Tested ecutwfc {optimal_ecutwfc} Ry is {margin:.1f}% ABOVE estimated value")
                    print(f"    → Could potentially use lower ecutwfc and still meet tolerance")
            
            print("="*80 + "\n")
        
        return recommendations
    
    def estimate_ecutwfc_for_tolerance(self, tolerance_meV: float) -> Optional[float]:
        """
        Estimate ecutwfc needed for a specific tolerance using Phase 1 exponential fit.
        
        ⭐ KEY BENEFIT: No new calculations needed! Reuses fit from Phase 1.
        
        This method allows you to estimate the ecutwfc for ANY tolerance value
        without refitting or retesting. Once you have a fit from Phase 1, you
        can extrapolate for any precision level.
        
        Args:
            tolerance_meV: Tolerance in meV/atom (e.g., 1.0 for 1 meV)
            
        Returns:
            Estimated ecutwfc for achieving this tolerance, or None if no fit available
            
        Raises:
            ValueError: If convergence study not run yet or fit failed
            
        Example:
            >>> wf = ConvergenceWorkflow(atoms, pseudopotentials, precision='low')
            >>> wf.run_convergence_study()  # ← Fit is created here
            >>> 
            >>> # Now estimate ecutwfc for different tolerances without new calculations!
            >>> ecut_1meV = wf.estimate_ecutwfc_for_tolerance(1.0)    # low
            >>> ecut_0p5meV = wf.estimate_ecutwfc_for_tolerance(0.5)  # medium
            >>> ecut_0p1meV = wf.estimate_ecutwfc_for_tolerance(0.1)  # high
            >>> ecut_0p01meV = wf.estimate_ecutwfc_for_tolerance(0.01) # ultra
        """
        if self.results is None or len(self.results) == 0:
            raise ValueError("No convergence results. Run convergence study first.")
        
        fit_result = getattr(self, 'phase1_fit_result', None)
        if not fit_result or not fit_result.get('success'):
            return None
        
        # Extract fit parameters
        A = fit_result['A']
        B = fit_result['B']
        
        # Convert tolerance from meV to eV
        tolerance_eV = tolerance_meV / 1000.0
        
        # Solve: |A * exp(-B * ecut)| = tolerance
        # exp(-B * ecut) = tolerance / |A|
        # -B * ecut = ln(tolerance / |A|)
        # ecut = -ln(tolerance / |A|) / B
        
        if A == 0 or B <= 0:
            return None
        
        if tolerance_eV >= abs(A):
            # Tolerance is larger than amplitude, convergence at ecut=0
            return 0.0
        
        try:
            ecut_for_tolerance = -np.log(tolerance_eV / abs(A)) / B
            return ecut_for_tolerance
        except (ValueError, ZeroDivisionError):
            return None
    
    def recommend_for_multiple_precisions(self, verbose: bool = True) -> Dict:
        """
        Estimate ecutwfc recommendations for ALL precision levels using Phase 1 fit.
        
        ⭐ KEY BENEFIT: Uses exponential fit from Phase 1 to estimate ecutwfc for
        'low', 'medium', 'high', 'ultra' without any additional calculations!
        
        This is powerful because:
        1. Run Phase 1 once with any precision (e.g., precision='low')
        2. Get exponential fit with E_inf, A, B
        3. Use fit to extrapolate ecutwfc for all other precisions
        4. Zero additional computational cost!
        
        Args:
            verbose: Print recommendations for each precision
            
        Returns:
            Dict mapping precision → ecutwfc recommendation
            
        Example:
            >>> wf = ConvergenceWorkflow(atoms, pseudopotentials, precision='low')
            >>> wf.run_convergence_study(ecut_max=100, ecut_step=10)
            >>>
            >>> # Single Phase 1 run, but get recommendations for all precisions!
            >>> multi_rec = wf.recommend_for_multiple_precisions(verbose=True)
            >>> 
            >>> for precision, ecut in multi_rec.items():
            ...     print(f"{precision:8s}: ecutwfc = {ecut:6.1f} Ry")
            >>>
            >>> # Output:
            >>> # low     : ecutwfc =   48.3 Ry
            >>> # medium  : ecutwfc =   76.5 Ry
            >>> # high    : ecutwfc =  115.2 Ry
            >>> # ultra   : ecutwfc =  192.7 Ry
        """
        if self.results is None or len(self.results) == 0:
            raise ValueError("No convergence results. Run convergence study first.")
        
        fit_result = getattr(self, 'phase1_fit_result', None)
        if not fit_result or not fit_result.get('success'):
            if verbose:
                print("\n⚠ No exponential fit available. Cannot extrapolate for multiple precisions.")
            return {}
        
        # Precision level → tolerance mapping (meV/atom)
        precision_tolerances = {
            'low': 1.0,
            'medium': 0.5,
            'high': 0.1,
            'ultra': 0.01,
        }
        
        recommendations = {}
        
        if verbose:
            print("\n" + "="*80)
            print("ECUTWFC ESTIMATES FOR MULTIPLE PRECISION LEVELS")
            print("(Using exponential fit - NO additional calculations needed!)")
            print("="*80)
            print(f"\n{'Precision':<12} {'Tolerance':<15} {'Estimated ecutwfc':<20} {'Safety factor':<15}")
            print("-"*80)
        
        current_precision = getattr(self, 'precision', 'low')
        current_ecut = getattr(self, 'optimal_ecutwfc', None)
        
        for precision, tolerance_meV in precision_tolerances.items():
            ecut_estimate = self.estimate_ecutwfc_for_tolerance(tolerance_meV)
            
            if ecut_estimate is not None:
                recommendations[precision] = ecut_estimate
                
                # Calculate safety factor (tested vs extrapolated)
                if current_ecut is not None and precision == current_precision:
                    safety_factor = current_ecut / ecut_estimate if ecut_estimate > 0 else 0
                else:
                    safety_factor = None
                
                if verbose:
                    ecut_str = f"{ecut_estimate:.1f} Ry"
                    if safety_factor is not None:
                        print(f"{precision:<12} {tolerance_meV:<15.2f} {ecut_str:<20} {safety_factor:<15.2f}x")
                    else:
                        print(f"{precision:<12} {tolerance_meV:<15.2f} {ecut_str:<20} {'(extrapolated)':<15}")
            else:
                recommendations[precision] = None
                if verbose:
                    print(f"{precision:<12} {tolerance_meV:<15.2f} {'Failed':<20}")
        
        if verbose:
            print("-"*80)
            if fit_result.get('success'):
                print(f"\nFit quality (R²): {fit_result.get('R_squared', 'N/A'):.6f}")
                if fit_result.get('R_squared', 0) > 0.99:
                    print("✅ Excellent fit - extrapolations are highly reliable")
                elif fit_result.get('R_squared', 0) > 0.95:
                    print("✓ Good fit - extrapolations are reasonably reliable")
                else:
                    print("⚠ Fair fit - use extrapolations with caution")
            print("\n💡 These are ESTIMATES based on Phase 1 exponential fit.")
            print("   For critical applications, consider validating with explicit testing.")
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
        for ecut in self.ecut_range:
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
    
    def _get_kpts_for_spacing(self, kspacing: float) -> Tuple[int, int, int]:
        """
        Calculate k-point mesh for a given k-spacing value.
        
        Uses ASE's kspacing_to_grid with automatic 2π normalization
        (same as kpts_from_spacing utility function).
        
        Args:
            kspacing: K-spacing in Angstrom^-1
            
        Returns:
            Tuple of (nk_x, nk_y, nk_z) k-point grid
        """
        return kpts_from_spacing(self.atoms, kspacing)
    
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
        valid_criteria = {'energy', 'forces', 'stress'}  # energy, forces, stress are implemented
        unsupported = set(convergence_criteria_list) - valid_criteria
        
        if unsupported:
            unsupported_str = ', '.join(sorted(unsupported))
            raise NotImplementedError(
                f"Convergence criteria not yet implemented: {unsupported_str}\n"
                f"Currently supported: 'energy', 'forces', 'stress'\n"
                f"Coming soon: 'geometry' (VC-RELAX), 'magnetic_moments' (nspin=2)"
            )
        
        # Add QE input flags for forces if needed
        if 'forces' in convergence_criteria_list:
            config['input_data_overrides']['tprnfor'] = True
        
        # Add QE input flags for stress if needed
        if 'stress' in convergence_criteria_list:
            config['input_data_overrides']['tstress'] = True
        
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
            # Extract hydrostatic pressure from stress tensor
            if 'stress' in completion:
                stress_array = np.array(completion['stress'])  # Shape: (3, 3) in kBar
                # Compute hydrostatic pressure: P = -(trace(σ) / 3)
                trace = np.trace(stress_array)
                hydrostatic_pressure = -(trace / 3.0)  # in kBar, negative for pressure
                # Convert kBar to GPa: 1 kBar = 0.1 GPa
                hydrostatic_pressure_gpa = hydrostatic_pressure * 0.1
                # Return absolute value of pressure in GPa
                return float(abs(hydrostatic_pressure_gpa))
            else:
                return np.nan
        
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
    
    def _detect_ecut_min_from_energy_delta(
        self,
        ecutwfc_vals: np.ndarray,
        energy_vals: np.ndarray,
        delta_threshold: float = 0.050,
        verbose: bool = False
    ) -> float:
        """
        Auto-detect ecut_min_for_fit using energy delta from E_inf.
        
        Strategy:
            1. Quick rough fit with ALL data to get initial E_inf
            2. Calculate ΔE = |E - E_inf| for each point
            3. Find FIRST point where ΔE < delta_threshold (50 meV fixed)
            4. Return that ecutwfc as ecut_min_for_fit
        
        Physics: The basis-incomplete region has large energy jumps (ΔE >> 50 meV).
                 The convergence region has small deviations (ΔE < 50 meV).
                 This separates basis-incomplete from exponential convergence regions.
        
        Args:
            ecutwfc_vals: Array of ecutwfc values (Ry)
            energy_vals: Array of corresponding energies (eV)
            delta_threshold: Energy delta threshold in eV (default 0.050 = 50 meV - FIXED)
            verbose: Print detection steps
            
        Returns:
            ecut_min value in Ry
        """
        from scipy.optimize import curve_fit
        
        ecutwfc_vals = np.array(ecutwfc_vals)
        energy_vals = np.array(energy_vals)
        
        # Sort by ecutwfc
        sort_idx = np.argsort(ecutwfc_vals)
        ecutwfc_vals = ecutwfc_vals[sort_idx]
        energy_vals = energy_vals[sort_idx]
        
        # STEP 1: Quick rough fit with ALL data
        def exponential_decay(x, E_inf, A, B):
            return E_inf + A * np.exp(-B * x)
        
        try:
            E_inf_guess = energy_vals[-1]
            A_guess = energy_vals[0] - E_inf_guess
            B_guess = 0.05
            
            popt, _ = curve_fit(
                exponential_decay,
                ecutwfc_vals,
                energy_vals,
                p0=[E_inf_guess, A_guess, B_guess],
                maxfev=10000
            )
            
            E_inf = popt[0]
            
            # STEP 2: Calculate ΔE for each point
            delta_E = np.abs(energy_vals - E_inf)
            
            # STEP 3: Find first point where ΔE < threshold (50 meV separates regions)
            candidates = ecutwfc_vals[delta_E < delta_threshold]
            
            if len(candidates) > 0:
                ecut_min = candidates[0]
            else:
                # Fallback: use point closest to threshold
                closest_idx = np.argmin(np.abs(delta_E - delta_threshold))
                ecut_min = ecutwfc_vals[closest_idx]
            
            if verbose:
                print(f"\n⚡ Energy Delta Method (region separation):")
                print(f"  Quick E_inf estimate: {E_inf:.8f} eV")
                print(f"  Delta threshold: {delta_threshold*1000:.1f} meV (separates basis-incomplete from convergence)")
                print(f"  First point with ΔE < {delta_threshold*1000:.1f} meV: ecut = {ecut_min:.1f} Ry")
                print(f"  (Data points and their ΔE):")
                for ecut, E, dE in zip(ecutwfc_vals, energy_vals, delta_E):
                    marker = " ← convergence region starts" if ecut >= ecut_min else ""
                    print(f"    ecutwfc={ecut:6.1f} Ry: ΔE = {dE*1000:8.2f} meV{marker}")
            
            return float(ecut_min)
            
        except Exception as e:
            if verbose:
                print(f"  ✗ Energy delta detection failed: {e}")
            return ecutwfc_vals[-1]  # Fallback to highest ecutwfc
    
    def _detect_ecut_min_from_pseudopotential_database(self, verbose: bool = False) -> Optional[float]:
        """
        Detect ecut_min_for_fit from pseudopotential database.
        
        Looks up standard values for known pseudopotential families.
        Returns None if pseudopotential not recognized (will use curvature method as fallback).
        
        Args:
            verbose: Print detected value
            
        Returns:
            ecut_min value in Ry, or None if not recognized
        """
        # Get first pseudopotential filename from our pseudopotentials dict
        if not self.pseudopotentials:
            return None
        
        # Pick any pseudo (for database lookup, they should all be from same source)
        pseudo_filename = list(self.pseudopotentials.values())[0]
        ecut_min = _get_ecut_min_from_pseudo_filename(pseudo_filename)
        
        if ecut_min is not None and verbose:
            print(f"\n📚 Pseudopotential Database Method:")
            print(f"  Recognized: {pseudo_filename}")
            print(f"  Detected ecut_min: {ecut_min:.0f} Ry")
        
        return ecut_min
    
    def _detect_ecut_min_from_curvature(
        self,
        ecutwfc_vals: np.ndarray,
        energy_vals: np.ndarray,
        verbose: bool = False
    ) -> float:
        """
        Auto-detect ecut_min_for_fit using second derivative (curvature) analysis.
        
        The basis-incomplete region has large curvature (rapidly changing energy slope).
        The convergence region has low curvature (smooth exponential behavior).
        This method finds the transition point (inflection point) automatically.
        
        Mathematical approach:
            1. Fit smooth spline to E(ecutwfc)
            2. Compute second derivative d²E/d(ecutwfc)²
            3. Find where curvature drops below mean (transition point)
            4. Return that ecutwfc as ecut_min_for_fit
        
        Args:
            ecutwfc_vals: Array of ecutwfc values (Ry)
            energy_vals: Array of corresponding energies (eV)
            verbose: Print analysis
            
        Returns:
            ecut_min value in Ry
        """
        from scipy.interpolate import UnivariateSpline
        
        ecutwfc_vals = np.array(ecutwfc_vals)
        energy_vals = np.array(energy_vals)
        
        # Sort by ecutwfc
        sort_idx = np.argsort(ecutwfc_vals)
        ecutwfc_vals = ecutwfc_vals[sort_idx]
        energy_vals = energy_vals[sort_idx]
        
        # Fit smooth spline (k=3 cubic spline)
        k = min(3, len(ecutwfc_vals) - 1)
        spline = UnivariateSpline(ecutwfc_vals, energy_vals, k=k, s=None)
        
        # Compute second derivative: d²E/d(ecutwfc)²
        d2E = spline.derivative(n=2)
        d2E_vals = np.abs(d2E(ecutwfc_vals))
        
        # Find threshold: mean curvature
        curvature_threshold = np.mean(d2E_vals)
        
        # Find first point where curvature < threshold
        candidates = ecutwfc_vals[d2E_vals < curvature_threshold]
        
        if len(candidates) > 0:
            ecut_min = candidates[0]
        else:
            # Fallback: use point with minimum curvature
            ecut_min = ecutwfc_vals[np.argmin(d2E_vals)]
        
        if verbose:
            print(f"\n📊 Curvature Method:")
            print(f"  Mean curvature: {curvature_threshold:.6f}")
            print(f"  Detected ecut_min: {ecut_min:.1f} Ry")
            print(f"  (First point where d²E/d(ecut)² < mean)")
        
        return float(ecut_min)
    
    def _detect_ecut_min_combined(
        self,
        ecutwfc_vals: np.ndarray,
        energy_vals: np.ndarray,
        verbose: bool = True
    ) -> float:
        """
        Auto-detect ecut_min_for_fit using combined methods.
        
        Strategy: Use energy delta as PRIMARY method (most physically meaningful).
        Use database and curvature as validation/fallback only.
        
        1. Energy delta method (PRIMARY - physics-based, 50 meV fixed threshold)
        2. Database lookup (validation, if pseudopotential recognized)
        3. Curvature analysis (validation, robust but can over-detect)
        
        Args:
            ecutwfc_vals: Array of ecutwfc values (Ry)
            energy_vals: Array of corresponding energies (eV)
            verbose: Print detection steps
            
        Returns:
            ecut_min value in Ry
        """
        results = {}
        
        # Method 1: Energy Delta (PRIMARY - physics-based, 50 meV threshold)
        results['energy_delta'] = self._detect_ecut_min_from_energy_delta(
            ecutwfc_vals, energy_vals, delta_threshold=0.050, verbose=verbose
        )
        
        # Method 2: Database (validation)
        results['database'] = self._detect_ecut_min_from_pseudopotential_database(verbose=verbose)
        
        # Method 3: Curvature (validation - often over-detects)
        results['curvature'] = self._detect_ecut_min_from_curvature(
            ecutwfc_vals, energy_vals, verbose=verbose
        )
        
        # Use energy delta as PRIMARY; validate with others
        ecut_min = results['energy_delta']
        
        if verbose:
            print(f"\n🔄 Combined Auto-Detection:")
            print(f"  Energy delta method (PRIMARY): {results['energy_delta']:.1f} Ry")
            if results['database'] is not None:
                print(f"  Database method (validation):  {results['database']:.0f} Ry")
            else:
                print(f"  Database method (validation):  Not recognized")
            print(f"  Curvature method (validation):  {results['curvature']:.1f} Ry")
            print(f"  ─────────────────────────────")
            print(f"  Selected (energy delta): {ecut_min:.1f} Ry ✓")
        
        return float(ecut_min)

    def _fit_exponential_decay_phase1(
        self,
        ecut_results: Dict[float, Dict[str, float]],
        criteria_tolerances: Dict[str, float],
        verbose: bool = True,
        ecut_min_for_fit: Optional[float] = None,
        auto_detect_ecut_min: bool = True
    ) -> Dict:
        """
        Fit exponential decay model to Phase 1 ecutwfc convergence data.
        
        Uses exponential fit: E(ecutwfc) = E_inf + A * exp(-B * ecutwfc)
        to extrapolate to asymptotic energy and determine minimum ecutwfc
        for achieving specified tolerance.
        
        **PHYSICS-BASED REGION SEPARATION** (NEW APPROACH):
        The convergence curve has TWO DISTINCT REGIONS with different physics:
        
        1. BASIS-INCOMPLETE REGION (ecut < ecut_min_for_fit):
           - Large energy jumps (100-2000 meV)
           - Non-exponential behavior
           - Pseudopotential basis set is too small
           - MUST BE EXCLUDED from exponential fit!
        
        2. CONVERGENCE REGION (ecut ≥ ecut_min_for_fit):
           - Smooth exponential behavior
           - Small energy changes (1-100 meV)
           - Basis set is sufficient
           - USE FOR EXPONENTIAL FIT
        
        This is NOT statistical outlier removal (which is wrong!).
        It's physics-based region separation: identifying where exponential
        convergence actually begins.
        
        Args:
            ecut_results: Dict mapping ecutwfc → Dict with 'energy' key
            criteria_tolerances: Dict with 'energy_tolerance' key
            verbose: Print fit results and diagnostics
            ecut_min_for_fit: Minimum ecutwfc for exponential fit.
                             If None: auto-detect using combined method (recommended)
                             If specified: use this value
            auto_detect_ecut_min: If ecut_min_for_fit is None, auto-detect it
                                 (default: True, recommended)
            
        Returns:
            Dict with keys:
                'E_inf': Extrapolated asymptotic energy (eV)
                'A': Exponential amplitude (eV)
                'B': Exponential decay constant (Ry⁻¹)
                'R_squared': Goodness of fit
                'min_ecutwfc_for_tolerance': Minimum ecutwfc to meet tolerance (RECOMMENDATION!)
                'success': Whether fit succeeded
                'fit_data': Tuple of (ecutwfc_vals, energy_vals) used in fit
                'ecut_min_for_fit': The ecut_min used to separate convergence regions
                'basis_incomplete_points': List of ecutwfc values excluded (basis-incomplete region)
                'convergence_region_points': List of ecutwfc values used in fit (convergence region)
        """
        from scipy.optimize import curve_fit
        
        # Extract energies and ecutwfc values
        ecutwfc_vals = np.array(sorted(ecut_results.keys()))
        energy_vals = np.array([ecut_results[ecut]['energy'] for ecut in ecutwfc_vals])
        
        if len(ecutwfc_vals) < 3:
            if verbose:
                print(f"⚠ Not enough data for exponential fit (need ≥3 points, got {len(ecutwfc_vals)})")
            return {
                'success': False,
                'fit_data': (ecutwfc_vals, energy_vals)
            }
        
        try:
            # STEP 1: AUTO-DETECT ecut_min_for_fit if not specified
            if ecut_min_for_fit is None and auto_detect_ecut_min:
                ecut_min_for_fit = self._detect_ecut_min_combined(
                    ecutwfc_vals, energy_vals, verbose=verbose
                )
            
            # Use conservative default if still None
            if ecut_min_for_fit is None:
                ecut_min_for_fit = 50.0
            
            # STEP 2: SEPARATE CONVERGENCE REGIONS
            # Physics: Only use data where basis is complete (ecut >= ecut_min_for_fit)
            basis_incomplete_mask = ecutwfc_vals < ecut_min_for_fit
            convergence_region_mask = ~basis_incomplete_mask
            
            basis_incomplete_ecut = ecutwfc_vals[basis_incomplete_mask].tolist()
            convergence_ecut = ecutwfc_vals[convergence_region_mask]
            convergence_energy = energy_vals[convergence_region_mask]
            
            if len(convergence_ecut) < 3:
                if verbose:
                    print(f"⚠ Not enough points in convergence region (ecut ≥ {ecut_min_for_fit:.1f})")
                    print(f"  Available: {ecutwfc_vals.tolist()}")
                return {
                    'success': False,
                    'fit_data': (ecutwfc_vals, energy_vals),
                    'ecut_min_for_fit': ecut_min_for_fit,
                    'basis_incomplete_points': basis_incomplete_ecut,
                    'convergence_region_points': convergence_ecut.tolist()
                }
            
            # STEP 3: Fit exponential to convergence region ONLY
            def exponential_decay(x, E_inf, A, B):
                return E_inf + A * np.exp(-B * x)
            
            E_inf_guess = convergence_energy[-1]
            A_guess = convergence_energy[0] - E_inf_guess
            B_guess = 0.05
            
            popt, _ = curve_fit(
                exponential_decay,
                convergence_ecut,
                convergence_energy,
                p0=[E_inf_guess, A_guess, B_guess],
                maxfev=10000
            )
            
            E_inf, A, B = popt
            
            # STEP 4: Calculate R² and validate fit quality
            residuals = convergence_energy - exponential_decay(convergence_ecut, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((convergence_energy - np.mean(convergence_energy))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # Check if fit is valid: R² should be reasonable AND A should be significant
            tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
            if r_squared < 0.8:
                if verbose:
                    print(f"\n⚠ Poor fit quality (R² = {r_squared:.4f} < 0.8)")
                    print(f"  This may indicate:")
                    print(f"  - Convergence region is too small (already converged)")
                    print(f"  - Data too noisy")
                    print(f"  - Model (exponential) is inadequate")
                    print(f"  - Try lowering ecut_min_for_fit to include earlier data")
                    
                # Still return but flag as lower confidence
                r_squared = max(r_squared, -999)  # Prevent impossibly low values
            
            if abs(A) < tolerance:
                if verbose:
                    print(f"\n⚠ Exponential amplitude A = {A:.8f} eV is smaller than tolerance {tolerance*1000:.2f} meV")
                    print(f"  Convergence cannot be evaluated - data range too small")
                    print(f"  Try lowering ecut_min_for_fit or testing larger ecutwfc range")
            
            # STEP 5: Calculate minimum ecutwfc for tolerance
            tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
            tolerance_meV = tolerance * 1000
            
            # Only extrapolate if fit is reasonable AND amplitude is significant
            if abs(A) > tolerance and r_squared > 0.7:
                # Fit is valid - use exponential extrapolation
                if A != 0 and B > 0:
                    ecut_min_for_tol = -np.log(tolerance / abs(A)) / B if tolerance < abs(A) else convergence_ecut[-1]
                else:
                    ecut_min_for_tol = convergence_ecut[-1]
            else:
                # Fit unreliable - use highest tested ecutwfc as recommendation
                ecut_min_for_tol = convergence_ecut[-1]
                if verbose:
                    print(f"  ⚠ Using maximum tested ecutwfc ({ecut_min_for_tol:.1f} Ry) as recommendation")
                    print(f"    (Fit not reliable for extrapolation)")
            
            fit_result = {
                'success': True,
                'E_inf': E_inf,
                'A': A,
                'B': B,
                'R_squared': r_squared,
                'min_ecutwfc_for_tolerance': ecut_min_for_tol,
                'tolerance_meV': tolerance_meV,
                'fit_data': (convergence_ecut, convergence_energy),
                'ecut_min_for_fit': ecut_min_for_fit,
                'basis_incomplete_points': basis_incomplete_ecut,
                'convergence_region_points': convergence_ecut.tolist(),
                'n_points_total': len(ecutwfc_vals),
                'n_points_used': len(convergence_ecut),
                'n_points_excluded': len(basis_incomplete_ecut)
            }
            
            if verbose:
                print(f"\n📊 EXPONENTIAL FIT RESULTS (Phase 1)")
                print(f"{'-'*70}")
                print(f"Function: E(ecutwfc) = E_inf + A * exp(-B * ecutwfc)")
                
                # Show convergence region separation
                print(f"\n🔬 CONVERGENCE REGION ANALYSIS:")
                print(f"  Separation point (ecut_min_for_fit): {ecut_min_for_fit:.1f} Ry")
                print(f"  Total points collected: {len(ecutwfc_vals)}")
                
                if len(basis_incomplete_ecut) > 0:
                    print(f"\n  BASIS-INCOMPLETE REGION (EXCLUDED FROM FIT):")
                    print(f"    Points: {basis_incomplete_ecut} (ecut < {ecut_min_for_fit:.1f})")
                    print(f"    Physics: Large energy jumps, non-exponential behavior")
                    print(f"    Action: Excluded to avoid fit distortion")
                
                print(f"\n  CONVERGENCE REGION (USED FOR FIT):")
                print(f"    Points: {convergence_ecut.tolist()} (ecut ≥ {ecut_min_for_fit:.1f})")
                print(f"    Physics: Smooth exponential behavior, basis complete")
                print(f"    Points used: {len(convergence_ecut)}")
                
                print(f"\nFitted parameters (convergence region):")
                print(f"  E_inf (asymptotic energy) = {E_inf:.8f} eV")
                print(f"  A (exponential amplitude) = {A:.8f} eV")
                print(f"  B (decay constant)        = {B:.6f} Ry⁻¹")
                print(f"  R² (goodness of fit)      = {r_squared:.6f}")
                
                # Quality checks
                if r_squared < 0.7:
                    print(f"\n  ⚠ WARNING: R² = {r_squared:.4f} is LOW (< 0.7)")
                    print(f"    Fit quality is poor. Possible causes:")
                    print(f"    - Too few points in convergence region")
                    print(f"    - Data already converged (A very small)")
                    print(f"    - ecut_min_for_fit set too high")
                    
                if abs(A) < tolerance:
                    print(f"\n  ⚠ WARNING: |A| = {abs(A):.8f} eV < tolerance {tolerance*1000:.2f} meV")
                    print(f"    Amplitude smaller than required tolerance!")
                    print(f"    Cannot reliably extrapolate")
                    
                if r_squared < 0.95:
                    print(f"\n  ⚠ Warning: R² = {r_squared:.4f} is below 0.95")
                    print(f"    Extrapolations should be validated with explicit testing.")
                
                # Recommendation with confidence level
                print(f"\n✅ RECOMMENDATION:")
                if abs(A) > tolerance and r_squared > 0.7:
                    print(f"  For ΔE < {tolerance_meV:.2f} meV: ecutwfc ≥ {ecut_min_for_tol:.1f} Ry (from fit)")
                else:
                    print(f"  For ΔE < {tolerance_meV:.2f} meV: ecutwfc ≥ {ecut_min_for_tol:.1f} Ry (maximum tested, fit unreliable)")
            
            return fit_result
            
        except Exception as e:
            if verbose:
                print(f"\n✗ Exponential fit failed: {e}")
            return {
                'success': False,
                'fit_data': (ecutwfc_vals, energy_vals),
                'error': str(e)
            }

    def run_convergence(
        self,
        label_prefix: str = 'convergence',
        ecut_min: Optional[float] = None,
        ecut_max: float = 200.0,
        ecut_step: float = 10.0,
        ecut_vals: Optional[List[float]] = None,
        ecut_range: Optional[Dict[str, float]] = None,
        n_ecut: Optional[int] = None,
        min_kspacing_allowed: float = 0.1,
        kspacing_step: float = 0.03,
        verbose: bool = True,
        batch_timeout: int = 3600,
        precision: Optional[str] = None,
        convergence_criteria_list_override: Optional[List[str]] = None,
        phases: str = 'both',
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
        
        CACHE & REUSE:
        - Results are cached automatically (ecut_results_cache, kspacing_results_cache)
        - Call multiple times with different precision levels to reuse calculations
        
        Benefits:
        - Converges to TRUE reference (not false convergence)
        - Expands ranges only as needed
        - Caches results (no redundant calculations)
        - Pseudo transferred only 2-3x
        - Can test multiple precision levels on same workflow instance
        
        Args:
            ecut_min: Minimum ecutwfc to test (ignored if ecut_range provided, default: uses __init__ value)
            ecut_max: Maximum ecutwfc to test (ignored if ecut_range provided, default: 200.0)
            ecut_step: Step size for ecutwfc increases (ignored if ecut_range provided, default: 10.0)
            ecut_vals: Optional explicit list of ecutwfc values to test (highest priority)
            ecut_range: Optional dict for range specification ({'min': ..., 'max': ..., 'step': ...} or {'min': ..., 'max': ..., 'n_points': ...})
            n_ecut: Shorthand for linspace mode - creates ecut_range={'min': ..., 'max': ..., 'n_points': n_ecut}
            phases: Which phases to run. Options:
                   'ecut' - Run only PHASE 1 (ecutwfc convergence)
                   'kpt'  - Run only PHASE 2 (kspacing convergence)
                   'both' - Run both phases (default)
            precision: Optional precision level override ('low', 'medium', 'high', 'ultra').
                      If None, uses self.precision from __init__
            convergence_criteria_list_override: Optional override for convergence criteria.
                      If None, uses self.convergence_criteria_list from __init__
        
        Returns:
            pandas.DataFrame with complete convergence results
        """
        # Validate phases parameter
        if phases not in ['ecut', 'kpt', 'both']:
            raise ValueError(
                f"phases must be 'ecut', 'kpt', or 'both', got '{phases}'"
            )
        
        # If only running kpt (PHASE 2), need a pre-computed ecutwfc
        if phases == 'kpt' and not hasattr(self, 'optimal_ecutwfc'):
            raise ValueError(
                "Cannot run PHASE 2 (kpt convergence) without a pre-computed ecutwfc. "
                "Either: 1) Run PHASE 1 first (phases='ecut' or 'both'), or "
                "2) Set self.optimal_ecutwfc manually, or "
                "3) Provide an ecutwfc value explicitly."
            )
        
        # Process ecutwfc range specification (priority order):
        # 1. ecut_vals (explicit list) - highest priority
        # 2. n_ecut (shorthand for linspace)
        # 3. ecut_range (dict-based with min/max/step or min/max/n_points)
        # 4. Defaults (ecut_min, ecut_max + ecut_step)
        if ecut_vals is not None:
            self.ecut_vals = ecut_vals
            self.ecut_range = None
        elif n_ecut is not None:
            # Convert n_points to linspace dict format
            self.ecut_vals = None
            range_min = ecut_min if ecut_min is not None else self.min_ecutwfc
            self.ecut_range = {
                'min': range_min,
                'max': ecut_max,
                'n_points': n_ecut
            }
        elif ecut_range is not None:
            # Use provided range dict
            self.ecut_vals = None
            self.ecut_range = ecut_range
        else:
            # Use existing instance values or update with new ones
            self.ecut_vals = None
            if ecut_min is not None:
                self.min_ecutwfc = ecut_min
            self.max_ecutwfc = ecut_max
        
        results_all = []
        
        # Use dynamic precision/criteria if provided, otherwise use instance defaults
        if precision is not None:
            criteria_tolerances = self._get_default_convergence_criteria(precision)
            convergence_criteria_list = convergence_criteria_list_override or self._get_default_convergence_criteria_list(precision)
        else:
            criteria_tolerances = self.convergence_criteria
            convergence_criteria_list = convergence_criteria_list_override or self.convergence_criteria_list
        
        # Validate calculation configuration for requested criteria
        calc_config = self._get_calculation_config(convergence_criteria_list)
        if calc_config['calc_type'] != 'scf':
            raise NotImplementedError(
                f"Calculation type '{calc_config['calc_type']}' not yet supported. "
                f"Only 'scf' (energy convergence) is currently implemented."
            )
        
        # ===== PHASE 1: DYNAMIC ECUTWFC CONVERGENCE =====
        if phases in ['ecut', 'both']:
            print("\n" + "="*80)
            print("PHASE 1: ECUTWFC CONVERGENCE (DYNAMIC)")
            print("="*80)
        
        fixed_kspacing_phase1 = self.initial_kspacing  # Use configured initial kspacing
        print(f"\nStructure: {self.atoms.get_chemical_formula()}")
        print(f"Fixed kspacing: {fixed_kspacing_phase1:.3f} Å⁻¹")
        print(f"Reference ecutwfc: {ecut_max:.1f} Ry (calculated separately)\n")
        
        # Initialize range: use step-based ARANGE mode by default
        # Reference (ecut_max) is calculated in first iteration for comparison, not in range
        min_ecutwfc = self.min_ecutwfc  # Use auto-detected or user-provided value
        ecut_range_list = self._build_ecutwfc_range(ecut_step=ecut_step)  # Build range using specified mode
        
        # Log which mode was used
        if self.ecut_vals is not None:
            if verbose:
                print(f"MODE 1 (Explicit list): Testing {len(ecut_range_list)} values")
        elif self.ecut_range is not None:
            mode_name = "LINSPACE" if 'n_points' in self.ecut_range else "ARANGE"
            if verbose:
                print(f"MODE 2 ({mode_name}): Testing {len(ecut_range_list)} values")
        else:
            if verbose:
                print(f"MODE 3 (Default ARANGE): Testing {len(ecut_range_list)} values with step={ecut_step} Ry")
        
        # Load cache from previous runs (if any) - REUTILIZA CALCULOS ANTERIORES
        ecut_results = self.ecut_results_cache.copy() if self.ecut_results_cache else {}
        reference_properties = self.reference_ecut_result  # Carrega referência se já foi calculada
        
        # PHASE 1: Submit all ecutwfc values in parallel
        # Since ecut_range_list is known (ARANGE mode), submit all at once
        # The batch scheduler (SLURM, PBS, etc) handles queue management
        to_calculate_all = [e for e in ecut_range_list if e not in ecut_results]
        
        # Always include reference (ecut_max) if not already calculated
        if ecut_max not in ecut_results:
            to_calculate_all.append(ecut_max)
        
        # Remove duplicates and sort
        to_calculate_all = sorted(set(to_calculate_all))
        
        if to_calculate_all:
            # Create workflow for PHASE 1
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
                wf_kwargs['pseudopotentials_base_path'] = self.pseudopotentials_base_path
            
            if self.queue is not None:
                wf_kwargs['queue'] = self.queue
            elif self.machine is not None:
                wf_kwargs['machine'] = self.machine
            
            if self.magnetic_config is not None:
                wf_kwargs['magnetic_config'] = self.magnetic_config
            
            if self.hubbard_config is not None:
                wf_kwargs['hubbard_config'] = self.hubbard_config
            
            wf_kwargs['enhance_nbands'] = self.enhance_nbands
            wf_kwargs.update(self.extra_kwargs)
            wf1 = CalculationWorkflow(**wf_kwargs)
            
            # Build all batch parameters for PHASE 1
            batch_params = []
            for ecutwfc in to_calculate_all:
                label = f"{label_prefix}/phase1/ecut{ecutwfc:.0f}_ksp{fixed_kspacing_phase1:.3f}"
                batch_params.append({
                    'label': label,
                    'ecutwfc': ecutwfc,
                    'ecutrho': ecutwfc * self.ecutrho_ratio,
                    'kspacing': fixed_kspacing_phase1,
                })
            
            has_reference = ecut_max in to_calculate_all
            if verbose:
                msg = f"Submitting all {len(batch_params)} ecutwfc tests in parallel"
                if has_reference:
                    msg += f" (including reference ecut={ecut_max})"
                print(f"{msg}...")
            
            # Submit ALL in parallel - scheduler handles queue
            batch_results = wf1.submit_scf_batch_multiple(batch_params, verbose=verbose)
            completion = wf1.wait_for_batch_jobs(batch_results, timeout=batch_timeout, verbose=verbose)
            
            # CRITICAL: Check if ALL calculations FAILED
            all_failed = all(not comp.get('success', False) for comp in completion)
            if all_failed:
                raise RuntimeError(
                    f"\n❌ CRITICAL ERROR: All PHASE 1 calculations FAILED\n"
                    f"   ALL {len(completion)} jobs failed\n"
                    f"   \n"
                    f"   Common causes:\n"
                    f"   1. Pseudopotential file not found or path is RELATIVE (must be ABSOLUTE)\n"
                    f"   2. Machine connection failed\n"
                    f"   3. Pseudopotential file is corrupted\n"
                    f"   \n"
                    f"   Check the error messages above for details."
                )
            
            # Build mapping from label to batch_params for proper matching
            label_to_param = {p['label']: p for p in batch_params}
            
            # Extract reference first (if not yet available)
            if reference_properties is None:
                for comp in completion:
                    if not comp['success']:
                        continue
                    # Find matching param by label
                    label = comp.get('label')
                    if not label or label not in label_to_param:
                        continue
                    param = label_to_param[label]
                    
                    if param['ecutwfc'] == ecut_max:
                        props = {}
                        for prop_name in convergence_criteria_list:
                            props[prop_name] = self._extract_property_from_result(
                                comp, len(self.atoms), prop_name
                            )
                        ecut_results[param['ecutwfc']] = props
                        reference_properties = props
                        
                        if verbose:
                            energy_str = f"{props.get('energy', np.nan):.6f}" if 'energy' in props else "N/A"
                            print(f"  ✓ [REFERENCE] ecutwfc={param['ecutwfc']:.1f}: E = {energy_str} eV/atom")
                        break
            
            # Store all successful results
            for comp in completion:
                if comp['success']:
                    # Find matching param by label (not index!)
                    label = comp.get('label')
                    if not label or label not in label_to_param:
                        if verbose:
                            print(f"  ⚠️  Cannot find matching params for label: {label}")
                        continue
                    param = label_to_param[label]
                    
                    # Skip if already stored
                    if param['ecutwfc'] in ecut_results:
                        continue
                    
                    props = {}
                    for prop_name in convergence_criteria_list:
                        props[prop_name] = self._extract_property_from_result(
                            comp, len(self.atoms), prop_name
                        )
                    ecut_results[param['ecutwfc']] = props
                    
                    is_reference = (param['ecutwfc'] == ecut_max)
                    result = {
                        'phase': 0 if is_reference else 1,
                        'ecutwfc': param['ecutwfc'],
                        'kspacing': param['kspacing'],
                        'energy_per_atom': props.get('energy', np.nan),
                        'label': param['label'],
                    }
                    results_all.append(result)
                    
                    # Print result
                    if not is_reference and verbose and reference_properties is not None:
                        energy = props.get('energy', np.nan)
                        ref_energy = reference_properties.get('energy', 0)
                        diff = abs(energy - ref_energy)
                        status = "✓" if diff < criteria_tolerances.get('energy_tolerance', 1e-3) else "✗"
                        print(f"  {status} ecutwfc={param['ecutwfc']:.1f}: E = {energy:.6f} eV/atom (ΔE = {diff:.6f})")
                else:
                    # Find matching param by label
                    label = comp.get('label')
                    if not label or label not in label_to_param:
                        if verbose:
                            print(f"  ⚠️  Cannot find matching params for failed label: {label}")
                        continue
                    param = label_to_param[label]
                    
                    # Critical: if reference failed, stop
                    if param['ecutwfc'] == ecut_max:
                        raise RuntimeError(
                            f"\n❌ CRITICAL ERROR: Reference calculation FAILED\n"
                            f"   ecutwfc={ecut_max:.1f} Ry: {comp.get('error', 'Failed')}\n"
                            f"   Cannot proceed without reference for convergence comparison.\n"
                            f"   \n"
                            f"   Possible causes:\n"
                            f"   1. ecutwfc={max_ecutwfc} is too high for this system\n"
                            f"   2. Pseudopotential issue or numerical instability\n"
                            f"   3. Machine/walltime limit reached\n"
                            f"   4. Memory/resource limit exceeded\n"
                            f"   \n"
                            f"   Try reducing max_ecutwfc and running again."
                        )
                    if verbose:
                        print(f"  ✗ ecutwfc={param['ecutwfc']}: {comp.get('error', 'Failed')}")
            
            if verbose:
                print(f"\n✓ PHASE 1: All {len(to_calculate_all)} ecutwfc calculations complete")
        
        # Ensure reference was calculated
        if reference_properties is None:
            raise RuntimeError("Could not obtain reference energy (ecutwfc=200)")
        
        # ⭐ Include ALL data in fit (including ecut_max) - fit will extract E_inf as reference
        # No need to exclude reference - the exponential fit determines the true E_inf
        test_ecut_results = ecut_results.copy()
        if not test_ecut_results:
            raise RuntimeError("PHASE 1 failed: no successful calculations")
        
        # ⭐ NEW: Fit exponential decay to extract E_inf (asymptotic energy)
        fit_result = self._fit_exponential_decay_phase1(
            test_ecut_results,
            criteria_tolerances,
            verbose=verbose
        )
        
        # Store fit results for later use in get_recommendations()
        self.phase1_fit_result = fit_result
        
        # Select ecutwfc using exponential fit (E_inf as reference)
        if fit_result.get('success'):
            # Get the extrapolated minimum ecutwfc for the required tolerance
            # Formula: A * exp(-B * ecut) = tolerance
            #          ecut = -ln(tolerance / A) / B
            ecut_min_for_tol = fit_result.get('min_ecutwfc_for_tolerance')
            
            if ecut_min_for_tol is not None:
                # ⭐ Check if this value was actually calculated or extrapolated
                if ecut_min_for_tol in ecut_results:
                    # Value was calculated - use it directly
                    optimal_ecutwfc = ecut_min_for_tol
                    converged_str = f"CONVERGED (fit-calculated for precision)"
                else:
                    # Value was extrapolated - find the smallest tested value that converges
                    E_inf = fit_result['E_inf']
                    tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
                    converged_ecutwfc = {}
                    for ecut, props_dict in ecut_results.items():
                        energy = props_dict.get('energy', np.nan)
                        delta_e = abs(energy - E_inf)
                        if delta_e < tolerance:
                            converged_ecutwfc[ecut] = delta_e
                    
                    if converged_ecutwfc:
                        optimal_ecutwfc = min(converged_ecutwfc.keys())
                        converged_str = f"CONVERGED (fit-extrapolated → tested minimum)"
                    else:
                        # Nothing converged - use closest to E_inf
                        best_ecut = min(ecut_results.keys(), 
                                       key=lambda e: abs(ecut_results[e].get('energy', np.nan) - E_inf))
                        optimal_ecutwfc = best_ecut
                        converged_str = f"NOT CONVERGED (best tested vs fit)"
            else:
                # Fallback if calculation failed
                E_inf = fit_result['E_inf']
                tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
                converged_ecutwfc = {}
                for ecut, props_dict in test_ecut_results.items():
                    energy = props_dict.get('energy', np.nan)
                    delta_e = abs(energy - E_inf)
                    if delta_e < tolerance:
                        converged_ecutwfc[ecut] = delta_e
                
                if converged_ecutwfc:
                    optimal_ecutwfc = min(converged_ecutwfc.keys())
                    converged_str = "CONVERGED (vs E_inf fit)"
                else:
                    best_ecut = min(test_ecut_results.keys(), 
                                   key=lambda e: abs(test_ecut_results[e].get('energy', np.nan) - E_inf))
                    optimal_ecutwfc = best_ecut
                    converged_str = "NOT CONVERGED (best available vs E_inf)"
        else:
            # Fallback to legacy method if fit fails: compare with max_ecutwfc reference
            if verbose:
                print(f"\n⚠ Using legacy reference method (ecutwfc={max_ecutwfc} Ry)\n")
            
            tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
            converged_ecutwfc = {}
            for ecut, props_dict in test_ecut_results.items():
                energy = props_dict.get('energy', np.nan)
                ref_energy = reference_properties.get('energy', 0)
                delta_e = abs(energy - ref_energy)
                if delta_e < tolerance:
                    converged_ecutwfc[ecut] = delta_e
            
            if converged_ecutwfc:
                # Select the MINIMUM (most efficient) converged ecutwfc
                optimal_ecutwfc = min(converged_ecutwfc.keys())
            else:
                # Fallback: if nothing converged, select the one with smallest error
                best_ecut = min(test_ecut_results.keys(), 
                               key=lambda e: abs(test_ecut_results[e].get('energy', np.nan) - 
                                               reference_properties.get('energy', 0)))
                optimal_ecutwfc = best_ecut
            
            converged_str = "CONVERGED (vs ecut=200)" if converged_ecutwfc else "NOT CONVERGED (best available)"
        
        self.optimal_ecutwfc = optimal_ecutwfc  # Store for later use in get_recommendations()
        optimal_props_phase1 = ecut_results[optimal_ecutwfc]
        
        print(f"\n✓ PHASE 1 COMPLETE: Selected ecutwfc = {optimal_ecutwfc:.1f} Ry ({converged_str})")
        
        # Persist cache for future runs with different precision levels
        self.ecut_results_cache = ecut_results.copy()
        self.reference_ecut_result = reference_properties
        print(f"  [Cache saved: {len(ecut_results)} ecutwfc values stored for reuse]")
        
        # If only running PHASE 1 (ecutwfc only), return results here
        if phases == 'ecut':
            print(f"\n" + "="*80)
            print("ECUTWFC CONVERGENCE COMPLETE (phases='ecut')")
            print("="*80 + "\n")
            self.results = pd.DataFrame(results_all)
            return self.results
        
        # ===== PHASE 2: DYNAMIC KSPACING CONVERGENCE =====
        # (only executed if phases='kpt' or 'both')
        print("\n" + "="*80)
        print("PHASE 2: KSPACING CONVERGENCE (DYNAMIC)")
        print("="*80)
        
        print(f"\nFixed ecutwfc: {optimal_ecutwfc:.1f} Ry (from PHASE 1)")
        print(f"Reference kspacing (fine): {min_kspacing_allowed:.3f} Å⁻¹\n")
        
        # Initialize range dynamically: [coarse, intermediate] 
        # Will decrease downward as needed by subtracting kspacing_step
        # Range: initial_kspacing → finer → ... → min_kspacing_allowed (fine reference)
        # Never go below min_kspacing_allowed (0.1 Å⁻¹ default)
        min_kspacing_allowed = min_kspacing_allowed  # Use parameter value, not DEFAULT_MIN_KSPACING
        kspacing_range = [self.initial_kspacing, self.initial_kspacing - kspacing_step]
        # Load cache from previous runs (if any) - REUTILIZA CALCULOS ANTERIORES
        ksp_results = self.kspacing_results_cache.copy() if self.kspacing_results_cache else {}
        
        # ⭐ SMART REUSE: Self.initial_kspacing was already calculated in Phase 1 with optimal_ecutwfc!
        # Copy Phase 1 result to Phase 2 cache to avoid recalculation
        # This value will be the INITIAL REFERENCE for sliding window
        reference_kspacing = None
        reference_properties_phase2 = None
        if optimal_ecutwfc in self.ecut_results_cache and self.initial_kspacing not in ksp_results:
            ksp_results[self.initial_kspacing] = self.ecut_results_cache[optimal_ecutwfc].copy()

            reference_kspacing = self.initial_kspacing
            reference_properties_phase2 = self.ecut_results_cache[optimal_ecutwfc].copy()
            if verbose:
                print(f"  [REUSE] Kspacing={self.initial_kspacing:.3f} Å⁻¹ from Phase 1 (optimal_ecutwfc={optimal_ecutwfc:.1f} Ry)")
                energy_str = f"{reference_properties_phase2.get('energy', np.nan):.6f}"
                print(f"  ✓ [REFERENCE] kspacing={self.initial_kspacing:.3f}: E = {energy_str} eV/atom")
        
        # Start testing from next finer (smaller) kspacing value
        ksp_to_test = self.initial_kspacing - kspacing_step
        iteration = 1
        tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
        converged_delta_e = None  # Store ΔE that caused convergence
        last_kpts = self._get_kpts_for_spacing(self.initial_kspacing)  # Track last mesh generated
        
        while ksp_to_test >= min_kspacing_allowed:
            
            print(f"\n--- Iteration {iteration} ---")
            print(f"Testing kspacing: {ksp_to_test:.3f} Å⁻¹")
            
            # FALLBACK: Check if kspacing would generate same k-mesh as previous value
            # If so, skip this value and try the next finer one automatically
            current_kpts = self._get_kpts_for_spacing(ksp_to_test)
            if current_kpts == last_kpts:
                if verbose:
                    print(f"  ⏭️  K-mesh unchanged: {current_kpts} (same as previous). Skipping to finer kspacing...")
                ksp_to_test = ksp_to_test - kspacing_step
                iteration += 1
                continue
            
            # K-mesh changed! Update for next comparison
            last_kpts = current_kpts
            
            # Prepare batch: only test values not yet in cache
            to_calculate = []
            if ksp_to_test not in ksp_results:
                to_calculate.append(ksp_to_test)
            
            # STEP 0: If value is already in cache, process it directly
            if not to_calculate and ksp_to_test in ksp_results:
                # Value is cached - process it without recalculation
                props = ksp_results[ksp_to_test]
                energy = props.get('energy', np.nan)
                ref_energy = reference_properties_phase2.get('energy', 0)
                tolerance = criteria_tolerances.get('energy_tolerance', 1e-3)
                diff = abs(energy - ref_energy)
                status = "✓" if diff < tolerance else "✗"
                
                if verbose:
                    print(f"  {status} kspacing={ksp_to_test:.3f}: E = {energy:.6f} eV/atom (ΔE = {diff:.6f}) [cached]")
                
                result = {
                    'phase': 2,
                    'ecutwfc': optimal_ecutwfc,
                    'kspacing': ksp_to_test,
                    'energy_per_atom': energy,
                    'label': f"{label_prefix}/phase2/ecut{optimal_ecutwfc:.0f}_ksp{ksp_to_test:.3f}",
                }
                results_all.append(result)
                
                # Check convergence
                if diff < tolerance:
                    if verbose:
                        print(f"\n✓ CONVERGED at kspacing={ksp_to_test:.3f} Å⁻¹")
                    break
                else:
                    # Not converged, try next finer value
                    current_ksp_index += 1
                    iteration += 1
                    continue
            
            if to_calculate:
                # Create workflow for this iteration
                wf_kwargs = {
                    'atoms': self.atoms,
                    'protocol': self.protocol,
                    'ecutwfc': optimal_ecutwfc,
                    'code_version': self.code_version,  # Will auto-convert to qe_version in CalculationWorkflow
                }
                
                if self._pseudo_config_name:
                    wf_kwargs['pseudopotentials_config'] = self._pseudo_config_name
                else:
                    wf_kwargs['pseudopotentials'] = self.pseudopotentials
                    # Pass base_path so CalculationWorkflow can set ESPRESSO_PSEUDO appropriately
                    wf_kwargs['pseudopotentials_base_path'] = self.pseudopotentials_base_path
                
                if self.queue is not None:
                    wf_kwargs['queue'] = self.queue
                elif self.machine is not None:
                    wf_kwargs['machine'] = self.machine
                
                if self.magnetic_config is not None:
                    wf_kwargs['magnetic_config'] = self.magnetic_config
                
                if self.hubbard_config is not None:
                    wf_kwargs['hubbard_config'] = self.hubbard_config
                
                # Pass extra kwargs (e.g., nbnd, conv_thr) to CalculationWorkflow
                wf_kwargs['enhance_nbands'] = self.enhance_nbands
                wf_kwargs.update(self.extra_kwargs)
                wf2 = CalculationWorkflow(**wf_kwargs)
                
                # Prepare batch for new values only
                batch_params = []
                for kspacing in to_calculate:
                    label = f"{label_prefix}/phase2/ecut{optimal_ecutwfc:.0f}_ksp{kspacing:.3f}"
                    batch_params.append({
                        'label': label,
                        'ecutwfc': optimal_ecutwfc,
                        'ecutrho': optimal_ecutwfc * self.ecutrho_ratio,
                        'kspacing': kspacing,
                    })
                
                if verbose:
                    is_first_batch = iteration == 1 and min_kspacing_allowed in to_calculate
                    msg = f"Submitting {len(batch_params)} kspacing tests"
                    if is_first_batch:
                        msg += f" (including reference kspacing={min_kspacing_allowed:.3f})"
                    print(f"{msg}...")
                
                # Submit batch
                batch_results = wf2.submit_scf_batch_multiple(batch_params, verbose=verbose)
                completion = wf2.wait_for_batch_jobs(batch_results, timeout=batch_timeout, verbose=verbose)
                
                # CRITICAL: Check if FIRST batch COMPLETELY FAILED (all jobs failed)
                if iteration == 1:
                    all_failed = all(not comp.get('success', False) for comp in completion)
                    if all_failed:
                        raise RuntimeError(
                            f"\n❌ CRITICAL ERROR: First batch of kspacing calculations FAILED (PHASE 2)\n"
                            f"   Iteration 1: ALL {len(completion)} jobs failed\n"
                            f"   \n"
                            f"   Common causes:\n"
                            f"   1. Pseudopotential file not found or path is RELATIVE (must be ABSOLUTE)\n"
                            f"   2. Machine connection failed\n"
                            f"   3. Pseudopotential file is corrupted\n"
                            f"   \n"
                            f"   Fix: Use absolute path for pseudopotential file, e.g.:\n"
                            f"   pseudopotentials={{\n"
                            f"       'Gd': '/home/vinicius/scratch/projects/spresso/pseudo/Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'\n"
                            f"   }}\n"
                            f"   \n"
                            f"   Check the error messages above for details."
                        )
                
                # STEP 1: Extract result and store
                for i, comp in enumerate(completion):
                    param = batch_params[i]
                    if comp['success']:
                        # Extract all properties for this result
                        props = {}
                        for prop_name in convergence_criteria_list:
                            props[prop_name] = self._extract_property_from_result(
                                comp, len(self.atoms), prop_name
                            )
                        ksp_results[param['kspacing']] = props
                
                # STEP 2: Check convergence for current kspacing value
                if ksp_to_test in ksp_results:
                    energy = ksp_results[ksp_to_test].get('energy', np.nan)
                    ref_energy = reference_properties_phase2.get('energy', 0)
                    diff = abs(energy - ref_energy)
                    status = "✓" if diff < tolerance else "✗"
                    
                    result = {
                        'phase': 2,
                        'ecutwfc': optimal_ecutwfc,
                        'kspacing': ksp_to_test,
                        'energy_per_atom': energy,
                        'label': f"{label_prefix}/phase2_iter{iteration}_ksp{int(ksp_to_test*1000)}",
                    }
                    results_all.append(result)
                    
                    if verbose:
                        print(f"  {status} kspacing={ksp_to_test:.3f}: E = {energy:.6f} eV/atom (ΔE = {diff*1000:.2f} meV/atom)")
                    
                    # Check if converged
                    if diff < tolerance:
                        if verbose:
                            print(f"\n✓ CONVERGED at kspacing={ksp_to_test:.3f} Å⁻¹")
                        reference_kspacing = ksp_to_test
                        converged_delta_e = diff  # Store the ΔE that caused convergence
                        break
                    else:
                        # Not converged, update reference and continue to finer kspacing
                        reference_properties_phase2 = ksp_results[ksp_to_test].copy()
                        reference_kspacing = ksp_to_test
            
            # Move to next finer (smaller) kspacing value
            ksp_to_test = ksp_to_test - kspacing_step
            iteration += 1
        
        # Ensure reference was set
        if reference_properties_phase2 is None:
            raise RuntimeError(f"Could not obtain reference energy from Phase 1")
        
        # PHASE 2: SELECT OPTIMAL KSPACING FROM CONVERGENCE RESULTS
        # If loop broke, it means reference_kspacing reached convergence
        # If loop completed naturally (while condition became false), no convergence was reached
        
        if not ksp_results:
            print("\n⚠️  PHASE 2: No successful kspacing tests. May need to adjust parameters.")
        else:
            # Check if the loop broke (converged) or completed naturally (no convergence)
            # If broke: reference_kspacing = the converged value (optimal)
            # If completed naturally: reference_kspacing = last tested value (did not converge)
            
            # The signal for break is: we exited with reference_kspacing != self.initial_kspacing
            # (since initial_kspacing is the starting reference and gets updated as we test)
            
            # Simpler approach: if reference_kspacing was updated to something finer than initial,
            # it means loop tested at least one value after setting up the reference
            
            if reference_kspacing is not None and reference_kspacing < self.initial_kspacing:
                # Loop tested values finer than initial, and reference_kspacing is the best one reached
                optimal_kspacing = reference_kspacing
                self.optimal_kspacing = optimal_kspacing
                
                if verbose:
                    print(f"\n✓ PHASE 2 CONVERGED")
                    print(f"  Optimal kspacing: {optimal_kspacing:.3f} Å⁻¹")
                    if converged_delta_e is not None:
                        print(f"  ΔE = {converged_delta_e*1000:.2f} meV/atom < tolerance = {tolerance*1000:.2f} meV/atom")
            else:
                # Loop completed without finding a converged value finer than initial
                if verbose:
                    print(f"\n✗ PHASE 2 NOT CONVERGED at any finer kspacing")
                    print(f"  Tolerance = {tolerance*1000:.2f} meV/atom (precision='{self.precision}')")
                    print(f"  Best result remains at initial kspacing: {self.initial_kspacing:.3f} Å⁻¹")
                
                # Use initial_kspacing as fallback
                self.optimal_kspacing = self.initial_kspacing
        
        # Persist kspacing cache for future runs with different precision levels
        self.kspacing_results_cache = ksp_results.copy()
        print(f"\n  [Cache saved: {len(ksp_results)} kspacing values stored for reuse]")
        
        print(f"\n" + "="*80)
        print("CONVERGENCE STUDY COMPLETE (INDEPENDENT WITH DYNAMIC RANGES)")
        print("="*80 + "\n")
        
        # Store results
        self.results = pd.DataFrame(results_all)
        return self.results
    
    def plot_phase1_fit(
        self,
        show: bool = True,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (14, 6)
    ):
        """
        Plot Phase 1 exponential fit with ALL points and CONVERGENCE REGION highlighted.
        
        Shows:
        - All tested ecutwfc points (blue dots)
        - Convergence region points used for fit (red dots - LARGER)
        - Basis-incomplete region points excluded (gray dots - smaller, crossed)
        - Fitted exponential curve (red line)
        - Energy tolerance band (green shaded area)
        - ecut_min_for_fit separation line (yellow dashed)
        
        Args:
            show: Display plot (default: True)
            save_path: Save plot to PDF (optional). Default: '<formula>_ecut_conv.pdf'
                      If None, uses compound formula automatically (e.g., 'Au_ecut_conv.pdf', 'Au2O3_ecut_conv.pdf')
            figsize: Figure size as (width, height) in inches
            
        Example:
            >>> wf = ConvergenceWorkflow(...)
            >>> wf.run_convergence_study()
            >>> wf.plot_phase1_fit()  # Saves as 'Au_ecut_conv.pdf' (or formula_ecut_conv.pdf)
            >>> wf.plot_phase1_fit(save_path='custom_name.pdf')
        """
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle
        except ImportError:
            logger.warning("matplotlib not available. Skipping plot.")
            return
        
        # Get fit results
        fit_result = getattr(self, 'phase1_fit_result', None)
        ecut_results = getattr(self, 'ecut_results_cache', {})
        
        if fit_result is None or not fit_result.get('success'):
            print("⚠ No Phase 1 fit available. Run convergence study first.")
            return
        
        if not ecut_results:
            print("⚠ No Phase 1 results. Run convergence study first.")
            return
        
        # Auto-generate filename if not provided
        if save_path is None:
            formula = self.atoms.get_chemical_formula()
            save_path = f"{formula}_ecut_conv.pdf"
        
        # Ensure PDF format
        if not save_path.lower().endswith('.pdf'):
            save_path = save_path.replace('.png', '').replace('.jpg', '') + '.pdf'
        
        # Get pseudopotential names
        pseudo_names = list(self.pseudopotentials.values()) if hasattr(self, 'pseudopotentials') else []
        pseudo_str = ', '.join(pseudo_names) if pseudo_names else 'Unknown'
        
        # Calculate recommendations for 3, 2, 1 meV
        rec_3meV = self.estimate_ecutwfc_for_tolerance(3.0)
        rec_2meV = self.estimate_ecutwfc_for_tolerance(2.0)
        rec_1meV = self.estimate_ecutwfc_for_tolerance(1.0)
        
        # Fit with ALL points (preliminary fit for comparison)
        from scipy.optimize import curve_fit
        def exponential_decay(x, E_inf, A, B):
            return E_inf + A * np.exp(-B * x)
        
        try:
            E_inf_guess = all_energies[-1]
            A_guess = all_energies[0] - E_inf_guess
            B_guess = 0.05
            
            popt_all, _ = curve_fit(
                exponential_decay,
                all_ecutwfc,
                all_energies,
                p0=[E_inf_guess, A_guess, B_guess],
                maxfev=10000
            )
            E_inf_all, A_all, B_all = popt_all
            fit_all_success = True
        except:
            fit_all_success = False
        
        # Extract fit data
        E_inf = fit_result['E_inf']
        A = fit_result['A']
        B = fit_result['B']
        R_squared = fit_result['R_squared']
        ecut_min_for_fit = fit_result['ecut_min_for_fit']
        
        basis_incomplete = fit_result.get('basis_incomplete_points', [])
        convergence_region = fit_result.get('convergence_region_points', [])
        
        # Separate all points
        all_ecutwfc = sorted(ecut_results.keys())
        all_energies = [ecut_results[ecut]['energy'] for ecut in all_ecutwfc]
        
        # Get tolerance
        tolerance = getattr(self, 'convergence_criteria', {}).get('energy_tolerance', 1e-3)
        tolerance_meV = tolerance * 1000
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        formula = self.atoms.get_chemical_formula()
        fig.suptitle('Ecutwfc Convergence', fontsize=14, fontweight='bold', y=0.98)
        
        # ============ LEFT PLOT: Convergence Region Points + Correct Fit ============
        delta_E = np.abs(np.array(all_energies) - E_inf) * 1000  # Convert to meV
        
        # Plot ONLY convergence region points
        for i, ecut in enumerate(all_ecutwfc):
            if ecut in convergence_region:
                dE = delta_E[i]
                ax1.plot(ecut, dE, 'o', color='#d62728', markersize=8, zorder=4)
        
        # Plot correct fit (convergence region only)
        ecut_range = np.linspace(min(all_ecutwfc), max(all_ecutwfc), 200)
        fitted_dE = np.abs(A) * np.exp(-B * ecut_range) * 1000  # Convert to meV
        ax1.plot(ecut_range, fitted_dE, '-', color='#d62728', linewidth=2.5, label='Fit', zorder=6)
        
        ax1.set_xlabel('ecutwfc (Ry)', fontsize=12, fontweight='bold')
        ax1.set_ylabel(r'$\Delta E = E - E_\infty$ (meV/ion)', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, linestyle=':')
        ax1.legend(loc='upper right', fontsize=10, frameon=False)
        
        # Add pseudopotential and fit quality info
        info_text = (
            f"Pseudopotentials: {pseudo_str}\n\n"
            f"Fit quality:\n"
            f"  R² = {R_squared:.6f}\n"
            f"  E∞ = {E_inf:.8f} eV/ion\n"
            f"  A = {abs(A)*1000:.2f} meV/ion\n"
            f"  B = {B:.6f} Ry⁻¹"
        )
        ax1.text(0.5, 0.98, info_text, transform=ax1.transAxes,
                fontsize=9, verticalalignment='top', horizontalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9), family='monospace')
        
        # ============ RIGHT PLOT: Convergence Region ONLY + Correct Fit ============
        # Plot only convergence region points
        for i, ecut in enumerate(all_ecutwfc):
            if ecut in convergence_region:
                dE = delta_E[i]
                ax2.plot(ecut, dE, 'o', color='#d62728', markersize=10, zorder=5)
        
        # Plot correct fit (convergence region only)
        fitted_dE = np.abs(A) * np.exp(-B * ecut_range) * 1000  # Convert to meV
        ax2.plot(ecut_range, fitted_dE, '-', color='#d62728', linewidth=2.5, label='Fit', zorder=6)
        
        # Mark ecut_min_for_fit separation line
        ax2.axvline(ecut_min_for_fit, color='#ff7f0e', linestyle='--', linewidth=2, alpha=0.7, label=rf'$E_{{\rm cut, min}}$ = {ecut_min_for_fit:.1f} Ry', zorder=2)
        
        ax2.set_xlabel('ecutwfc (Ry)', fontsize=12, fontweight='bold')
        ax2.set_ylabel(r'$\Delta E = E - E_\infty$ (meV/ion)', fontsize=12, fontweight='bold')
        ax2.set_yscale('log')
        ax2.grid(True, alpha=0.3, which='both', linestyle=':')
        ax2.legend(loc='best', fontsize=10, frameon=False)
        
        # Add recommendations to right plot
        rec_text = f"RECOMMENDATIONS:\n"
        
        # Add tolerance recommendations
        if rec_3meV is not None:
            rec_text += f"3 meV/ion: {rec_3meV:.1f} Ry\n"
        if rec_2meV is not None:
            rec_text += f"2 meV/ion: {rec_2meV:.1f} Ry\n"
        if rec_1meV is not None:
            rec_text += f"1 meV/ion: {rec_1meV:.1f} Ry"
        
        ax2.text(0.35, 0.98, rec_text, transform=ax2.transAxes,
                fontsize=9, verticalalignment='top', horizontalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9), family='monospace')
        
        plt.tight_layout()
        
        # Save as PDF
        plt.savefig(save_path, dpi=150, bbox_inches='tight', format='pdf')
        logger.info(f"Phase 1 fit plot saved to {save_path}")
        print(f"✓ Plot saved: {save_path}")
        
        if show:
            plt.show()

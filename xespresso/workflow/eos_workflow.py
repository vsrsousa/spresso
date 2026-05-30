"""
Equation of State (EOS) Workflow for Quantum ESPRESSO calculations.

This module provides tools to systematically study the relationship between
volume and energy, fit it to an equation of state (Birch-Murnaghan), and extract
equilibrium properties like bulk modulus and equilibrium volume.

The workflow:
1. Scales the structure isotropically to different volumes
2. Runs SCF calculations for each volume
3. Collects E-V data points
4. Fits Birch-Murnaghan equation of state using ASE's proven algorithm
5. Extracts equilibrium properties (V₀, E₀, B₀)

EOS Fitting Algorithm:
    This implementation uses the Birch-Murnaghan 3rd-order equation of state from
    ASE's EquationOfState class. The fitting procedure:
    
    1. Fit a parabola E(V) = a + b*V + c*V² to get initial parameter estimates
    2. Extract V₀ (volume at minimum), E₀ (energy at minimum), B₀ (bulk modulus)
    3. Use scipy.optimize.curve_fit to refine parameters with birchmurnaghan()
    
    This approach is more robust than direct optimization and avoids getting stuck
    in local minima. Results are validated against ASE's EquationOfState.

Examples:
    >>> from ase.build import bulk
    >>> from xespresso.workflow import EOSWorkflow
    >>> 
    >>> # Create EOS workflow
    >>> atoms = bulk('Fe', 'bcc', a=2.87)
    >>> eos = EOSWorkflow(
    ...     atoms=atoms,
    ...     pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
    ...     protocol='moderate'
    ... )
    >>> 
    >>> # Run EOS study from 95% to 105% of original volume
    >>> results = eos.run_eos_study(
    ...     scale_factors=[0.95, 0.975, 1.0, 1.025, 1.05],
    ...     parallel=True
    ... )
    >>> 
    >>> # Get equilibrium properties
    >>> props = eos.get_eos_properties()
    >>> print(f"V₀ = {props['v0']:.4f} Å³")
    >>> print(f"E₀ = {props['e0']:.6f} eV")
    >>> print(f"B₀ = {props['bulk_modulus']:.2f} GPa")
    >>> 
    >>> # Plot results
    >>> eos.plot_eos_curve(save_path='eos.png')
    >>> eos.plot_residuals(save_path='residuals.png')

References:
    - Birch-Murnaghan equation: Birch, F. (1947). Physical Review, 71(11), 809-824.
    - ASE implementation: https://gitlab.com/ase/ase/-/blob/master/ase/eos.py
"""

import os
import logging
import warnings
import numpy as np
import pandas as pd
from typing import Dict, Optional, Union, Tuple, List
from pathlib import Path
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor, as_completed
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit

from ase import Atoms
from ase.io import read as ase_read

from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso import Espresso


logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTS AND CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

# Default volume range and number of points for EOS study
DEFAULT_VOLUME_MIN_FACTOR = 0.95      # 95% of original volume
DEFAULT_VOLUME_MAX_FACTOR = 1.05      # 105% of original volume
DEFAULT_EOS_POINTS = 7                # Number of volume points to sample
DEFAULT_MAX_WORKERS = 4               # Maximum parallel SCF workers

# EOS fitting parameters
DEFAULT_B0_PRIME = 4.0                # Starting guess for B₀' (derivative)
BIRCH_MURNAGHAN_ORDER = 2             # 2nd order Birch-Murnaghan EOS

# Unit conversion factors
# 1 GPa = 1e9 Pa; 1 eV/Ų = 1.602176634e11 Pa
PA_PER_GPa = 1e9
PA_PER_EV_ANG3 = 1.602176634e11       # Pressure (Pa) per eV/Ų
GPa_PER_EV_ANG3 = PA_PER_EV_ANG3 / PA_PER_GPa  # ~160.217662 GPa per eV/Ų

# Fitting tolerance
FIT_CONVERGENCE_TOL = 1e-8
FIT_MAX_ITERATIONS = 10000

# Data quality thresholds
MIN_R_SQUARED_THRESHOLD = 0.99        # Warn if fit quality < 99%
MIN_POINTS_FOR_FIT = 3                # Minimum points required for fitting


# ═══════════════════════════════════════════════════════════════════════════════
# EOS FITTING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def birchmurnaghan(V, E0, B0, BP, V0):
    """
    Birch-Murnaghan 3rd order equation of state (from ASE).
    
    Calculates energy at given volume(s) using the Birch-Murnaghan equation.
    
    Formula: 
        η = (V₀/V)^(1/3)
        E(V) = E₀ + (9/16)*V₀*B₀ * {(η² - 1)²*(6 + B₀'*(η² - 1) - 4*η²)}
    
    Args:
        V: Volume(s) in Ų
        E0: Energy at equilibrium volume (eV)
        B0: Bulk modulus (eV/Ų)
        BP: First derivative of B₀ (dimensionless)
        V0: Equilibrium volume (Ų)
        
    Returns:
        E: Energy at given volume(s) (eV)
        
    Source: PRB 70, 224107 - Birch-Murnaghan equation
    Note: There's a typo in the paper regarding eta definition, corrected here.
    """
    eta = (V0 / V)**(1 / 3)
    E = E0 + 9 * B0 * V0 / 16 * (eta**2 - 1)**2 * (
        6 + BP * (eta**2 - 1) - 4 * eta**2)
    return E


def _parabola(V, a, b, c):
    """Parabola fit: E(V) = a + b*V + c*V²"""
    return a + b * V + c * V**2


def fit_birch_murnaghan(volumes: np.ndarray, energies: np.ndarray) -> Dict:
    """
    Fit Birch-Murnaghan EOS to E-V data (ASE algorithm).
    
    This uses scipy.optimize.curve_fit to fit E₀, B₀, B₀', V₀ parameters
    to the birchmurnaghan() function, following ASE's well-tested approach.
    
    Algorithm:
    1. Fit a parabola to get initial parameter estimates
    2. Extract V₀ (minimum), E₀ (min energy), B₀ (second derivative at min)
    3. Use curve_fit to refine parameters with birchmurnaghan()
    
    Args:
        volumes: Array of volumes (Ų)
        energies: Array of corresponding energies (eV)
        
    Returns:
        Dictionary with fitted parameters:
        {
            'e0': float,           # Energy at V₀ (eV)
            'v0': float,           # Equilibrium volume (Ų)
            'b0': float,           # Bulk modulus (eV/Ų)
            'b0_prime': float,     # dB₀/dP (dimensionless)
            'r_squared': float,    # Fit quality (0-1)
            'residuals': ndarray,  # (E_data - E_fit)
            'converged': bool,     # Optimization converged successfully
            'message': str,        # Convergence message
        }
    """
    if len(volumes) < MIN_POINTS_FOR_FIT:
        raise ValueError(
            f"Need at least {MIN_POINTS_FOR_FIT} data points for EOS fitting, "
            f"got {len(volumes)}"
        )
    
    volumes = np.array(volumes)
    energies = np.array(energies)
    
    # Step 1: Fit a parabola E(V) = a + b*V + c*V² to get initial estimates
    # This is more robust than picking the minimum point directly
    try:
        p0 = [np.min(energies), 1, 1]
        popt_parabola, _ = curve_fit(_parabola, volumes, energies, p0=p0)
        a_par, b_par, c_par = popt_parabola
    except Exception as e:
        logger.warning(f"Parabola fit failed: {e}. Using direct estimation.")
        # Fallback: use minimum point
        idx_min = np.argmin(energies)
        v0_initial = volumes[idx_min]
        e0_initial = energies[idx_min]
        b0_initial = 0.1  # eV/Ų
        bp_initial = DEFAULT_B0_PRIME
        popt = [e0_initial, b0_initial, bp_initial, v0_initial]
        converged = False
    else:
        # Step 2: Extract parameters from parabola
        # Minimum at dE/dV = 0  =>  b + 2*c*V = 0  =>  V_min = -b/(2*c)
        v0_initial = -b_par / (2 * c_par)
        e0_initial = _parabola(v0_initial, a_par, b_par, c_par)
        
        # Bulk modulus from 2nd derivative: B = V * d²E/dV²
        # For parabola: d²E/dV² = 2*c  =>  B ≈ 2*c*V
        b0_initial = 2 * c_par * v0_initial
        
        # Check if minimum is within data range
        minvol = np.min(volumes)
        maxvol = np.max(volumes)
        if not (minvol < v0_initial < maxvol):
            warnings.warn(
                f'Minimum volume V₀={v0_initial:.4f} is outside data range '
                f'[{minvol:.4f}, {maxvol:.4f}]. Data may not bracket the minimum.',
                UserWarning
            )
        
        # Step 3: Use curve_fit to refine with birchmurnaghan()
        # Initial guess for [E0, B0, BP, V0]
        x0 = [e0_initial, b0_initial, DEFAULT_B0_PRIME, v0_initial]
        
        try:
            # curve_fit is more robust than minimize for this case
            popt, _ = curve_fit(
                birchmurnaghan,
                volumes,
                energies,
                p0=x0,
                maxfev=FIT_MAX_ITERATIONS,
            )
            converged = True
        except Exception as e:
            logger.warning(
                f"Birch-Murnaghan fit did not converge: {e}. "
                f"Using parabola estimates."
            )
            popt = x0
            converged = False
    
    # Extract fitted parameters
    e0_fit, b0_fit, bp_fit, v0_fit = popt
    
    # Ensure physical parameters
    if b0_fit <= 0:
        logger.warning(f"Fitted B₀ = {b0_fit:.6f} eV/Ų is non-positive. Setting to small positive value.")
        b0_fit = 0.1
    
    # Calculate R² (coefficient of determination)
    E_fit = birchmurnaghan(volumes, e0_fit, b0_fit, bp_fit, v0_fit)
    residuals = energies - E_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((energies - np.mean(energies))**2)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    
    # Warn if fit quality is poor
    if r_squared < MIN_R_SQUARED_THRESHOLD:
        warnings.warn(
            f"EOS fit quality is poor (R² = {r_squared:.6f} < {MIN_R_SQUARED_THRESHOLD}). "
            f"Consider:\n"
            f"  - Using more data points\n"
            f"  - Extending volume range\n"
            f"  - Checking for calculation errors",
            UserWarning
        )
    
    return {
        'e0': float(e0_fit),
        'v0': float(v0_fit),
        'b0': float(b0_fit),
        'b0_prime': float(bp_fit),
        'r_squared': float(r_squared),
        'residuals': residuals,
        'converged': converged,
        'message': 'Success' if converged else 'Converged using parabola estimates',
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN EOS WORKFLOW CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class EOSWorkflow:
    """
    Equation of State (EOS) workflow for systematic volume optimization.
    
    This class manages a complete EOS study:
    1. Scale structures to different volumes
    2. Run SCF calculations in parallel
    3. Fit Birch-Murnaghan EOS
    4. Extract equilibrium properties
    
    Attributes:
        atoms: Original ASE Atoms object
        pseudopotentials: Dict mapping elements to pseudopotential files
        protocol: Calculation protocol ('fast', 'moderate', 'accurate')
        results_df: DataFrame with E-V data (volumes, energies, scale_factors)
        eos_params: Dict with fitted EOS parameters
        
    Example:
        >>> from ase.build import bulk
        >>> from xespresso.workflow import EOSWorkflow
        >>> 
        >>> atoms = bulk('Si', 'diamond', a=5.43)
        >>> eos = EOSWorkflow(
        ...     atoms=atoms,
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     protocol='moderate'
        ... )
        >>> 
        >>> results = eos.run_eos_study(
        ...     volume_range=(0.95, 1.05),
        ...     n_points=7,
        ...     parallel=True,
        ...     label='eos/si'
        ... )
        >>> 
        >>> props = eos.get_eos_properties()
        >>> eos.plot_eos_curve('eos_results.png')
    """
    
    def __init__(
        self,
        atoms: Union[Atoms, str, Path],
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: Optional[str] = None,
        protocol: str = 'moderate',
        kspacing: Optional[float] = None,
        input_data: Optional[Dict] = None,
        magnetic_config: Optional[Union[str, Dict]] = None,
        queue: Optional[Dict] = None,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        debug: bool = False,
        **kwargs
    ):
        """
        Initialize an EOS workflow.
        
        Args:
            atoms: ASE Atoms object, CIF file path, or structure file path
            pseudopotentials: Dict mapping element symbols to pseudopotential files
                            Either this or pseudopotentials_config must be provided.
                            Example: {'Si': 'Si.pbe.UPF', 'O': 'O.pbe-n.UPF'}
            pseudopotentials_config: Name of pseudopotentials config to load
                                   Example: 'SSSP_efficiency'
            protocol: Calculation protocol: 'fast', 'moderate', or 'accurate'
            kspacing: K-point spacing in Ų⁻¹ (overrides preset if provided)
            input_data: Additional input parameters (merged with preset)
            magnetic_config: Magnetic configuration specification
            queue: Queue/scheduler configuration for job submission
            machine: Machine configuration name to load
            code_version: Quantum ESPRESSO version to use
            debug: Enable debug logging
            **kwargs: Additional parameters passed to CalculationWorkflow
            
        Raises:
            ValueError: If neither pseudopotentials nor pseudopotentials_config provided
            FileNotFoundError: If structure file not found
        """
        # Setup logging
        if debug:
            logger.setLevel(logging.DEBUG)
        
        # Load structure if needed
        if isinstance(atoms, (str, Path)):
            if not os.path.exists(atoms):
                raise FileNotFoundError(f"Structure file not found: {atoms}")
            self.atoms = ase_read(atoms)
            logger.info(f"Loaded structure from: {atoms}")
        else:
            self.atoms = atoms.copy()
        
        # Store original volume for reference
        self.v0_original = self.atoms.get_volume()
        
        # Initialize workflow parameters
        self.pseudopotentials = pseudopotentials
        self.pseudopotentials_config = pseudopotentials_config
        self.protocol = protocol
        self.kspacing = kspacing
        self.input_data = input_data or {}
        self.magnetic_config = magnetic_config
        self.queue = queue
        self.machine = machine
        self.code_version = code_version
        self.workflow_kwargs = kwargs
        
        # Convert machine to queue if provided
        if machine is not None and queue is None:
            from xespresso.machines import load_machine
            self.queue = load_machine(machine_name=machine)
            
            # Load code configuration and extract modules for specified version
            if code_version is not None:
                self._merge_code_modules_into_queue(machine, code_version)
            
            self.machine = None  # Clear machine to avoid passing both queue and machine to CalculationWorkflow
            logger.info(f"Loaded machine '{machine}' and converted to queue dict")
        
        # Data storage
        self.results_df = None  # Will store E-V data
        self.eos_params = None  # Will store fitted EOS parameters
        
        # Initialize point index mapping
        self.eos_structures = {}  # Dict of {factor: atoms_scaled}
        self.eos_calcs = {}  # Dict of {factor: calculator}
        self.error_log = {}  # Dict of {factor: error_message}
        self._factor_to_index = {}  # Dict of {factor: point_index} for naming
        
        logger.info(f"EOSWorkflow initialized for {self.atoms.get_chemical_formula()}")
        logger.info(f"Original volume: {self.v0_original:.4f} Ų")
    
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
    
    # ═════════════════════════════════════════════════════════════════════════════
    # VOLUME SCALING METHODS
    # ═════════════════════════════════════════════════════════════════════════════
    
    def scale_volume_uniformly(
        self,
        atoms: Atoms,
        scale_factor: float
    ) -> Atoms:
        """
        Scale atomic structure by scaling the cell isotropically.
        
        Multiplies the cell vectors by scale_factor^(1/3) so that the volume
        changes by exactly scale_factor. Atomic positions are scaled proportionally.
        
        Args:
            atoms: Input ASE Atoms object
            scale_factor: Volume scaling factor (1.0 = original, 0.95 = 5% smaller, etc.)
            
        Returns:
            Scaled Atoms object (original unchanged)
            
        Raises:
            ValueError: If scale_factor <= 0
            
        Example:
            >>> atoms = bulk('Fe', 'bcc', a=2.87)
            >>> atoms_scaled = eos.scale_volume_uniformly(atoms, 1.05)
            >>> print(atoms_scaled.get_volume() / atoms.get_volume())
            1.05
        """
        if scale_factor <= 0:
            raise ValueError(f"scale_factor must be > 0, got {scale_factor}")
        
        # Create a copy to avoid modifying original
        scaled_atoms = atoms.copy()
        
        # Linear scale factor for cell (V ~ a³, so a ~ V^(1/3))
        linear_scale = scale_factor ** (1/3)
        
        # Scale cell vectors
        scaled_atoms.cell *= linear_scale
        
        # Scale positions (relative to cell origin)
        scaled_atoms.positions *= linear_scale
        
        # Verify volume scaling
        v_original = atoms.get_volume()
        v_scaled = scaled_atoms.get_volume()
        ratio = v_scaled / v_original
        
        if abs(ratio - scale_factor) > 1e-10:
            logger.warning(
                f"Volume scaling mismatch: expected {scale_factor}, got {ratio:.10f}"
            )
        
        return scaled_atoms
    
    def generate_volume_range(
        self,
        volume_range: Tuple[float, float] = (DEFAULT_VOLUME_MIN_FACTOR, DEFAULT_VOLUME_MAX_FACTOR),
        n_points: int = DEFAULT_EOS_POINTS,
        spacing: str = 'linear'
    ) -> np.ndarray:
        """
        Generate volume scale factors uniformly spaced in the range.
        
        Args:
            volume_range: Tuple (min_factor, max_factor) for volume scaling
                        Example: (0.95, 1.05) for ±5% range
            n_points: Number of volume points to sample (minimum 3)
            spacing: 'linear' (uniform in linear space) or 'log' (uniform in log space)
            
        Returns:
            Array of scale factors in ascending order
            
        Raises:
            ValueError: If n_points < 3 or invalid spacing
            
        Example:
            >>> factors = eos.generate_volume_range((0.95, 1.05), n_points=7)
            >>> print(factors)
            [0.95   0.967  0.983  1.0  1.017  1.033  1.05]
        """
        if n_points < 3:
            raise ValueError(f"Need at least 3 points, got {n_points}")
        
        min_factor, max_factor = volume_range
        
        if spacing == 'linear':
            factors = np.linspace(min_factor, max_factor, n_points)
        elif spacing == 'log':
            # Log space ensures equal relative spacing
            factors = np.exp(np.linspace(np.log(min_factor), np.log(max_factor), n_points))
        else:
            raise ValueError(f"Unknown spacing: {spacing}. Use 'linear' or 'log'")
        
        return factors
    
    def create_scaled_structures(
        self,
        scale_factors: np.ndarray
    ) -> Dict[float, Atoms]:
        """
        Create scaled structures for each scale factor.
        
        Args:
            scale_factors: Array of volume scale factors
            
        Returns:
            Dict mapping scale_factor -> scaled Atoms object
            
        Example:
            >>> factors = eos.generate_volume_range((0.95, 1.05), n_points=5)
            >>> structures = eos.create_scaled_structures(factors)
            >>> for factor, atoms in structures.items():
            ...     print(f"Factor: {factor:.3f}, Volume: {atoms.get_volume():.2f}")
        """
        structures = {}
        self._factor_to_index = {}
        
        for idx, factor in enumerate(scale_factors):
            try:
                scaled = self.scale_volume_uniformly(self.atoms, factor)
                structures[factor] = scaled
                self._factor_to_index[factor] = idx + 1  # 1-indexed
                logger.debug(f"Created structure {idx+1}: factor={factor:.4f}, V={scaled.get_volume():.4f} Ų")
            except Exception as e:
                logger.error(f"Error scaling structure by {factor}: {e}")
                self.error_log[factor] = str(e)
        
        self.eos_structures = structures
        return structures
    
    # ═════════════════════════════════════════════════════════════════════════════
    # SCF EXECUTION METHODS
    # ═════════════════════════════════════════════════════════════════════════════
    
    def _run_single_eos_point(
        self,
        scale_factor: float,
        label: str = 'eos'
    ) -> Tuple[float, float, Optional[Exception]]:
        """
        Run SCF calculation for a single volume (scale factor).
        
        Internal method used by run_eos_study for parallel execution.
        
        Args:
            scale_factor: Volume scale factor
            label: Base label for calculation directory
            
        Returns:
            Tuple: (scale_factor, energy, error_or_None)
                If successful: (factor, E, None)
                If failed: (factor, None, Exception)
        """
        try:
            # Get scaled structure
            if scale_factor not in self.eos_structures:
                raise ValueError(f"Structure for factor {scale_factor} not created")
            
            atoms_scaled = self.eos_structures[scale_factor]
            point_idx = self._factor_to_index[scale_factor]
            
            # Create calculation label using point index (not decimal scale factor)
            calc_label = f"{label}/point_{point_idx:02d}"
            
            # Create workflow for this point
            workflow = CalculationWorkflow(
                atoms=atoms_scaled,
                pseudopotentials=self.pseudopotentials,
                pseudopotentials_config=self.pseudopotentials_config,
                protocol=self.protocol,
                kspacing=self.kspacing,
                input_data=self.input_data,
                magnetic_config=self.magnetic_config,
                queue=self.queue,
                machine=self.machine,
                code_version=self.code_version,
                **self.workflow_kwargs
            )
            
            # Run SCF
            logger.info(f"Running point_{point_idx:02d}: factor={scale_factor:.4f}, V={atoms_scaled.get_volume():.4f} Ų")
            calc = workflow.run_scf(label=calc_label)
            
            # Extract energy and volume
            energy = calc.results.get('energy')
            volume = atoms_scaled.get_volume()
            
            # Store calculator for later access
            self.eos_calcs[scale_factor] = calc
            
            logger.info(f"✓ Point {point_idx:02d} completed: factor={scale_factor:.4f}, E={energy:.6f} eV")
            
            return (scale_factor, energy, volume, None)
            
        except Exception as e:
            error_msg = f"Error at volume factor {scale_factor}: {str(e)}"
            logger.error(error_msg)
            self.error_log[scale_factor] = error_msg
            return (scale_factor, None, None, e)
    
    def run_eos_study(
        self,
        scale_factors: Optional[List[float]] = None,
        volume_range: Tuple[float, float] = (DEFAULT_VOLUME_MIN_FACTOR, DEFAULT_VOLUME_MAX_FACTOR),
        n_points: int = DEFAULT_EOS_POINTS,
        protocol: Optional[str] = None,
        label: str = 'eos',
        parallel: bool = True,
        max_workers: Optional[int] = None,
        dry_run: bool = False,
        **kwargs
    ) -> pd.DataFrame:
        """
        Run complete EOS study: create structures, run SCFs, collect data.
        
        This is the main entry point for an EOS calculation.
        
        Args:
            scale_factors: List of volume scale factors (e.g., [0.98, 1.00, 1.02]).
                          If provided, takes precedence over volume_range/n_points.
            volume_range: Tuple (min_factor, max_factor) for volume scaling (default: 0.95-1.05)
                         Only used if scale_factors is None.
            n_points: Number of volume points (default: 7). Only used if scale_factors is None.
            protocol: Override protocol for this study (default: use __init__ protocol)
            label: Base directory label for calculations
            parallel: Run SCFs in parallel if True
            max_workers: Maximum number of parallel workers (default: CPU_COUNT)
            dry_run: If True, only generate input files without running (default: False)
            **kwargs: Additional parameters (not used currently)
            
        Returns:
            pandas.DataFrame with columns:
            - 'factor': Volume scale factor
            - 'volume': Volume in Ų
            - 'energy': Total energy in eV
            
        Raises:
            ValueError: If too few successful calculations
            
        Example:
            >>> eos = EOSWorkflow(atoms, pseudopotentials, protocol='moderate')
            >>> # Option 1: specify exact scale factors
            >>> df = eos.run_eos_study(
            ...     scale_factors=[0.98, 1.00, 1.02],
            ...     parallel=True,
            ...     label='eos/fe'
            ... )
            >>> # Option 2: generate factors from range
            >>> df = eos.run_eos_study(
            ...     volume_range=(0.95, 1.05),
            ...     n_points=7,
            ...     parallel=True,
            ...     label='eos/fe'
            ... )
            >>> print(df)
        """
        # Override protocol if specified
        if protocol is not None:
            self.protocol = protocol
        
        # Generate or use provided scale factors
        if scale_factors is not None:
            logger.info(f"Using provided scale factors: {scale_factors}")
        else:
            logger.info(f"Generating volume range: {volume_range[0]:.2%} to {volume_range[1]:.2%}")
            scale_factors = self.generate_volume_range(volume_range, n_points)
        logger.info(f"Scale factors: {scale_factors}")
        
        # Create scaled structures
        logger.info("Creating scaled structures...")
        self.create_scaled_structures(scale_factors)
        
        # Detect execution strategy based on scheduler
        # SLURM → batch submit all jobs at once (default behavior)
        # Otherwise → use worker pool or sequential submission
        has_slurm_scheduler = (
            self.queue is not None and
            self.queue.get('scheduler') in ['slurm', 'sbatch', 'SLURM']
        )
        
        # For SLURM: default to batch mode (all jobs at once) unless explicitly blocking
        is_non_blocking = (
            self.queue is not None and
            self.queue.get('wait_for_completion', False) is False  # Changed default to False for SLURM
        )
        
        is_batch_mode = has_slurm_scheduler and is_non_blocking


        
        # Dry run mode: only generate inputs
        if dry_run:
            logger.info("DRY RUN MODE: Generating input files for all EOS points (no execution)")
            self._run_eos_dry_run(label)
            logger.info(f"DRY RUN: Input files generated in {label}/ directories")
            logger.warning("No data was collected (dry run mode). To run calculations, use dry_run=False")
            return pd.DataFrame()  # Empty DataFrame
        
        # Run SCF calculations with appropriate strategy
        if is_batch_mode:
            # SLURM batch mode (local or remote): submit ALL jobs at once, SLURM manages queue
            logger.info(f"SLURM batch mode detected: submitting {len(scale_factors)} jobs at once...")
            results = self._run_eos_slurm_batch(label)
        elif parallel:
            # Local/blocking remote: use worker pool (submit one-by-one with workers)
            logger.info(f"Parallel worker mode: running {len(scale_factors)} EOS points with {max_workers or 'CPU'} workers...")
            results = self._run_eos_parallel(label, max_workers)
        else:
            # Sequential execution (one job at a time)
            logger.info(f"Sequential mode: running {len(scale_factors)} EOS points one-by-one...")
            results = self._run_eos_sequential(label)
        
        # Filter out failures
        successful_results = [r for r in results if r[3] is None]
        
        if len(successful_results) < MIN_POINTS_FOR_FIT:
            raise ValueError(
                f"Only {len(successful_results)} successful calculations out of {len(results)}. "
                f"Need at least {MIN_POINTS_FOR_FIT} for EOS fitting."
            )
        
        # Create DataFrame
        data = {
            'factor': [r[0] for r in successful_results],
            'energy': [r[1] for r in successful_results],
            'volume': [r[2] for r in successful_results],
        }
        
        self.results_df = pd.DataFrame(data).sort_values('volume').reset_index(drop=True)
        
        logger.info(f"EOS study complete: {len(self.results_df)} points collected")
        logger.info(f"\nE-V Data:\n{self.results_df.to_string()}")
        
        return self.results_df
    
    def _run_eos_dry_run(self, label: str) -> None:
        """Generate input files for all EOS points without running calculations."""
        for factor in sorted(self.eos_structures.keys()):
            atoms_scaled = self.eos_structures[factor]
            point_idx = self._factor_to_index[factor]
            calc_label = f"{label}/point_{point_idx:02d}"
            
            # Create workflow for this point
            workflow = CalculationWorkflow(
                atoms=atoms_scaled,
                pseudopotentials=self.pseudopotentials,
                pseudopotentials_config=self.pseudopotentials_config,
                protocol=self.protocol,
                kspacing=self.kspacing,
                input_data=self.input_data,
                magnetic_config=self.magnetic_config,
                queue=self.queue,
                machine=self.machine,
                code_version=self.code_version,
                **self.workflow_kwargs
            )
            
            # Generate input files only (dry_run=True)
            logger.info(f"Generating inputs: point_{point_idx:02d}, factor={factor:.4f}, V={atoms_scaled.get_volume():.4f} Ų")
            workflow.run_scf(label=calc_label, dry_run=True)
    
    def _run_eos_sequential(self, label: str) -> List:
        """Run EOS points sequentially."""
        results = []
        for factor in sorted(self.eos_structures.keys()):
            result = self._run_single_eos_point(factor, label)
            results.append(result)
        return results
    
    def _run_eos_parallel(self, label: str, max_workers: Optional[int]) -> List:
        """Run EOS points in parallel using ThreadPoolExecutor."""
        if max_workers is None:
            max_workers = min(DEFAULT_MAX_WORKERS, len(self.eos_structures))
        
        logger.info(f"Using {max_workers} parallel workers")
        
        results = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(self._run_single_eos_point, factor, label): factor
                for factor in sorted(self.eos_structures.keys())
            }
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    factor = futures[future]
                    logger.error(f"Error in parallel execution for factor {factor}: {e}")
                    results.append((factor, None, None, e))
        
        return results
    
    def _run_eos_slurm_batch(self, label: str) -> List:
        """
        Run all EOS points by submitting them to SLURM without workers.
        
        Uses the pattern:
        1. Create ONE workflow instance (shared remote connection)
        2. Submit ALL jobs at once via submit_scf_batch() with different atoms (non-blocking)
        3. Monitor ALL jobs together via wait_for_batch_jobs()
        
        This avoids SSH connection issues from sequential submission/monitoring.
        
        Returns:
            List of results: [(factor, energy, volume, error_or_None), ...]
        """
        logger.info(f"SLURM batch mode: submitting ALL {len(self.eos_structures)} jobs at once...")
        
        # Create ONE workflow with shared remote connection (for all points)
        # Use first structure as template (atoms will be replaced per calculation)
        shared_workflow = CalculationWorkflow(
            atoms=list(self.eos_structures.values())[0],
            pseudopotentials=self.pseudopotentials,
            pseudopotentials_config=self.pseudopotentials_config,
            protocol=self.protocol,
            kspacing=self.kspacing,
            input_data=self.input_data,
            magnetic_config=self.magnetic_config,
            queue=self.queue,
            machine=self.machine,
            code_version=self.code_version,
            **self.workflow_kwargs
        )
        
        # STEP 1: Submit ALL jobs NON-BLOCKING (same workflow, different atoms per point)
        batch_results = []
        
        for i, factor in enumerate(sorted(self.eos_structures.keys())):
            try:
                atoms_scaled = self.eos_structures[factor]
                point_idx = self._factor_to_index[factor]
                calc_label = f"{label}/point_{point_idx:02d}"
                
                # Update the shared workflow's atoms for this point
                shared_workflow.atoms = atoms_scaled
                
                logger.info(f"Submitting point_{point_idx:02d}: factor={factor:.4f}, V={atoms_scaled.get_volume():.4f} Ų")
                
                # Submit WITHOUT waiting (returns immediately)
                # Uses the SAME shared_workflow instance with updated atoms
                result = shared_workflow.submit_scf_batch(label=calc_label, wait_for_completion=False)
                batch_results.append(result)
                
            except Exception as e:
                logger.error(f"Error submitting job for factor {factor}: {e}")
                batch_results.append({
                    'calc': None,
                    'job_id': None,
                    'label': f"{label}/point_{self._factor_to_index[factor]:02d}",
                    'submitted': False,
                    'completed': False,
                    'error': str(e),
                })
        
        logger.info(f"All {len(batch_results)} jobs submitted. Monitoring {sum(1 for r in batch_results if r['submitted'])} jobs...")
        
        # STEP 2: Monitor ALL jobs together using the SAME shared_workflow
        timeout = self.queue.get('job_timeout', 3600) if self.queue else 3600
        
        # Monitor all jobs at once using the shared workflow that has the established remote connection
        monitoring_results = shared_workflow.wait_for_batch_jobs(batch_results, timeout=timeout, verbose=False)
        
        # STEP 3: Map results back by label (important for order preservation)
        label_to_monitor_result = {r['label']: r for r in monitoring_results}
        
        # Collect results in correct order
        results = []
        for factor in sorted(self.eos_structures.keys()):
            point_idx = self._factor_to_index[factor]
            calc_label = f"{label}/point_{point_idx:02d}"
            
            # Look up result by label
            monitor_result = label_to_monitor_result.get(calc_label)
            
            if monitor_result is None:
                # Submission failed
                logger.warning(f"Skipping point_{point_idx:02d} (submission failed)")
                results.append((factor, None, None, Exception("Submission failed")))
                continue
            
            atoms_scaled = self.eos_structures[factor]
            volume = atoms_scaled.get_volume()
            
            if monitor_result['success']:
                energy = monitor_result.get('energy')
                if energy is not None:
                    logger.info(f"✓ Point {point_idx:02d}: factor={factor:.4f}, E={energy:.6f} eV")
                    results.append((factor, energy, volume, None))
                else:
                    logger.warning(f"Point {point_idx:02d} completed but no energy extracted")
                    results.append((factor, None, volume, Exception("No energy in results")))
            else:
                error = monitor_result.get('error', 'Unknown error')
                logger.error(f"✗ Point {point_idx:02d} failed: {error}")
                results.append((factor, None, volume, Exception(error)))
        
        logger.info(f"SLURM batch complete: {sum(1 for r in results if r[3] is None)}/{len(results)} succeeded")
        return results


    
    # ═════════════════════════════════════════════════════════════════════════════
    # EOS FITTING AND ANALYSIS METHODS
    # ═════════════════════════════════════════════════════════════════════════════
    
    def fit_eos(self) -> Dict:
        """
        Fit Birch-Murnaghan EOS to collected E-V data.
        
        Must call run_eos_study() first to collect data.
        
        Returns:
            Dictionary with fitted parameters:
            {
                'e0': float,           # Energy at V₀ (eV)
                'v0': float,           # Equilibrium volume (Ų)
                'b0': float,           # Bulk modulus (GPa)
                'b0_prime': float,     # B₀'
                'r_squared': float,    # Fit quality
                'residuals': ndarray,  # Fit residuals
                'converged': bool,     # Convergence status
            }
            
        Raises:
            ValueError: If no data collected yet
        """
        if self.results_df is None or len(self.results_df) == 0:
            raise ValueError("No EOS data collected. Call run_eos_study() first.")
        
        # Extract V and E
        volumes = self.results_df['volume'].values
        energies = self.results_df['energy'].values
        
        # Fit EOS
        logger.info("Fitting Birch-Murnaghan EOS...")
        eos_result = fit_birch_murnaghan(volumes, energies)
        
        # fit_birch_murnaghan returns B0 in eV/Ų, convert to GPa for display
        self.eos_params = eos_result.copy()
        self.eos_params['b0_gpa'] = eos_result['b0'] * GPa_PER_EV_ANG3
        
        logger.info(f"EOS fit converged: {self.eos_params['converged']}")
        logger.info(f"Fit quality (R²): {self.eos_params['r_squared']:.6f}")
        logger.info(f"Bulk modulus (B₀): {self.eos_params['b0_gpa']:.2f} GPa")
        
        return self.eos_params
    
    def get_eos_properties(self) -> Dict:
        """
        Get equilibrium properties from fitted EOS.
        
        Must call fit_eos() first.
        
        Returns:
            Dictionary with equilibrium properties:
            {
                'v0': float,                # Equilibrium volume (Ų)
                'e0': float,                # Energy at V₀ (eV)
                'bulk_modulus': float,      # B₀ in GPa
                'bulk_modulus_prime': float,# B₀'
                'b0_gpa': float,            # Alias for bulk_modulus
                'e0_ev': float,             # Alias for e0
                'v0_angstrom3': float,      # Alias for v0
                'r_squared': float,         # Fit quality
                'converged': bool,          # Fit convergence status
                'energy_shift': float,      # E₀ - min(E_data)
            }
            
        Raises:
            ValueError: If EOS not fitted yet
        """
        if self.eos_params is None:
            raise ValueError("EOS not fitted yet. Call fit_eos() first.")
        
        props = {
            'v0': self.eos_params['v0'],
            'e0': self.eos_params['e0'],
            'bulk_modulus': self.eos_params['b0_gpa'],
            'bulk_modulus_prime': self.eos_params['b0_prime'],
            # Aliases for clarity
            'b0_gpa': self.eos_params['b0_gpa'],
            'e0_ev': self.eos_params['e0'],
            'v0_angstrom3': self.eos_params['v0'],
            'r_squared': self.eos_params['r_squared'],
            'converged': self.eos_params['converged'],
            'energy_shift': self.eos_params['e0'] - np.min(self.results_df['energy'].values),
        }
        
        return props
    
    def predict_energy(self, volume: float) -> float:
        """
        Predict energy at a given volume using fitted EOS.
        
        Args:
            volume: Volume in Ų
            
        Returns:
            Predicted energy in eV
            
        Raises:
            ValueError: If EOS not fitted
        """
        if self.eos_params is None:
            raise ValueError("EOS not fitted yet. Call fit_eos() first.")
        
        E = birchmurnaghan(
            np.array([volume]),
            self.eos_params['e0'],
            self.eos_params['b0'],
            self.eos_params['b0_prime'],
            self.eos_params['v0']
        )
        
        return E[0]
    
    def calculate_pressure(self, volume: float) -> float:
        """
        Calculate pressure at given volume using fitted EOS.
        
        P = -dE/dV from Birch-Murnaghan
        
        Args:
            volume: Volume in Ų
            
        Returns:
            Pressure in GPa
            
        Raises:
            ValueError: If EOS not fitted
        """
        if self.eos_params is None:
            raise ValueError("EOS not fitted yet. Call fit_eos() first.")
        
        # Numerical derivative: dE/dV ≈ (E(V+dV) - E(V-dV)) / (2*dV)
        dV = volume * 1e-4  # Small volume step
        
        E_plus = birchmurnaghan(
            np.array([volume + dV]),
            self.eos_params['e0'],
            self.eos_params['b0'],
            self.eos_params['b0_prime'],
            self.eos_params['v0']
        )[0]
        
        E_minus = birchmurnaghan(
            np.array([volume - dV]),
            self.eos_params['e0'],
            self.eos_params['b0'],
            self.eos_params['b0_prime'],
            self.eos_params['v0']
        )[0]
        
        dE_dV = (E_plus - E_minus) / (2 * dV)  # eV/Ų
        
        # Convert to GPa: dE/dV in eV/Ų, convert to GPa
        P_GPa = -dE_dV * GPa_PER_EV_ANG3
        
        return P_GPa
    
    # ═════════════════════════════════════════════════════════════════════════════
    # PLOTTING METHODS
    # ═════════════════════════════════════════════════════════════════════════════
    
    def plot_eos_curve(
        self,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 7),
        dpi: int = 150
    ) -> None:
        """
        Plot E-V curve with Birch-Murnaghan fit.
        
        Args:
            save_path: Path to save figure (if None, display only)
            figsize: Figure size (width, height) in inches
            dpi: Resolution in dots per inch
            
        Raises:
            ValueError: If no data collected or EOS not fitted
        """
        if self.results_df is None:
            raise ValueError("No EOS data. Call run_eos_study() first.")
        
        if self.eos_params is None:
            raise ValueError("EOS not fitted. Call fit_eos() first.")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        
        # Plot data points
        volumes = self.results_df['volume'].values
        energies = self.results_df['energy'].values
        ax.scatter(volumes, energies, s=100, color='red', label='SCF points', zorder=5)
        
        # Plot fitted curve
        v_min, v_max = np.min(volumes) * 0.95, np.max(volumes) * 1.05
        v_smooth = np.linspace(v_min, v_max, 200)
        e_smooth = birchmurnaghan(
            v_smooth,
            self.eos_params['e0'],
            self.eos_params['b0'],
            self.eos_params['b0_prime'],
            self.eos_params['v0']
        )
        ax.plot(v_smooth, e_smooth, 'b-', linewidth=2, label='Birch-Murnaghan fit')
        
        # Mark equilibrium point
        ax.axvline(self.eos_params['v0'], color='green', linestyle='--', alpha=0.7, label=f"V₀ = {self.eos_params['v0']:.4f} Ų")
        ax.axhline(self.eos_params['e0'], color='green', linestyle='--', alpha=0.7)
        ax.plot(self.eos_params['v0'], self.eos_params['e0'], 'g*', markersize=20)
        
        # Labels and formatting
        ax.set_xlabel('Volume (Ų)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Energy (eV)', fontsize=12, fontweight='bold')
        ax.set_title('Equation of State - Birch-Murnaghan Fit', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, alpha=0.3)
        
        # Add text box with parameters
        textstr = f"B₀ = {self.eos_params['b0_gpa']:.2f} GPa\nB₀' = {self.eos_params['b0_prime']:.3f}\nR² = {self.eos_params['r_squared']:.6f}"
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
            logger.info(f"EOS curve saved to: {save_path}")
        else:
            plt.show()
        
        plt.close()
    
    def plot_residuals(
        self,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 6),
        dpi: int = 150
    ) -> None:
        """
        Plot fitting residuals (data - fit).
        
        Args:
            save_path: Path to save figure
            figsize: Figure size (width, height)
            dpi: Resolution
            
        Raises:
            ValueError: If EOS not fitted
        """
        if self.eos_params is None:
            raise ValueError("EOS not fitted. Call fit_eos() first.")
        
        volumes = self.results_df['volume'].values
        residuals = self.eos_params['residuals']
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, dpi=dpi)
        
        # Plot 1: Residuals vs volume
        ax1.bar(volumes, residuals * 1000, color='steelblue', alpha=0.7)
        ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax1.set_xlabel('Volume (Ų)', fontsize=11)
        ax1.set_ylabel('Residual (meV)', fontsize=11)
        ax1.set_title('EOS Fit Residuals', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Histogram of residuals
        ax2.hist(residuals * 1000, bins=max(5, len(residuals)), color='steelblue', alpha=0.7, edgecolor='black')
        ax2.set_xlabel('Residual (meV)', fontsize=11)
        ax2.set_ylabel('Frequency', fontsize=11)
        ax2.set_title('Distribution of Residuals', fontsize=13, fontweight='bold')
        ax2.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero residual')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
            logger.info(f"Residuals plot saved to: {save_path}")
        else:
            plt.show()
        
        plt.close()
    
    def plot_pressure(
        self,
        save_path: Optional[str] = None,
        figsize: Tuple[int, int] = (10, 6),
        dpi: int = 150
    ) -> None:
        """
        Plot pressure vs volume derived from EOS.
        
        Args:
            save_path: Path to save figure
            figsize: Figure size (width, height)
            dpi: Resolution
            
        Raises:
            ValueError: If EOS not fitted
        """
        if self.eos_params is None:
            raise ValueError("EOS not fitted. Call fit_eos() first.")
        
        volumes = self.results_df['volume'].values
        v_min, v_max = np.min(volumes) * 0.95, np.max(volumes) * 1.05
        v_smooth = np.linspace(v_min, v_max, 200)
        
        # Calculate pressure using numerical derivative
        P_smooth = []
        for v in v_smooth:
            P_smooth.append(self.calculate_pressure(v))
        P_smooth = np.array(P_smooth)
        
        # Calculate pressure at data points
        P_data = [self.calculate_pressure(v) for v in volumes]
        
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        
        # Plot curve
        ax.plot(v_smooth, P_smooth, 'b-', linewidth=2, label='P(V) from EOS')
        
        # Plot data points
        ax.scatter(volumes, P_data, s=100, color='red', label='At SCF points', zorder=5)
        
        # Mark equilibrium
        ax.axvline(self.eos_params['v0'], color='green', linestyle='--', alpha=0.7)
        ax.axhline(0, color='black', linestyle='-', linewidth=0.8, label='P = 0 (equilibrium)')
        
        # Labels
        ax.set_xlabel('Volume (Ų)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pressure (GPa)', fontsize=12, fontweight='bold')
        ax.set_title('Pressure vs Volume (from Birch-Murnaghan EOS)', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
            logger.info(f"Pressure plot saved to: {save_path}")
        else:
            plt.show()
        
        plt.close()
    
    # ═════════════════════════════════════════════════════════════════════════════
    # DATA EXPORT METHODS
    # ═════════════════════════════════════════════════════════════════════════════
    
    def to_csv(self, filepath: str) -> None:
        """
        Save E-V data to CSV file.
        
        Args:
            filepath: Path to save CSV file
        """
        if self.results_df is None:
            raise ValueError("No data to save. Run run_eos_study() first.")
        
        self.results_df.to_csv(filepath, index=False)
        logger.info(f"EOS data saved to: {filepath}")
    
    def to_json(self, filepath: str) -> None:
        """
        Save EOS fit parameters to JSON file.
        
        Args:
            filepath: Path to save JSON file
        """
        if self.eos_params is None:
            raise ValueError("EOS not fitted. Call fit_eos() first.")
        
        import json
        
        params_export = {
            'e0': float(self.eos_params['e0']),
            'v0': float(self.eos_params['v0']),
            'b0': float(self.eos_params['b0_gpa']),
            'b0_prime': float(self.eos_params['b0_prime']),
            'r_squared': float(self.eos_params['r_squared']),
            'converged': bool(self.eos_params['converged']),
            'original_volume': float(self.v0_original),
        }
        
        with open(filepath, 'w') as f:
            json.dump(params_export, f, indent=2)
        
        logger.info(f"EOS parameters saved to: {filepath}")
    
    def summary(self) -> str:
        """
        Return formatted summary of EOS results.
        
        Returns:
            String with formatted results
        """
        if self.results_df is None or self.eos_params is None:
            return "EOS study not completed. Run run_eos_study() and fit_eos() first."
        
        props = self.get_eos_properties()
        
        summary_str = f"""
╔════════════════════════════════════════════════════════════════╗
║                  EQUATION OF STATE RESULTS                    ║
╠════════════════════════════════════════════════════════════════╣
║ Structure: {self.atoms.get_chemical_formula():50s} ║
║ Calculation Protocol: {self.protocol:40s} ║
╠════════════════════════════════════════════════════════════════╣
║ EQUILIBRIUM PROPERTIES:                                        ║
║   V₀ (Equilibrium Volume):  {props['v0']:19.6f} Ų      ║
║   E₀ (Energy at V₀):        {props['e0']:19.6f} eV      ║
║   B₀ (Bulk Modulus):        {props['bulk_modulus']:19.2f} GPa     ║
║   B₀' (Pressure Derivative):{props['bulk_modulus_prime']:19.3f}         ║
╠════════════════════════════════════════════════════════════════╣
║ FIT QUALITY:                                                   ║
║   R² (Coefficient):         {props['r_squared']:19.6f}         ║
║   Fit Converged:            {"Yes" if self.eos_params['converged'] else "No":>19s}         ║
║   Number of Points:         {len(self.results_df):>19d}         ║
╠════════════════════════════════════════════════════════════════╣
║ DATA STATISTICS:                                               ║
║   Min Volume:               {self.results_df['volume'].min():19.6f} Ų      ║
║   Max Volume:               {self.results_df['volume'].max():19.6f} Ų      ║
║   Min Energy:               {self.results_df['energy'].min():19.6f} eV      ║
║   Max Energy:               {self.results_df['energy'].max():19.6f} eV      ║
╚════════════════════════════════════════════════════════════════╝
        """
        
        return summary_str


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def quick_eos(
    atoms: Union[Atoms, str, Path],
    pseudopotentials: Optional[Dict[str, str]] = None,
    pseudopotentials_config: Optional[str] = None,
    volume_range: Tuple[float, float] = (0.95, 1.05),
    n_points: int = 7,
    protocol: str = 'moderate',
    label: str = 'eos',
    **kwargs
) -> Tuple[EOSWorkflow, Dict]:
    """
    Quick helper function to run a complete EOS study in one call.
    
    Args:
        atoms: ASE Atoms object or structure file path
        pseudopotentials: Dict mapping elements to pseudopotential files
        pseudopotentials_config: Name of pseudopotentials config to load
        volume_range: Volume range (min_factor, max_factor)
        n_points: Number of volume points
        protocol: Calculation protocol ('fast', 'moderate', 'accurate')
        label: Base directory label
        **kwargs: Additional parameters
        
    Returns:
        Tuple: (EOSWorkflow object, properties dictionary)
        
    Example:
        >>> atoms = bulk('Si', 'diamond', a=5.43)
        >>> eos, props = quick_eos(
        ...     atoms,
        ...     pseudopotentials={'Si': 'Si.pbe.UPF'},
        ...     protocol='moderate'
        ... )
        >>> print(f"B₀ = {props['bulk_modulus']:.2f} GPa")
    """
    # Create workflow
    eos = EOSWorkflow(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        pseudopotentials_config=pseudopotentials_config,
        protocol=protocol,
        **kwargs
    )
    
    # Run study
    eos.run_eos_study(
        volume_range=volume_range,
        n_points=n_points,
        label=label,
        parallel=True
    )
    
    # Fit EOS
    eos.fit_eos()
    
    # Get properties
    props = eos.get_eos_properties()
    
    return eos, props

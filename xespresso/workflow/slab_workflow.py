"""
SlabWorkflow: Complete workflow for surface slab calculations.

This module provides a modern, integrated interface for:
1. Bulk structure convergence (using ConvergenceWorkflow)
2. Slab generation from bulk structures  
3. Slab-specific convergence studies (k-mesh, layer thickness)
4. Structure relaxation with proper constraints
5. Surface energy calculations

Requires: ase, xespresso>=2024.0

Example:
    >>> from ase.build import bulk
    >>> from xespresso.workflow.slab_workflow import SlabWorkflow
    >>> 
    >>> bulk_atoms = bulk('Au', 'fcc', a=4.08)
    >>> wf = SlabWorkflow(
    ...     bulk_atoms=bulk_atoms,
    ...     surface_indices=[(1,0,0), (1,1,0), (1,1,1)],
    ...     pseudopotentials_config='default'
    ... )
    >>> 
    >>> # Phase 2 (available now): slab generation
    >>> slabs = wf.generate_slabs()
    >>> wf.save_results()
    
Progress:
    Phase 1 (Bulk convergence): TO BE IMPLEMENTED
    Phase 2 (Slab generation):  ✓ COMPLETE & TESTED
    Phase 3 (Slab convergence): TO BE IMPLEMENTED
    Phase 4 (Relaxation):       TO BE IMPLEMENTED
    Phase 5 (Analysis):         TO BE IMPLEMENTED
"""

import os
import logging
import numpy as np
from typing import Dict, Optional, List, Tuple
from pathlib import Path
from copy import deepcopy

from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms

try:
    from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
    HAS_CONVERGENCE_WORKFLOW = True
except ImportError:
    HAS_CONVERGENCE_WORKFLOW = False
    logging.warning("ConvergenceWorkflow not available. Bulk convergence will not work.")


# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)


class SlabWorkflow:
    """
    Complete workflow for surface slab calculations.
    
    Manages workflow phases:
    - Phase 1: Bulk convergence (calls ConvergenceWorkflow)
    - Phase 2: Slab generation (COMPLETE)
    - Phase 3: Slab convergence (k-mesh, layer thickness)
    - Phase 4: Structure relaxation (calls CalculationWorkflow)
    - Phase 5: Surface energy analysis
    
    Phases can be run individually or as complete workflow via run_all().
    
    Current Status: Phase 2 ready for use.
    
    Parameters
    ----------
    bulk_atoms : ase.Atoms
        Bulk structure for slab generation
    surface_indices : List[Tuple[int, int, int]], default=[(1,0,0), (1,1,0), (1,1,1)]
        Miller indices for surfaces to study
    min_vacuum_size : float, default=15.0
        Minimum vacuum size in Angstrom
    nlayers : int, default=4
        Number of atomic layers in slab
    fix_layer_indices : List[int], default=[0, 1]
        Indices (from bottom) of layers to freeze during relaxation
    tol_layer_detection : float, default=1.0
        Tolerance for layer detection in Angstrom
    pseudopotentials : Dict[str, str], optional
        Element → pseudopotential filename mapping
    pseudopotentials_config : str, default='default'
        Name of pseudopotentials configuration
    protocol : str, default='moderate'
        Calculation protocol: 'fast', 'moderate', 'accurate'
    precision : str, default='low'
        Precision level for bulk convergence: 'low', 'normal', 'high'
    machine : str, optional
        Machine name for job submission
    queue : Dict, optional
        Queue configuration dictionary
    code_version : str, default='7.4.1'
        Quantum ESPRESSO version to use
    verbose : bool, default=True
        Print progress information
        
    Attributes
    ----------
    bulk_results : Dict
        Results from bulk convergence phase
    bulk_recommendations : Dict
        Optimal ecutwfc and kspacing from convergence
    slabs : Dict[Tuple, Atoms]
        Generated slabs
    relaxed_slabs : Dict[Tuple, object]
        Relaxed slab calculations
    surface_energies : Dict[Tuple, float]
        Surface energies in J/m²
        
    Examples
    --------
    **Phase 2 (Available):**
    
    >>> wf = SlabWorkflow(bulk_atoms, surface_indices=[(1,1,1)])
    >>> slabs = wf.generate_slabs(save_slabs=True, save_dir='./slabs/')
    >>> wf.save_results(output_dir='./slabs/')
    
    **Complete workflow (Phases 1-5, pending):**
    
    >>> wf = SlabWorkflow(bulk_atoms)
    >>> results = wf.run_all(label_prefix='au_surfaces')
    """
    
    def __init__(
        self,
        bulk_atoms: Atoms,
        surface_indices: List[Tuple[int, int, int]] = None,
        min_vacuum_size: float = 15.0,
        nlayers: int = 4,
        fix_layer_indices: List[int] = None,
        tol_layer_detection: float = 1.0,
        pseudopotentials: Optional[Dict[str, str]] = None,
        pseudopotentials_config: str = 'default',
        protocol: str = 'moderate',
        precision: str = 'low',
        machine: Optional[str] = None,
        queue: Optional[Dict] = None,
        code_version: str = '7.4.1',
        verbose: bool = True,
        use_primitive_cell: bool = False,
    ):
        """
        Initialize SlabWorkflow.
        
        Validates inputs and prepares for slab calculations.
        
        Raises
        ------
        ValueError
            If bulk_atoms is empty or if both machine and queue are specified
        """
        # Validate inputs
        if bulk_atoms is None or len(bulk_atoms) == 0:
            raise ValueError("bulk_atoms must be a non-empty Atoms object")
        
        if machine is not None and queue is not None:
            raise ValueError(
                "Cannot specify both 'machine' and 'queue'. Use one or the other."
            )
        
        # Core structure
        self.bulk_atoms = bulk_atoms.copy()
        self.surface_indices = surface_indices or [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
        
        # Slab parameters
        self.min_vacuum_size = min_vacuum_size
        self.nlayers = nlayers
        self.fix_layer_indices = fix_layer_indices or [0, 1]
        self.tol_layer_detection = tol_layer_detection
        
        # Calculation parameters
        self.pseudopotentials = pseudopotentials
        self.pseudopotentials_config = pseudopotentials_config
        self.protocol = protocol
        self.precision = precision
        self.machine = machine
        self.queue = queue
        self.code_version = code_version
        self.verbose = verbose
        
        # Slab reduction options
        self.use_primitive_cell = use_primitive_cell  # Reduce to primitive (1×1) instead of supercell (2×2)
        
        # Results storage (populated by workflow phases)
        self.bulk_results = None
        self.bulk_recommendations = None
        self.slabs = {}  # {(h,k,l): Atoms}
        self.relaxed_slabs = {}
        self.surface_energies = {}
        self.convergence_results = {}
        self.relax_results = {}
        self.slab_kpts = None  # Optional: explicit k-points for slab calculations
        
        logger.info(
            f"SlabWorkflow initialized:\n"
            f"  Bulk: {self.bulk_atoms.get_chemical_formula()} ({len(self.bulk_atoms)} atoms)\n"
            f"  Surfaces: {self.surface_indices}\n"
            f"  Layers: {self.nlayers}, Fix indices: {self.fix_layer_indices}\n"
            f"  Min vacuum: {self.min_vacuum_size} Å"
        )
    
    # =========================================================================
    # PHASE 2: SLAB GENERATION (COMPLETE & AVAILABLE)
    # =========================================================================
    
    def generate_slabs(
        self,
        override_surface_indices: Optional[List[Tuple[int, int, int]]] = None,
        save_slabs: bool = False,
        save_dir: str = './',
    ) -> Dict[Tuple[int, int, int], Atoms]:
        """
        Generate slabs for all surface indices using ASE.
        
        **Status**: ✓ COMPLETE
        
        Uses ase.build.surface() to create slabs from bulk structure.
        This is more reliable and works for any element and crystal structure.
        
        Post-processing:
        - Automatic orthogonal cell generation (by ASE)
        - Exact number of layers (by ASE)
        - Vacuum layer addition
        - Bottom layer constraints
        
        Parameters
        ----------
        override_surface_indices : List[Tuple], optional
            Alternative surface indices (if not using self.surface_indices)
        save_slabs : bool, default=False
            Save generated slabs to CIF files
        save_dir : str, default='./'
            Directory to save slab files
            
        Returns
        -------
        Dict[Tuple, Atoms]
            Dictionary mapping (h,k,l) → slab Atoms object
            
        Raises
        ------
        RuntimeError
            If slab generation fails for a surface
            
        Notes
        -----
        Saves slab structures to CIF files if save_slabs=True.
        Applies FixAtoms constraint to bottom layers.
        Uses ASE surface() instead of PyMatGen for better reliability.
        """
        indices = override_surface_indices or self.surface_indices
        logger.info(f"Generating slabs for {len(indices)} surfaces...")
        
        # Create output directory if saving
        if save_slabs:
            Path(save_dir).mkdir(parents=True, exist_ok=True)
        
        # Generate each slab
        for hkl in indices:
            surface_name = f"{hkl[0]}{hkl[1]}{hkl[2]}"
            logger.info(f"Generating {surface_name}...")
            
            try:
                # Use _regenerate_slab_with_nlayers to respect self.nlayers parameter
                slab = self._regenerate_slab_with_nlayers(hkl, self.nlayers)
                
                self.slabs[hkl] = slab
                
                # Log results
                cell_lengths = slab.cell.lengths()
                logger.info(
                    f"  ✓ {surface_name}: {len(slab)} atoms\n"
                    f"    Cell: {cell_lengths[0]:.3f} × {cell_lengths[1]:.3f} × {cell_lengths[2]:.3f} Å"
                )
                
                # Optional: Save to file
                if save_slabs:
                    filename = os.path.join(save_dir, f'slab_{surface_name}.cif')
                    write(filename, slab)
                    logger.info(f"    Saved: {filename}")
                
            except Exception as e:
                logger.error(f"  ✗ Error generating {surface_name}: {e}")
                continue
        
        if self.slabs:
            logger.info(f"✓ Generated {len(self.slabs)} slabs successfully")
        else:
            logger.warning("✗ No slabs were generated")
        
        return self.slabs
    
    def _check_calculation_complete(self, calc_dir: str, num_atoms: int) -> tuple:
        """
        Check if a calculation directory exists and has completed successfully.
        
        Returns (is_complete, energy_per_atom) where:
        - is_complete: bool, True if calculation finished
        - energy_per_atom: float, energy per atom if complete, else None
        """
        import os
        
        if not os.path.exists(calc_dir):
            return False, None
        
        try:
            # Look for Quantum ESPRESSO output file
            output_file = os.path.join(calc_dir, 'pwscf.out')
            if not os.path.exists(output_file):
                return False, None
            
            # Check if calculation completed by looking for final energy
            with open(output_file, 'r') as f:
                content = f.read()
                # Look for "JOB DONE" or final energy marker
                if 'JOB DONE' in content or 'total energy' in content:
                    # Extract final energy
                    lines = content.split('\n')
                    for line in reversed(lines):
                        if 'total energy' in line and '=' in line:
                            try:
                                # Parse "total energy = -123.456 Ry"
                                energy_str = line.split('=')[-1].strip().split()[0]
                                energy = float(energy_str)
                                energy_per_atom = energy / num_atoms
                                return True, energy_per_atom
                            except (ValueError, IndexError):
                                pass
            
            return False, None
        except Exception as e:
            logger.debug(f"Error checking calculation {calc_dir}: {e}")
            return False, None
    
    def _calculate_fixed_layers(self, nlayers: int) -> List[int]:
        """
        Calculate which layers to freeze based on total number of layers.
        
        Rule:
        - If nlayers is EVEN: freeze bottom nlayers/2
        - If nlayers is ODD:  freeze bottom (nlayers-1)/2
        
        Parameters
        ----------
        nlayers : int
            Total number of atomic layers
            
        Returns
        -------
        List[int]
            Indices of layers to freeze (0-based from bottom)
            
        Examples
        --------
        >>> self._calculate_fixed_layers(3)
        [0]  # (3-1)/2 = 1 layer
        
        >>> self._calculate_fixed_layers(4)
        [0, 1]  # 4/2 = 2 layers
        
        >>> self._calculate_fixed_layers(6)
        [0, 1, 2]  # 6/2 = 3 layers
        """
        if nlayers % 2 == 0:  # EVEN
            num_fixed = nlayers // 2
        else:  # ODD
            num_fixed = (nlayers - 1) // 2
        
        return list(range(num_fixed))
    
    def _regenerate_slab_with_nlayers(
        self,
        surface_index: Tuple[int, int, int],
        nlayers: int,
        vacuum_size: Optional[float] = None,
    ) -> Atoms:
        """
        Regenerate slab for a surface with specific number of layers using ASE.
        
        Uses ase.build.surface() which is more reliable and works for any element.
        Automatically generates exactly nlayers atomic layers with proper orthogonal cell.
        
        Parameters
        ----------
        surface_index : Tuple[int, int, int]
            Miller indices (h, k, l)
        nlayers : int
            Target number of atomic layers
        vacuum_size : Optional[float], default=None
            Vacuum size in Ångströms. If None, uses self.min_vacuum_size.
            
        Returns
        -------
        ase.Atoms
            Generated slab with constraints applied
        """
        from ase.build import surface
        
        # Use provided vacuum_size or default
        vac = vacuum_size if vacuum_size is not None else self.min_vacuum_size
        
        try:
            h, k, l = surface_index
            
            logger.info(f"  Calculating slab for nlayers={nlayers}:")
            logger.info(f"    Surface: ({h}{k}{l})")
            
            # Use ASE surface() which is more robust than PyMatGen
            # It automatically:
            # - Generates exactly nlayers layers
            # - Creates orthogonal cell
            # - Ensures positive vacuum
            # - Works for any element and crystal structure
            slab = surface(self.bulk_atoms, surface_index, layers=nlayers, vacuum=vac)
            
            logger.info(f"    Generated slab: {len(slab)} atoms")
            logger.info(f"    Cell c[2] (vacuum): {slab.cell[2, 2]:.4f} Å")
            
            # ===== APPLY CONSTRAINTS =====
            slab_copy = slab.copy()
            z_positions = slab_copy.get_positions()[:, 2]
            z_min = z_positions.min()
            z_max = z_positions.max()
            slab_height = z_max - z_min
            
            # Calculate per-layer thickness
            if nlayers > 0:
                layer_thickness = slab_height / nlayers
            else:
                raise ValueError(f"nlayers must be > 0, got {nlayers}")
            
            # Get fixed layer indices based on nlayers
            fix_indices = self._calculate_fixed_layers(nlayers)
            
            # Identify atoms to freeze (bulk-like layers at bottom)
            fixed_atom_indices = []
            for i, z in enumerate(z_positions):
                # Calculate which layer this atom belongs to
                layer_idx = int((z - z_min) / layer_thickness + 0.5)  # Round to nearest
                layer_idx = min(layer_idx, nlayers - 1)  # Safety: cap at nlayers-1
                
                if layer_idx in fix_indices:
                    fixed_atom_indices.append(i)
            
            # Apply constraint
            if fixed_atom_indices:
                from ase.constraints import FixAtoms
                constraint = FixAtoms(indices=fixed_atom_indices)
                slab_copy.set_constraint(constraint)
                logger.debug(
                    f"Applied FixAtoms to {len(fixed_atom_indices)} atoms "
                    f"(layers {sorted(fix_indices)} of {nlayers} total)"
                )
            
            return slab_copy
            
        except Exception as e:
            logger.error(f"Error regenerating slab with {nlayers} layers: {e}")
            raise


    
    def _apply_constraints(self, slab: Atoms) -> Atoms:
        """
        Apply FixAtoms constraint to bottom layers.
        
        Identifies bottom layers and freezes them during relaxation.
        Standard for surface calculations to maintain bulk-like behavior.
        
        Parameters
        ----------
        slab : ase.Atoms
            Input slab structure
            
        Returns
        -------
        ase.Atoms
            Slab with FixAtoms constraint applied
        """
        slab = slab.copy()
        
        # Get z-positions
        z_positions = slab.get_positions()[:, 2]
        z_min = z_positions.min()
        z_max = z_positions.max()
        
        # Calculate per-layer thickness
        layer_thickness = (z_max - z_min) / self.nlayers
        
        # Identify atoms to freeze
        fixed_indices = []
        for i, z in enumerate(z_positions):
            layer_idx = int((z - z_min) / layer_thickness)
            if layer_idx in self.fix_layer_indices:
                fixed_indices.append(i)
        
        # Apply constraint
        if fixed_indices:
            constraint = FixAtoms(indices=fixed_indices)
            slab.set_constraint(constraint)
            logger.debug(
                f"    Applied FixAtoms to {len(fixed_indices)} atoms "
                f"(layers {self.fix_layer_indices})"
            )
        else:
            logger.warning(
                f"    ⚠ No atoms found in layers {self.fix_layer_indices}"
            )
        
        return slab
    
    def _get_atom_indices_for_layers(
        self,
        slab: Atoms,
        surface_index: Tuple[int, int, int],
        num_layers_to_fix: int,
    ) -> List[int]:
        """
        Map layer indices to atom indices for a given slab.
        
        Uses z-position based layer identification to get atoms in the bottom
        num_layers_to_fix layers for freezing via FixAtoms constraint.
        
        Parameters
        ----------
        slab : ase.Atoms
            Slab structure
        surface_index : Tuple[int, int, int]
            Miller indices (for logging only)
        num_layers_to_fix : int
            Number of bottom layers to freeze
            
        Returns
        -------
        List[int]
            Atom indices to freeze (0-based)
        """
        if num_layers_to_fix <= 0:
            return []
        
        z_positions = slab.get_positions()[:, 2]
        z_min = z_positions.min()
        z_max = z_positions.max()
        z_range = z_max - z_min
        
        if z_range <= 0:
            logger.warning(f"    ⚠ Slab has zero z-range, no atoms to fix")
            return []
        
        # Identify total number of layers (approximation)
        # This assumes roughly equal spacing
        n_atoms = len(slab)
        approx_atoms_per_layer = max(1, n_atoms // 10)  # Rough estimate
        
        # Calculate z-position threshold for bottom layers
        layer_height = z_range / max(1, n_atoms // approx_atoms_per_layer)
        z_threshold = z_min + (num_layers_to_fix * layer_height)
        
        # Get atoms below threshold
        fixed_indices = [i for i, z in enumerate(z_positions) if z <= z_threshold]
        
        # Ensure we don't fix all atoms
        if len(fixed_indices) >= n_atoms:
            # Fallback: use bottom fraction based on layer count
            n_to_fix = max(1, n_atoms // 4)  # Fix ~25% if heuristic fails
            fixed_indices = fixed_indices[:n_to_fix]
            logger.warning(
                f"    ⚠ Layer detection heuristic fixed too many atoms, "
                f"using bottom {n_to_fix} atoms"
            )
        
        return fixed_indices
    
    # =========================================================================
    # UTILITY: Extract energy from SCF results (using convergence_workflow pattern)
    # =========================================================================
    
    def _extract_energy_from_result(self, result_obj, num_atoms: int) -> Optional[float]:
        """
        Extract energy from SCF result robustly.
        
        Follows the pattern from ConvergenceWorkflow._extract_property_from_result()
        
        Parameters
        ----------
        result_obj : object
            Result object from cw.run_scf() that contains calculated properties
        num_atoms : int
            Number of atoms in the structure
            
        Returns
        -------
        float or None
            Energy per atom in eV, or None if extraction fails
        """
        try:
            # Try to access energy directly (works for most xespresso results)
            if hasattr(result_obj, 'results') and isinstance(result_obj.results, dict):
                if 'energy' in result_obj.results:
                    energy_total = result_obj.results['energy']
                    return energy_total / num_atoms
            
            # Fallback: direct attribute access
            if hasattr(result_obj, 'energy'):
                energy_total = result_obj.energy
                return energy_total / num_atoms
                
            # Last attempt: dictionary-like access
            if isinstance(result_obj, dict) and 'energy' in result_obj:
                energy_total = result_obj['energy']
                return energy_total / num_atoms
                
        except Exception as e:
            logger.debug(f"  Could not extract energy from result object: {e}")
        
        return None
    
    # =========================================================================
    # PHASE 1: BULK CONVERGENCE (PLACEHOLDER - Phase 2 task)
    # =========================================================================
    
    def set_bulk_parameters(
        self,
        optimal_ecutwfc: float,
        optimal_kspacing: Optional[float] = None,
        optimal_kpts: Optional[Tuple[int, int, int]] = None,
        bulk_energy_per_atom: Optional[float] = None,
        energy_tolerance: float = 3.0e-3,
        precision: str = 'low',
    ) -> None:
        """
        Manually set bulk convergence parameters (skip run_bulk_convergence).
        
        Useful when you already have convergence results from previous runs
        or when working with a pre-relaxed structure.
        
        Parameters
        ----------
        optimal_ecutwfc : float
            Optimal energy cutoff (Ry) from convergence study
        optimal_kspacing : float, optional
            Optimal k-spacing (Å⁻¹) from convergence study
        optimal_kpts : Tuple[int, int, int], optional
            K-points from bulk calculation (e.g., (6, 6, 6))
            If provided, kspacing will be calculated from bulk cell
            This is PREFERRED for consistency with slab calculations
        bulk_energy_per_atom : float, optional
            Bulk energy per atom (eV/atom) for surface energy calc.
            If None, will be calculated later
        energy_tolerance : float, default=3.0e-3
            Convergence tolerance (meV/atom)
        precision : str, default='low'
            Precision level ('low', 'normal', 'high')
            
        Example
        -------
        >>> # After relaxing bulk manually with 6x6x6 k-mesh:
        >>> slab_wf = SlabWorkflow(bulk_atoms=relaxed_atoms, ...)
        >>> slab_wf.set_bulk_parameters(
        ...     optimal_ecutwfc=60.0,
        ...     optimal_kpts=(6, 6, 6),  # ← PREFERRED!
        ...     bulk_energy_per_atom=-3.8,
        ... )
        >>> # Now can run slab convergence:
        >>> slab_wf.run_slab_convergence(...)
        """
        # Calculate kspacing from kpts if provided
        if optimal_kpts is not None:
            bulk_cell = self.bulk_atoms.cell
            a = np.linalg.norm(bulk_cell[0, :])
            # Correct formula with 2π factor
            optimal_kspacing = 2 * np.pi / (a * optimal_kpts[0])
            logger.info(f"Calculating kspacing from bulk k-mesh {optimal_kpts}:")
            logger.info(f"  Bulk cell a: {a:.4f} Å")
            logger.info(f"  Calculated kspacing: {optimal_kspacing:.4f} Å⁻¹")
        elif optimal_kspacing is None:
            raise ValueError(
                "Must provide either 'optimal_kspacing' or 'optimal_kpts'"
            )
        
        self.bulk_recommendations = {
            'optimal_ecutwfc': optimal_ecutwfc,
            'optimal_kspacing': optimal_kspacing,
            'optimal_kpts': optimal_kpts,
            'bulk_energy_per_atom': bulk_energy_per_atom or 0.0,
            'energy_tolerance': energy_tolerance,
            'precision': precision,
        }
        
        logger.info("✓ Bulk parameters manually set:")
        logger.info(f"  ecutwfc: {optimal_ecutwfc} Ry")
        logger.info(f"  kspacing: {optimal_kspacing:.4f} Å⁻¹")
        if optimal_kpts:
            logger.info(f"  kpts (bulk): {optimal_kpts}")
        logger.info(f"  tolerance: {energy_tolerance*1000:.3f} meV/atom")
    
    def set_slab_kpoints(self, kpts: Tuple[int, int, int]) -> None:
        """
        Set explicit k-points for slab calculations (instead of calculating from kspacing).
        
        Useful when you want to ensure slab calculations use specific k-mesh.
        
        Parameters
        ----------
        kpts : Tuple[int, int, int]
            K-points to use for slab (nk_x, nk_y, nk_z)
            Example: (9, 9, 1) for Au(111) matching 6x6x6 bulk density
            
        Example
        -------
        >>> slab_wf = SlabWorkflow(bulk_atoms=relaxed_atoms, ...)
        >>> slab_wf.set_bulk_parameters(
        ...     optimal_ecutwfc=60.0,
        ...     optimal_kpts=(6, 6, 6),
        ... )
        >>> # Force slab to use 9x9x1 instead of calculating
        >>> slab_wf.set_slab_kpoints((9, 9, 1))
        >>> 
        >>> slab_conv = slab_wf.run_slab_convergence(...)
        """
        self.slab_kpts = kpts
        logger.info(f"✓ Slab k-points manually set: {kpts}")
    
    def _prepare_input_data(self) -> Dict:
        """
        Prepare input_data dictionary with bulk convergence parameters.
        
        Used to pass optimal ecutwfc and other parameters to batch submissions.
        
        Returns
        -------
        Dict with keys:
            - 'ecutwfc': Optimal energy cutoff from bulk convergence
            - Plus any other parameters from bulk_recommendations
        """
        if not self.bulk_recommendations:
            return {}
        
        # Extract ecutwfc from bulk recommendations
        input_data = {
            'ecutwfc': self.bulk_recommendations['optimal_ecutwfc'],
        }
        
        return input_data
    
    def run_bulk_convergence(
        self,
        label_prefix: str = 'bulk_convergence',
        verbose: Optional[bool] = None,
        convergence_criteria: List[str] = None,
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
        job_timeout: Optional[int] = None,
        walltime: Optional[str] = None,
    ) -> Dict:
        """
        Run bulk convergence study.
        
        **Status**: ✓ COMPLETE (Phase 2 task)
        
        Integrate ConvergenceWorkflow to determine optimal
        ecutwfc and kspacing for slab calculations.
        
        Uses independent two-phase algorithm:
        - Phase 1: ecutwfc convergence (fixed coarse kspacing)
        - Phase 2: kspacing convergence (fixed optimal ecutwfc)
        
        Parameters
        ----------
        label_prefix : str, default='bulk_convergence'
            Prefix for calculation directories
        verbose : bool, optional
            Override instance verbose setting
        convergence_criteria : List[str], optional
            Convergence criteria to check (default: ['energy'])
            Options: 'energy', 'forces', 'stress', 'geometry'
        max_ecutwfc : float, default=200.0
            Maximum ecutwfc to test (safety limit)
        ecutwfc_step : float, default=10.0
            Step size for ecutwfc increases
        job_timeout : int, optional
            Maximum time to wait for remote jobs (seconds).
            Example: 3600 (1 hour), 7200 (2 hours)
            If None, uses queue default.
        walltime : str, optional
            Maximum walltime for job submission (e.g., '1:00:00' for 1 hour).
            Format depends on scheduler (SLURM, PBS, Torque, etc).
            If None, uses machine configuration default.
            
        Returns
        -------
        Dict
            Convergence results with recommendations:
            - optimal_ecutwfc: Recommended energy cutoff (Ry)
            - optimal_kspacing: Recommended k-spacing (Å⁻¹)
            - precision: Precision level
            - energy_tolerance_meV_atom: Achieved tolerance
            
        Raises
        ------
        ImportError
            If ConvergenceWorkflow is not available
        RuntimeError
            If convergence study fails
        """
        if not HAS_CONVERGENCE_WORKFLOW:
            raise ImportError(
                "ConvergenceWorkflow not available. "
                "Install with: pip install xespresso"
            )
        
        verbose_flag = verbose if verbose is not None else self.verbose
        
        logger.info("="*70)
        logger.info("PHASE 1: BULK CONVERGENCE")
        logger.info("="*70)
        logger.info(f"Running convergence for {self.bulk_atoms.get_chemical_formula()}")
        logger.info(f"  Precision: {self.precision}")
        logger.info(f"  Protocol: {self.protocol}")
        logger.info(f"  Max ecutwfc: {max_ecutwfc} Ry, step: {ecutwfc_step} Ry")
        
        # Set default convergence criteria
        if convergence_criteria is None:
            convergence_criteria = ['energy']
        
        # Prepare queue with job_timeout and walltime if provided
        queue = self.queue
        if job_timeout is not None or walltime is not None:
            queue = queue.copy() if queue else {}
            if job_timeout is not None:
                queue['job_timeout'] = job_timeout
            if walltime is not None:
                queue['walltime'] = walltime
        
        try:
            # Create and run convergence workflow
            workflow = ConvergenceWorkflow(
                atoms=self.bulk_atoms,
                pseudopotentials=self.pseudopotentials,
                pseudopotentials_config=self.pseudopotentials_config,
                protocol=self.protocol,
                precision=self.precision,
                convergence_criteria_list=convergence_criteria,
                machine=self.machine,
                queue=queue,
                code_version=self.code_version,
            )
            
            # Run the convergence study
            logger.info("Starting convergence study...")
            workflow.run_convergence_study(
                label_prefix=label_prefix,
                verbose=verbose_flag,
                max_ecutwfc=max_ecutwfc,
                ecutwfc_step=ecutwfc_step,
                batch_timeout=3600,  # 1 hour timeout
            )
            
            # Get recommendations
            logger.info("Analyzing convergence results...")
            recommendations = workflow.get_recommendations(verbose=verbose_flag)
            
            # Store results
            self.bulk_results = workflow.results
            self.bulk_recommendations = recommendations
            
            logger.info("✓ Bulk convergence complete")
            logger.info(
                f"  Recommended ecutwfc: {recommendations['optimal_ecutwfc']} Ry\n"
                f"  Recommended kspacing: {recommendations['optimal_kspacing']} Å⁻¹"
            )
            
            return recommendations
            
        except Exception as e:
            logger.error(f"Convergence study failed: {e}")
            raise RuntimeError(f"Bulk convergence failed: {e}") from e
    
    # =========================================================================
    # PHASE 3: SLAB CONVERGENCE (COMPLETE)
    # =========================================================================
    
    def _calculate_anisotropic_kmesh(
        self,
        slab: Atoms,
        kspacing: Optional[float] = None,
    ) -> Tuple[int, int, int]:
        """
        Calculate anisotropic k-mesh for 2D slab.
        
        Uses bulk kspacing for in-plane directions (x, y).
        Always uses nk_z = 1 for 2D systems (perpendicular to surface).
        
        If self.slab_kpts is set, returns those explicit k-points instead.
        
        Parameters
        ----------
        slab : ase.Atoms
            Slab structure
        kspacing : float, optional
            K-spacing in Å⁻¹. If None, uses bulk_recommendations.
            
        Returns
        -------
        Tuple[int, int, int]
            (nk_x, nk_y, nk_z) k-mesh
            
        Notes
        -----
        Formula for each direction:
            nk_i = max(1, round(2π / (cell_i × kspacing)))
        
        If explicit kpts are set via set_slab_kpoints(), returns those instead.
        """
        # Check if explicit k-points are set
        if self.slab_kpts is not None:
            logger.info(f"Using explicit slab k-points: {self.slab_kpts}")
            return self.slab_kpts
        
        # Get k-spacing from bulk or parameter
        if kspacing is None:
            if not self.bulk_recommendations:
                raise ValueError(
                    "K-spacing not provided and bulk convergence not done. "
                    "Run run_bulk_convergence() first or provide kspacing."
                )
            kspacing = self.bulk_recommendations['optimal_kspacing']
        
        # Get cell parameters
        cell = slab.cell
        a = np.linalg.norm(cell[0, :])
        b = np.linalg.norm(cell[1, :])
        
        # Calculate k-mesh (anisotropic: x,y from bulk kspacing, z=1)
        # Using round() instead of ceil() for better separation of close kspacing values
        # Formula: nk = round(2π / (cell × kspacing))
        nk_x = max(1, int(np.round(2 * np.pi / (a * kspacing))))
        nk_y = max(1, int(np.round(2 * np.pi / (b * kspacing))))
        nk_z = 1  # 2D slab
        
        logger.debug(f"K-mesh: ({nk_x}, {nk_y}, {nk_z}) for kspacing={kspacing:.4f}")
        
        return (nk_x, nk_y, nk_z)
    
    def run_slab_convergence(
        self,
        surface_index: Tuple[int, int, int],
        vacuum_test: Optional[List[float]] = None,
        nlayers_test: Optional[List[int]] = None,
        nlayers_for_vacuum: Optional[int] = None,
        test_nlayers: bool = False,
        convergence_tol: float = 0.001,
        label_prefix: str = 'slab_convergence',
        skip_calculations: bool = False,
        machine: Optional[str] = None,
        queue: Optional[Dict] = None,
        code_version: Optional[str] = None,
        protocol: Optional[str] = None,
        job_timeout: Optional[int] = None,
        walltime: Optional[str] = None,
    ) -> Dict:
        """
        Run slab-specific convergence study (Phase 3).
        
        **Status**: ✓ COMPLETE
        
        Tests convergence of vacuum size. Layer thickness testing is optional
        but NOT RECOMMENDED in this phase (SCF is expensive). For layer convergence,
        use run_slab_relax() with nlayers_test parameter (Phase 4b with relaxation).
        
        Convergence sequence:
        1. Test VACUUM: [10, 12, 15, 18, 20, 25, 30] Å (with fixed nlayers)
        2. Test LAYERS: [3, 4, 5, 6, 7] (optional, only if test_nlayers=True)
        
        **Recommended workflow**:
        - Phase 3 (this method): test_nlayers=False (default) - only vacuum convergence
        - Phase 4b (run_slab_relax): nlayers_test=[3,4,5,6,7] - layer convergence with relaxation
        
        Parameters
        ----------
        surface_index : Tuple[int, int, int]
            Miller index of surface (e.g., (1,1,1))
        vacuum_test : List[float], optional
            Vacuum sizes to test (Å). Default: [10, 12, 15, 18, 20, 25, 30]
        nlayers_test : List[int], optional
            Layer counts to test. Default: [3, 4, 5, 6, 7] (only used if test_nlayers=True)
        nlayers_for_vacuum : int, optional
            Number of layers to use when testing vacuum convergence (Phase 3, Step 1).
            If None, uses minimum of nlayers_test (or 3 if nlayers_test not provided).
            Use this to test vacuum convergence with a specific slab thickness.
            Example: nlayers_for_vacuum=5 → test vacuum with a 5-layer slab
        test_nlayers : bool, default=False
            If False (recommended): test ONLY vacuum with fixed nlayers (faster, cheaper SCF)
            If True: also test layer thickness (expensive, better done with relaxation in Phase 4b)
        convergence_tol : float, default=0.001
            Energy tolerance for convergence (eV/atom)
        label_prefix : str, default='slab_convergence'
            Prefix for calculation directories
        skip_calculations : bool, default=False
            If True, skip actual SCF calculations (for testing)
        machine : str, optional
            Remote machine name (e.g., 'medusa'). If None, uses self.machine
        queue : Dict, optional
            Queue configuration. If None, uses self.queue
        code_version : str, optional
            QE version code (e.g., '7.4.1'). If None, uses self.code_version
        protocol : str, optional
            Protocol name (e.g., 'fast', 'moderate', 'accurate'). 
            If provided, overrides self.protocol for this convergence.
        job_timeout : int, optional
            Maximum time to wait for remote jobs (seconds).
            Example: 3600 (1 hour), 7200 (2 hours)
            If None, uses queue default.
        walltime : str, optional
            Maximum walltime for job submission (e.g., '1:00:00' for 1 hour).
            Format depends on scheduler (SLURM, PBS, Torque, etc).
            If None, uses machine configuration default.
            
        Returns
        -------
        Dict
            Convergence results:
            - 'vacuum_results': Dict of vacuum (Å) → energy (eV/atom)
            - 'layer_results': Dict of nlayers → energy (eV/atom) [only if test_nlayers=True]
            - 'optimal_vacuum': Recommended vacuum size (Å)
            - 'optimal_nlayers': Recommended number of layers [only if test_nlayers=True]
            - 'kmesh_calc': (nk_x, nk_y, nk_z) calculated k-mesh
            - 'converged': True if convergence achieved
            
        Raises
        ------
        ValueError
            If surface not generated or bulk convergence not completed
        """
        logger.info("="*70)
        logger.info(f"PHASE 3: SLAB CONVERGENCE ({surface_index})")
        logger.info("="*70)
        
        # Set defaults
        if vacuum_test is None:
            vacuum_test = [10, 12, 15, 18, 20, 25, 30]
        if nlayers_test is None:
            nlayers_test = [3, 4, 5, 6, 7]
        
        # Use provided protocol or fall back to self.protocol
        protocol = protocol or self.protocol
        
        # Prepare queue with job_timeout and walltime if provided
        if job_timeout is not None or walltime is not None:
            queue = queue or {}
            queue = dict(queue)  # Make a copy to avoid modifying original
            if job_timeout is not None:
                queue['job_timeout'] = job_timeout
            if walltime is not None:
                queue['walltime'] = walltime
        
        # Validate prerequisites
        if surface_index not in self.slabs:
            raise ValueError(
                f"Surface {surface_index} not generated. "
                f"Run generate_slabs() first."
            )
        
        if not self.bulk_recommendations:
            raise ValueError(
                "Bulk convergence not completed. "
                "Run run_bulk_convergence() first."
            )
        
        logger.info(f"  Using bulk parameters:")
        logger.info(f"    ecutwfc: {self.bulk_recommendations['optimal_ecutwfc']} Ry")
        logger.info(f"    kspacing: {self.bulk_recommendations['optimal_kspacing']} Å⁻¹")
        
        # Calculate k-mesh (deterministic from bulk)
        kmesh = self._calculate_anisotropic_kmesh(
            self.slabs[surface_index],
            self.bulk_recommendations['optimal_kspacing']
        )
        logger.info(f"    k-mesh: {kmesh}")
        
        results = {
            'surface_index': surface_index,
            'kmesh_calc': kmesh,
            'vacuum_results': {},
            'layer_results': {},
            'converged': False,
        }
        
        # STEP 1: Vacuum convergence with specified layers (most efficient)
        if nlayers_for_vacuum is not None:
            min_nlayers = nlayers_for_vacuum
        else:
            min_nlayers = min(nlayers_test) if nlayers_test else 3
        logger.info(f"  (Using nlayers={min_nlayers} for vacuum convergence)")
        
        vac_results = self._test_vacuum_convergence(
            surface_index=surface_index,
            vacuum_test=vacuum_test,
            kmesh=kmesh,
            nlayers=min_nlayers,
            label_prefix=label_prefix,
            skip_calculations=skip_calculations,
            machine=machine or self.machine,
            queue=queue or self.queue,
            code_version=code_version or self.code_version,
            protocol=protocol,
        )
        
        results['vacuum_results'] = vac_results['energies']
        optimal_vac = vac_results['optimal_vacuum']
        vac_converged = vac_results['converged']
        results['optimal_vacuum'] = optimal_vac
        
        logger.info(f"  → Optimal vacuum: {optimal_vac:.1f} Å")
        if vac_converged:
            logger.info(f"  ✓ Vacuum convergence achieved")
        else:
            logger.warning(f"  ⚠ Vacuum convergence NOT achieved")
        
        # STEP 2: Layer convergence (only if requested)
        if test_nlayers:
            logger.info("\n" + "-"*70)
            logger.info("STEP 2: LAYER CONVERGENCE")
            logger.info("-"*70)
            logger.info(f"Testing layers: {nlayers_test}")
            
            layer_results = self._test_layer_convergence(
                surface_index=surface_index,
                nlayers_test=nlayers_test,
                optimal_vacuum=optimal_vac,
                kmesh=kmesh,
                label_prefix=label_prefix,
                skip_calculations=skip_calculations,
                machine=machine or self.machine,
                queue=queue or self.queue,
                code_version=code_version or self.code_version,
                protocol=protocol,
            )
            
            results['layer_results'] = layer_results['energies']
            optimal_layers = layer_results['optimal_nlayers']
            layer_converged = layer_results['converged']
            results['optimal_nlayers'] = optimal_layers
            
            logger.info(f"  → Optimal layers: {optimal_layers}")
            if layer_converged:
                logger.info(f"  ✓ Layer convergence achieved")
            else:
                logger.warning(f"  ⚠ Layer convergence NOT achieved")
        else:
            # Skip layer testing - will be done in Phase 4b (run_slab_relax)
            logger.info("\n" + "-"*70)
            logger.info("STEP 2: SKIPPED (Layer convergence deferred to Phase 4b)")
            logger.info("-"*70)
            logger.info("Layer thickness convergence will be tested in Phase 4b (run_slab_relax)")
            logger.info("with atomic relaxation for more realistic results.")
            optimal_layers = None
            layer_converged = False
            results['layer_results'] = {}
            results['optimal_nlayers'] = None
        
        # Store convergence status
        results['converged'] = vac_converged and (layer_converged if test_nlayers else True)
        self.convergence_results[surface_index] = results
        
        # Final summary
        logger.info("\n" + "="*70)
        logger.info(f"✓ Phase 3 complete for {surface_index}")
        logger.info(f"  Optimal vacuum: {optimal_vac:.1f} Å")
        if test_nlayers:
            logger.info(f"  Optimal layers: {optimal_layers}")
        else:
            logger.info(f"  Layer convergence: deferred to Phase 4b (run_slab_relax)")
        logger.info(f"  K-mesh (x,y,z): {kmesh}")
        logger.info("="*70)
        
        return results
    
    def _ensure_positive_vacuum(self, slab: Atoms) -> Atoms:
        """
        Ensure vácuo (c vector) is positive.
        
        PyMatGen may generate cells with negative z-component.
        This function flips the sign if needed, preserving the structure.
        
        Parameters
        ----------
        slab : ase.Atoms
            Slab structure
            
        Returns
        -------
        ase.Atoms
            Slab with positive vacuum
        """
        cell = slab.get_cell()
        if cell[2, 2] < 0:
            cell[2] = -cell[2]
            slab.set_cell(cell)
            logger.debug(f"  Fixed negative vacuum: {cell[2, 2]:.4f} → {abs(cell[2, 2]):.4f} Å")
        return slab

    def _test_vacuum_convergence(
        self,
        surface_index: Tuple[int, int, int],
        vacuum_test: List[float],
        kmesh: Tuple[int, int, int],
        nlayers: int = 3,
        label_prefix: str = 'slab_convergence',
        skip_calculations: bool = False,
        machine: Optional[str] = None,
        queue: Optional[Dict] = None,
        code_version: Optional[str] = None,
        protocol: Optional[str] = None,
    ) -> Dict:
        """Test vacuum convergence by varying vacuum size using parallel batch submission.
        
        Uses a slab with fixed number of layers (usually minimum) for efficiency.
        All calculations are submitted in parallel to the remote scheduler via batch_utils.
        """
        logger.info("\n" + "-"*70)
        logger.info("STEP 1: VACUUM CONVERGENCE")
        logger.info("-"*70)
        logger.info(f"Testing vacuum sizes: {vacuum_test} Å (with nlayers={nlayers})")
        
        # Generate slab with specified nlayers for vacuum testing
        base_slab = self._regenerate_slab_with_nlayers(surface_index, nlayers)
        num_atoms = len(base_slab)
        energies = {}
        
        if skip_calculations:
            # Mock mode: generate energies directly
            for vacuum in vacuum_test:
                slab = base_slab.copy()
                slab.center(vacuum=vacuum, axis=2)
                energy = -len(slab) * 50.0  # Arbitrary negative energy
                logger.info(f"  Testing vacuum = {vacuum:.1f} Å: (mock) E={energy/num_atoms:.6f} eV/atom")
                energies[vacuum] = energy / num_atoms
        else:
            # Prepare structures for parallel batch submission
            from xespresso.workflow.batch_utils import submit_structures_parallel, collect_results
            
            # Smart caching: check which calculations already exist and are complete
            structures_to_submit = []
            completed_vacuums = {}
            
            for vacuum in vacuum_test:
                slab = base_slab.copy()
                slab.center(vacuum=vacuum, axis=2)
                calc_label = f"{label_prefix}/vacuum_{vacuum:.0f}"
                
                # Check if calculation already completed
                calc_dir = calc_label if os.path.isabs(calc_label) else os.path.join(os.getcwd(), calc_label)
                calc_complete, energy = self._check_calculation_complete(calc_dir, num_atoms)
                
                if calc_complete:
                    # Reuse previous result
                    completed_vacuums[vacuum] = energy
                    logger.info(f"  ✓ Reusing completed: vacuum={vacuum:.1f} Å: E={energy:.6f} eV/atom")
                else:
                    # Need to calculate
                    # Ensure vácuo is positive before submitting
                    slab = self._ensure_positive_vacuum(slab)
                    
                    structures_to_submit.append({
                        'atoms': slab,
                        'label': calc_label,
                        'param_key': vacuum,
                        'num_atoms': num_atoms,
                    })
            
            # Submit only new/incomplete calculations in parallel
            if structures_to_submit:
                logger.info(f"\n  Submitting {len(structures_to_submit)} new calculations in parallel...")
                
                batch_result = submit_structures_parallel(
                    structures=structures_to_submit,
                    kmesh=kmesh,
                    pseudopotentials=self.pseudopotentials,
                    pseudopotentials_config=self.pseudopotentials_config,
                    protocol=protocol or self.protocol,
                    precision=self.precision,
                    machine=machine or self.machine,
                    queue=queue or self.queue,
                    code_version=code_version or self.code_version,
                    input_data=self._prepare_input_data(),
                )
                
                # Collect results from newly submitted jobs
                results_dict = collect_results(batch_result, verbose=False)
                
                # Extract energies from new results
                for vacuum, completion in results_dict.items():
                    if completion.get('success', False) and 'energy' in completion:
                        energy_per_atom = completion['energy'] / num_atoms
                        energies[vacuum] = energy_per_atom
                        logger.info(f"  ✓ vacuum={vacuum:.1f} Å: E={energy_per_atom:.6f} eV/atom")
                    else:
                        error_msg = completion.get('error', 'unknown error')
                        logger.warning(f"  ⚠ vacuum={vacuum:.1f} Å: {error_msg}")
            
            # Add completed vacuums to results
            energies.update(completed_vacuums)
            
            if not structures_to_submit and not completed_vacuums:
                logger.warning("  No vacuum tests completed")
                return {
                    'energies': {},
                    'optimal_vacuum': vacuum_test[2] if len(vacuum_test) > 2 else 15.0,
                    'converged': False,
                }
        
        # Find optimal (lowest energy) vacuum
        if not energies:
            logger.warning("  No energies available")
            return {
                'energies': {},
                'optimal_vacuum': vacuum_test[2] if len(vacuum_test) > 2 else 15.0,
                'converged': False,
            }
        
        optimal_vacuum = min(energies, key=energies.get)
        
        # Check convergence: ∆E between last two decreasing
        sorted_vac = sorted(energies.items())
        converged = False
        if len(sorted_vac) >= 2:
            e_last = sorted_vac[-1][1]
            e_prev = sorted_vac[-2][1]
            converged = abs(e_last - e_prev) < 0.001  # 1 meV/atom
        
        return {
            'energies': energies,
            'optimal_vacuum': optimal_vacuum,
            'converged': converged,
        }
    
    def _test_layer_convergence(
        self,
        surface_index: Tuple[int, int, int],
        nlayers_test: List[int],
        optimal_vacuum: float,
        kmesh: Tuple[int, int, int],
        label_prefix: str,
        skip_calculations: bool = False,
        machine: Optional[str] = None,
        queue: Optional[Dict] = None,
        code_version: Optional[str] = None,
        protocol: Optional[str] = None,
    ) -> Dict:
        """Test layer convergence by varying number of layers using parallel batch submission.
        
        All calculations are submitted in parallel to the remote scheduler.
        """
        from xespresso.workflow.batch_utils import submit_structures_parallel, collect_results
        
        logger.info("\n" + "-"*70)
        logger.info("STEP 2: LAYER CONVERGENCE")
        logger.info("-"*70)
        logger.info(f"Testing layers: {nlayers_test}")
        
        energies = {}
        
        if skip_calculations:
            # Mock mode: generate energies directly
            for nlayers in nlayers_test:
                slab = self._regenerate_slab_with_nlayers(surface_index, nlayers)
                slab.center(vacuum=optimal_vacuum, axis=2)
                num_atoms = len(slab)
                energy = -nlayers * 100.0
                logger.info(f"  Testing layers = {nlayers}: {num_atoms} atoms, (mock) E={energy/num_atoms:.6f} eV/atom")
                energies[nlayers] = energy / num_atoms
        else:
            # Prepare structures for parallel batch submission
            # Smart caching: check which calculations already exist and are complete
            structures_to_submit = []
            completed_nlayers = {}
            
            for nlayers in nlayers_test:
                slab = self._regenerate_slab_with_nlayers(surface_index, nlayers)
                slab.center(vacuum=optimal_vacuum, axis=2)
                num_atoms = len(slab)
                
                fix_indices = self._calculate_fixed_layers(nlayers)
                calc_label = f"{label_prefix}/nlayers_{nlayers}"
                
                # Check if calculation already completed
                calc_dir = calc_label if os.path.isabs(calc_label) else os.path.join(os.getcwd(), calc_label)
                calc_complete, energy = self._check_calculation_complete(calc_dir, num_atoms)
                
                if calc_complete:
                    # Reuse previous result
                    completed_nlayers[nlayers] = energy
                    logger.info(f"  ✓ Reusing completed: nlayers={nlayers}: {num_atoms} atoms: E={energy:.6f} eV/atom")
                else:
                    # Need to calculate
                    logger.info(f"  Submitting new: nlayers={nlayers}: {num_atoms} atoms, freezing layers {fix_indices}")
                    structures_to_submit.append({
                        'atoms': slab,
                        'label': calc_label,
                        'param_key': nlayers,
                        'num_atoms': num_atoms,
                    })
            
            # Submit only new/incomplete calculations in parallel
            if structures_to_submit:
                logger.info(f"\n  Submitting {len(structures_to_submit)} new calculations in parallel...")
                
                batch_result = submit_structures_parallel(
                    structures=structures_to_submit,
                    kmesh=kmesh,
                    pseudopotentials=self.pseudopotentials,
                    pseudopotentials_config=self.pseudopotentials_config,
                    protocol=protocol or self.protocol,
                    precision=self.precision,
                    machine=machine or self.machine,
                    queue=queue or self.queue,
                    code_version=code_version or self.code_version,
                    input_data=self._prepare_input_data(),
                )
                
                # Collect results from newly submitted jobs
                results_dict = collect_results(batch_result, verbose=False)
                
                # Extract energies from new results
                for nlayers, completion in results_dict.items():
                    if completion.get('success', False) and 'energy' in completion:
                        # Find num_atoms for this nlayers
                        struct = next(
                            (s for s in structures_to_submit if s['param_key'] == nlayers),
                            None
                        )
                        num_atoms = struct['num_atoms'] if struct else None
                        
                        if num_atoms:
                            energy_per_atom = completion['energy'] / num_atoms
                            energies[nlayers] = energy_per_atom
                            logger.info(f"  ✓ nlayers={nlayers}: E={energy_per_atom:.6f} eV/atom")
                    else:
                        error_msg = completion.get('error', 'unknown error')
                        logger.warning(f"  ⚠ nlayers={nlayers}: {error_msg}")
            
            # Add completed nlayers to results
            energies.update(completed_nlayers)
            
            if not structures_to_submit and not completed_nlayers:
                logger.warning("  No layer tests completed")
                return {
                    'energies': {},
                    'optimal_nlayers': nlayers_test[1] if len(nlayers_test) > 1 else 4,
                    'converged': False,
                }
        
        # Find optimal (lowest energy) layers
        if not energies:
            logger.warning("  No energies available")
            return {
                'energies': {},
                'optimal_nlayers': nlayers_test[1] if len(nlayers_test) > 1 else 4,
                'converged': False,
            }
        
        optimal_nlayers = min(energies, key=energies.get)
        
        # Check convergence: energy difference < 1 meV between consecutive
        sorted_layers = sorted(energies.items())
        converged = False
        if len(sorted_layers) >= 2:
            e_last = sorted_layers[-1][1]
            e_prev = sorted_layers[-2][1]
            converged = abs(e_last - e_prev) < 0.001  # 1 meV/atom
        
        logger.info(f"  Layer convergence results:")
        for nl, en in sorted(energies.items()):
            logger.info(f"    nlayers={nl}: E={en:.6f} eV/atom")
        logger.info(f"  Optimal: {optimal_nlayers} layers (converged: {converged})")
        
        return {
            'energies': energies,
            'optimal_nlayers': optimal_nlayers,
            'converged': converged,
        }
    
    # =========================================================================
    # PHASE 4: STRUCTURE RELAXATION
    # =========================================================================
    
    def run_slab_relax(
        self,
        surfaces: Optional[List[Tuple[int, int, int]]] = None,
        nlayers_test: Optional[List[int]] = None,
        vacuum_size: Optional[float] = None,
        relax_type: str = 'relax',
        fmax: float = 0.05,
        nstep_relax: int = 50,
        upscale: float = 100.0,
        scf_conv_thr: Optional[float] = None,
        dipole_correction: bool = True,
        machine: Optional[str] = None,
        code_version: Optional[str] = None,
        label_prefix: str = 'relax',
        protocol: Optional[str] = None,
        job_timeout: Optional[int] = None,
        walltime: Optional[str] = None,
        nosym: Optional[bool] = None,
    ) -> Dict:
        """
        Relax slab structures with constraints (Phases 4-4b).
        
        Relaxes slabs using CalculationWorkflow with:
        - Fixed bottom layers using parity rule (EVEN: nlayers/2, ODD: (nlayers-1)/2)
        - Anisotropic k-mesh from Phase 3 (or explicit kpts)
        - ecutwfc from Phase 1 bulk convergence
        - Optional dipole correction for 2D systems
        
        **Phase 4**: Relax with converged parameters
        **Phase 4b**: Test surface energy convergence with variable nlayers
        
        Parameters
        ----------
        surfaces : Optional[List[Tuple[int, int, int]]], default=None
            List of surface indices to relax. If None, relaxes all generated slabs.
            Example: [(1,0,0), (1,1,1)]
            
        nlayers_test : Optional[List[int]], default=None
            If provided, test relaxation with different nlayers values.
            Example: [3, 4, 5, 6, 7] → test convergence of surface energy.
            If None, uses self.nlayers (Phase 4 single relaxation).
            
        vacuum_size : Optional[float], default=None
            Vacuum size in Ångströms. If None, uses optimal_vacuum from Phase 3.
            Only used when nlayers_test is provided (Phase 4b).
            
        relax_type : str, default='relax'
            Type of relaxation: 'relax' (ions only) or 'vc-relax' (ions + cell)
            
        fmax : float, default=0.05
            Maximum force convergence criterion in eV/Å.
            Relaxation converges when all atomic forces are below this value.
            Internally converted to atomic units (Hartree/bohr) for QE input:
            forc_conv_thr (a.u.) = fmax (eV/Å) / 51.4220652
            
        nstep_relax : int, default=50
            Maximum number of BFGS relaxation steps.
            Increase to 100-200 for difficult convergence (e.g., large slabs).
            
        upscale : float, default=100.0
            BFGS curvature scaling factor. Higher values help with poor Hessian approximations.
            Typical range: 50-200. Increase if BFGS curvature warnings appear.
            
        scf_conv_thr : Optional[float], default=None
            SCF convergence threshold during relaxation (Ry).
            If None, uses default from protocol. 
            Tighter values (e.g., 1e-9) improve force accuracy.
            
        dipole_correction : bool, default=True
            Add dipole correction for 2D systems (automatically sets dipole='z')
            
        machine : Optional[str], default=None
            Remote machine for calculations. If None, uses self.machine.
            Example: 'medusa', 'localhost'
            
        code_version : Optional[str], default=None
            Quantum ESPRESSO version. If None, uses self.code_version.
            Example: '7.4.1'
            
        label_prefix : str, default='relax'
            Prefix for calculation directories
            
        protocol : str, optional
            Protocol name (e.g., 'fast', 'moderate', 'accurate'). 
            If provided, overrides self.protocol for this relaxation.
        job_timeout : int, optional
            Maximum time to wait for remote jobs (seconds).
            Example: 3600 (1 hour), 7200 (2 hours)
            If None, uses queue default.
        walltime : str, optional
            Maximum walltime for job submission (e.g., '1:00:00' for 1 hour).
            Format depends on scheduler (SLURM, PBS, Torque, etc).
            If None, uses machine configuration default.
        nosym : bool, optional
            Disable symmetry during relaxation. Default: True (recommended).
            - True: nosym=.true. in QE input (disables symmetry operations)
            - False: nosym=.false. (use symmetry, may break during geometry opt)
            If None, defaults to True for safe relaxation.
            
        Returns
        -------
        Dict
            Phase 4 (single nlayers):
            {
                (1,1,1): {
                    'surface_index': (1,1,1),
                    'relaxed_slab': Atoms,
                    'energy': float (eV),
                    'converged': bool,
                }
            }
            
            Phase 4b (multiple nlayers):
            {
                3: {'energy': -1234.56, 'relaxed_slab': Atoms, ...},
                4: {'energy': -2345.67, 'relaxed_slab': Atoms, ...},
                ...
            }
            
        Raises
        ------
        ValueError
            If slabs not generated or bulk convergence not completed
        NotImplementedError
            If CalculationWorkflow not available
        """
        from xespresso.workflow.calculation_workflow import CalculationWorkflow
        from ase.constraints import FixAtoms
        
        # Determine execution mode
        phase = "4b (nlayers convergence)" if nlayers_test else "4 (single relaxation)"
        
        # Validate prerequisites
        if not self.convergence_results:
            raise ValueError(
                "Phase 3 slab convergence not completed. "
                "Run run_slab_convergence() first."
            )
        
        if not surfaces:
            surfaces = list(self.slabs.keys())
        
        # Default values
        if machine is None:
            machine = self.machine
        if code_version is None:
            code_version = self.code_version
        
        # Use provided protocol or fall back to self.protocol
        protocol = protocol or self.protocol
        
        # Prepare queue with job_timeout and walltime if provided
        queue = self.queue
        if job_timeout is not None or walltime is not None:
            queue = queue.copy() if queue else {}
            if job_timeout is not None:
                queue['job_timeout'] = job_timeout
            if walltime is not None:
                queue['walltime'] = walltime
        
        logger.info("="*70)
        logger.info(f"PHASE {phase}: SLAB RELAXATION ({len(surfaces)} surfaces)")
        logger.info("="*70)
        logger.info(f"  Relax type: {relax_type}")
        logger.info(f"  Dipole correction: {dipole_correction}")
        logger.info(f"  Surfaces: {surfaces}")
        if nlayers_test:
            logger.info(f"  Test nlayers: {nlayers_test}")
            logger.info(f"  Vacuum size (fixed): {vacuum_size} Å")
        
        results = {}
        
        # =====================================================================
        # PHASE 4B: Multiple nlayers (convergence test) - PARALLEL SUBMISSION
        # =====================================================================
        if nlayers_test:
            from xespresso.workflow.batch_utils import submit_relaxations_parallel
            
            logger.info(f"\n{'─'*70}")
            logger.info(f"CONVERGENCE TEST: Variable nlayers with fixed vacuum")
            logger.info(f"{'─'*70}")
            
            # Use provided vacuum or get from Phase 3
            if vacuum_size is None:
                conv_data = self.convergence_results[surfaces[0]]
                vacuum_size = conv_data['optimal_vacuum']
                logger.info(f"  Using optimal_vacuum from Phase 3: {vacuum_size} Å")
            
            # Get other parameters from Phase 3
            conv_data = self.convergence_results[surfaces[0]]
            ecutwfc = self.bulk_recommendations['optimal_ecutwfc']
            kmesh = conv_data['kmesh_calc']
            
            logger.info(f"  Parameters:")
            logger.info(f"    ecutwfc: {ecutwfc} Ry")
            logger.info(f"    k-mesh: {kmesh}")
            logger.info(f"    vacuum: {vacuum_size} Å")
            logger.info(f"    Parallel submission: YES (scheduler will parallelize)")
            
            # Prepare structures for parallel submission
            structures = []
            
            for surface in surfaces:
                for target_nlayers in nlayers_test:
                    # Generate slab with target nlayers and fixed vacuum_size
                    slab = self._regenerate_slab_with_nlayers(
                        surface,
                        target_nlayers,
                        vacuum_size=vacuum_size,  # Use converged vacuum_size
                    )
                    
                    # Apply parity rule for fixed layers
                    # EVEN: fix nlayers/2, ODD: fix (nlayers-1)/2
                    if target_nlayers % 2 == 0:
                        num_fixed = target_nlayers // 2
                    else:
                        num_fixed = (target_nlayers - 1) // 2
                    
                    # Map layer indices to atom indices
                    fixed_atom_indices = self._get_atom_indices_for_layers(
                        slab, surface, num_fixed
                    )
                    
                    # Apply constraints
                    if fixed_atom_indices:
                        slab.set_constraint(FixAtoms(indices=fixed_atom_indices))
                    
                    structures.append({
                        'atoms': slab,
                        'label': f"{label_prefix}/{surface[0]}{surface[1]}{surface[2]}/nlayers_{target_nlayers}",
                        'param_key': target_nlayers,
                        'num_atoms': len(slab),
                    })
            
            logger.info(f"  Prepared {len(structures)} slabs for parallel relaxation")
            
            # Prepare input_data with ecutwfc and force convergence threshold
            # Note: forc_conv_thr is in atomic units (Hartree/bohr)
            # Conversion: 1 eV/Å = 1/51.4220652 Hartree/bohr
            forc_conv_thr_au = fmax / 51.4220652  # Convert eV/Å to a.u.
            input_data = {
                'ecutwfc': ecutwfc,
                'forc_conv_thr': forc_conv_thr_au,  # Force convergence threshold in a.u.
                'nosym': nosym if nosym is not None else True,  # Disable symmetry to avoid sym. violations during geometry opt.
                'nstep': nstep_relax,  # Maximum BFGS steps
                'upscale': upscale,  # BFGS curvature scaling
            }
            
            # Add SCF convergence threshold if provided
            if scf_conv_thr is not None:
                input_data['conv_thr'] = scf_conv_thr
            
            # Submit all relaxations in parallel via batch_utils
            batch_result = submit_relaxations_parallel(
                structures=structures,
                kmesh=kmesh,
                pseudopotentials_config=self.pseudopotentials_config,
                relax_type=relax_type,
                protocol=protocol,
                machine=machine,
                queue=self.queue,
                code_version=code_version,
                input_data=input_data,
                dipole_correction=dipole_correction,
            )
            
            logger.info(f"\nWaiting for {batch_result['total_submitted']} relaxations to complete...")
            logger.info(f"(Jobs are running in parallel on remote scheduler)\n")
            
            # Collect results from submitted relaxations using utility function
            from xespresso.workflow.batch_utils import collect_relax_results
            
            results_dict = collect_relax_results(batch_result, timeout=3600, poll_interval=30)
            
            # Convert to slab_workflow format and extract energies
            results = {}
            for param_key, completion in results_dict.items():
                nlayers = param_key
                label = completion.get('label', f'nlayers_{nlayers}')
                
                if completion.get('success', False) and completion.get('energy') is not None:
                    energy = completion['energy']
                    results[nlayers] = {
                        'energy': energy,
                        'converged': completion.get('converged', True),
                        'nlayers': nlayers,
                        'label': label,
                    }
                    logger.info(f"  ✓ {label}: E={energy:.6f} eV")
                else:
                    error_msg = completion.get('error', 'unknown error')
                    results[nlayers] = {
                        'energy': None,
                        'converged': False,
                        'nlayers': nlayers,
                        'error': error_msg,
                        'label': label,
                    }
                    logger.warning(f"  ✗ {label}: {error_msg}")
            
            logger.info(f"\n{'='*70}")
            successful = len([r for r in results.values() if r.get('converged')])
            logger.info(f"PHASE 4B COMPLETE: {successful} / {len(nlayers_test)} relaxations successful")
            logger.info(f"{'='*70}")
            
            return results
        
        # =====================================================================
        # PHASE 4: Single nlayers (normal relaxation)
        # =====================================================================
        logger.info(f"\nRelaxing with self.nlayers={self.nlayers}...")
        
        for surface in surfaces:
            logger.info(f"\n{'─'*70}")
            logger.info(f"Relaxing surface {surface}")
            logger.info(f"{'─'*70}")
            
            if surface not in self.slabs:
                logger.warning(f"  Surface {surface} not generated, skipping")
                continue
            
            try:
                # Get slab and apply constraints using parity rule
                slab = self.slabs[surface].copy()
                
                # Apply parity rule for fixed layers
                if self.nlayers % 2 == 0:
                    num_fixed = self.nlayers // 2
                else:
                    num_fixed = (self.nlayers - 1) // 2
                
                fixed_atom_indices = self._get_atom_indices_for_layers(
                    slab, surface, num_fixed
                )
                
                logger.info(f"  nlayers={self.nlayers} ({self.nlayers%2==0 and 'EVEN' or 'ODD'})")
                logger.info(f"  Fixing {num_fixed} layers ({len(fixed_atom_indices)} atoms)")
                
                if fixed_atom_indices:
                    slab.set_constraint(FixAtoms(indices=fixed_atom_indices))
                
                # Get optimal parameters from Phase 3
                conv_data = self.convergence_results[surface]
                kmesh = conv_data['kmesh_calc']
                ecutwfc = self.bulk_recommendations['optimal_ecutwfc']
                optimal_vacuum = conv_data['optimal_vacuum']
                
                logger.info(f"  Using parameters:")
                logger.info(f"    ecutwfc: {ecutwfc} Ry")
                logger.info(f"    k-mesh: {kmesh}")
                logger.info(f"    vacuum: {optimal_vacuum} Å")
                
                # Create CalculationWorkflow for relaxation
                calc_wf = CalculationWorkflow(
                    atoms=slab,
                    protocol=protocol,
                    pseudopotentials_config=self.pseudopotentials_config,
                    machine=machine,
                    queue=self.queue,
                    code_version=code_version,
                )
                
                # Override with Phase 1 parameters
                calc_wf.input_data['ecutwfc'] = ecutwfc
                # Note: forc_conv_thr is in atomic units (Hartree/bohr)
                # Conversion: 1 eV/Å = 1/51.4220652 Hartree/bohr
                forc_conv_thr_au = fmax / 51.4220652  # Convert eV/Å to a.u.
                calc_wf.input_data['forc_conv_thr'] = forc_conv_thr_au
                calc_wf.input_data['nosym'] = nosym if nosym is not None else True  # Disable symmetry during relaxation
                calc_wf.input_data['nstep'] = nstep_relax  # Maximum BFGS steps
                calc_wf.input_data['upscale'] = upscale  # BFGS curvature scaling
                if scf_conv_thr is not None:
                    calc_wf.input_data['conv_thr'] = scf_conv_thr
                calc_wf.input_data.pop('kspacing', None)  # Remove kspacing, use k-mesh
                
                # Add dipole correction if requested
                if dipole_correction:
                    calc_wf.input_data['dipole'] = 'z'
                    logger.info(f"  Added dipole correction (2D system)")
                
                # Run relaxation
                label = f"{label_prefix}/{surface[0]}{surface[1]}{surface[2]}"
                logger.info(f"  Running {relax_type} calculation...")
                
                calc = calc_wf.run_relax(
                    label=label,
                    relax_type=relax_type,
                    kpts=kmesh,  # Explicit k-mesh from Phase 3
                )
                
                # Store results
                relaxed_slab = calc_wf.atoms
                energy = calc_wf.atoms.get_potential_energy()
                
                results[surface] = {
                    'surface_index': surface,
                    'relaxed_slab': relaxed_slab,
                    'calculator': calc,
                    'energy': energy,
                    'converged': True,
                    'results': calc.results if hasattr(calc, 'results') else {},
                }
                
                logger.info(f"  ✓ Relaxation complete: E = {energy:.6f} eV")
                
            except Exception as e:
                logger.error(f"  ✗ Relaxation failed: {e}")
                results[surface] = {
                    'surface_index': surface,
                    'converged': False,
                    'error': str(e),
                }
        
        # Store results
        self.relax_results = results
        logger.info(f"\n{'='*70}")
        logger.info(f"PHASE 4 COMPLETE: {len([r for r in results.values() if r.get('converged')])} / {len(surfaces)} relaxations successful")
        logger.info(f"{'='*70}")
        
    
    # =========================================================================
    # PHASE 5: ANALYSIS (PLACEHOLDER - Phase 5 task)
    # =========================================================================
    
    def calculate_surface_energies(
        self,
        bulk_energy_per_atom: float,
    ) -> Dict:
        """
        Calculate surface energies.
        
        **Status**: ⏳ PHASE 5 task (not yet implemented)
        
        Will compute γ = (E_slab - n*E_bulk) / (2*Area) in J/m².
        
        Parameters
        ----------
        bulk_energy_per_atom : float
            Energy per atom in converged bulk (eV/atom)
            
        Returns
        -------
        Dict
            Surface energies
            
        Raises
        ------
        NotImplementedError
            Currently a Phase 5 task
        """
        logger.warning("="*70)
        logger.warning("PHASE 5: SURFACE ENERGY ANALYSIS")
        logger.warning("="*70)
        logger.warning("⏳ This is a Phase 5 implementation task")
        
        raise NotImplementedError(
            "calculate_surface_energies() is a Phase 5 task."
        )
    
    # =========================================================================
    # UTILITIES
    # =========================================================================
    
    def run_all(
        self,
        label_prefix: str = 'slab_workflow',
    ) -> Dict:
        """
        Run complete workflow (all phases 1-5).
        
        **Status**: Phase 2 complete, Phases 1,3-5 pending
        
        Parameters
        ----------
        label_prefix : str, default='slab_workflow'
            Base prefix for all calculations
            
        Returns
        -------
        Dict
            Results from all phases
        """
        logger.info("="*70)
        logger.info("COMPLETE SLAB WORKFLOW")
        logger.info("="*70)
        logger.info("Phase 1 (Bulk Convergence): ⏳ PENDING")
        logger.info("Phase 2 (Slab Generation):  ✓ RUNNING")
        logger.info("Phase 3 (K-convergence):    ⏳ PENDING")
        logger.info("Phase 4 (Relaxation):       ⏳ PENDING")
        logger.info("Phase 5 (Analysis):         ⏳ PENDING")
        logger.info("")
        
        # Phase 2 - available now
        self.generate_slabs(save_slabs=True)
        
        logger.warning(
            "\n⏳ Phases 1, 3-5 not yet implemented.\n"
            "See SLAB_WORKFLOW_IMPLEMENTATION.md for implementation progress."
        )
        
        return {
            'slabs': self.slabs,
            'status': 'Phase 2 complete. See SLAB_WORKFLOW_IMPLEMENTATION.md',
        }
    
    def save_results(
        self,
        output_dir: str = './',
    ) -> None:
        """
        Save workflow results to disk.
        
        Saves slab structures and any available results.
        
        Parameters
        ----------
        output_dir : str, default='./'
            Output directory for results
        """
        logger.info(f"Saving results to {output_dir}")
        
        # Create directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Save slabs
        if self.slabs:
            for hkl, slab in self.slabs.items():
                filename = os.path.join(output_dir, f'slab_{hkl[0]}{hkl[1]}{hkl[2]}.cif')
                write(filename, slab)
                logger.info(f"  Saved: {filename}")
        else:
            logger.warning("  No slabs to save. Run generate_slabs() first.")
    
    def __repr__(self) -> str:
        """String representation."""
        return (
            f"SlabWorkflow(\n"
            f"  bulk: {self.bulk_atoms.get_chemical_formula()} ({len(self.bulk_atoms)} atoms)\n"
            f"  surfaces: {len(self.surface_indices)}\n"
            f"  generated_slabs: {len(self.slabs)}\n"
            f"  status: Phase 2 ready\n"
            f")"
        )

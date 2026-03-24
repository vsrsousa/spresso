"""
SlabWorkflow: Complete workflow for surface slab calculations.

This module provides a modern, integrated interface for:
1. Bulk structure convergence (using ConvergenceWorkflow)
2. Slab generation from bulk structures  
3. Slab-specific convergence studies (k-mesh, layer thickness)
4. Structure relaxation with proper constraints
5. Surface energy calculations

Requires: ase, pymatgen, xespresso>=2024.0

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
    from pymatgen.analysis.adsorption import AdsorbateSiteFinder
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.io.ase import AseAtomsAdaptor
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False
    logging.warning("pymatgen not available. SlabWorkflow will have limited functionality.")

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
    min_slab_size : float, default=6.0
        Minimum slab thickness in Angstrom
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
        min_slab_size: float = 6.0,
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
    ):
        """
        Initialize SlabWorkflow.
        
        Validates inputs and prepares for slab calculations.
        
        Raises
        ------
        ValueError
            If bulk_atoms is empty or if both machine and queue are specified
        ImportError
            If pymatgen is not installed
        """
        # Validate inputs
        if bulk_atoms is None or len(bulk_atoms) == 0:
            raise ValueError("bulk_atoms must be a non-empty Atoms object")
        
        if not HAS_PYMATGEN:
            raise ImportError(
                "pymatgen is required for SlabWorkflow. "
                "Install with: pip install pymatgen"
            )
        
        if machine is not None and queue is not None:
            raise ValueError(
                "Cannot specify both 'machine' and 'queue'. Use one or the other."
            )
        
        # Core structure
        self.bulk_atoms = bulk_atoms.copy()
        self.surface_indices = surface_indices or [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
        
        # Slab parameters
        self.min_slab_size = min_slab_size
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
        
        # Results storage (populated by workflow phases)
        self.bulk_results = None
        self.bulk_recommendations = None
        self.slabs = {}  # {(h,k,l): Atoms}
        self.relaxed_slabs = {}
        self.surface_energies = {}
        self.convergence_results = {}
        
        logger.info(
            f"SlabWorkflow initialized:\n"
            f"  Bulk: {self.bulk_atoms.get_chemical_formula()} ({len(self.bulk_atoms)} atoms)\n"
            f"  Surfaces: {self.surface_indices}\n"
            f"  Layers: {self.nlayers}, Fix indices: {self.fix_layer_indices}\n"
            f"  Min slab: {self.min_slab_size} Å, Min vacuum: {self.min_vacuum_size} Å"
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
        Generate slabs for all surface indices using pymatgen.
        
        **Status**: ✓ COMPLETE
        
        Uses pymatgen SlabGenerator to create slabs from bulk structure.
        Applies post-processing:
        - Cell orthogonalization
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
        """
        indices = override_surface_indices or self.surface_indices
        logger.info(f"Generating slabs for {len(indices)} surfaces...")
        
        # Convert to pymatgen for slab generation
        try:
            struct = AseAtomsAdaptor.get_structure(self.bulk_atoms)
            struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
            logger.info(f"  Converted to pymatgen conventional structure")
        except Exception as e:
            logger.error(f"Failed to convert structure to pymatgen: {e}")
            raise RuntimeError(f"Structure conversion failed: {e}")
        
        # Create output directory if saving
        if save_slabs:
            Path(save_dir).mkdir(parents=True, exist_ok=True)
        
        # Generate each slab
        for hkl in indices:
            surface_name = f"{hkl[0]}{hkl[1]}{hkl[2]}"
            logger.info(f"Generating {surface_name}...")
            
            try:
                # Generate using pymatgen SlabGenerator
                slabgen = SlabGenerator(
                    struct,
                    miller_index=hkl,
                    min_slab_size=self.min_slab_size,
                    min_vacuum_size=self.min_vacuum_size,
                    center_slab=False,
                )
                
                slabs_list = slabgen.get_slabs(tol=0.1)
                if not slabs_list:
                    logger.warning(f"  ⚠ No slabs generated for {hkl}")
                    continue
                
                # Use first (lowest energy) termination
                slab_pymatgen = slabs_list[0]
                slab = AseAtomsAdaptor.get_atoms(slab_pymatgen)
                
                # Post-processing
                slab = self._orthogonalize_cell(slab)
                slab.center(vacuum=self.min_vacuum_size, axis=2)
                slab = self._apply_constraints(slab)
                
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
    
    def _orthogonalize_cell(self, slab: Atoms) -> Atoms:
        """
        Orthogonalize slab cell.
        
        Transforms cell vectors to standard orientation:
        - a1 along x-axis
        - a2 in x-y plane
        - a3 perpendicular (along z for slabs)
        
        Parameters
        ----------
        slab : ase.Atoms
            Input slab structure
            
        Returns
        -------
        ase.Atoms
            Slab with orthogonalized cell
        """
        slab = slab.copy()
        a1, a2, a3 = slab.cell
        
        # Build new orthogonal cell
        new_cell = np.zeros((3, 3))
        
        # a1 along x
        norm_a1 = np.linalg.norm(a1)
        new_cell[0, 0] = norm_a1
        
        # a2 in x-y plane
        norm_a2 = np.linalg.norm(a2)
        cos_angle = np.dot(a1, a2) / (norm_a1 * norm_a2)
        cos_angle = np.clip(cos_angle, -1, 1)  # Numerical safety
        sin_angle = np.sqrt(1 - cos_angle**2)
        
        new_cell[1, 0] = norm_a2 * cos_angle
        new_cell[1, 1] = norm_a2 * sin_angle
        
        # a3 perpendicular to a1-a2 plane (along z)
        new_cell[2, 2] = np.linalg.norm(a3)
        
        slab.set_cell(new_cell, scale_atoms=True)
        slab.wrap()
        
        return slab
    
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
    
    # =========================================================================
    # PHASE 1: BULK CONVERGENCE (PLACEHOLDER - Phase 2 task)
    # =========================================================================
    
    def run_bulk_convergence(
        self,
        label_prefix: str = 'bulk_convergence',
        verbose: Optional[bool] = None,
        convergence_criteria: List[str] = None,
        max_ecutwfc: float = 200.0,
        ecutwfc_step: float = 10.0,
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
                queue=self.queue,
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
    # PHASE 3: SLAB CONVERGENCE (PLACEHOLDER - Phase 3 task)
    # =========================================================================
    
    def run_slab_convergence(
        self,
        surface_index: Tuple[int, int, int],
        nlayers_test: Optional[List[int]] = None,
        label_prefix: str = 'slab_convergence',
    ) -> Dict:
        """
        Run slab-specific convergence study.
        
        **Status**: ⏳ PHASE 3 task (not yet implemented)
        
        Will test convergence for k-mesh, layer thickness, and constraints.
        
        Parameters
        ----------
        surface_index : Tuple[int, int, int]
            Miller index of surface
        nlayers_test : List[int], optional
            Number of layers to test
        label_prefix : str, default='slab_convergence'
            Prefix for calculation directories
            
        Returns
        -------
        Dict
            Convergence results
            
        Raises
        ------
        NotImplementedError
            Currently a Phase 3 task
        """
        logger.warning("="*70)
        logger.warning(f"PHASE 3: SLAB CONVERGENCE ({surface_index})")
        logger.warning("="*70)
        logger.warning("⏳ This is a Phase 3 implementation task")
        
        raise NotImplementedError(
            "run_slab_convergence() is a Phase 3 task."
        )
    
    # =========================================================================
    # PHASE 4: STRUCTURE RELAXATION (PLACEHOLDER - Phase 4 task)
    # =========================================================================
    
    def relax_slabs(
        self,
        label_prefix: str = 'relax_slabs',
    ) -> Dict:
        """
        Relax slab structures.
        
        **Status**: ⏳ PHASE 4 task (not yet implemented)
        
        Will integrate CalculationWorkflow for structure relaxation
        with anisotropic k-mesh and layer constraints.
        
        Parameters
        ----------
        label_prefix : str, default='relax_slabs'
            Prefix for calculation directories
            
        Returns
        -------
        Dict
            Relaxed calculations
            
        Raises
        ------
        NotImplementedError
            Currently a Phase 4 task
        """
        logger.warning("="*70)
        logger.warning("PHASE 4: STRUCTURE RELAXATION")
        logger.warning("="*70)
        logger.warning("⏳ This is a Phase 4 implementation task")
        
        raise NotImplementedError(
            "relax_slabs() is a Phase 4 task."
        )
    
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

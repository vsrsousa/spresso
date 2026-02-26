"""
DatabaseWorkflow: Integration of CalculationWorkflow with ASE Database and Provenance.

Provides automatic caching of calculations, result retrieval, and provenance tracking.
"""

import logging
import hashlib
import json
from pathlib import Path
from typing import Dict, Optional, Union

import ase.db
from ase import Atoms

from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.db.provenance import ProvenanceDB

logger = logging.getLogger(__name__)


class DatabaseWorkflow:
    """
    Wrapper around CalculationWorkflow that integrates with ASE Database and Provenance.
    
    Features:
    - Automatic caching of calculations based on structure + parameters
    - Retrieval of previous results without recalculation
    - Provenance tracking of all calculations and derivations
    - Convergence parameter inheritance for related calculations
    """
    
    def __init__(self, 
                 db_path: Union[str, Path] = '~/.xespresso/structures.db',
                 provenance_path: Union[str, Path] = '~/.xespresso/provenance.db'):
        """
        Initialize DatabaseWorkflow.
        
        Args:
            db_path: Path to ASE database file (default: ~/.xespresso/structures.db)
            provenance_path: Path to provenance SQLite database (default: ~/.xespresso/provenance.db)
        """
        self.db_path = Path(db_path).expanduser()
        self.provenance_path = Path(provenance_path).expanduser()
        
        # Connect to databases
        self.db = ase.db.connect(str(self.db_path))
        self.provenance = ProvenanceDB(str(self.provenance_path))
        
        logger.info(f"DatabaseWorkflow initialized")
        logger.info(f"  ASE Database: {self.db_path}")
        logger.info(f"  Provenance DB: {self.provenance_path}")
    
    def _compute_calculation_hash(self, atoms: Atoms, calc_params: Dict, calculation_method: str = 'scf') -> str:
        """
        Compute hash of calculation inputs.
        
        Hash includes:
        - Atomic structure (composition + positions)
        - Code type (pw, ph, dos, bands, etc) and calculation method (scf, relax, etc)
        - Pseudopotential configuration
        - All relevant QE parameters (ecutwfc, kspacing, conv_thr, etc)
        - Hubbard parameters (U values, magnetic config)
        - Complete input_data with all QE settings
        - NOT machine configuration or code_version
        
        Args:
            atoms: ASE Atoms object
            calc_params: Calculation parameters dict
            calculation_method: Type of calculation (scf, relax, vc-relax, etc)
        
        Returns:
            str: Hex hash
        """
        # Extract relevant calculation parameters (exclude machine and code_version)
        # Include complete input_data which contains all QE settings
        qe_params = {
            'code': calc_params.get('code'),  # Which QE code (pw, ph, dos, bands, wannier, etc)
            'calculation_method': calculation_method,  # Type of calculation (scf, relax, vc-relax, phonon, etc)
            'protocol': calc_params.get('protocol'),
            'pseudopotentials_config': calc_params.get('pseudopotentials_config'),  # Pseudopotential configuration
            'pseudopotentials': calc_params.get('pseudopotentials'),
            'ecutwfc': calc_params.get('ecutwfc'),
            'ecutrho': calc_params.get('ecutrho'),
            'kspacing': calc_params.get('kspacing'),
            'conv_thr': calc_params.get('conv_thr'),
            'mixing_beta': calc_params.get('mixing_beta'),
            'electron_maxstep': calc_params.get('electron_maxstep'),
            'magnetic_config': calc_params.get('magnetic_config'),
            # Relaxation parameters
            'cell_dofree': calc_params.get('cell_dofree'),  # Cell degrees of freedom (for vc-relax)
            # Hubbard parameters (DFT+U settings)
            'hubbard': calc_params.get('hubbard'),
            'hubbard_v': calc_params.get('hubbard_v'),
            'lda_plus_u': calc_params.get('lda_plus_u'),
            # Complete input_data contains all QE parameters (ibrav, volume, cell options, etc)
            'input_data': json.dumps(calc_params.get('input_data', {}), sort_keys=True),
        }
        
        # Include atomic structure
        structure_str = f"{atoms.get_chemical_formula()}:{atoms.get_positions().tobytes()}"
        
        # Combine and hash
        combined = json.dumps(qe_params, sort_keys=True) + structure_str
        calc_hash = hashlib.sha256(combined.encode()).hexdigest()
        
        return calc_hash
    
    def get_structure_with_params(self, structure_id: int):
        """
        Retrieve structure and its convergence parameters.
        
        Args:
            structure_id: ASE database ID
        
        Returns:
            Tuple of (Atoms, convergence_params_dict or None)
        """
        row = self.db.get_atoms(structure_id)
        atoms = row.toatoms()
        
        # Get convergence parameters if stored
        conv_params = row.get('convergence_params')
        if isinstance(conv_params, str):
            conv_params = json.loads(conv_params)
        
        return atoms, conv_params
    
    def get_or_calculate(self, 
                        atoms: Atoms,
                        calculation_params: Dict,
                        calculation_method: str = 'scf',
                        input_structure_id: Optional[int] = None,
                        force_recalculate: bool = False) -> tuple:
        """
        Get calculation result from cache or execute if not found.
        
        Args:
            atoms: ASE Atoms object
            calculation_params: Dict with protocol, machine, etc
            calculation_method: Type of calculation ('scf', 'relax', 'phonon', etc)
            input_structure_id: ASE DB ID of input structure (for provenance)
            force_recalculate: If True, ignore cache and recalculate
        
        Returns:
            Tuple of (result, from_cache)
            - result: Atoms object (with calc results) or dict with results
            - from_cache: Boolean indicating if result was retrieved from cache
        """
        # Compute hash of this calculation
        calc_hash = self._compute_calculation_hash(atoms, calculation_params, calculation_method)
        logger.info(f"Calculation hash: {calc_hash[:8]}...")
        
        # Check if calculation exists in provenance
        if not force_recalculate:
            existing = self.provenance.query_by_hash(calc_hash)
            if existing:
                logger.info(f"Found existing calculation in provenance")
                output_structure_id = existing['output_structure_id']
                
                if output_structure_id is not None:
                    # Retrieve from ASE database
                    result = self.db.get_atoms(output_structure_id)
                    
                    logger.info(f"Retrieved result from database (structure ID: {output_structure_id})")
                    return result, True
        
        # Calculation not found, need to compute
        logger.info(f"Calculation not in cache, executing...")
        
        # Create workflow
        workflow = CalculationWorkflow(atoms, **calculation_params)
        
        result = None
        success = False
        error_message = None
        
        try:
            if calculation_method == 'scf':
                result = workflow.run_scf()
            elif calculation_method in ['relax', 'vc-relax']:
                # Automatically detect relax_type from calculation_method
                relax_type = calculation_method if calculation_method == 'vc-relax' else 'relax'
                label = calculation_method  # Use different labels: 'relax' vs 'vc-relax'
                result = workflow.run_relax(label=label, relax_type=relax_type)
            elif calculation_method == 'phonon':
                result = workflow.run_phonon()
            elif calculation_method == 'band':
                result = workflow.run_band()
            else:
                raise ValueError(f"Unknown calculation method: {calculation_method}")
            
            success = True
            
        except Exception as e:
            success = False
            error_message = str(e)
            logger.error(f"Calculation failed: {error_message}")
            
            # Log failed calculation to provenance
            calc_id = self.provenance.log_calculation(
                calc_hash=calc_hash,
                input_structure_id=input_structure_id,
                output_structure_id=None,  # No output for failed calc
                energy=None,
                convergence_steps=None,
                convergence_params=calculation_params,
                calculation_method=calculation_method,
                success=False,
                error_message=error_message
            )
            
            # Log execution details for failed calc
            machine = calculation_params.get('machine')
            self.provenance.log_execution(
                calculation_id=calc_id,
                machine=machine,
                xespresso_version='1.0.0',
                success=False,
                error_message=error_message
            )
            
            logger.info(f"Logged failed calculation to provenance (calc ID: {calc_id})")
            # Re-raise the exception
            raise
        
        # Only reaches here on success
        # Store result in ASE database
        # Detach calculator before writing - ASE database doesn't need it
        # and the Espresso calc.check_state() returns bool instead of list
        result_atoms = workflow.atoms.copy()
        result_atoms.calc = None  # Remove calculator reference
        
        output_structure_id = self.db.write(result_atoms,
                                           calculation_hash=calc_hash,
                                           calculation_method=calculation_method,
                                           convergence_params=json.dumps(calculation_params))
        
        # Log to provenance
        # result is the Espresso calculator object with results from read_results()
        energy = result.results.get('energy') if hasattr(result, 'results') and result.results else None
        convergence_steps = result.results.get('convergence_steps') if hasattr(result, 'results') and result.results else None
        
        calc_id = self.provenance.log_calculation(
            calc_hash=calc_hash,
            input_structure_id=input_structure_id,
            output_structure_id=output_structure_id,
            energy=energy,
            convergence_steps=convergence_steps,
            convergence_params=calculation_params,
            calculation_method=calculation_method,
            success=True
        )
        
        # Log execution details
        machine = calculation_params.get('machine')
        self.provenance.log_execution(
            calculation_id=calc_id,
            machine=machine,
            xespresso_version='1.0.0'  # TODO: get from xespresso.__version__
        )
        
        # Log derivation if input structure was provided
        if input_structure_id is not None:
            self.provenance.log_derivation(
                source_structure_id=input_structure_id,
                derived_structure_id=output_structure_id,
                derivation_method=calculation_method,
                energy_change=energy
            )
        
        logger.info(f"Stored result in database (structure ID: {output_structure_id})")
        return result_atoms, False
    
    def run_properties_on_structure(self,
                                   structure_id: int,
                                   calc_type: str,
                                   override_params: Optional[Dict] = None) -> tuple:
        """
        Calculate properties on a structure using inherited convergence parameters.
        
        Args:
            structure_id: ASE database ID of structure
            calc_type: Type of property calculation (phonon, band, dos, etc)
            override_params: Optional dict to override inherited parameters
        
        Returns:
            Tuple of (result, from_cache)
        """
        # Get structure and its convergence parameters
        atoms, inherited_params = self.get_structure_with_params(structure_id)
        
        if inherited_params is None:
            raise ValueError(
                f"Structure {structure_id} has no convergence parameters. "
                "Use get_or_calculate() first on this structure."
            )
        
        # Merge with overrides
        calc_params = inherited_params.copy()
        if override_params:
            calc_params.update(override_params)
        
        logger.info(f"Running {calc_type} on structure {structure_id}")
        logger.info(f"Using inherited convergence params: {inherited_params}")
        
        # Mark that this depends on the parent structure calculation
        result, from_cache = self.get_or_calculate(
            atoms,
            calc_params,
            calculation_method=calc_type,
            input_structure_id=structure_id
        )
        
        # Log dependency in provenance
        if not from_cache:
            parent_calc = self.provenance.query_by_structure(structure_id)
            if parent_calc:
                # TODO: Log dependency relationship
                pass
        
        return result, from_cache
    
    def log_structure_derivation(self,
                                source_structure_id: int,
                                derived_structure_id: int,
                                derivation_method: str,
                                energy_change: Optional[float] = None):
        """
        Log that one structure was derived from another (e.g., relaxation).
        
        Args:
            source_structure_id: ASE DB ID of parent
            derived_structure_id: ASE DB ID of derived structure
            derivation_method: Type of derivation (vc-relax, relax, etc)
            energy_change: Energy difference (optional)
        """
        self.provenance.log_derivation(
            source_structure_id=source_structure_id,
            derived_structure_id=derived_structure_id,
            derivation_method=derivation_method,
            energy_change=energy_change
        )
        logger.info(f"Logged derivation: {source_structure_id} → {derived_structure_id}")
    
    def get_structure_lineage(self, structure_id: int) -> list:
        """
        Get complete evolution of a structure.
        
        Args:
            structure_id: ASE DB ID
        
        Returns:
            List of dicts showing structure evolution
        """
        return self.provenance.get_calculation_history(structure_id)
    
    def validate_consistency(self, structure_id: int) -> bool:
        """
        Verify all calculations on a structure used consistent parameters.
        
        Args:
            structure_id: ASE DB ID
        
        Returns:
            bool: True if consistent, False otherwise
        """
        calcs = self.provenance.query_by_structure(structure_id)
        
        if not calcs:
            return True
        
        first_params = json.loads(calcs[0]['convergence_params']) if calcs[0]['convergence_params'] else None
        
        for calc in calcs[1:]:
            calc_params = json.loads(calc['convergence_params']) if calc['convergence_params'] else None
            if calc_params != first_params:
                logger.warning(
                    f"Inconsistency found: Calculation {calc['id']} uses different parameters"
                )
                return False
        
        logger.info(f"✅ All calculations on structure {structure_id} are consistent")
        return True
    
    def close(self):
        """Close all database connections."""
        self.provenance.close()
        logger.info("DatabaseWorkflow closed")

"""
Workflow module for xespresso GUI.

This module provides workflow orchestration for the GUI, coordinating
multiple calculations and using xespresso's logic and definitions.
"""

from typing import Dict, List, Tuple, Optional
from ase import Atoms
from xespresso import Espresso
from gui.calculations import prepare_calculation_from_gui
import logging

logger = logging.getLogger(__name__)


class GUIWorkflow:
    """
    Workflow orchestration for GUI calculations.
    
    This class manages multiple calculation steps, using the calculation
    preparation modules to create atoms and Espresso objects, then
    coordinating their execution.
    
    Following xespresso patterns: workflows orchestrate, calculations prepare,
    job submission executes.
    """
    
    def __init__(
        self,
        atoms: Atoms,
        config: Dict,
        base_label: str = "workflow"
    ):
        """
        Initialize workflow.
        
        Args:
            atoms: Initial ASE Atoms object
            config: GUI workflow configuration
            base_label: Base directory for workflow calculations
        """
        self.atoms = atoms.copy()
        self.config = config
        self.base_label = base_label
        self.calculations = {}
        self.results = {}
        
        logger.info(f"Initialized workflow for {atoms.get_chemical_formula()}")
    
    def add_calculation(
        self,
        name: str,
        calc_config: Dict,
        atoms: Optional[Atoms] = None
    ) -> Tuple[Atoms, Espresso]:
        """
        Add a calculation step to the workflow.
        
        This method uses the calculation preparation module to create
        the atoms and Espresso objects for this calculation step.
        
        Args:
            name: Name for this calculation step
            calc_config: Configuration for this calculation
            atoms: Optional atoms object (uses workflow atoms if None)
            
        Returns:
            tuple: (atoms, calculator) prepared for this step
        """
        if atoms is None:
            atoms = self.atoms
        
        label = f"{self.base_label}/{name}"
        
        # Use calculation module to prepare objects
        calc_atoms, calculator = prepare_calculation_from_gui(
            atoms,
            calc_config,
            label=label
        )
        
        self.calculations[name] = {
            'atoms': calc_atoms,
            'calculator': calculator,
            'config': calc_config
        }
        
        logger.info(f"Added calculation step: {name}")
        
        return calc_atoms, calculator
    
    def get_calculation(self, name: str) -> Tuple[Atoms, Espresso]:
        """
        Get prepared calculation by name.
        
        Args:
            name: Calculation step name
            
        Returns:
            tuple: (atoms, calculator) for this step
        """
        if name not in self.calculations:
            raise KeyError(f"Calculation '{name}' not found in workflow")
        
        calc = self.calculations[name]
        return calc['atoms'], calc['calculator']
    
    def run_step(self, name: str) -> Dict:
        """
        Run a specific calculation step.
        
        Args:
            name: Calculation step name
            
        Returns:
            dict: Results from this calculation
        """
        if name not in self.calculations:
            raise KeyError(f"Calculation '{name}' not found in workflow")
        
        calc_info = self.calculations[name]
        atoms = calc_info['atoms']
        calculator = calc_info['calculator']
        
        logger.info(f"Running workflow step: {name}")
        
        # Execute the calculation using xespresso
        calculator.run(atoms=atoms)
        
        # Store results
        self.results[name] = calculator.results.copy()
        
        return self.results[name]
    
    def get_results(self, name: Optional[str] = None) -> Dict:
        """
        Get results from workflow.
        
        Args:
            name: Optional calculation name. If None, returns all results.
            
        Returns:
            dict: Results dictionary
        """
        if name is None:
            return self.results
        
        if name not in self.results:
            raise KeyError(f"No results for calculation '{name}'")
        
        return self.results[name]


def create_scf_relax_workflow(
    atoms: Atoms,
    pseudopotentials: Dict[str, str],
    base_label: str = "workflow",
    **kwargs
) -> GUIWorkflow:
    """
    Create a workflow with SCF followed by relaxation.
    
    This is an example of using the workflow module to coordinate
    multiple calculations, following xespresso patterns.
    
    Args:
        atoms: Initial structure
        pseudopotentials: Pseudopotential mapping
        base_label: Base directory
        **kwargs: Additional parameters for calculations
        
    Returns:
        GUIWorkflow: Configured workflow ready to run
        
    Example:
        >>> workflow = create_scf_relax_workflow(
        ...     atoms,
        ...     {'Fe': 'Fe.pbe.UPF'},
        ...     base_label='fe_workflow'
        ... )
        >>> # Get prepared calculations
        >>> scf_atoms, scf_calc = workflow.get_calculation('scf')
        >>> # Execute step
        >>> workflow.run_step('scf')
    """
    # Base configuration
    base_config = {
        'pseudopotentials': pseudopotentials,
        'ecutwfc': kwargs.get('ecutwfc', 50.0),
        'ecutrho': kwargs.get('ecutrho', 400.0),
        'kpts': kwargs.get('kpts', (4, 4, 4)),
        'occupations': kwargs.get('occupations', 'smearing'),
        'degauss': kwargs.get('degauss', 0.02),
    }
    
    # Create workflow
    workflow = GUIWorkflow(atoms, base_config, base_label)
    
    # Add SCF step
    scf_config = base_config.copy()
    scf_config['calc_type'] = 'scf'
    workflow.add_calculation('scf', scf_config)
    
    # Add relaxation step
    relax_config = base_config.copy()
    relax_config['calc_type'] = 'relax'
    relax_config['forc_conv_thr'] = kwargs.get('forc_conv_thr', 1.0e-3)
    workflow.add_calculation('relax', relax_config)
    
    return workflow

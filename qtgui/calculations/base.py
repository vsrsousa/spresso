"""
Base calculation preparation module for xespresso GUI.

This module provides the base class for preparing calculations in the GUI.
Following the principle that calculation modules should create and prepare
atoms and Espresso calculator objects, not the main GUI code.
"""

from typing import Dict, Optional, Tuple
from ase import Atoms
from xespresso import Espresso
import logging

logger = logging.getLogger(__name__)


class BaseCalculationPreparation:
    """
    Base class for preparing calculations in the GUI.
    
    This class encapsulates the logic for creating atoms and Espresso
    calculator objects from GUI configuration, following xespresso's
    design patterns.
    
    The key principle is: calculation modules prepare the objects,
    job submission modules execute them.
    """
    
    def __init__(
        self,
        atoms: Atoms,
        config: Dict,
        label: str = "calculation"
    ):
        """
        Initialize calculation preparation.
        
        Args:
            atoms: ASE Atoms object representing the structure
            config: Configuration dictionary from GUI (workflow_config)
            label: Directory/prefix for calculation files
        """
        self.atoms = atoms.copy()  # Work with a copy
        self.config = config
        self.label = label
        self.calculator = None
        
        logger.info(f"Initialized {self.__class__.__name__} for {atoms.get_chemical_formula()}")
    
    def prepare(self) -> Tuple[Atoms, Espresso]:
        """
        Prepare atoms and Espresso calculator for the calculation.
        
        This method creates the Espresso calculator with appropriate
        parameters based on the GUI configuration, following xespresso patterns.
        
        Returns:
            tuple: (atoms, calculator) ready for execution
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement prepare() method"
        )
    
    def get_calculator(self) -> Espresso:
        """
        Get the prepared calculator.
        
        If prepare() hasn't been called, calls it first.
        
        Returns:
            Espresso: The prepared calculator
        """
        if self.calculator is None:
            self.prepare()
        return self.calculator
    
    def get_atoms(self) -> Atoms:
        """
        Get the atoms object.
        
        Returns:
            Atoms: The atoms for this calculation
        """
        return self.atoms

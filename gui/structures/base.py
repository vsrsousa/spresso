"""
Base structure handling module for xespresso GUI.

This module provides the base class for structure loading and management,
following the modular design pattern where dedicated modules handle
specific functionality.
"""

from typing import Optional, List, Tuple
from pathlib import Path
from ase import Atoms
import logging

logger = logging.getLogger(__name__)


class BaseStructureHandler:
    """
    Base class for structure handling in the GUI.
    
    This class encapsulates structure loading, validation, and management
    operations, following the principle that dedicated modules should handle
    specific functionality rather than having all logic in GUI pages.
    
    Attributes:
        atoms (Atoms): The loaded ASE Atoms object
        source_path (str): Path to the source file
        format (str): Structure file format
    """
    
    def __init__(self, atoms: Optional[Atoms] = None, source_path: Optional[str] = None):
        """
        Initialize structure handler.
        
        Args:
            atoms: Optional ASE Atoms object
            source_path: Optional path to source file
        """
        self.atoms = atoms
        self.source_path = source_path
        self.format = None
        
        if atoms and source_path:
            self.format = self._detect_format(source_path)
            logger.info(f"Initialized structure handler with {len(atoms)} atoms from {source_path}")
    
    def _detect_format(self, filepath: str) -> str:
        """
        Detect structure file format from extension.
        
        Args:
            filepath: Path to structure file
            
        Returns:
            str: Detected format
        """
        ext = Path(filepath).suffix.lower()
        format_map = {
            '.cif': 'cif',
            '.xyz': 'xyz',
            '.pdb': 'pdb',
            '.vasp': 'vasp',
            '.poscar': 'vasp',
            '.traj': 'traj',
            '.json': 'json'
        }
        return format_map.get(ext, 'unknown')
    
    def get_atoms(self) -> Optional[Atoms]:
        """
        Get the atoms object.
        
        Returns:
            Atoms: The loaded atoms object
        """
        return self.atoms
    
    def get_info(self) -> dict:
        """
        Get structure information.
        
        Returns:
            dict: Structure information dictionary
        """
        if not self.atoms:
            return {}
        
        return {
            'formula': self.atoms.get_chemical_formula(),
            'n_atoms': len(self.atoms),
            'cell': self.atoms.get_cell().tolist(),
            'pbc': self.atoms.get_pbc().tolist(),
            'source_path': self.source_path,
            'format': self.format
        }

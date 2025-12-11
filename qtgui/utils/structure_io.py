"""
Utility functions for reading and writing atomic structures.

This module provides functions for reading structure files with proper
handling of symmetry information, particularly for CIF files.
"""

from ase import io as ase_io
from ase import Atoms
from typing import Union
from pathlib import Path


def read_structure(file_path: Union[str, Path], **kwargs) -> Atoms:
    """
    Read a structure file with proper handling of CIF symmetry.
    
    For CIF files, this function reads the primitive cell by default to preserve
    symmetry information. For other file formats, it uses the standard ASE reader.
    
    Args:
        file_path: Path to the structure file
        **kwargs: Additional keyword arguments passed to ase.io.read()
        
    Returns:
        Atoms: The atomic structure
        
    Notes:
        For CIF files, the function automatically sets primitive_cell=True to read
        the primitive cell instead of expanding all symmetry operations.
        
    Examples:
        >>> atoms = read_structure('structure.cif')  # Reads primitive cell
        >>> atoms = read_structure('structure.xyz')  # Standard reading
        >>> atoms = read_structure('structure.cif', primitive_cell=False)  # Force expanded cell
    """
    file_path = str(file_path)
    
    # Check if it's a CIF file
    if file_path.lower().endswith('.cif'):
        # Set default parameters for CIF files to read primitive cell
        # unless explicitly overridden by user
        if 'primitive_cell' not in kwargs:
            kwargs['primitive_cell'] = True
        # When using primitive_cell=True, we need to set subtrans_included=False.
        # The subtrans_included parameter controls whether sublattice translations
        # (e.g., centering operations like face-centering, body-centering) are included
        # in the symmetry operations. When these are included, ASE cannot determine
        # the primitive cell and will raise a RuntimeError. Setting it to False
        # allows ASE to extract the primitive unit cell from the symmetry information.
        if kwargs.get('primitive_cell', False) and 'subtrans_included' not in kwargs:
            kwargs['subtrans_included'] = False
    
    return ase_io.read(file_path, **kwargs)

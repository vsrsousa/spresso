"""
Band path generation using standardized k-paths from seekpath.

This module provides a function to generate band structure k-paths by:
1. Using seekpath to identify the crystal's Bravais lattice and space group
2. Looking up the standardized band path from hardcoded seekpath data
3. Mapping GAMMA → G for compatibility with Wannier90

The module avoids runtime seekpath path generation by using hardcoded
definitions extracted from seekpath's GitHub repository.
"""

from ase.dft.kpoints import BandPath
try:
    import seekpath
except ImportError:
    raise ImportError("seekpath is required. Install with: pip install seekpath")

from .spresso_seekpath_data import SEEKPATH_SPACE_GROUP_DATA


def get_bandpath(atoms, with_time_reversal=True):
    """
    Generate band structure k-path using seekpath conventions.
    
    Uses the crystal's space group (via seekpath) to select the appropriate
    standardized band path from a pre-computed hardcoded database.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atomic structure for which to generate band path
    with_time_reversal : bool, optional
        Whether to use time reversal symmetry (default: True)
    
    Returns
    -------
    ase.dft.kpoints.BandPath
        Band path object with path string and special points
        
    Raises
    ------
    ImportError
        If seekpath is not installed
    ValueError
        If the crystal system is not supported
        
    Examples
    --------
    >>> from ase.build import bulk
    >>> import numpy as np
    >>> from xespresso.utils.bandpath import get_bandpath
    >>> 
    >>> # Silicon (FCC, Fd-3m)
    >>> si = bulk('Si', 'diamond', a=5.4)
    >>> bp = get_bandpath(si)
    >>> print(bp.path)
    'GXU, KGLWX'
    
    >>> # Iron (BCC)
    >>> fe = bulk('Fe', 'bcc', a=2.87)
    >>> bp = get_bandpath(fe)
    >>> print(bp.path)
    'GHNGPH, PN'
    """
    
    # Convert to seekpath format (cell, positions, atomic_numbers)
    cell_tuple = (atoms.cell[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    
    # Use seekpath to identify Bravais lattice and space group
    sp_result = seekpath.get_path(cell_tuple, with_time_reversal=with_time_reversal)
    
    # Use bravais_lattice_extended which gives the specific space group type
    # e.g., 'cF1', 'cF2', 'cI1', 'hP1', etc.
    bravais_lattice_extended = sp_result['bravais_lattice_extended']
    
    # Look up the band path in our hardcoded data
    if bravais_lattice_extended not in SEEKPATH_SPACE_GROUP_DATA:
        raise ValueError(
            f"Unsupported space group type: {bravais_lattice_extended}. "
            f"Supported types: {list(SEEKPATH_SPACE_GROUP_DATA.keys())}"
        )
    
    sg_data = SEEKPATH_SPACE_GROUP_DATA[bravais_lattice_extended]
    
    # Extract path and special points
    path_str = sg_data['path']
    special_points = sg_data['special_points'].copy()
    
    # Map GAMMA → G for Wannier90 compatibility
    path_str = path_str.replace('GAMMA', 'G')
    special_points = {
        k.replace('GAMMA', 'G'): v 
        for k, v in special_points.items()
    }
    
    # Create and return ASE BandPath object
    return BandPath(path=path_str, special_points=special_points, cell=atoms.cell)

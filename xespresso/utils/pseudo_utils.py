"""
Pseudopotential utility functions.

Provides utilities for working with pseudopotential configurations,
including ratio calculations for ecutrho based on pseudopotential types.
"""

import logging
from typing import Optional, Set
from xespresso.pseudopotentials.config import PseudopotentialsConfig

logger = logging.getLogger(__name__)


def get_ecutrho_ratio(
    elements: Set[str],
    pseudo_config: Optional[PseudopotentialsConfig] = None
) -> float:
    """
    Calculate the appropriate ecutrho/ecutwfc ratio based on pseudopotential types.
    
    - Norm-Conserving (NC): ratio = 4.0
    - Ultrasoft (US): ratio = 8.0
    - PAW: ratio = 8.0
    
    When mixing different types, returns the MAXIMUM ratio needed to ensure
    compatibility with all pseudopotentials in the structure.
    
    Args:
        elements: Set of element symbols in the structure (e.g., {'Si', 'O'})
        pseudo_config: PseudopotentialsConfig object loaded from JSON config.
                      If None, returns default 4.0 (assumes Norm-Conserving)
    
    Returns:
        float: Ratio (ecutrho/ecutwfc) - default 4.0 if no pseudo_config
        
    Example:
        >>> from xespresso.pseudopotentials import load_pseudopotentials_config
        >>> config = load_pseudopotentials_config('SSSP_efficiency')
        >>> ratio = get_ecutrho_ratio({'Si'}, config)  # Returns 8.0 (Ultrasoft)
        >>> ecutrho = 50 * ratio  # = 400.0
    """
    if pseudo_config is None:
        logger.debug("No pseudo_config available, using default ratio 4.0 (NC assume)")
        return 4.0  # Default ratio for NC pseudos
    
    max_ratio = 4.0  # Start with NC ratio
    
    try:
        logger.debug(f"Elements in structure: {elements}")
        
        for element in elements:
            pseudo_obj = pseudo_config.get_pseudopotential(element)
            
            if pseudo_obj is None:
                logger.warning(f"Element {element}: pseudopotential not found in config")
                continue
            
            logger.debug(f"Element {element}: pseudo_obj type={type(pseudo_obj)}, has 'type' attr={hasattr(pseudo_obj, 'type')}")
            
            if hasattr(pseudo_obj, 'type'):
                logger.debug(f"Element {element}: pseudo_obj.type = {pseudo_obj.type} (raw value)")
            
            if pseudo_obj and hasattr(pseudo_obj, 'type') and pseudo_obj.type:
                pseudo_type = pseudo_obj.type.upper()
                logger.debug(f"Element {element}: pseudo_type uppercased = {pseudo_type}")
                
                # Determine ratio based on pseudopotential type
                if 'PAW' in pseudo_type or 'PROJECTOR' in pseudo_type:
                    ratio = 8.0
                    type_label = 'PAW'
                elif 'ULTRASOFT' in pseudo_type or 'US' in pseudo_type:
                    ratio = 8.0
                    type_label = 'Ultrasoft'
                elif 'NORM-CONSERVING' in pseudo_type or 'NC' in pseudo_type or 'ONCV' in pseudo_type:
                    ratio = 4.0
                    type_label = 'Norm-Conserving'
                else:
                    # Unknown type, use conservative ratio (8.0 is safer)
                    ratio = 8.0
                    type_label = f"Unknown ({pseudo_obj.type})"
                
                max_ratio = max(max_ratio, ratio)
                logger.info(f"Element {element}: type={type_label}, ratio={ratio}")
            else:
                logger.warning(f"Element {element}: no type info available (pseudo_obj.type = {getattr(pseudo_obj, 'type', 'NO ATTR')})")
    except Exception as e:
        logger.warning(f"Could not determine ecutrho ratio from pseudos: {e}", exc_info=True)
    
    if max_ratio > 4.0:
        logger.info(f"Using ecutrho/ecutwfc ratio: {max_ratio} (Ultrasoft/PAW pseudopotentials detected)")
    else:
        logger.info(f"Using ecutrho/ecutwfc ratio: {max_ratio} (Norm-Conserving pseudopotentials)")
    
    return max_ratio


def discover_pseudopotential_directory(
    pseudopotentials: dict,
) -> tuple[dict, Optional[str]]:
    """
    Discover the base directory for pseudopotential files and return resolved paths.
    
    If pseudopotentials are passed as filenames (e.g., 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'),
    this function searches for them in:
    1. Current working directory (os.getcwd())
    2. ESPRESSO_PSEUDO environment variable
    3. Common locations
    
    Returns the full paths and the base directory where they were found.
    
    Args:
        pseudopotentials: Dictionary mapping element symbols to UPF files.
                         Can be:
                         - Absolute paths: '/path/to/Gd.UPF' (returned as-is)
                         - Filenames only: 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF' (searched for)
                         - Relative paths: 'pseudo/Gd.UPF' (resolved from cwd)
    
    Returns:
        Tuple of (resolved_dict, base_dir):
        - resolved_dict: Dictionary with full absolute paths for each pseudopotential
        - base_dir: Directory where pseudopotentials were found (or None if all were absolute)
    
    Raises:
        FileNotFoundError: If any pseudopotential file cannot be found
    
    Example:
        >>> pseudos = {'Gd': 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}
        >>> resolved, base_dir = discover_pseudopotential_directory(pseudos)
        >>> # resolved = {'Gd': '/home/user/spresso/Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}
        >>> # base_dir = '/home/user/spresso'
    """
    import os
    
    resolved = {}
    base_dir = None
    
    # Search locations in order of preference
    search_paths = [
        os.getcwd(),  # Current directory (highest priority)
        os.environ.get('ESPRESSO_PSEUDO'),  # ESPRESSO_PSEUDO env var
        os.path.expandvars('$HOME/.espresso/pseudo'),  # Home xespresso pseudo dir
    ]
    # Filter out None values
    search_paths = [p for p in search_paths if p]
    
    for element, path in pseudopotentials.items():
        # If absolute path, use as-is but extract base directory
        if os.path.isabs(path):
            resolved[element] = path
            # Extract directory for base_dir if not yet set
            if base_dir is None:
                base_dir = os.path.dirname(path)
            continue
        
        # If it's a relative path or just a filename
        found = False
        found_path = None
        found_in_dir = None
        
        # First try joining with cwd
        test_path = os.path.join(os.getcwd(), path)
        if os.path.exists(test_path):
            resolved[element] = os.path.abspath(test_path)
            # Extract the directory where the file was found (not cwd)
            found_in_dir = os.path.dirname(resolved[element])
            found = True
        else:
            # Try to find just the filename in search paths
            filename = os.path.basename(path)
            for search_dir in search_paths:
                test_path = os.path.join(search_dir, filename)
                if os.path.exists(test_path):
                    resolved[element] = os.path.abspath(test_path)
                    found_in_dir = search_dir
                    found = True
                    break
        
        if not found:
            # List where we looked
            looked_in = [os.path.join(p, os.path.basename(path)) for p in search_paths]
            raise FileNotFoundError(
                f"\n❌ Pseudopotential file not found: {element} -> {path}\n"
                f"\n   Searched in:\n"
                f"   - {os.path.join(os.getcwd(), path)}\n"
                + "\n".join([f"   - {p}" for p in looked_in]) +
                f"\n\n   Solutions:\n"
                f"   1. Place pseudopotential in current directory: {os.getcwd()}/\n"
                f"   2. Set ESPRESSO_PSEUDO env var: export ESPRESSO_PSEUDO=/path/to/pseudo\n"
                f"   3. Use absolute path in pseudopotentials dict\n"
            )
        
        # Track the base directory (use the first one found)
        if base_dir is None:
            base_dir = found_in_dir
    
    # Double-check: all paths in resolved dict should be absolute
    for element, path in resolved.items():
        if not os.path.isabs(path):
            # Should not happen, but ensure it's absolute
            resolved[element] = os.path.abspath(path)
    
    return resolved, base_dir


def calculate_nbnd_from_structure(
    atoms,
    pseudopotentials: dict,
    pseudopotentials_base_path: Optional[str] = None,
    buffer: int = 0
) -> int:
    """
    Calculate number of bands (nbnd) directly from total valence electrons in structure.
    
    Reads z_valence from pseudopotential UPF files and calculates total valence electrons
    by counting atoms of each type in the structure:
    
    nbnd = (N_X × valence_X) + (N_Y × valence_Y) + ... + buffer
    
    For example, X₂Y₃ structure:
    - X: 2 atoms with valence_X electrons each
    - Y: 3 atoms with valence_Y electrons each
    - Total valence = 2×valence_X + 3×valence_Y
    - nbnd = total_valence + buffer
    
    Args:
        atoms: ASE Atoms object representing the structure
        pseudopotentials: Dict mapping element symbols to pseudopotential filenames or paths
        pseudopotentials_base_path: Base directory where pseudopotentials are located.
                                   If provided, used to find pseudos not found as absolute paths.
        buffer: Number of additional bands beyond minimum (default: 0)
    
    Returns:
        int: Number of bands (nbnd) equal to total valence electrons
        
    Examples:
        >>> from ase.build import bulk
        >>> from xespresso.utils.pseudo_utils import calculate_nbnd_from_structure
        >>> 
        >>> # Si diamond (2 Si atoms × 4 valence = 8 electrons)
        >>> atoms = bulk('Si', 'diamond', a=5.43)
        >>> pseudos = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}
        >>> nbnd = calculate_nbnd_from_structure(atoms, pseudos)
        >>> print(f"Si₂: nbnd = {nbnd}")  # nbnd = 8
        >>> 
        >>> # X₂Y₃ structure (2 atoms X with 4e⁻ + 3 atoms Y with 6e⁻)
        >>> nbnd = calculate_nbnd_from_structure(atoms_x2y3, pseudos)
        >>> print(f"X₂Y₃: nbnd = {nbnd}")  # nbnd = 2×4 + 3×6 = 26
    """
    import os
    from collections import Counter
    import re
    
    # Count atoms of each element
    element_counts = Counter(atoms.get_chemical_symbols())
    logger.debug(f"Structure composition: {dict(element_counts)}")
    
    total_valence = 0
    
    for element, count in element_counts.items():
        if element not in pseudopotentials:
            logger.warning(f"Element {element} not found in pseudopotentials dict")
            # Assume a default valence (conservative estimate)
            z_valence = 8
        else:
            pseudo_file = pseudopotentials[element]
            
            # Try to find the pseudopotential file
            upf_path = None
            if os.path.isabs(pseudo_file):
                upf_path = pseudo_file
            else:
                # Try relative to current directory first
                if os.path.exists(pseudo_file):
                    upf_path = os.path.abspath(pseudo_file)
                # Try in pseudopotentials_base_path if provided
                elif pseudopotentials_base_path:
                    candidate = os.path.join(pseudopotentials_base_path, os.path.basename(pseudo_file))
                    if os.path.exists(candidate):
                        upf_path = candidate
            
            z_valence = 8  # Default fallback
            
            if upf_path and os.path.exists(upf_path):
                try:
                    with open(upf_path, "r", encoding="utf-8", errors="ignore") as f:
                        text = f.read()
                    
                    # Look for z_valence in UPF file (multiple formats)
                    # Format 1: z_valence = "4"
                    m = re.search(r'z_valence\s*=\s*"(\d+(?:\.\d+)?)"', text, re.IGNORECASE)
                    if not m:
                        # Format 2: z_valence: 4
                        m = re.search(r'z_valence\s*:\s*(\d+(?:\.\d+)?)', text, re.IGNORECASE)
                    if not m:
                        # Format 3: Z_valence = "4"
                        m = re.search(r'Z_valence\s*=\s*"(\d+(?:\.\d+)?)"', text)
                    if not m:
                        # Format 4: <z_valence>4</z_valence> (XML)
                        m = re.search(r'<z_valence>(\d+(?:\.\d+)?)</z_valence>', text, re.IGNORECASE)
                    if not m:
                        # Format 5: <Z_valence>4</Z_valence> (XML uppercase)
                        m = re.search(r'<Z_valence>(\d+(?:\.\d+)?)</Z_valence>', text)
                    
                    if m:
                        z_valence = int(m.group(1))
                        logger.debug(f"  {element}: z_valence = {z_valence} (from {os.path.basename(upf_path)})")
                    else:
                        logger.debug(f"  {element}: z_valence not found in {os.path.basename(upf_path)}, using default {z_valence}")
                
                except Exception as e:
                    logger.debug(f"  {element}: Error reading pseudo {upf_path}: {e}, using default {z_valence}")
            else:
                logger.debug(f"  {element}: Pseudo file not found, using default z_valence = {z_valence}")
        
        element_total = z_valence * count
        logger.debug(f"  {element}: {count} atoms × {z_valence} electrons = {element_total} valence electrons")
        total_valence += element_total
    
    logger.info(f"Total valence electrons in structure: {total_valence}")
    
    # nbnd = total valence electrons + buffer
    nbnd = total_valence + buffer
    
    logger.info(f"Calculated nbnd: {total_valence} + {buffer} = {nbnd}")
    
    return nbnd


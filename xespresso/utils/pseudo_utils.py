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
            found_in_dir = os.getcwd()
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

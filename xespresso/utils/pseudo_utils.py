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

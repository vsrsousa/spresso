#!/usr/bin/env python3
"""
Snippet to add to xespresso/utils/pseudo_utils.py

This adds the calculate_nbnd_from_structure function for enhancing
nbands based on actual structure valence electron counts.
"""

import os
import re
from typing import Dict, Optional
from ase import Atoms
import logging

logger = logging.getLogger(__name__)


def calculate_nbnd_from_structure(
    atoms: Atoms,
    pseudopotentials: Dict[str, str],
    pseudopotentials_base_path: Optional[str] = None,
    buffer: int = 10,
    spin_factor: float = 2.0
) -> int:
    """
    Calculate recommended number of bands (nbnd) based on structure and pseudopotentials.
    
    Reads z_valence from pseudopotential UPF files and calculates total valence electrons
    by counting atoms of each type in the structure. This ensures nbnd is sufficient
    to hold all valence electrons.
    
    For spin-unpolarized calculations: nbnd = (total_valence / 2) + buffer
    For spin-polarized: nbnd = total_valence + buffer
    
    Args:
        atoms: ASE Atoms object representing the structure
        pseudopotentials: Dict mapping element symbols to pseudopotential filenames or paths
        pseudopotentials_base_path: Base directory where pseudopotentials are located.
                                   If provided, used to find pseudos not found as absolute paths.
        buffer: Number of additional bands beyond minimum for convergence (default: 10)
        spin_factor: Divisor for converting electrons to bands. Default 2.0 for spin-degenerate.
                    Use 1.0 for spin-polarized calculations.
    
    Returns:
        int: Recommended nbnd value
        
    Examples:
        >>> from ase.build import bulk
        >>> from xespresso.utils.pseudo_utils import calculate_nbnd_from_structure
        >>> atoms = bulk('Si', 'diamond', a=5.43)
        >>> pseudos = {'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'}
        >>> nbnd = calculate_nbnd_from_structure(atoms, pseudos, buffer=15)
        >>> print(f"Recommended nbnd for Si diamond: {nbnd}")  # 8 electrons / 2 + 15 = 19
    """
    from collections import Counter
    
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
                    
                    # Look for z_valence in UPF file
                    m = re.search(r'z_valence\s*=\s*"(\d+)"', text, re.IGNORECASE)
                    if not m:
                        m = re.search(r'z_valence\s*:\s*(\d+)', text, re.IGNORECASE)
                    if not m:
                        m = re.search(r'Z_valence\s*=\s*"(\d+)"', text)
                    
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
    
    # Convert electrons to bands
    # For spin-degenerate (non-magnetic): divide by 2
    # For spin-polarized: divide by 1 (or multiply spin_factor = 1.0)
    nbnd = max(64, int(total_valence / spin_factor) + buffer)
    
    logger.info(f"Calculated nbnd: {total_valence} / {spin_factor} + {buffer} = {nbnd}")
    
    return nbnd

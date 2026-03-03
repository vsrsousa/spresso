"""
detector.py

Utilities for detecting and extracting information from pseudopotential files.
"""

import os
import re
from typing import Dict, List, Optional, Tuple


def detect_upf_files(directory: str, recursive: bool = True) -> List[str]:
    """
    Detect all UPF (pseudopotential) files in a directory.
    
    Handles both uppercase and lowercase extensions (.UPF, .upf, .Upf, etc.)
    
    Args:
        directory: Directory path to search
        recursive: Whether to search recursively in subdirectories
    
    Returns:
        List of full paths to UPF files
    """
    upf_files = []
    
    if not os.path.exists(directory):
        return upf_files
    
    if recursive:
        for root, dirs, files in os.walk(directory):
            for file in files:
                # Case-insensitive extension check
                if file.lower().endswith('.upf'):
                    upf_files.append(os.path.join(root, file))
    else:
        for file in os.listdir(directory):
            full_path = os.path.join(directory, file)
            # Case-insensitive extension check
            if os.path.isfile(full_path) and file.lower().endswith('.upf'):
                upf_files.append(full_path)
    
    return sorted(upf_files)


def extract_element_from_filename(filename: str) -> Optional[str]:
    """
    Extract element symbol from pseudopotential filename.
    
    Handles various case combinations:
    - Fe.pbe-spn-kjpaw_psl.0.2.1.UPF -> Fe
    - si_pbe_v1.4.uspp.F.UPF -> Si
    - C.pbe-n-kjpaw_psl.1.0.0.UPF -> C
    - FE.pbe-spn.UPF -> Fe
    - fe.pbe.upf -> Fe
    - SI_pbe.UPF -> Si
    
    Args:
        filename: Pseudopotential filename
    
    Returns:
        Element symbol in proper case (e.g., 'Fe', 'Si', 'C') or None if not detected
    """
    basename = os.path.basename(filename)
    
    # Try various patterns (case-insensitive)
    patterns = [
        r'^([A-Za-z]{1,2})[\._-]',  # Element at start followed by . _ or -
        r'^([A-Za-z]{1,2})_',        # Element_something
    ]
    
    for pattern in patterns:
        match = re.match(pattern, basename, re.IGNORECASE)
        if match:
            element = match.group(1)
            # Normalize to proper case (e.g., 'fe' -> 'Fe', 'FE' -> 'Fe', 'SI' -> 'Si')
            return element[0].upper() + element[1:].lower() if len(element) > 1 else element.upper()
    
    return None


def extract_functional_from_filename(filename: str) -> Optional[str]:
    """
    Extract functional type from pseudopotential filename.
    
    Args:
        filename: Pseudopotential filename
    
    Returns:
        Functional name (e.g., 'PBE', 'LDA', 'PBEsol') or None
    """
    basename = os.path.basename(filename).lower()
    
    if 'pbe' in basename and 'pbesol' not in basename:
        return 'PBE'
    elif 'pbesol' in basename:
        return 'PBEsol'
    elif 'lda' in basename:
        return 'LDA'
    elif 'blyp' in basename:
        return 'BLYP'
    elif 'pz' in basename:
        return 'PZ'
    
    return None


def extract_type_from_filename(filename: str) -> Optional[str]:
    """
    Extract pseudopotential type from filename.
    
    Args:
        filename: Pseudopotential filename
    
    Returns:
        Pseudopotential type (e.g., 'ultrasoft', 'paw', 'norm-conserving')
    """
    basename = os.path.basename(filename).lower()

    # Prefer explicit ONCV/ONCVPSP markers first (norm-conserving)
    if re.search(r'oncv|oncvpsp|oncv_psp|oncvp', basename):
        return 'Norm-conserving'

    # Common norm-conserving markers
    if re.search(r'\b(nc|norm-conserving|normconserving)\b', basename):
        return 'Norm-conserving'

    # PAW markers (includes KJPaw/psl.kjpaw variants)
    if re.search(r'\b(paw|kjpaw|atompaw|psl\.kjpaw)\b', basename):
        return 'PAW'

    # Ultrasoft / USPP markers
    if re.search(r'\b(uspp|rrkjus|uspp\.f|uspp|uspp_)\b', basename) or 'us.' in basename:
        return 'Ultrasoft'

    # Fallback checks
    if 'ultrasoft' in basename or 'rrkj' in basename:
        return 'Ultrasoft'

    return None


def parse_upf_header(filepath: str) -> Dict[str, any]:
    """
    Parse UPF file header to extract metadata.
    
    Handles case variations in element symbols (Fe, FE, fe -> Fe).
    Extracts suggested cutoff energy if available.
    
    Args:
        filepath: Path to UPF file
    
    Returns:
        Dictionary with extracted information (element, z_valence, functional, 
        suggested_ecutwfc, etc.)
    """
    info = {}
    
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            # Read first 150 lines (header is usually at the top)
            lines = [f.readline() for _ in range(150)]
            content = ''.join(lines)
            
            # Extract element symbol (case-insensitive)
            element_match = re.search(r'<PP_INFO>.*?Element:\s*([A-Za-z]{1,2})', content, re.DOTALL | re.IGNORECASE)
            if not element_match:
                element_match = re.search(r'element\s*=\s*["\']?([A-Za-z]{1,2})["\']?', content, re.IGNORECASE)
            if element_match:
                element = element_match.group(1)
                # Normalize to proper case (Fe, Si, C, etc.)
                info['element'] = element[0].upper() + element[1:].lower() if len(element) > 1 else element.upper()
            
            # Extract valence charge
            z_match = re.search(r'z_valence\s*=\s*["\']?([\d.]+)["\']?', content, re.IGNORECASE)
            if z_match:
                info['z_valence'] = float(z_match.group(1))
            
            # Extract suggested cutoff energy (in Ry)
            # Try different formats:
            # UPF v2 format: suggested_ecutwfc="40.0"
            # SSSP format: Suggested minimum cutoff for wavefunctions:  29. Ry
            ecut_match = re.search(r'suggested_ecutwfc\s*=\s*["\']?([\d.]+)["\']?', content, re.IGNORECASE)
            if not ecut_match:
                # Try SSSP format
                ecut_match = re.search(r'Suggested minimum cutoff for wavefunctions:\s*([\d.]+)', content, re.IGNORECASE)
            
            if ecut_match:
                try:
                    value = float(ecut_match.group(1))
                    # Only store if value > 0 (some files have 0.0 meaning no suggestion)
                    if value > 0:
                        info['suggested_ecutwfc'] = value
                except (ValueError, IndexError):
                    pass
            
            # Extract functional
            functional_match = re.search(r'functional\s*=\s*["\']?([A-Za-z0-9-]+)["\']?', content, re.IGNORECASE)
            if functional_match:
                functional = functional_match.group(1).upper()
                if 'PBE' in functional:
                    info['functional'] = 'PBE'
                elif 'LDA' in functional:
                    info['functional'] = 'LDA'
            
            # Extract pseudopotential type from header
            # UPF files use XML-like headers with type attribute
            # Look for patterns like: type="PAW", type="ultrasoft", type="NC"
            
            # First try XML-style type attribute (most reliable)
            type_match = re.search(r'type\s*=\s*["\']([^"\']+)["\']', content, re.IGNORECASE)
            if type_match:
                pseudo_type = type_match.group(1).strip().upper()
                if 'PAW' in pseudo_type or 'PROJECTOR' in pseudo_type:
                    info['type'] = 'PAW'
                elif 'US' in pseudo_type or 'ULTRASOFT' in pseudo_type:
                    info['type'] = 'Ultrasoft'
                elif 'NC' in pseudo_type or 'NORM' in pseudo_type:
                    info['type'] = 'Norm-conserving'
            
            # If not found in type attribute, search descriptive text markers
            if 'type' not in info:
                content_upper = content.upper()
                # Order matters: check PAW first (most specific), then ultrasoft, then NC
                # XML-like or descriptive markers
                if re.search(r'PROJECTOR\s+AUGMENTED|\bPAW\b', content_upper) and 'PAWXML' not in content_upper:
                    info['type'] = 'PAW'
                elif re.search(r'ULTRASOFT|RRKJUS|USPP', content_upper):
                    info['type'] = 'Ultrasoft'
                elif re.search(r'NORM-?CONSERVING|ONCV|ONCVPSP|ONCV_PSP|ONCVPSP', content_upper):
                    info['type'] = 'Norm-conserving'

            # Additional lookups: some UPF variants include 'pseudo_type' or similar attributes
            if 'type' not in info:
                extra_match = re.search(r'(pseudo[-_ ]?type|pseudotype|pseudo_type)\s*=\s*["\']?([^"\'\s>]+)', content, re.IGNORECASE)
                if extra_match:
                    pseudo_type = extra_match.group(2).strip().upper()
                    if 'PAW' in pseudo_type:
                        info['type'] = 'PAW'
                    elif 'US' in pseudo_type or 'ULTRASOFT' in pseudo_type or 'USPP' in pseudo_type:
                        info['type'] = 'Ultrasoft'
                    elif 'NC' in pseudo_type or 'NORM' in pseudo_type or 'ONCV' in pseudo_type:
                        info['type'] = 'Norm-conserving'
    
    except Exception as e:
        # If parsing fails, return empty dict (gracefully handle errors)
        pass
    
    return info


def detect_pseudopotentials(directory: str, recursive: bool = True) -> Dict[str, Dict]:
    """
    Detect all pseudopotentials in a directory and extract their information.
    
    Args:
        directory: Directory path to search
        recursive: Whether to search recursively
    
    Returns:
        Dictionary mapping element to pseudopotential info dict
    """
    upf_files = detect_upf_files(directory, recursive)
    pseudopotentials = {}
    
    for upf_path in upf_files:
        filename = os.path.basename(upf_path)
        
        # Extract element from filename first
        element = extract_element_from_filename(filename)
        
        # Try to parse UPF file for more accurate info
        header_info = parse_upf_header(upf_path)
        
        # Use header info if available, otherwise use filename-based detection
        if 'element' in header_info:
            element = header_info['element']
        
        if not element:
            # Skip if we can't determine the element
            continue
        
        # Build pseudopotential info
        pseudo_info = {
            'filename': filename,
            'path': upf_path,
            'element': element
        }
        
        # Add functional
        functional = header_info.get('functional') or extract_functional_from_filename(filename)
        if functional:
            pseudo_info['functional'] = functional
        
        # Add type (only from UPF header parsing, don't use filename-based detection)
        if 'type' in header_info:
            pseudo_info['type'] = header_info['type']
        
        # Add z_valence
        if 'z_valence' in header_info:
            pseudo_info['z_valence'] = header_info['z_valence']
        
        # If element already exists, keep the first one found
        # (or we could implement logic to prefer certain types)
        if element not in pseudopotentials:
            pseudopotentials[element] = pseudo_info
    
    return pseudopotentials


def detect_pseudopotentials_remote(directory: str, 
                                   ssh_connection: Dict,
                                   recursive: bool = True) -> Dict[str, Dict]:
    """
    Detect pseudopotentials on a remote machine via SSH.
    
    Handles case-insensitive file extensions (.UPF, .upf, .Upf, etc.).
    
    Args:
        directory: Directory path on remote machine
        ssh_connection: SSH connection info with keys: 'host', 'username', 'port'
        recursive: Whether to search recursively
    
    Returns:
        Dictionary mapping element to pseudopotential info dict
    """
    import subprocess
    
    host = ssh_connection.get('host')
    username = ssh_connection.get('username', os.environ.get('USER'))
    port = ssh_connection.get('port', 22)
    
    try:
        # Find UPF files remotely (case-insensitive)
        if recursive:
            find_cmd = f"find {directory} -iname '*.upf'"
        else:
            find_cmd = f"find {directory} -maxdepth 1 -iname '*.upf'"
        
        ssh_cmd = f"ssh -p {port} {username}@{host} '{find_cmd}'"
        result = subprocess.run(ssh_cmd, shell=True, capture_output=True, text=True, timeout=30)
        
        if result.returncode != 0:
            return {}
        
        upf_files = [f.strip() for f in result.stdout.strip().split('\n') if f.strip()]
        
        pseudopotentials = {}
        for upf_path in upf_files:
            filename = os.path.basename(upf_path)
            element = extract_element_from_filename(filename)
            
            if not element:
                continue
            
            pseudo_info = {
                'filename': filename,
                'path': upf_path,
                'element': element
            }
            
            functional = extract_functional_from_filename(filename)
            if functional:
                pseudo_info['functional'] = functional
            
            if element not in pseudopotentials:
                pseudopotentials[element] = pseudo_info
        
        return pseudopotentials
    
    except Exception as e:
        print(f"Error detecting remote pseudopotentials: {e}")
        return {}

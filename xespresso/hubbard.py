"""
hubbard.py

Utilities for handling Hubbard parameters in both old and new QE formats.

This module supports:
- Old format (QE < 7.0): Hubbard_U(i), Hubbard_V(na,nb,k) in SYSTEM namelist
- New format (QE >= 7.0): HUBBARD card with explicit orbital specification

The new format is more flexible and recommended for QE 7.0+.
"""

from typing import Dict, List, Optional, Tuple, Any


class HubbardConfig:
    """
    Configuration for Hubbard parameters supporting both old and new QE formats.
    
    Attributes:
        use_new_format: If True, use HUBBARD card (QE >= 7.0). If False, use old format.
        projector: Projector type for new format ('atomic', 'ortho-atomic', 'norm-atomic', 'wf', 'pseudo')
        u_params: Dictionary of U parameters {species-orbital: value}
        v_params: List of V parameters [(species1-orbital1, species2-orbital2, i, j, value)]
        j_params: Dictionary of J parameters (for old format)
        alpha_params: Dictionary of alpha parameters (for old format)
        beta_params: Dictionary of beta parameters (for old format)
    """
    
    def __init__(self, use_new_format: bool = None, projector: str = 'atomic'):
        """
        Initialize HubbardConfig.
        
        Args:
            use_new_format: If None, auto-detect based on parameters. 
                          If True, force new format. If False, force old format.
            projector: Projector type for new format (default: 'atomic')
        """
        self.use_new_format = use_new_format
        self.projector = projector
        
        # Parameters storage
        self.u_params = {}  # {species-orbital: value} or {species: value} for old format
        self.v_params = []  # [(species1-orbital1, species2-orbital2, i, j, value)] or old format
        self.j_params = {}
        self.alpha_params = {}
        self.beta_params = {}
    
    def add_u(self, species: str, value: float, orbital: Optional[str] = None):
        """
        Add a U parameter.
        
        Args:
            species: Species name (e.g., 'Fe1', 'Mn')
            value: U value in eV
            orbital: Orbital specification for new format (e.g., '3d', '4f')
        """
        if orbital and self.use_new_format is not False:
            key = f"{species}-{orbital}"
        else:
            key = species
        self.u_params[key] = value
    
    def add_v(self, species1: str, species2: str, value: float,
             orbital1: Optional[str] = None, orbital2: Optional[str] = None,
             i: int = 1, j: int = 1):
        """
        Add a V parameter (inter-site interaction).
        
        Args:
            species1: First species name
            species2: Second species name
            value: V value in eV
            orbital1: First orbital for new format (e.g., '3d')
            orbital2: Second orbital for new format (e.g., '3d')
            i: Site index 1 (for new format, typically 1)
            j: Site index 2 (for new format, typically 1)
        """
        if orbital1 and orbital2 and self.use_new_format is not False:
            # New format
            self.v_params.append((f"{species1}-{orbital1}", f"{species2}-{orbital2}", 
                                i, j, value))
        else:
            # Old format - store as tuple
            self.v_params.append((species1, species2, i, j, value))
    
    def add_j(self, species: str, value: float, ityp: Optional[int] = None):
        """Add J parameter (old format only)."""
        if ityp:
            key = (species, ityp)
        else:
            key = species
        self.j_params[key] = value
    
    def add_alpha(self, species: str, value: float):
        """Add alpha parameter (old format only)."""
        self.alpha_params[species] = value
    
    def add_beta(self, species: str, value: float):
        """Add beta parameter (old format only)."""
        self.beta_params[species] = value
    
    def should_use_new_format(self) -> bool:
        """
        Determine if new format should be used.
        
        Returns:
            True if new format should be used, False otherwise
        """
        if self.use_new_format is not None:
            return self.use_new_format
        
        # Auto-detect: use new format if any orbital specifications are present
        for key in self.u_params.keys():
            if '-' in key:  # Has orbital specification like 'Fe1-3d'
                return True
        
        for v_param in self.v_params:
            if '-' in v_param[0]:  # Has orbital specification
                return True
        
        return False
    
    def to_old_format_dict(self, species_info: Dict) -> Dict[str, Any]:
        """
        Convert to old format dictionary for SYSTEM namelist.
        
        Args:
            species_info: Dictionary with species information including indices
        
        Returns:
            Dictionary of parameters to add to SYSTEM namelist
        """
        params = {}
        
        # U parameters
        for species, value in self.u_params.items():
            # Remove orbital suffix if present
            species_base = species.split('-')[0]
            if species_base in species_info:
                idx = species_info[species_base]['index']
                params[f'Hubbard_U({idx})'] = value
        
        # V parameters (old format with indices)
        for v_param in self.v_params:
            if len(v_param) == 5:
                species1, species2, i, j, value = v_param
                species1_base = species1.split('-')[0]
                species2_base = species2.split('-')[0]
                
                # Get atom indices (this is a simplified version)
                # In practice, you'd need to map to actual atom indices
                params[f'Hubbard_V({i},{j},1)'] = value
        
        # J parameters
        for (species, ityp), value in self.j_params.items():
            if species in species_info:
                idx = species_info[species]['index']
                params[f'Hubbard_J({ityp},{idx})'] = value
        
        # Alpha parameters
        for species, value in self.alpha_params.items():
            if species in species_info:
                idx = species_info[species]['index']
                params[f'Hubbard_alpha({idx})'] = value
        
        # Beta parameters
        for species, value in self.beta_params.items():
            if species in species_info:
                idx = species_info[species]['index']
                params[f'Hubbard_beta({idx})'] = value
        
        return params
    
    def to_new_format_card(self) -> List[str]:
        """
        Convert to new format HUBBARD card.
        
        Returns:
            List of strings forming the HUBBARD card
        """
        lines = []
        lines.append(f"HUBBARD {{{self.projector}}}\n")
        
        # U parameters
        for species_orbital, value in self.u_params.items():
            lines.append(f"  U {species_orbital} {value}\n")
        
        # V parameters
        for v_param in self.v_params:
            if len(v_param) == 5:
                spec1_orb, spec2_orb, i, j, value = v_param
                lines.append(f"  V {spec1_orb} {spec2_orb} {i} {j} {value}\n")
        
        return lines
    
    @classmethod
    def from_input_data(cls, input_data: Dict, qe_version: Optional[str] = None) -> 'HubbardConfig':
        """
        Create HubbardConfig from input_data dictionary.
        
        Args:
            input_data: Dictionary with input parameters
            qe_version: QE version string (e.g., '7.2', '6.8')
        
        Returns:
            HubbardConfig object
        """
        # Determine format based on QE version if provided
        use_new_format = None
        if qe_version:
            try:
                major, minor = map(int, qe_version.split('.')[:2])
                use_new_format = (major >= 7)
            except (ValueError, AttributeError):
                pass
        
        # Check for explicit format specification
        if 'hubbard_format' in input_data:
            use_new_format = (input_data['hubbard_format'] == 'card')
        
        # Check for new format data
        if 'hubbard' in input_data and isinstance(input_data['hubbard'], dict):
            use_new_format = True
        
        config = cls(use_new_format=use_new_format)
        
        # Parse new format
        if 'hubbard' in input_data and isinstance(input_data['hubbard'], dict):
            hubbard_data = input_data['hubbard']
            
            if 'projector' in hubbard_data:
                config.projector = hubbard_data['projector']
            
            # U parameters
            if 'u' in hubbard_data:
                for spec_orb, value in hubbard_data['u'].items():
                    parts = spec_orb.split('-')
                    if len(parts) == 2:
                        config.add_u(parts[0], value, parts[1])
                    else:
                        config.add_u(spec_orb, value)
            
            # V parameters
            if 'v' in hubbard_data:
                for v_spec in hubbard_data['v']:
                    config.add_v(
                        v_spec.get('species1', ''),
                        v_spec.get('species2', ''),
                        v_spec.get('value', 0.0),
                        v_spec.get('orbital1'),
                        v_spec.get('orbital2'),
                        v_spec.get('i', 1),
                        v_spec.get('j', 1)
                    )
        
        # Parse old format from INPUT_NTYP
        if 'input_ntyp' in input_data or 'INPUT_NTYP' in input_data:
            input_ntyp = input_data.get('input_ntyp') or input_data.get('INPUT_NTYP')
            
            if 'Hubbard_U' in input_ntyp:
                for species, value in input_ntyp['Hubbard_U'].items():
                    config.add_u(species, value)
            
            if 'Hubbard_J' in input_ntyp:
                for species, value in input_ntyp['Hubbard_J'].items():
                    config.add_j(species, value)
            
            if 'Hubbard_alpha' in input_ntyp:
                for species, value in input_ntyp['Hubbard_alpha'].items():
                    config.add_alpha(species, value)
            
            if 'Hubbard_beta' in input_ntyp:
                for species, value in input_ntyp['Hubbard_beta'].items():
                    config.add_beta(species, value)
        
        # Parse old format hubbard_v
        if 'hubbard_v' in input_data:
            for key, value in input_data['hubbard_v'].items():
                # Parse key like '(1,1,1)' or '(3,3,1)'
                try:
                    parts = key.strip('()').split(',')
                    if len(parts) >= 3:
                        na, nb, k = map(int, parts[:3])
                        config.v_params.append(('', '', na, nb, value))
                except (ValueError, AttributeError):
                    pass
        
        return config


def build_hubbard_str(input_data: Dict, species_info: Dict, 
                     qe_version: Optional[str] = None) -> List[str]:
    """
    Build Hubbard parameters section for QE input file.
    
    Args:
        input_data: Dictionary with input parameters
        species_info: Dictionary with species information
        qe_version: QE version string (e.g., '7.2')
    
    Returns:
        List of strings for the Hubbard section
    """
    # Check if there are any Hubbard parameters
    has_hubbard = any([
        'hubbard' in input_data,
        'input_ntyp' in input_data and 'Hubbard_U' in input_data.get('input_ntyp', {}),
        'INPUT_NTYP' in input_data and 'Hubbard_U' in input_data.get('INPUT_NTYP', {}),
        'hubbard_v' in input_data,
    ])
    
    if not has_hubbard:
        return []
    
    # Create HubbardConfig from input data
    config = HubbardConfig.from_input_data(input_data, qe_version)
    
    # Return appropriate format
    if config.should_use_new_format():
        return config.to_new_format_card()
    else:
        # Old format is handled in build_section_str by adding to SYSTEM namelist
        # So we return empty here
        return []


def apply_hubbard_to_system(input_parameters: Dict, input_data: Dict, 
                           species_info: Dict, qe_version: Optional[str] = None) -> Dict:
    """
    Apply Hubbard parameters to SYSTEM namelist (for old format).
    
    This function modifies input_parameters['system'] in place.
    
    Args:
        input_parameters: Dictionary with namelist parameters
        input_data: Dictionary with input data
        species_info: Dictionary with species information
        qe_version: QE version string
    
    Returns:
        Modified input_parameters
    """
    config = HubbardConfig.from_input_data(input_data, qe_version)
    
    # Only apply if using old format
    if not config.should_use_new_format():
        old_params = config.to_old_format_dict(species_info)
        input_parameters['system'].update(old_params)
    
    return input_parameters

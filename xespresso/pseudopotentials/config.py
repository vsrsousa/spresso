"""
config.py

Data structures for pseudopotential configuration management.
"""

import json
from typing import Dict, Optional, List
from dataclasses import dataclass, asdict, field


@dataclass
class Pseudopotential:
    """
    Represents a single pseudopotential file for an element.
    
    Attributes:
        element: Chemical symbol (e.g., 'Fe', 'Si', 'O')
        filename: Pseudopotential filename (e.g., 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF')
        path: Full path to pseudopotential file
        functional: Exchange-correlation functional (e.g., 'PBE', 'LDA', 'PBEsol')
        type: Pseudopotential type (e.g., 'ultrasoft', 'paw', 'norm-conserving')
        z_valence: Number of valence electrons
    """
    element: str
    filename: str
    path: str
    functional: Optional[str] = None
    type: Optional[str] = None
    z_valence: Optional[float] = None
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {k: v for k, v in asdict(self).items() if v is not None}
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'Pseudopotential':
        """Create Pseudopotential from dictionary."""
        return cls(**data)


@dataclass
class PseudopotentialsConfig:
    """
    Configuration for pseudopotentials.
    
    Pseudopotentials are primarily stored on the local user computer after downloading.
    The machine_name is optional and only needed if pseudopotentials are stored on a remote machine.
    
    Supports multiple pseudopotential libraries/sets with version tracking.
    
    Attributes:
        name: Configuration name/label (e.g., 'SSSP_efficiency', 'pbe_standard')
        base_path: Base directory containing pseudopotential files (local or remote path)
        pseudopotentials: Dictionary mapping element to Pseudopotential object
        machine_name: Optional name of remote machine (if pseudopotentials are on remote system)
        description: Optional description of this pseudopotential set
        functional: Primary functional used (e.g., 'PBE', 'LDA')
        library: Pseudopotential library name (e.g., 'SSSP', 'PSLibrary', 'Pseudo Dojo')
        version: Library version (e.g., '1.1.2', '1.0.0', '0.4')
    """
    name: str
    base_path: str
    pseudopotentials: Dict[str, Pseudopotential] = field(default_factory=dict)
    machine_name: Optional[str] = None
    description: Optional[str] = None
    functional: Optional[str] = None
    library: Optional[str] = None
    version: Optional[str] = None
    
    def add_pseudopotential(self, pseudo: Pseudopotential):
        """
        Add a pseudopotential to the configuration.
        
        Args:
            pseudo: Pseudopotential object to add
        """
        self.pseudopotentials[pseudo.element] = pseudo
    
    def get_pseudopotential(self, element: str) -> Optional[Pseudopotential]:
        """
        Get a pseudopotential by element symbol.
        
        Args:
            element: Chemical symbol (e.g., 'Fe', 'Si')
        
        Returns:
            Pseudopotential object or None if not found
        """
        return self.pseudopotentials.get(element)
    
    def list_elements(self) -> List[str]:
        """
        Get list of available elements.
        
        Returns:
            List of element symbols
        """
        return sorted(self.pseudopotentials.keys())
    
    def get_pseudo_group_dict(self) -> Dict[str, str]:
        """
        Get pseudopotentials in xespresso pseudo_group format.
        
        Returns a dictionary mapping UPPERCASE element symbols to filenames,
        compatible with xespresso's pseudo_gropus dictionary format.
        This allows the configuration to be used with xespresso's pseudo_group parameter.
        
        Returns:
            Dictionary mapping element symbols (uppercase) to pseudopotential filenames
            
        Example:
            >>> config = load_pseudopotentials_config("SSSP_efficiency")
            >>> pseudo_dict = config.get_pseudo_group_dict()
            >>> # Returns: {'FE': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF', 'O': 'O.pbe-n-kjpaw_psl.0.1.UPF', ...}
        """
        result = {}
        for element, pseudo in self.pseudopotentials.items():
            if isinstance(pseudo, Pseudopotential):
                # Use uppercase for consistency with xespresso's pseudo_gropus format
                result[element.upper()] = pseudo.filename
            else:
                # Handle case where pseudo might be a dict (from deserialization)
                filename = pseudo.get('filename') if isinstance(pseudo, dict) else str(pseudo)
                result[element.upper()] = filename
        return result
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        result = {
            "name": self.name,
            "base_path": self.base_path,
            "pseudopotentials": {}
        }
        
        if self.machine_name:
            result["machine_name"] = self.machine_name
        if self.description:
            result["description"] = self.description
        if self.functional:
            result["functional"] = self.functional
        if self.library:
            result["library"] = self.library
        if self.version:
            result["version"] = self.version
        
        # Convert pseudopotentials to dict
        for element, pseudo in self.pseudopotentials.items():
            if isinstance(pseudo, Pseudopotential):
                result["pseudopotentials"][element] = pseudo.to_dict()
            else:
                result["pseudopotentials"][element] = pseudo
        
        return result
    
    def to_json(self, filepath: str):
        """
        Save configuration to JSON file.
        
        Args:
            filepath: Path to save the configuration
        """
        import os
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'PseudopotentialsConfig':
        """
        Create PseudopotentialsConfig from dictionary.
        
        Args:
            data: Dictionary containing configuration data
        
        Returns:
            PseudopotentialsConfig object
        """
        # Extract pseudopotentials and convert them
        pseudos_data = data.pop('pseudopotentials', {})
        
        config = cls(**data)
        
        # Add pseudopotentials
        for element, pseudo_data in pseudos_data.items():
            if isinstance(pseudo_data, Pseudopotential):
                config.pseudopotentials[element] = pseudo_data
            else:
                config.pseudopotentials[element] = Pseudopotential.from_dict(pseudo_data)
        
        return config
    
    @classmethod
    def from_json(cls, filepath: str) -> 'PseudopotentialsConfig':
        """
        Load configuration from JSON file.
        
        Args:
            filepath: Path to JSON configuration file
        
        Returns:
            PseudopotentialsConfig object
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        return cls.from_dict(data)
    
    def get_pseudopotentials_dict(self) -> Dict[str, str]:
        """
        Get pseudopotentials as a simple element->filename dictionary.
        
        This format is compatible with xespresso calculation inputs.
        
        Returns:
            Dictionary mapping element to filename
        """
        return {element: pseudo.filename for element, pseudo in self.pseudopotentials.items()}

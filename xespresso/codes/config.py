"""
config.py

Data structures for Quantum ESPRESSO code configuration.
"""

import json
from typing import Dict, Optional, List
from dataclasses import dataclass, asdict, field


@dataclass
class Code:
    """
    Represents a single Quantum ESPRESSO executable.
    
    Attributes:
        name: Code name (e.g., 'pw', 'hp', 'dos', 'bands')
        path: Full path to executable (e.g., '/usr/bin/pw.x')
        version: QE version string (e.g., '7.2', '6.8')
        parallel_command: Optional MPI command prefix (e.g., 'mpirun -np {nprocs}')
        default_parallel: Default parallelization options (e.g., '-nk 4')
    """
    name: str
    path: str
    version: Optional[str] = None
    parallel_command: Optional[str] = None
    default_parallel: Optional[str] = None
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {k: v for k, v in asdict(self).items() if v is not None}
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'Code':
        """Create Code from dictionary."""
        return cls(**data)


@dataclass
class CodesConfig:
    """
    Configuration for all QE codes on a machine.
    
    Supports multiple QE versions on the same machine. Codes can be organized by version,
    allowing users to choose which version to use for each calculation.
    
    Attributes:
        machine_name: Name of the machine
        qe_prefix: Common prefix for QE executables (e.g., '/opt/qe-7.2/bin')
        qe_version: Default QE version on this machine
        label: Optional custom label for this version (e.g., 'qe-7.2', 'qe-dev')
        codes: Dictionary mapping code name to Code object (for single version configs)
        versions: Dictionary mapping version string to version-specific configuration
                  (for multi-version configs)
        modules: Optional list of modules to load (e.g., ['quantum-espresso/7.2'])
        environment: Optional environment variables to set
    
    Multi-version structure:
        versions = {
            "7.2": {
                "qe_prefix": "/opt/qe-7.2/bin",
                "modules": ["quantum-espresso/7.2"],
                "codes": {...}
            },
            "6.8": {
                "qe_prefix": "/opt/qe-6.8/bin",
                "modules": ["quantum-espresso/6.8"],
                "codes": {...}
            }
        }
    """
    machine_name: str
    codes: Dict[str, Code] = field(default_factory=dict)
    qe_prefix: Optional[str] = None
    qe_version: Optional[str] = None
    label: Optional[str] = None
    modules: Optional[List[str]] = None
    environment: Optional[Dict[str, str]] = None
    versions: Optional[Dict[str, Dict]] = None
    
    def add_code(self, code: Code, version: Optional[str] = None):
        """
        Add a code to the configuration.
        
        Args:
            code: Code object to add
            version: Optional version string. If provided, adds to version-specific codes.
        """
        if version:
            # Add to version-specific configuration
            if not self.versions:
                self.versions = {}
            if version not in self.versions:
                self.versions[version] = {"codes": {}}
            if "codes" not in self.versions[version]:
                self.versions[version]["codes"] = {}
            self.versions[version]["codes"][code.name] = code
        else:
            # Add to main codes dictionary
            self.codes[code.name] = code
    
    def get_code(self, name: str, version: Optional[str] = None) -> Optional[Code]:
        """
        Get a code by name, optionally from a specific version.
        
        Args:
            name: Code name (e.g., 'pw', 'hp')
            version: Optional version string. If None, uses default version or main codes.
        
        Returns:
            Code object or None if not found
        """
        if version and self.versions and version in self.versions:
            # Get from version-specific codes
            codes_data = self.versions[version].get("codes", {})
            if name in codes_data:
                code_data = codes_data[name]
                if isinstance(code_data, Code):
                    return code_data
                else:
                    return Code.from_dict(code_data)
        
        # Fall back to main codes dictionary
        return self.codes.get(name)
    
    def get_all_codes(self, version: Optional[str] = None) -> Dict[str, Code]:
        """
        Get all codes as a dictionary, from either main codes or a specific version.
        
        Args:
            version: Optional version string. If provided, gets codes from that version.
                    If None, gets codes from main codes dict or qe_version if available.
        
        Returns:
            Dictionary mapping code name to Code object
        """
        # If no version specified, try to use qe_version
        if version is None and self.qe_version and self.versions and self.qe_version in self.versions:
            version = self.qe_version
        
        if version and self.versions and version in self.versions:
            # Get from version-specific codes
            codes_data = self.versions[version].get("codes", {})
            result = {}
            for name, code_data in codes_data.items():
                if isinstance(code_data, Code):
                    result[name] = code_data
                else:
                    result[name] = Code.from_dict(code_data)
            return result
        
        # Return main codes dictionary
        return self.codes.copy()
    
    def has_any_codes(self) -> bool:
        """
        Check if there are any codes configured, in either main codes or versions.
        
        Returns:
            True if any codes are configured
        """
        if self.codes:
            return True
        
        if self.versions:
            for version_config in self.versions.values():
                if version_config.get('codes'):
                    return True
        
        return False
    
    def list_versions(self) -> List[str]:
        """List all available QE versions."""
        if self.versions:
            return list(self.versions.keys())
        elif self.qe_version:
            return [self.qe_version]
        return []
    
    def get_version_config(self, version: str) -> Optional[Dict]:
        """
        Get version-specific configuration.
        
        Args:
            version: Version string
        
        Returns:
            Dictionary with version-specific config or None
        """
        if self.versions and version in self.versions:
            return self.versions[version]
        return None
    
    def to_dict(self) -> Dict:
        """
        Convert CodesConfig to dictionary for JSON serialization.
        
        This method implements smart redundancy reduction:
        - Top-level 'modules' are excluded when versions structure exists
        - Top-level 'qe_version' is kept to indicate default/active version
        - Top-level 'label' excluded when versions exist (stored per-version)
        
        Returns:
            Dictionary representation suitable for JSON serialization
        """
        result = {
            'machine_name': self.machine_name,
            'codes': {name: code.to_dict() for name, code in self.codes.items()},
        }
        if self.qe_prefix:
            result['qe_prefix'] = self.qe_prefix
        
        # qe_version can be kept at top level as it indicates the default/active version
        # even when multiple versions exist
        if self.qe_version:
            result['qe_version'] = self.qe_version
        
        # Only include top-level label and modules if no versions structure exists
        # (backward compatibility for single-version configs)
        has_versions = self.versions and len(self.versions) > 0
        
        if self.label and not has_versions:
            result['label'] = self.label
        if self.modules and not has_versions:
            result['modules'] = self.modules
            
        if self.environment:
            result['environment'] = self.environment
        if self.versions:
            # Serialize versions
            result['versions'] = {}
            for ver, ver_config in self.versions.items():
                result['versions'][ver] = {}
                for key, value in ver_config.items():
                    if key == 'codes':
                        # Serialize codes in version
                        result['versions'][ver]['codes'] = {
                            name: (code.to_dict() if isinstance(code, Code) else code)
                            for name, code in value.items()
                        }
                    else:
                        result['versions'][ver][key] = value
        return result
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'CodesConfig':
        """Create CodesConfig from dictionary."""
        # Make a copy to avoid modifying the original
        data = data.copy()
        
        # Handle main codes
        codes_data = data.pop('codes', {})
        codes = {name: Code.from_dict(code_data) if isinstance(code_data, dict) else code_data 
                 for name, code_data in codes_data.items()}
        
        # Handle versions
        versions_data = data.pop('versions', None)
        versions = None
        if versions_data:
            versions = {}
            for ver, ver_config in versions_data.items():
                versions[ver] = ver_config.copy()
                if 'codes' in ver_config:
                    versions[ver]['codes'] = {
                        name: Code.from_dict(code_data) if isinstance(code_data, dict) else code_data
                        for name, code_data in ver_config['codes'].items()
                    }
        
        return cls(codes=codes, versions=versions, **data)
    
    def to_json(self, filepath: str, indent: int = 2):
        """Save to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=indent)
    
    @classmethod
    def from_json(cls, filepath: str) -> 'CodesConfig':
        """Load from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)
    
    def list_codes(self, version: Optional[str] = None) -> List[str]:
        """
        List all available code names.
        
        Args:
            version: Optional version string. If provided, lists codes for that version.
        
        Returns:
            List of code names
        """
        if version and self.versions and version in self.versions:
            codes = self.versions[version].get('codes', {})
            return list(codes.keys())
        return list(self.codes.keys())
    
    def has_code(self, name: str, version: Optional[str] = None) -> bool:
        """
        Check if a code is available.
        
        Args:
            name: Code name
            version: Optional version string. If None, checks main codes and all versions.
        
        Returns:
            True if code is available
        """
        if version and self.versions and version in self.versions:
            codes = self.versions[version].get('codes', {})
            return name in codes
        
        # Check main codes
        if name in self.codes:
            return True
        
        # If no version specified, also check if code exists in any version
        if not version and self.versions:
            for version_config in self.versions.values():
                if 'codes' in version_config and name in version_config['codes']:
                    return True
        
        return False

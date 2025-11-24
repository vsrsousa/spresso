"""
manager.py

Utilities for managing pseudopotential configurations.
"""

import os
from typing import Dict, Optional, List
from pathlib import Path

from .config import Pseudopotential, PseudopotentialsConfig
from .detector import detect_pseudopotentials, detect_pseudopotentials_remote


DEFAULT_PSEUDOPOTENTIALS_DIR = os.path.expanduser("~/.xespresso/pseudopotentials")


class PseudopotentialsManager:
    """
    Manager for pseudopotential configurations.
    
    This class provides utilities to:
    - Detect available pseudopotentials on local or remote machines
    - Create and manage pseudopotential configurations
    - Load and save configurations
    """
    
    def __init__(self, config: Optional[PseudopotentialsConfig] = None):
        """
        Initialize PseudopotentialsManager.
        
        Args:
            config: Optional PseudopotentialsConfig object
        """
        self.config = config
    
    @classmethod
    def detect_pseudopotentials(cls,
                                base_path: str,
                                ssh_connection: Optional[Dict] = None,
                                recursive: bool = True) -> Dict[str, Dict]:
        """
        Detect available pseudopotentials in a directory.
        
        Args:
            base_path: Directory containing pseudopotential files
            ssh_connection: SSH connection info for remote detection
                           {'host': 'hostname', 'username': 'user', 'port': 22}
            recursive: Whether to search recursively in subdirectories
        
        Returns:
            Dictionary mapping element to pseudopotential info
        """
        if ssh_connection:
            return detect_pseudopotentials_remote(base_path, ssh_connection, recursive)
        else:
            return detect_pseudopotentials(base_path, recursive)
    
    @classmethod
    def create_config(cls,
                     name: str,
                     base_path: str,
                     detected_pseudos: Dict[str, Dict],
                     machine_name: Optional[str] = None,
                     description: Optional[str] = None,
                     functional: Optional[str] = None,
                     library: Optional[str] = None,
                     version: Optional[str] = None) -> PseudopotentialsConfig:
        """
        Create a PseudopotentialsConfig from detected pseudopotentials.
        
        Args:
            name: Configuration name
            base_path: Base directory containing pseudopotential files
            detected_pseudos: Dictionary from detect_pseudopotentials
            machine_name: Optional name of remote machine (if pseudopotentials are on remote system)
            description: Optional description
            functional: Primary functional (e.g., 'PBE', 'LDA')
            library: Library name (e.g., 'SSSP', 'PSLibrary', 'Pseudo Dojo')
            version: Library version (e.g., '1.1.2', '1.0.0', '0.4')
        
        Returns:
            PseudopotentialsConfig object
        """
        config = PseudopotentialsConfig(
            name=name,
            base_path=base_path,
            machine_name=machine_name,
            description=description,
            functional=functional,
            library=library,
            version=version
        )
        
        for element, pseudo_info in detected_pseudos.items():
            pseudo = Pseudopotential(
                element=pseudo_info['element'],
                filename=pseudo_info['filename'],
                path=pseudo_info['path'],
                functional=pseudo_info.get('functional'),
                type=pseudo_info.get('type'),
                z_valence=pseudo_info.get('z_valence')
            )
            config.add_pseudopotential(pseudo)
        
        return config
    
    @staticmethod
    def save_config(config: PseudopotentialsConfig,
                   output_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR,
                   filename: Optional[str] = None,
                   overwrite: bool = False) -> str:
        """
        Save a PseudopotentialsConfig to JSON file.
        
        Args:
            config: PseudopotentialsConfig object to save
            output_dir: Directory to save the file
            filename: Optional filename (default: <name>.json)
            overwrite: If True, overwrites existing file without asking
        
        Returns:
            Path to saved file
        
        Raises:
            FileExistsError: If file exists and overwrite is False
        """
        os.makedirs(output_dir, exist_ok=True)
        
        if not filename:
            filename = f"{config.name}.json"
        
        filepath = os.path.join(output_dir, filename)
        
        # Check if file already exists
        if os.path.exists(filepath) and not overwrite:
            raise FileExistsError(f"Configuration file already exists: {filepath}")
        
        config.to_json(filepath)
        
        return filepath
    
    @staticmethod
    def load_config(name: str,
                   pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR) -> Optional[PseudopotentialsConfig]:
        """
        Load a PseudopotentialsConfig from JSON file.
        
        Args:
            name: Name of the configuration
            pseudopotentials_dir: Directory containing pseudopotentials configurations
        
        Returns:
            PseudopotentialsConfig object or None if not found
        """
        filepath = os.path.join(pseudopotentials_dir, f"{name}.json")
        
        if not os.path.exists(filepath):
            return None
        
        try:
            return PseudopotentialsConfig.from_json(filepath)
        except Exception as e:
            print(f"Error loading pseudopotentials config: {e}")
            return None
    
    @staticmethod
    def list_configs(pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR) -> List[str]:
        """
        List all available pseudopotential configurations.
        
        Args:
            pseudopotentials_dir: Directory containing configurations
        
        Returns:
            List of configuration names
        """
        if not os.path.exists(pseudopotentials_dir):
            return []
        
        configs = []
        for file in os.listdir(pseudopotentials_dir):
            if file.endswith('.json'):
                configs.append(file[:-5])  # Remove .json extension
        
        return sorted(configs)
    
    @staticmethod
    def delete_config(name: str,
                     pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR):
        """
        Delete a pseudopotential configuration.
        
        Args:
            name: Name of the configuration to delete
            pseudopotentials_dir: Directory containing configurations
        
        Raises:
            FileNotFoundError: If configuration doesn't exist
        """
        filepath = os.path.join(pseudopotentials_dir, f"{name}.json")
        
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Configuration '{name}' not found")
        
        os.remove(filepath)
    
    @staticmethod
    def set_default_config(config_name: str,
                          pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR) -> str:
        """
        Set a configuration as the default by copying it to default.json.
        
        Args:
            config_name: Name of the configuration to set as default
            pseudopotentials_dir: Directory containing configurations
        
        Returns:
            Path to the default.json file
        
        Raises:
            FileNotFoundError: If the specified configuration doesn't exist
        """
        source_filepath = os.path.join(pseudopotentials_dir, f"{config_name}.json")
        default_filepath = os.path.join(pseudopotentials_dir, "default.json")
        
        if not os.path.exists(source_filepath):
            raise FileNotFoundError(f"Configuration '{config_name}' not found")
        
        # Copy the configuration to default.json
        import shutil
        shutil.copy2(source_filepath, default_filepath)
        
        return default_filepath
    
    @staticmethod
    def get_default_config(pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR) -> Optional[PseudopotentialsConfig]:
        """
        Load the default pseudopotential configuration.
        
        Args:
            pseudopotentials_dir: Directory containing configurations
        
        Returns:
            PseudopotentialsConfig object or None if no default exists
        """
        return PseudopotentialsManager.load_config("default", pseudopotentials_dir)
    
    @staticmethod
    def has_default_config(pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR) -> bool:
        """
        Check if a default configuration exists.
        
        Args:
            pseudopotentials_dir: Directory containing configurations
        
        Returns:
            True if default.json exists, False otherwise
        """
        default_filepath = os.path.join(pseudopotentials_dir, "default.json")
        return os.path.exists(default_filepath)
    
    @staticmethod
    def clear_default_config(pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR):
        """
        Remove the default configuration.
        
        Args:
            pseudopotentials_dir: Directory containing configurations
        """
        default_filepath = os.path.join(pseudopotentials_dir, "default.json")
        if os.path.exists(default_filepath):
            os.remove(default_filepath)


def create_pseudopotentials_config(name: str,
                                  base_path: str,
                                  machine_name: Optional[str] = None,
                                  ssh_connection: Optional[Dict] = None,
                                  recursive: bool = True,
                                  description: Optional[str] = None,
                                  functional: Optional[str] = None,
                                  library: Optional[str] = None,
                                  version: Optional[str] = None,
                                  save: bool = True,
                                  output_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR,
                                  overwrite: bool = False) -> PseudopotentialsConfig:
    """
    Create a pseudopotentials configuration (with optional auto-save).
    
    Pseudopotentials are primarily stored on the local user computer.
    The machine_name is optional and only needed if pseudopotentials are on a remote machine.
    
    Args:
        name: Configuration name
        base_path: Directory containing pseudopotential files (local or remote path)
        machine_name: Optional name of remote machine (only if pseudopotentials are on remote system)
        ssh_connection: SSH connection info for remote detection
        recursive: Whether to search recursively
        description: Optional description
        functional: Primary functional (e.g., 'PBE', 'LDA')
        library: Library name (e.g., 'SSSP', 'PSLibrary', 'Pseudo Dojo')
        version: Library version (e.g., '1.1.2', '1.0.0', '0.4')
        save: Whether to save the configuration to file
        output_dir: Directory to save the configuration
        overwrite: If True, overwrites existing file
    
    Returns:
        PseudopotentialsConfig object
    """
    print(f"🔍 Detecting pseudopotentials in '{base_path}'...")
    
    detected_pseudos = PseudopotentialsManager.detect_pseudopotentials(
        base_path=base_path,
        ssh_connection=ssh_connection,
        recursive=recursive
    )
    
    if not detected_pseudos:
        print("⚠️  No pseudopotentials detected")
        return PseudopotentialsConfig(
            name=name,
            base_path=base_path,
            machine_name=machine_name,
            description=description,
            functional=functional,
            library=library,
            version=version
        )
    
    location = f"on remote machine '{machine_name}'" if machine_name else "locally"
    lib_info = f" ({library} v{version})" if library and version else (f" ({library})" if library else "")
    print(f"✅ Found pseudopotentials for {len(detected_pseudos)} elements {location}{lib_info}: {', '.join(sorted(detected_pseudos.keys()))}")
    
    config = PseudopotentialsManager.create_config(
        name=name,
        base_path=base_path,
        detected_pseudos=detected_pseudos,
        machine_name=machine_name,
        description=description,
        functional=functional,
        library=library,
        version=version
    )
    
    if save and config.pseudopotentials:
        try:
            filepath = PseudopotentialsManager.save_config(
                config,
                output_dir=output_dir,
                overwrite=overwrite
            )
            print(f"💾 Configuration saved to: {filepath}")
        except FileExistsError as e:
            print(f"⚠️  {e}")
    
    return config


def load_pseudopotentials_config(name: str,
                                 pseudopotentials_dir: str = DEFAULT_PSEUDOPOTENTIALS_DIR) -> Optional[PseudopotentialsConfig]:
    """
    Load a pseudopotentials configuration from file.
    
    Args:
        name: Name of the configuration
        pseudopotentials_dir: Directory containing configurations
    
    Returns:
        PseudopotentialsConfig object or None if not found
    """
    config = PseudopotentialsManager.load_config(name, pseudopotentials_dir)
    
    if config:
        location = f"on remote machine '{config.machine_name}'" if config.machine_name else "locally"
        print(f"✅ Loaded pseudopotentials configuration '{name}' ({location})")
        print(f"   Base path: {config.base_path}")
        print(f"   Elements: {len(config.pseudopotentials)} ({', '.join(config.list_elements())})")
        if config.functional:
            print(f"   Functional: {config.functional}")
        if config.library:
            lib_str = f"   Library: {config.library}"
            if config.version:
                lib_str += f" v{config.version}"
            print(lib_str)
    else:
        print(f"⚠️  No pseudopotentials configuration found for '{name}'")
    
    return config

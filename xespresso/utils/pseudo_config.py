"""
Pseudopotential configuration management for xespresso.

This module provides utilities to manage pseudopotential configurations
in JSON format stored in ~/.xespresso directory.
"""

import os
import json
from pathlib import Path
from typing import Dict, Optional


def get_config_dir() -> Path:
    """Get the xespresso configuration directory path.
    
    Returns:
        Path: Path to ~/.xespresso directory
    """
    config_dir = Path.home() / ".xespresso"
    return config_dir


def ensure_config_dir() -> Path:
    """Ensure the configuration directory exists.
    
    Returns:
        Path: Path to the configuration directory
    """
    config_dir = get_config_dir()
    config_dir.mkdir(exist_ok=True)
    return config_dir


def save_pseudo_config(config_name: str, config_data: Dict, overwrite: bool = False) -> None:
    """Save a pseudopotential configuration to JSON file.
    
    Args:
        config_name: Name of the configuration (will be used as filename)
        config_data: Dictionary containing the pseudopotential configuration
        overwrite: Whether to overwrite existing configuration
        
    Raises:
        FileExistsError: If configuration exists and overwrite is False
    """
    config_dir = ensure_config_dir()
    config_file = config_dir / f"{config_name}.json"
    
    if config_file.exists() and not overwrite:
        raise FileExistsError(
            f"Configuration '{config_name}' already exists. "
            "Set overwrite=True to replace it."
        )
    
    with open(config_file, 'w') as f:
        json.dump(config_data, f, indent=2)


def load_pseudo_config(config_name: str) -> Dict:
    """Load a pseudopotential configuration from JSON file.
    
    Args:
        config_name: Name of the configuration to load
        
    Returns:
        Dict: The loaded configuration
        
    Raises:
        FileNotFoundError: If configuration does not exist
    """
    config_dir = get_config_dir()
    config_file = config_dir / f"{config_name}.json"
    
    if not config_file.exists():
        raise FileNotFoundError(
            f"Configuration '{config_name}' not found in {config_dir}"
        )
    
    with open(config_file, 'r') as f:
        return json.load(f)


def list_pseudo_configs() -> list:
    """List all available pseudopotential configurations.
    
    Returns:
        list: List of configuration names (without .json extension)
    """
    config_dir = get_config_dir()
    if not config_dir.exists():
        return []
    
    return [f.stem for f in config_dir.glob("*.json")]


def delete_pseudo_config(config_name: str) -> None:
    """Delete a pseudopotential configuration.
    
    Args:
        config_name: Name of the configuration to delete
        
    Raises:
        FileNotFoundError: If configuration does not exist
    """
    config_dir = get_config_dir()
    config_file = config_dir / f"{config_name}.json"
    
    if not config_file.exists():
        raise FileNotFoundError(
            f"Configuration '{config_name}' not found in {config_dir}"
        )
    
    config_file.unlink()


def get_pseudo_info(config_name: str, element: str) -> Optional[str]:
    """Get pseudopotential filename for a specific element from a configuration.
    
    Args:
        config_name: Name of the configuration
        element: Chemical symbol of the element
        
    Returns:
        Optional[str]: Pseudopotential filename, or None if not found
    """
    try:
        config = load_pseudo_config(config_name)
        return config.get('pseudopotentials', {}).get(element)
    except FileNotFoundError:
        return None

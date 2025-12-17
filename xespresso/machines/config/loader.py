"""
loader.py

Utility for loading machine configurations for xespresso workflows.

This module supports:
- Parsing machine profiles from a JSON config file
- Loading individual machine JSON files from a machines directory
- Machine class-based configuration with validation
- Local and remote execution modes
- Key-based SSH authentication only (password-based authentication is no longer supported)
- Optional job resources, environment setup, and cleanup commands
- Normalization of script blocks (prepend/postpend) to ensure compatibility
- Interactive fallback if machine name is invalid
- Default machine configuration
- Logging and warnings for traceability

Default config path: ~/.xespresso/machines.json
Default machines directory: ~/.xespresso/machines/
Default machine name: "local_desktop"

Example usage:
from xespresso.machines.config.loader import load_machine
queue = load_machine()  # Loads default machine
machine_obj = load_machine(return_object=True)  # Returns Machine object
"""

import os
import json
from typing import Union, Optional
try:
    from xespresso.utils import warnings as warnings
    from xespresso.utils.logging import get_logger
    from xespresso.machines.machine import Machine
    logger = get_logger()
    warnings.apply_custom_format()
except ImportError:
    # Fallback to standard Python warnings and logging
    import warnings
    import logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)
    from xespresso.machines.machine import Machine

DEFAULT_CONFIG_PATH = os.path.expanduser("~/.xespresso/machines.json")
DEFAULT_MACHINES_DIR = os.path.expanduser("~/.xespresso/machines")
DEFAULT_MACHINE_NAME = "local_desktop"

def normalize_script_block(block):
    """
    Ensures that script blocks (prepend/postpend) are returned as strings.
    Accepts either a string or a list of strings.
    """
    if isinstance(block, list):
        return "\n".join(block)
    return block or ""


def sanitize_string_value(value):
    """
    Remove leading and trailing quotes from string values that may have been
    inadvertently added during manual JSON editing or serialization.
    
    This fixes issues where commands like 'export OMP_NUM_THREADS=1' were
    stored as '"export OMP_NUM_THREADS=1"' in JSON, causing quotes to appear
    in the generated job scripts.
    
    Parameters:
        value: String value to sanitize, or None
        
    Returns:
        Sanitized string with outer quotes removed, or original value if not a string
    """
    if not isinstance(value, str):
        return value
    
    # Strip leading and trailing whitespace first
    value = value.strip()
    
    # Remove outer quotes if present (both single and double)
    if (value.startswith('"') and value.endswith('"')) or \
       (value.startswith("'") and value.endswith("'")):
        value = value[1:-1]
    
    return value


def sanitize_list_values(items):
    """
    Sanitize all string values in a list by removing quotes.
    
    Parameters:
        items: List of items to sanitize
        
    Returns:
        List with sanitized string values
    """
    if not isinstance(items, list):
        return items
    
    return [sanitize_string_value(item) if isinstance(item, str) else item 
            for item in items]


def _load_from_individual_file(machine_name: str, machines_dir: str = DEFAULT_MACHINES_DIR) -> Optional[Machine]:
    """
    Load a machine from an individual JSON file in the machines directory.
    
    Parameters:
        machine_name (str): Name of the machine
        machines_dir (str): Directory containing individual machine files
        
    Returns:
        Machine: Machine object if found, None otherwise
    """
    machine_file = os.path.join(machines_dir, f"{machine_name}.json")
    if os.path.exists(machine_file):
        try:
            machine = Machine.from_file(machine_file)
            logger.info(f"Loaded machine '{machine_name}' from individual file: {machine_file}")
            return machine
        except Exception as e:
            logger.error(f"Failed to load machine from {machine_file}: {e}")
            return None
    return None


def sanitize_machine_config(config):
    """
    Sanitize a machine configuration dictionary by removing unwanted quotes
    from string values and list items.
    
    This fixes issues where JSON values like '"export VAR=1"' contain embedded
    quotes that end up in generated job scripts.
    
    Parameters:
        config (dict): Machine configuration dictionary
        
    Returns:
        dict: Sanitized configuration dictionary
    """
    if not isinstance(config, dict):
        return config
    
    # Fields that should be sanitized
    string_fields = ['launcher', 'workdir', 'host', 'username', 'env_setup']
    list_fields = ['prepend', 'postpend', 'modules']
    
    sanitized = config.copy()
    
    # Sanitize string fields
    for field in string_fields:
        if field in sanitized:
            sanitized[field] = sanitize_string_value(sanitized[field])
    
    # Sanitize list fields
    for field in list_fields:
        if field in sanitized:
            sanitized[field] = sanitize_list_values(sanitized[field])
    
    return sanitized


def _load_from_machines_json(machine_name: str, config_path: str = DEFAULT_CONFIG_PATH) -> Optional[Machine]:
    """
    Load a machine from the traditional machines.json file.
    
    Parameters:
        machine_name (str): Name of the machine
        config_path (str): Path to machines.json
        
    Returns:
        Machine: Machine object if found, None otherwise
    """
    if not os.path.exists(config_path):
        return None
    
    try:
        with open(config_path) as f:
            config = json.load(f)
        
        machines = config.get("machines", {})
        if machine_name in machines:
            machine_config = machines[machine_name]
            # Sanitize the configuration to remove unwanted quotes
            machine_config = sanitize_machine_config(machine_config)
            machine = Machine.from_dict(machine_name, machine_config)
            logger.info(f"Loaded machine '{machine_name}' from {config_path}")
            return machine
    except Exception as e:
        logger.error(f"Failed to load machine from {config_path}: {e}")
    
    return None


def _get_default_machine_name(config_path: str = DEFAULT_CONFIG_PATH, machines_dir: str = DEFAULT_MACHINES_DIR) -> Optional[str]:
    """
    Get the default machine name from configuration.
    
    Checks:
    1. machines.json for a "default" key
    2. default.json file in machines directory
    
    Parameters:
        config_path (str): Path to machines.json
        machines_dir (str): Directory containing individual machine files
        
    Returns:
        str: Default machine name, or None if not found
    """
    # Check machines.json for default key
    if os.path.exists(config_path):
        try:
            with open(config_path) as f:
                config = json.load(f)
            default = config.get("default")
            if default:
                logger.info(f"Found default machine in {config_path}: {default}")
                return default
        except Exception as e:
            logger.debug(f"Could not read default from {config_path}: {e}")
    
    # Check for default.json in machines directory
    default_file = os.path.join(machines_dir, "default.json")
    if os.path.exists(default_file):
        try:
            with open(default_file) as f:
                config = json.load(f)
            default = config.get("default")
            if default:
                logger.info(f"Found default machine in {default_file}: {default}")
                return default
        except Exception as e:
            logger.debug(f"Could not read default from {default_file}: {e}")
    
    return None


def load_machine(
    config_path: str = DEFAULT_CONFIG_PATH, 
    machine_name: str = DEFAULT_MACHINE_NAME,
    machines_dir: str = DEFAULT_MACHINES_DIR,
    return_object: bool = False
) -> Union[dict, Machine, None]:
    """
    Loads and parses a machine configuration into a queue dictionary or Machine object.
    
    Supports two configuration formats:
    1. Traditional: All machines in a single machines.json file
    2. Modular: Individual JSON files in machines/ directory with optional default
    
    The loader tries to find the machine in this order:
    1. Individual file in machines/ directory (e.g., ~/.xespresso/machines/cluster1.json)
    2. Entry in machines.json file (e.g., ~/.xespresso/machines.json)
    
    If machine_name is DEFAULT_MACHINE_NAME, also checks for a configured default machine.
    If the machine name is invalid, suggests available options interactively.

    Parameters:
        config_path (str): Path to the JSON config file
        machine_name (str): Name of the machine to load
        machines_dir (str): Directory containing individual machine files
        return_object (bool): If True, returns Machine object; if False, returns queue dict

    Returns:
        dict | Machine | None: Parsed configuration for scheduler and remote execution,
                                or Machine object if return_object=True,
                                or None if config not found
    """
    # If requesting default machine, check if a different default is configured
    if machine_name == DEFAULT_MACHINE_NAME:
        configured_default = _get_default_machine_name(config_path, machines_dir)
        if configured_default:
            machine_name = configured_default
            logger.info(f"Using configured default machine: {machine_name}")
    
    # Try loading from individual file first
    machine = _load_from_individual_file(machine_name, machines_dir)
    
    # Fall back to machines.json
    if machine is None:
        machine = _load_from_machines_json(machine_name, config_path)
    
    # If still not found, offer alternatives
    if machine is None:
        available_machines = list_machines(config_path, machines_dir)
        
        if not available_machines:
            warnings.warn(
                f"No machine configurations found.\n"
                f"Checked:\n"
                f" - {config_path}\n"
                f" - {machines_dir}\n"
                f"To create one, run:\n"
                f" from xespresso.machines.config.creator import create_machine\n"
                f" create_machine(path='{config_path}')\n"
                f"Returning None."
            )
            logger.error(f"No machine configurations found")
            return None
        
        logger.warning(f"Machine '{machine_name}' not found in config.")
        print(f"❌ Machine '{machine_name}' not found.")
        print("🧭 Available machines:")
        for name in available_machines:
            print(f" - {name}")
        retry = input("Enter a valid machine name or press Enter to cancel: ").strip()
        if retry and retry in available_machines:
            machine_name = retry
            logger.info(f"Retrying with machine: {machine_name}")
            # Recursive call with user's choice
            return load_machine(config_path, machine_name, machines_dir, return_object)
        else:
            print("❌ Aborted. No machine loaded.")
            logger.info("User aborted machine selection.")
            return None
    
    # Return Machine object or queue dict based on preference
    if return_object:
        return machine
    else:
        return machine.to_queue()

def save_machine(
    machine: Machine,
    config_path: str = DEFAULT_CONFIG_PATH,
    machines_dir: str = DEFAULT_MACHINES_DIR,
    use_individual_file: bool = True
) -> str:
    """
    Save a Machine configuration to disk.
    
    Supports two save formats:
    1. Individual file: Save to ~/.xespresso/machines/{machine_name}.json (default)
    2. Machines.json: Add to ~/.xespresso/machines.json
    
    Parameters:
        machine (Machine): Machine object to save
        config_path (str): Path to machines.json file
        machines_dir (str): Directory for individual machine files
        use_individual_file (bool): If True, save as individual file (recommended)
        
    Returns:
        str: Path to the saved file
        
    Raises:
        ValueError: If machine is not a Machine object
        IOError: If unable to save the file
    """
    if not isinstance(machine, Machine):
        raise ValueError("Expected Machine object")
    
    if use_individual_file:
        # Save as individual file (recommended)
        os.makedirs(machines_dir, exist_ok=True)
        filepath = os.path.join(machines_dir, f"{machine.name}.json")
        
        try:
            machine.to_file(filepath)
            logger.info(f"Machine '{machine.name}' saved to individual file: {filepath}")
            return filepath
        except Exception as e:
            logger.error(f"Failed to save machine to {filepath}: {e}")
            raise IOError(f"Failed to save machine: {e}")
    else:
        # Save to machines.json (traditional format)
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        
        # Load existing config or create new one
        config = {"machines": {}}
        if os.path.exists(config_path):
            try:
                with open(config_path) as f:
                    config = json.load(f)
                if "machines" not in config:
                    config["machines"] = {}
            except Exception as e:
                logger.warning(f"Could not load existing config from {config_path}: {e}")
                config = {"machines": {}}
        
        # Add/update machine in config
        config["machines"][machine.name] = machine.to_dict()
        
        try:
            with open(config_path, "w") as f:
                json.dump(config, f, indent=2)
            logger.info(f"Machine '{machine.name}' saved to {config_path}")
            return config_path
        except Exception as e:
            logger.error(f"Failed to save machine to {config_path}: {e}")
            raise IOError(f"Failed to save machine: {e}")


def list_machines(config_path: str = DEFAULT_CONFIG_PATH, machines_dir: str = DEFAULT_MACHINES_DIR) -> list:
    """
    Returns a list of machine names defined in config files.
    
    Scans both:
    1. Individual machine files in machines/ directory
    2. machines.json file

    Parameters:
        config_path (str): Path to the JSON config file
        machines_dir (str): Directory containing individual machine files

    Returns:
        list[str]: List of machine names
    """
    machines = set()
    
    # Check individual machine files
    if os.path.exists(machines_dir):
        try:
            for filename in os.listdir(machines_dir):
                if filename.endswith('.json') and filename != 'default.json':
                    machine_name = filename[:-5]  # Remove .json extension
                    machines.add(machine_name)
                    logger.debug(f"Found machine file: {filename}")
        except Exception as e:
            logger.error(f"Failed to list machines from {machines_dir}: {e}")
    
    # Check machines.json
    if os.path.exists(config_path):
        try:
            with open(config_path) as f:
                config = json.load(f)
            machines.update(config.get("machines", {}).keys())
            logger.debug(f"Found machines in {config_path}: {list(config.get('machines', {}).keys())}")
        except Exception as e:
            logger.error(f"Failed to read machine list from {config_path}: {e}")
    
    return sorted(list(machines))

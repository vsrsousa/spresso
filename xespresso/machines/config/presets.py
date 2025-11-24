"""
presets.py

Utility module for managing machine configuration presets in xespresso.

Presets are reusable templates stored as JSON files in the 'templates' directory.
They define default machine settings such as execution mode, scheduler, resources,
launcher command, and module environment. Presets can be used to pre-fill values
during machine creation or exported from existing machine configurations.

Functions provided:
- list_presets(): Lists available preset names (without .json extension)
- load_preset(name): Loads a preset and returns its configuration as a dict
- preset_exists(name): Checks if a preset exists by name
- create_preset_from_machine(machine, name): Saves a machine config as a new preset

This module is intended to be used by creator.py and other tools that need
standardized machine configurations. It does not validate presets internally;
validation should be handled by the consumer module.

Logging is applied for traceability. Warnings are not used here, as this module
does not interact directly with the user.
"""

import os
import json
try:
    from xespresso.utils import warnings as warnings
    from xespresso.utils.logging import get_logger
    logger = get_logger()
    warnings.apply_custom_format()
except ImportError:
    # Fallback to standard Python warnings and logging
    import warnings
    import logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "..", "templates")

def list_presets():
    """Returns a list of available preset names (without .json extension)."""
    if not os.path.isdir(TEMPLATE_DIR):
        logger.warning(f"Preset directory not found: {TEMPLATE_DIR}")
        return []
    return [f.replace(".json", "") for f in os.listdir(TEMPLATE_DIR) if f.endswith(".json")]

def load_preset(name: str) -> dict:
    """Loads a preset by name and returns its dictionary content."""
    path = os.path.join(TEMPLATE_DIR, name + ".json")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Preset '{name}' not found.")
    with open(path) as f:
        logger.info(f"Loaded preset: {name}")
        return json.load(f)

def preset_exists(name: str) -> bool:
    """Checks if a preset exists."""
    return os.path.isfile(os.path.join(TEMPLATE_DIR, name + ".json"))

def create_preset_from_machine(machine: dict, name: str):
    """Saves a machine configuration as a new preset."""
    os.makedirs(TEMPLATE_DIR, exist_ok=True)
    path = os.path.join(TEMPLATE_DIR, name + ".json")
    if os.path.exists(path):
        warnings.warn(f"Preset '{name}' already exists and will be overwritten.")
    with open(path, "w") as f:
        json.dump(machine, f, indent=2)
    logger.info(f"Preset '{name}' created at {path}")

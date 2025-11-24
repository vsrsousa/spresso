"""
migrate.py

Utility for migrating machine configurations from machines.json to individual JSON files.

This module helps transition from the traditional single-file format (machines.json)
to the newer modular format where each machine has its own JSON file.

Usage:
    from xespresso.machines.config.migrate import migrate_machines
    
    # Migrate all machines
    migrate_machines()
    
    # Migrate specific machines
    migrate_machines(machine_names=["cluster1", "gpu_node"])
    
    # Specify custom paths
    migrate_machines(
        machines_json_path="~/.xespresso/machines.json",
        output_dir="~/.xespresso/machines"
    )
"""

import os
import json
from typing import List, Optional
try:
    from xespresso.utils.logging import get_logger
    logger = get_logger()
except ImportError:
    import logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)


DEFAULT_MACHINES_JSON = os.path.expanduser("~/.xespresso/machines.json")
DEFAULT_OUTPUT_DIR = os.path.expanduser("~/.xespresso/machines")


def migrate_machines(
    machines_json_path: str = DEFAULT_MACHINES_JSON,
    output_dir: str = DEFAULT_OUTPUT_DIR,
    machine_names: Optional[List[str]] = None,
    overwrite: bool = False,
    preserve_original: bool = True
) -> dict:
    """
    Migrate machine configurations from machines.json to individual JSON files.
    
    Each machine in the machines.json file will be saved to a separate file
    in the output directory (e.g., machine_name.json).
    
    Parameters:
        machines_json_path (str): Path to the machines.json file
        output_dir (str): Directory to save individual machine files
        machine_names (List[str], optional): List of specific machines to migrate.
                                             If None, migrates all machines.
        overwrite (bool): Whether to overwrite existing individual machine files
        preserve_original (bool): Whether to keep the original machines.json file
    
    Returns:
        dict: Migration results with keys:
            - success (bool): Overall migration success
            - migrated (List[str]): List of successfully migrated machines
            - skipped (List[str]): List of skipped machines (already exist)
            - failed (List[str]): List of machines that failed to migrate
            - errors (dict): Errors encountered during migration
            - default (str, optional): Default machine if specified in config
    """
    results = {
        "success": False,
        "migrated": [],
        "skipped": [],
        "failed": [],
        "errors": {},
        "default": None
    }
    
    # Expand paths
    machines_json_path = os.path.expanduser(machines_json_path)
    output_dir = os.path.expanduser(output_dir)
    
    # Check if machines.json exists
    if not os.path.exists(machines_json_path):
        error_msg = f"machines.json not found at {machines_json_path}"
        logger.error(error_msg)
        print(f"❌ {error_msg}")
        results["errors"]["file_not_found"] = error_msg
        return results
    
    # Load machines.json
    try:
        with open(machines_json_path, 'r') as f:
            config = json.load(f)
        logger.info(f"Loaded machines.json from {machines_json_path}")
    except Exception as e:
        error_msg = f"Failed to load machines.json: {e}"
        logger.error(error_msg)
        print(f"❌ {error_msg}")
        results["errors"]["load_error"] = str(e)
        return results
    
    # Extract machines and default
    machines = config.get("machines", {})
    default_machine = config.get("default")
    
    if default_machine:
        results["default"] = default_machine
        logger.info(f"Found default machine: {default_machine}")
    
    if not machines:
        error_msg = "No machines found in machines.json"
        logger.warning(error_msg)
        print(f"⚠️ {error_msg}")
        results["errors"]["no_machines"] = error_msg
        return results
    
    # Filter machines if specific names provided
    if machine_names:
        machines_to_migrate = {name: machines[name] for name in machine_names if name in machines}
        missing = [name for name in machine_names if name not in machines]
        if missing:
            logger.warning(f"Machines not found in config: {missing}")
            print(f"⚠️ Machines not found: {', '.join(missing)}")
    else:
        machines_to_migrate = machines
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Migrate each machine
    print(f"\n🚀 Migrating {len(machines_to_migrate)} machine(s) from {machines_json_path}")
    print(f"   Output directory: {output_dir}\n")
    
    for machine_name, machine_config in machines_to_migrate.items():
        output_path = os.path.join(output_dir, f"{machine_name}.json")
        
        # Check if file already exists
        if os.path.exists(output_path) and not overwrite:
            logger.info(f"Skipping {machine_name}: file already exists")
            print(f"⏭️  Skipped '{machine_name}' (file exists, use overwrite=True to replace)")
            results["skipped"].append(machine_name)
            continue
        
        # Save machine to individual file
        try:
            with open(output_path, 'w') as f:
                json.dump(machine_config, f, indent=2)
            logger.info(f"Migrated {machine_name} to {output_path}")
            print(f"✅ Migrated '{machine_name}' → {output_path}")
            results["migrated"].append(machine_name)
        except Exception as e:
            error_msg = f"Failed to save {machine_name}: {e}"
            logger.error(error_msg)
            print(f"❌ Failed to migrate '{machine_name}': {e}")
            results["failed"].append(machine_name)
            results["errors"][machine_name] = str(e)
    
    # Create default.json if default machine is specified
    if default_machine:
        default_path = os.path.join(output_dir, "default.json")
        try:
            with open(default_path, 'w') as f:
                json.dump({"default": default_machine}, f, indent=2)
            logger.info(f"Created default.json pointing to {default_machine}")
            print(f"✅ Created default.json (points to '{default_machine}')")
        except Exception as e:
            logger.error(f"Failed to create default.json: {e}")
            print(f"⚠️ Failed to create default.json: {e}")
    
    # Print summary
    print(f"\n📊 Migration Summary:")
    print(f"   Migrated: {len(results['migrated'])} machine(s)")
    print(f"   Skipped:  {len(results['skipped'])} machine(s)")
    print(f"   Failed:   {len(results['failed'])} machine(s)")
    
    if results["migrated"] or (results["skipped"] and not results["failed"]):
        results["success"] = True
        
        # Offer to backup or remove original
        if preserve_original and results["migrated"]:
            print(f"\n💡 Original file preserved at: {machines_json_path}")
            print(f"   You can safely delete it after verifying the migration.")
    
    return results


def rollback_migration(
    output_dir: str = DEFAULT_OUTPUT_DIR,
    machines_json_path: str = DEFAULT_MACHINES_JSON,
    machine_names: Optional[List[str]] = None
) -> dict:
    """
    Rollback migration by deleting individual machine files.
    
    WARNING: This will delete individual machine JSON files. Use with caution.
    
    Parameters:
        output_dir (str): Directory containing individual machine files
        machines_json_path (str): Path to machines.json (for reference)
        machine_names (List[str], optional): Specific machines to rollback.
                                             If None, removes all.
    
    Returns:
        dict: Rollback results
    """
    results = {
        "success": False,
        "removed": [],
        "failed": [],
        "errors": {}
    }
    
    output_dir = os.path.expanduser(output_dir)
    
    if not os.path.exists(output_dir):
        logger.warning(f"Output directory not found: {output_dir}")
        print(f"⚠️ Directory not found: {output_dir}")
        return results
    
    # Get list of files to remove
    if machine_names:
        files_to_remove = [f"{name}.json" for name in machine_names]
    else:
        files_to_remove = [f for f in os.listdir(output_dir) if f.endswith('.json')]
    
    print(f"\n🔄 Rolling back {len(files_to_remove)} machine file(s)\n")
    
    for filename in files_to_remove:
        filepath = os.path.join(output_dir, filename)
        
        if not os.path.exists(filepath):
            logger.warning(f"File not found: {filepath}")
            print(f"⏭️  Skipped '{filename}' (not found)")
            continue
        
        try:
            os.remove(filepath)
            machine_name = filename[:-5]  # Remove .json
            logger.info(f"Removed {filepath}")
            print(f"✅ Removed '{filename}'")
            results["removed"].append(machine_name)
        except Exception as e:
            logger.error(f"Failed to remove {filepath}: {e}")
            print(f"❌ Failed to remove '{filename}': {e}")
            results["failed"].append(filename)
            results["errors"][filename] = str(e)
    
    print(f"\n📊 Rollback Summary:")
    print(f"   Removed: {len(results['removed'])} file(s)")
    print(f"   Failed:  {len(results['failed'])} file(s)")
    
    if results["removed"]:
        results["success"] = True
    
    return results

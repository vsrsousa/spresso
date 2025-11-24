#!/usr/bin/env python3
"""
Example: Migrating machines.json to individual machine files

This script demonstrates how to use the migrate_machines function to convert
a traditional machines.json file into individual JSON files for each machine.

Usage:
    python migrate_example.py
"""

from xespresso.machines import migrate_machines

# Example 1: Migrate all machines with default paths
print("=" * 60)
print("Example 1: Migrate all machines from ~/.xespresso/machines.json")
print("=" * 60)

# This will:
# - Read all machines from ~/.xespresso/machines.json
# - Create individual files in ~/.xespresso/machines/
# - Preserve the original machines.json file
result = migrate_machines()

if result["success"]:
    print(f"\n✅ Successfully migrated {len(result['migrated'])} machines!")
    for machine in result["migrated"]:
        print(f"   - {machine}")
else:
    print("\n❌ Migration failed or no machines found")

# Example 2: Migrate specific machines
print("\n" + "=" * 60)
print("Example 2: Migrate only specific machines")
print("=" * 60)

# This will migrate only the specified machines
result = migrate_machines(
    machine_names=["local_desktop", "slurm_cluster"],
    overwrite=False  # Don't overwrite existing files
)

# Example 3: Migrate with custom paths
print("\n" + "=" * 60)
print("Example 3: Migrate with custom paths")
print("=" * 60)

# This allows you to specify custom input/output paths
result = migrate_machines(
    machines_json_path="./examples/machines.json",
    output_dir="./examples/machines_migrated",
    overwrite=True  # Overwrite existing files
)

if result["success"]:
    print(f"\n✅ Migration successful!")
    print(f"   Migrated: {result['migrated']}")
    if result["default"]:
        print(f"   Default machine: {result['default']}")

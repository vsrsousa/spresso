# Machine Configuration Migration Guide

This guide explains how to migrate from the traditional `machines.json` format to individual machine JSON files.

## Overview

XEspresso now supports two formats for machine configurations:

1. **Traditional format**: All machines in a single `machines.json` file
2. **Modular format**: Individual JSON files for each machine in a `machines/` directory

The modular format is recommended as it:
- Makes version control easier (smaller diffs)
- Allows sharing individual machine configs
- Simplifies management of many machines
- Enables better organization

## Migration Process

### Automatic Migration

Use the `migrate_machines()` function to automatically convert your existing `machines.json`:

```python
from xespresso.machines import migrate_machines

# Migrate all machines with default paths
result = migrate_machines()

# Migrate with custom paths
result = migrate_machines(
    machines_json_path="~/.xespresso/machines.json",
    output_dir="~/.xespresso/machines"
)

# Migrate only specific machines
result = migrate_machines(
    machine_names=["cluster1", "cluster2"]
)

# Overwrite existing files
result = migrate_machines(overwrite=True)
```

### What Gets Migrated

1. **Individual machine files**: Each machine becomes `<machine_name>.json`
2. **Default machine**: If specified in `machines.json`, a `default.json` is created
3. **Original file**: Preserved by default (can be deleted after verification)

### Example

**Before** (`~/.xespresso/machines.json`):
```json
{
  "default": "local_desktop",
  "machines": {
    "local_desktop": {
      "execution": "local",
      "scheduler": "direct",
      "workdir": "./calculations",
      "nprocs": 4
    },
    "cluster1": {
      "execution": "remote",
      "scheduler": "slurm",
      "workdir": "/home/user/calc",
      "host": "cluster.edu",
      "username": "user"
    }
  }
}
```

**After** (`~/.xespresso/machines/`):
- `local_desktop.json`
- `cluster1.json`
- `default.json` (contains `{"default": "local_desktop"}`)

## Creating New Machines

When creating new machines with `create_machine()`, you now have the option to choose the save format:

```python
from xespresso.machines import create_machine

create_machine()  # Interactive prompt will ask for save format
```

During the interactive setup, you'll see:

```
💾 Choose how to save the machine configuration:
 [1] Add to machines.json (traditional, all machines in one file)
 [2] Save as individual JSON file (recommended, one file per machine)
Choose save format [1/2] [2]:
```

## Loading Machines

The `load_machine()` function automatically handles both formats:

```python
from xespresso.machines import load_machine

# Load from either format (automatically detected)
queue = load_machine(machine_name="cluster1")

# Load as Machine object
machine = load_machine(machine_name="cluster1", return_object=True)
```

The loader checks in this order:
1. Individual file: `~/.xespresso/machines/<machine_name>.json`
2. Entry in: `~/.xespresso/machines.json`

## Rollback

If you need to revert the migration:

```python
from xespresso.machines import rollback_migration

# Remove all individual machine files
result = rollback_migration()

# Remove specific machines
result = rollback_migration(machine_names=["cluster1", "cluster2"])
```

**Note**: This only removes the individual files. Your original `machines.json` remains intact.

## Best Practices

1. **Migrate existing configurations**: If you have a `machines.json`, migrate it to individual files
2. **Use version control**: Add `~/.xespresso/machines/*.json` to your dotfiles repository
3. **Share configs**: Individual files are easier to share with colleagues
4. **Set default**: Create a `default.json` to specify your default machine
5. **Keep organized**: Use descriptive machine names (e.g., `lab_cluster`, `home_workstation`)

## Command-line Usage

You can also use the migration function from the command line:

```bash
# Create a simple migration script
python -c "from xespresso.machines import migrate_machines; migrate_machines()"
```

Or use the example script:

```bash
python examples/migrate_example.py
```

## Troubleshooting

### Migration fails with "file not found"
- Check that `~/.xespresso/machines.json` exists
- Verify the path with: `ls -la ~/.xespresso/machines.json`

### Machines are skipped during migration
- Individual files already exist
- Use `overwrite=True` to replace them

### Can't load machine after migration
- Check that the file exists: `ls ~/.xespresso/machines/`
- Verify JSON syntax: `python -m json.tool ~/.xespresso/machines/machine_name.json`

## API Reference

### `migrate_machines()`

Migrate machines from `machines.json` to individual files.

**Parameters:**
- `machines_json_path` (str): Path to machines.json (default: `~/.xespresso/machines.json`)
- `output_dir` (str): Directory for individual files (default: `~/.xespresso/machines`)
- `machine_names` (List[str], optional): Specific machines to migrate (default: all)
- `overwrite` (bool): Overwrite existing files (default: False)
- `preserve_original` (bool): Keep machines.json (default: True)

**Returns:**
- dict: Migration results with keys `success`, `migrated`, `skipped`, `failed`, `errors`, `default`

### `rollback_migration()`

Remove individual machine files created by migration.

**Parameters:**
- `output_dir` (str): Directory containing machine files (default: `~/.xespresso/machines`)
- `machines_json_path` (str): Path to machines.json for reference (default: `~/.xespresso/machines.json`)
- `machine_names` (List[str], optional): Specific machines to remove (default: all)

**Returns:**
- dict: Rollback results with keys `success`, `removed`, `failed`, `errors`

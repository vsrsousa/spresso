# Quick Reference: Machine Configuration Migration

## Quick Start

### Migrate All Machines
```python
from xespresso.machines import migrate_machines

result = migrate_machines()
```

### Migrate with Command Line
```bash
python scripts/migrate_machines_cli.py
```

### Create New Machine with Format Choice
```python
from xespresso.machines import create_machine

create_machine()  # Interactive prompt will ask for save format
```

## Common Commands

### Migrate specific machines
```python
migrate_machines(machine_names=["cluster1", "cluster2"])
```

### Migrate with custom paths
```python
migrate_machines(
    machines_json_path="./my_machines.json",
    output_dir="./machines/"
)
```

### Overwrite existing files
```python
migrate_machines(overwrite=True)
```

### Rollback migration
```python
from xespresso.machines import rollback_migration

rollback_migration()  # Remove all individual machine files
```

## File Locations

- **Traditional format**: `~/.xespresso/machines.json`
- **Individual files**: `~/.xespresso/machines/<machine_name>.json`
- **Default machine**: `~/.xespresso/machines/default.json`

## Import Statements

```python
# Main functions
from xespresso.machines import (
    Machine,              # Machine class
    load_machine,         # Load a machine configuration
    create_machine,       # Create new machine (interactive)
    migrate_machines,     # Migrate machines.json to individual files
    rollback_migration,   # Rollback migration
    list_machines        # List all available machines
)

# Direct imports
from xespresso.machines.config.migrate import migrate_machines, rollback_migration
from xespresso.machines.config.loader import load_machine, list_machines
from xespresso.machines.config.creator import create_machine
from xespresso.machines.machine import Machine
```

## Machine Properties

```python
machine = Machine(name='test', scheduler='slurm', ...)

# Access properties
machine.name          # Machine name
machine.scheduler     # Scheduler type ('direct', 'slurm', etc.)
machine.execution     # Execution mode ('local' or 'remote')
machine.workdir       # Working directory
machine.nprocs        # Number of processes
machine.is_remote     # Boolean: is this remote?
machine.is_local      # Boolean: is this local?

# Serialization
machine.to_dict()     # Convert to dictionary
machine.to_file(path) # Save to JSON file
machine.to_queue()    # Convert to queue format (for schedulers)

# Load from file
machine = Machine.from_file('path/to/machine.json')
```

## Migration Result

```python
result = migrate_machines()

# Check result
result['success']      # bool: Overall success
result['migrated']     # list: Successfully migrated machines
result['skipped']      # list: Skipped machines (already exist)
result['failed']       # list: Failed migrations
result['errors']       # dict: Error details
result['default']      # str: Default machine name (if set)
```

## Troubleshooting

### Problem: "Machine config file not found"
**Solution**: Create machines.json first or use `create_machine()`

### Problem: "Machines are skipped during migration"
**Solution**: Use `overwrite=True` to replace existing files

### Problem: "Can't import migrate_machines"
**Solution**: 
```python
# Use full import path
from xespresso.machines.config.migrate import migrate_machines
# Or use package import
from xespresso.machines import migrate_machines
```

### Problem: "Migration successful but can't load machine"
**Solution**: Check file exists and has valid JSON
```bash
ls ~/.xespresso/machines/
python -m json.tool ~/.xespresso/machines/machine_name.json
```

## Examples in Repository

- `examples/migrate_example.py` - Python migration examples
- `scripts/migrate_machines_cli.py` - Command-line migration tool
- `docs/MACHINE_MIGRATION.md` - Complete migration guide
- `docs/MACHINE_CONFIGURATION.md` - Configuration reference

## Testing

Run tests to verify everything works:
```bash
python -m pytest tests/test_machine.py tests/test_migrate.py -v
```

## Support

- See `docs/MACHINE_MIGRATION.md` for detailed guide
- See `IMPLEMENTATION_DETAILS.md` for technical details
- Check `tests/test_migrate.py` for usage examples

# Remote Job Execution Support in Simplified Workflow

## Summary

This implementation adds comprehensive remote job execution support to the xespresso simplified workflow system. Users can now easily run Quantum ESPRESSO calculations on remote HPC clusters using the same simple API as local execution.

## Problem Statement

The user asked: "I forgot to ask you how the workflow you implemented is handling remote job executions"

The original workflow implementation (`CalculationWorkflow`, `quick_scf`, `quick_relax`) did not support remote execution. While the underlying scheduler system had full remote execution capabilities via `RemoteExecutionMixin`, the workflow API didn't expose this functionality.

## Solution

### 1. Added Remote Execution Parameters

**CalculationWorkflow class:**
- `queue` parameter: Accept direct queue configuration dictionary
- `machine` parameter: Load machine configuration from `~/.xespresso/machines/`
- Validation: Prevents using both parameters simultaneously

**Helper functions:**
- Updated `quick_scf()` and `quick_relax()` to accept `queue` and `machine` parameters

### 2. Integration with Existing Infrastructure

The implementation leverages existing remote execution infrastructure:
- **RemoteExecutionMixin**: Handles SSH connections and file transfers
- **Connection caching**: Automatic reuse of SSH connections
- **SLURM support**: Job submission and monitoring
- **File management**: Automatic transfer of inputs, pseudopotentials, and outputs

No changes were needed to the existing scheduler system.

### 3. Usage Examples

**Using machine configuration (recommended):**
```python
from xespresso import quick_scf

calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    machine='cluster1'  # Loads from ~/.xespresso/machines/cluster1.json
)
```

**Using direct queue configuration:**
```python
queue = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "cluster.edu",
    "remote_user": "username",
    "remote_dir": "/scratch/calculations",
    "remote_auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"},
    "nodes": 1,
    "ntasks-per-node": 16,
    "time": "02:00:00"
}

calc = quick_scf(
    'structure.cif',
    {'Fe': 'Fe.pbe-spn.UPF'},
    protocol='moderate',
    queue=queue
)
```

**Combining all features:**
```python
# Remote execution + protocol presets + magnetic config + Hubbard U + k-spacing
calc = quick_scf(
    atoms,
    {'Fe': 'Fe.pbe-spn.UPF', 'O': 'O.pbe.UPF'},
    protocol='accurate',
    magnetic_config={'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}},
    kspacing=0.15,
    machine='cluster1'
)
```

## Files Changed

1. **xespresso/workflow/simple_workflow.py**
   - Added `queue` and `machine` parameters to `CalculationWorkflow.__init__`
   - Added `queue` and `machine` parameters to `from_cif()` classmethod
   - Modified `run_scf()` and `run_relax()` to pass queue to calculator
   - Updated `quick_scf()` and `quick_relax()` helper functions
   - Added `load_machine` import at module level

2. **tests/test_workflow.py**
   - Added `test_workflow_with_queue()` - test queue configuration
   - Added `test_workflow_queue_and_machine_conflict()` - test validation
   - Added `test_quick_scf_with_queue()` - test quick function with queue
   - Added `test_quick_relax_with_queue()` - test quick function with queue
   - Added `test_workflow_none_queue()` - test backward compatibility

3. **WORKFLOW_DOCUMENTATION.md**
   - Added extensive "Remote Job Execution" section
   - Examples of machine and queue configurations
   - Integration examples with other features
   - Best practices and troubleshooting

4. **examples/workflow_remote_execution.py**
   - 6 complete examples demonstrating remote execution
   - Machine configuration examples
   - Direct queue configuration examples
   - Combination with other features (magnetic, Hubbard, k-spacing)

## Testing

**All tests pass (21/21):**
- Backward compatibility maintained (workflow works without queue/machine)
- Queue configuration works correctly
- Machine configuration works correctly
- Validation prevents using both queue and machine
- All existing workflow tests still pass
- All magnetic configuration tests still pass

**Security:**
- CodeQL security scan: No vulnerabilities found
- Code review completed and feedback addressed

## Key Features

1. **Backward Compatible**: Existing workflow code continues to work without changes
2. **Simple API**: Same easy-to-use interface for local and remote execution
3. **Connection Persistence**: SSH connections automatically cached and reused
4. **Automatic File Transfer**: Inputs, pseudopotentials, and outputs handled automatically
5. **SLURM Integration**: Job submission and monitoring built-in
6. **Works with All Features**: Protocol presets, k-spacing, magnetic configs, Hubbard parameters all work with remote execution
7. **Validation**: Clear error messages for configuration issues
8. **Documentation**: Comprehensive documentation and examples

## Answer to Problem Statement

**How does the workflow handle remote job executions?**

The workflow system now fully supports remote job execution through two approaches:

1. **Machine configurations** (recommended): Users create and save machine configurations in `~/.xespresso/machines/` and reference them by name in workflows
2. **Direct queue configuration**: Users can pass queue configuration dictionaries directly

The workflow leverages the existing `RemoteExecutionMixin` infrastructure which provides:
- Automatic SSH connection caching and reuse (no redundant connections)
- Automatic file transfer (inputs, pseudopotentials, outputs)
- SLURM job submission and monitoring
- Error handling and retries

Users can run calculations on remote clusters with the exact same simple API as local execution, just adding a `machine='cluster1'` or `queue={...}` parameter.

## Examples

See:
- `examples/workflow_remote_execution.py` - Complete examples
- `WORKFLOW_DOCUMENTATION.md` - Full documentation
- `docs/MACHINE_CONFIGURATION.md` - Machine configuration guide
- `docs/REMOTE_CONNECTION_PERSISTENCE.md` - Connection management details

# Machine/Code/Version Persistence and Adaptive Resources Implementation

## Overview

This document describes the implementation of cross-page persistence for machine, code, and version selections, and the addition of adaptive resources configuration to the Calculation Setup page.

## Problem Statement

Users reported two issues:

1. **Machine/Code/Version Not Persisting**: When users selected a machine, code, or version in Calculation Setup and then navigated to Workflow Builder (or vice versa), their selections were not preserved. Each page maintained its own separate state.

2. **Missing Adaptive Resources in Calculation Setup**: The Workflow Builder page had adaptive resources that showed only relevant options based on machine type (direct execution vs scheduler), but Calculation Setup did not have this feature.

## Solution

### 1. Shared Session State Variables

**Before:**
- Calculation Setup used: `selected_machine_for_calc`, `calc_selected_version`, `calc_selected_code`
- Workflow Builder used: `selected_machine_for_workflow`, `workflow_selected_version`, `workflow_selected_code`

**After:**
Both pages now use shared variables:
- `selected_machine` (shared)
- `selected_version` (shared)
- `selected_code` (shared)

This ensures that selections made in one page are immediately available in the other page.

### 2. Default Index Handling

Added logic to restore previously selected values when the page loads:

```python
# Determine default index for machine selector
default_machine_idx = 0
if (
    st.session_state.get("selected_machine")
    and st.session_state.selected_machine in available_machines
):
    default_machine_idx = available_machines.index(
        st.session_state.selected_machine
    )

selected_machine_name = st.selectbox(
    "Select Machine:",
    options=available_machines,
    index=default_machine_idx,  # Restores previous selection
    key="calc_machine_selector",
)
```

### 3. Adaptive Resources Configuration

The Resources Configuration section now adapts based on machine type:

#### Direct Execution Mode
For machines with `scheduler="direct"`:
- Shows only **Number of Processors (nprocs)**
- Shows the launcher command that will be used
- Automatically resolves `{nprocs}` placeholders in launcher templates
- Example: `mpirun -np {nprocs}` becomes `mpirun -np 8` when nprocs=8

#### Scheduler Mode  
For machines with scheduler (SLURM, PBS, SGE):
- Shows full resource configuration:
  - Nodes
  - Tasks per Node
  - Memory
  - Time Limit
  - Partition/Queue
  - Account (optional)

## User Experience Improvements

### Before
1. User selects Machine 3 in Calculation Setup
2. User navigates to Workflow Builder
3. Machine dropdown shows Machine 1 (default) ❌
4. User must re-select Machine 3

### After
1. User selects Machine 3 in Calculation Setup
2. User navigates to Workflow Builder  
3. Machine dropdown shows Machine 3 (preserved) ✓
4. No need to re-select

## Code Changes

### Files Modified

1. **xespresso/gui/pages/calculation_setup.py**
   - Updated to use shared session state variables
   - Added default index handling for machine selector
   - Implemented adaptive resources configuration
   - Added scheduler type detection
   - Added nprocs configuration for direct execution
   - Added launcher placeholder resolution

2. **xespresso/gui/pages/workflow_builder.py**
   - Updated to use shared session state variables
   - Added default index handling for machine selector
   - No changes to resources section (already had adaptive resources)

## Testing

Created comprehensive tests to verify:

1. **Shared State Variables**: Both pages use the same session state variables
2. **Default Index Handling**: Machine selectors properly restore selections
3. **Adaptive Resources**: Calculation Setup adapts to machine type
4. **Consistency**: Both pages have consistent behavior

All tests passed:
```
✓ Shared State Variables
✓ Default Index Handling  
✓ Adaptive Resources in Calculation Setup
✓ Consistency Between Pages
```

## Security

CodeQL security analysis found 0 alerts. The changes do not introduce any security vulnerabilities.

## Backward Compatibility

- Existing machine configurations work without modification
- Session state additions are backward compatible
- UI enhancements are additive, not replacing existing functionality
- No breaking changes to the API or file formats

## Benefits

1. **Seamless Navigation**: Users can switch between pages without losing their selections
2. **Better UX**: Relevant options shown based on machine configuration
3. **Consistency**: Both pages now have identical behavior
4. **Clarity**: Users see only options that apply to their machine type
5. **Less Confusion**: Direct execution users don't see scheduler-specific options

## Implementation Details

### Session State Flow

```
User in Calculation Setup:
1. Selects Machine 3
2. st.session_state.selected_machine = "Machine 3"

User navigates to Workflow Builder:
3. Page reads st.session_state.selected_machine
4. Finds "Machine 3" in available_machines list
5. Sets default_machine_idx to index of "Machine 3"
6. Selectbox shows "Machine 3" as selected ✓

User changes to Machine 5:
7. st.session_state.selected_machine = "Machine 5"

User navigates back to Calculation Setup:
8. Page reads st.session_state.selected_machine
9. Finds "Machine 5" in list
10. Selectbox shows "Machine 5" as selected ✓
```

### Adaptive Resources Logic

```python
# Get machine and determine scheduler type
machine = st.session_state.calc_machine
scheduler_type = getattr(machine, "scheduler", "direct")

if scheduler_type == "direct":
    # Show only nprocs
    nprocs = st.number_input("Number of Processors (nprocs):", ...)
    config["nprocs"] = nprocs
    
    # Resolve launcher template
    if "{nprocs}" in launcher:
        resolved = launcher.replace("{nprocs}", str(nprocs))
    
else:
    # Show full resources (nodes, tasks, memory, time, partition)
    nodes = st.number_input("Nodes:", ...)
    ntasks_per_node = st.number_input("Tasks per Node:", ...)
    # ... etc
```

## Future Enhancements

Possible future improvements:
- Add session persistence to disk (save/load sessions)
- Extend persistence to other configuration options
- Add validation for resource values
- Add resource presets for common job types

## References

- Original implementation: `IMPLEMENTATION_ADAPTIVE_RESOURCES_AND_PERSISTENCE.md`
- Session management: `xespresso/gui/utils/session_manager.py`
- Machine configuration: `xespresso/machines/config/loader.py`
- Codes configuration: `xespresso/codes/manager.py`

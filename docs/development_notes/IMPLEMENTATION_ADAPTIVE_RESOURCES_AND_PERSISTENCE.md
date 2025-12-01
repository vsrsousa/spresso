# Implementation Summary: Adaptive Resources and Persistent Structure Viewer

## Changes Made

This implementation addresses two main requirements from the problem statement:

### 1. Adaptive Resources Configuration

**Problem:** The resource adjustment should adapt to the type of machine. For direct execution, only the number of processor cores (nprocs) makes sense. Other parameters (nodes, ntasks-per-node, memory, time, partition) only make sense for a scheduler. Also, if the launcher is defined as `"mpirun -np {nprocs}"`, it should detect and use the nprocs value.

**Solution:**

#### File: `xespresso/gui/pages/workflow_builder.py`

Modified the resources configuration section (lines 397-488) to:

1. **Detect scheduler type** from the selected machine
2. **For direct scheduler (`scheduler="direct"`):**
   - Only show the `nprocs` (Number of Processors) input field
   - Display an info message explaining that scheduler resources don't apply to direct execution
   - Detect if launcher contains `{nprocs}` placeholder and show resolved command
   - Example: `mpirun -np {nprocs}` becomes `mpirun -np 8` when nprocs=8

3. **For schedulers (slurm, pbs, sge):**
   - Show full resource configuration:
     - Nodes
     - Tasks per Node (ntasks-per-node)
     - Memory (mem)
     - Time Limit
     - Partition/Queue
     - Account (optional)
   - Display info message explaining scheduler mode

4. **When resources are not adjusted:**
   - For direct execution: Show nprocs and launcher
   - For schedulers: Show full resource configuration from machine defaults

### 2. Persistent Structure Viewer

**Problem:** The structure viewer should have a persistent structure that was defined when the user chose the structure. The structure should persist across page navigation, and there should be a dedicated section to view the currently selected structure.

**Solution:**

#### File: `xespresso/gui/pages/structure_viewer.py`

1. **Added persistent structure display** (lines 18-47):
   - Shows "Currently Selected Structure" section at the top of the page
   - Displays key metrics: Formula, Number of Atoms, Source
   - Includes expandable "View Current Structure" section with full visualization
   - Only appears when a structure is loaded

2. **Enhanced structure loading** with source tracking:
   - `render_upload_tab()`: Stores source as `"Upload: {filename}"`
   - `render_browse_tab()`: Stores source as `"File: {basename}"`
   - `render_build_structure_tab()`: Stores source as `"Built: {element} {crystal_type}"` or `"Built: {molecule_name} molecule"`
   - `render_ase_database_tab()`: Stores source as `"Database: ID {id}"`

3. **Session state management:**
   - `current_structure`: Stores the ASE Atoms object
   - `structure_source`: Stores the source information for display
   - These persist across page navigation in Streamlit's session state

## Testing

Created comprehensive test scripts to verify functionality:

### Test 1: Adaptive Resources (`/tmp/test_adaptive_resources.py`)
- ✓ Direct scheduler correctly configured with only nprocs
- ✓ SLURM scheduler correctly configured with full resources
- ✓ PBS scheduler correctly configured with full resources
- ✓ Launcher with `{nprocs}` placeholder detected and resolved correctly
- ✓ `to_queue()` serialization works for both scheduler types

### Test 2: Structure Persistence (`/tmp/test_structure_persistence.py`)
- ✓ Structure loading from build (bulk/molecule) works correctly
- ✓ Structure loading from upload (simulated) works correctly
- ✓ Structure loading from database (simulated) works correctly
- ✓ Structure source information is stored and persisted
- ✓ Structure persists across page navigation (simulated)

### Test 3: Existing Tests
- ✓ 77 GUI tests run
- ✓ 61 tests passed
- ✓ Structure persistence test specifically passed
- Failed tests are unrelated to our changes (missing pseudopotentials, environment issues)

## User Benefits

### For Resources Configuration:
1. **Clarity**: Users now see only relevant options based on their machine type
2. **Simplicity**: Direct execution users don't get confused by scheduler-specific options
3. **Transparency**: The resolved launcher command is shown, making it clear what will be executed
4. **Consistency**: Machine configuration page already follows this pattern; now workflow builder does too

### For Structure Viewer:
1. **Persistence**: Structure remains accessible across all pages
2. **Visibility**: Users can always see what structure is currently selected
3. **Context**: Source information helps users remember where the structure came from
4. **Quick Access**: Expandable viewer allows quick checks without leaving the page
5. **Confidence**: Clear indication that structure is "ready for calculations"

## Backward Compatibility

All changes are backward compatible:
- Existing machine configurations work without modification
- Session state additions don't break existing workflows
- The UI enhancements are additive, not replacing existing functionality

## Files Modified

1. `xespresso/gui/pages/workflow_builder.py`
   - Lines 397-488: Adaptive resources configuration

2. `xespresso/gui/pages/structure_viewer.py`
   - Lines 11-47: Persistent structure display
   - Lines 45-76: Enhanced upload tab with source tracking
   - Lines 235-286: Enhanced browse tab with source tracking
   - Lines 329-346: Enhanced build crystal with source tracking
   - Lines 362-375: Enhanced build molecule with source tracking
   - Lines 499-510: Enhanced database load with source tracking

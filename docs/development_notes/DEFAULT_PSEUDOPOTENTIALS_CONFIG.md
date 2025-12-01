# Default Pseudopotentials Configuration

## Overview

The pseudopotentials module now supports setting a **default configuration** that will be automatically used in Calculation Setup and Workflow Builder pages. This makes it easier for users to have a preferred pseudopotential set that's always ready to use.

## How It Works

### Setting a Default

1. **Configure a pseudopotential library** (one-time):
   - Go to "🧪 Pseudopotentials Configuration" page
   - Auto-detect and save your preferred library (e.g., "SSSP_1.1.2_PBE_efficiency")

2. **Set it as default**:
   - In the existing configurations list, load your preferred configuration
   - Click the **⭐ Set as Default** button
   - A copy is created as `default.json` in `~/.xespresso/pseudopotentials/`

3. **Use in calculations**:
   - Navigate to Calculation Setup or Workflow Builder
   - The default configuration is automatically selected
   - It appears first in the list with a ⭐ indicator
   - You can still choose other configurations if needed

### Clearing a Default

- Click **🚫 Clear Default** button to remove the default setting
- The original configuration remains intact, only `default.json` is removed
- After clearing, users must manually select a configuration

## User Interface

### Pseudopotentials Configuration Page

When viewing a configuration:

```
┌────────────────────────────────────────────────────────────┐
│ Configuration Details                                      │
├────────────────────────────────────────────────────────────┤
│ Name: SSSP_1.1.2_PBE_efficiency                           │
│ Library: SSSP v1.1.2                                      │
│ Functional: PBE                                           │
│ Elements: 118                                             │
├────────────────────────────────────────────────────────────┤
│                                                            │
│ [⭐ Set as Default]    [🚫 Clear Default]                 │
│                                                            │
│ ───────────────────────────────────────────────────────── │
│                                                            │
│ [🗑️ Delete This Configuration]                            │
│                                                            │
└────────────────────────────────────────────────────────────┘
```

**Button Behavior:**
- **Set as Default**: Only shown for non-default configurations
- **Clear Default**: Only shown when a default exists (regardless of which config you're viewing)
- **Delete**: Disabled for default configuration (must clear default first)

### Calculation Setup / Workflow Builder

When selecting pseudopotentials:

```
┌────────────────────────────────────────────────────────────┐
│ 🧪 Pseudopotentials                                        │
├────────────────────────────────────────────────────────────┤
│                                                            │
│ ℹ️ Using default pseudopotentials configuration           │
│                                                            │
│ Filter by Library: [All ▼]  Filter by Functional: [All ▼]│
│                                                            │
│ Select Pseudopotential Configuration:                     │
│ ┌────────────────────────────────────────────────────────┐│
│ │ ⭐ Default (SSSP 1.1.2, PBE, 118 elements) ▼          ││
│ ├────────────────────────────────────────────────────────┤│
│ │ • ⭐ Default (SSSP 1.1.2, PBE, 118 elements)          ││
│ │ • SSSP_1.1.2_PBE_precision (SSSP 1.1.2, PBE, ...)    ││
│ │ • PSLibrary_PBE (PSLibrary 1.0.0, PBE, ...)          ││
│ └────────────────────────────────────────────────────────┘│
│                                                            │
│ ✅ All required elements available                         │
│                                                            │
│ Selected Pseudopotentials:                                │
│ │ Element │ File                                         ││
│ │ Fe      │ Fe.pbe-spn-kjpaw_psl.0.2.1.UPF              ││
│ │ Si      │ Si.pbe-n-rrkjus_psl.1.0.0.UPF               ││
│                                                            │
└────────────────────────────────────────────────────────────┘
```

## Technical Details

### File Structure

```
~/.xespresso/pseudopotentials/
├── SSSP_1.1.2_PBE_efficiency.json    # Original configuration
├── SSSP_1.1.2_PBE_precision.json     # Another configuration
├── PSLibrary_PBE.json                 # Yet another configuration
└── default.json                       # Copy of the default (e.g., copy of SSSP_1.1.2_PBE_efficiency.json)
```

### API Methods

**PseudopotentialsManager class methods:**

```python
# Set a configuration as default
PseudopotentialsManager.set_default_config(
    config_name="SSSP_1.1.2_PBE_efficiency",
    pseudopotentials_dir=DEFAULT_PSEUDOPOTENTIALS_DIR
)
# Creates ~/.xespresso/pseudopotentials/default.json

# Get the default configuration
default_config = PseudopotentialsManager.get_default_config(
    pseudopotentials_dir=DEFAULT_PSEUDOPOTENTIALS_DIR
)
# Returns PseudopotentialsConfig object or None

# Check if default exists
has_default = PseudopotentialsManager.has_default_config(
    pseudopotentials_dir=DEFAULT_PSEUDOPOTENTIALS_DIR
)
# Returns True or False

# Clear the default
PseudopotentialsManager.clear_default_config(
    pseudopotentials_dir=DEFAULT_PSEUDOPOTENTIALS_DIR
)
# Removes default.json
```

### Behavior in GUI Selector

The `render_pseudopotentials_selector()` component:

1. **On first load**: Checks for `default.json`
   - If exists: Automatically selects it and shows info message
   - If not exists: User must select from list

2. **In configuration list**:
   - Default appears first with ⭐ indicator
   - Named as "⭐ Default (library version, functional, elements)"
   - User can select other configurations if desired

3. **Persistence**:
   - Selected configuration stored in session state
   - Separate for Calculation Setup and Workflow Builder
   - Default is only auto-selected on first visit to page

## Use Cases

### Use Case 1: Single Preferred Library

**Scenario**: User always uses SSSP efficiency for PBE calculations

**Setup**:
1. Configure "SSSP_1.1.2_PBE_efficiency"
2. Set as default
3. Done!

**Usage**:
- Every time user opens Calculation Setup or Workflow Builder
- SSSP efficiency is automatically selected
- No need to choose from list each time

### Use Case 2: Multiple Libraries with Preferred Default

**Scenario**: User has SSSP efficiency (general), SSSP precision (accurate), and PSLibrary (special cases)

**Setup**:
1. Configure all three libraries
2. Set SSSP efficiency as default (most common)
3. When needed, user can switch to precision or PSLibrary

**Usage**:
- Most calculations use default (SSSP efficiency)
- For high-precision: User selects SSSP precision
- For special cases: User selects PSLibrary

### Use Case 3: Team Standard

**Scenario**: Research group has a standard pseudopotential set

**Setup**:
1. Administrator configures the standard set
2. Sets it as default
3. Shares `~/.xespresso/pseudopotentials/` with team

**Usage**:
- All team members use the same default
- Ensures consistency across calculations
- Individual users can add their own configs if needed

## Benefits

1. **Convenience**: No need to select pseudopotentials every time
2. **Consistency**: Same pseudopotentials used across calculations
3. **Flexibility**: Can still choose other configurations when needed
4. **Safety**: Original configurations are preserved
5. **Discovery**: New users see a sensible default immediately

## Example Workflow

```python
from xespresso.pseudopotentials import (
    create_pseudopotentials_config,
    PseudopotentialsManager,
    DEFAULT_PSEUDOPOTENTIALS_DIR
)

# Step 1: Create configuration
config = create_pseudopotentials_config(
    name="SSSP_1.1.2_PBE_efficiency",
    base_path="/home/user/pseudopotentials/SSSP_1.1.2_PBE_efficiency",
    library="SSSP",
    version="1.1.2",
    functional="PBE",
    save=True
)

# Step 2: Set as default
PseudopotentialsManager.set_default_config(
    "SSSP_1.1.2_PBE_efficiency",
    DEFAULT_PSEUDOPOTENTIALS_DIR
)

print("✅ Default configuration set!")

# Later, in calculation code:
# The GUI will automatically use this default
# Or in script:
default = PseudopotentialsManager.get_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
if default:
    pseudopotentials = default.get_pseudopotentials_dict()
    pseudo_dir = default.base_path
else:
    # Fallback to manual specification
    pseudopotentials = {"Fe": "Fe.UPF", "O": "O.UPF"}
    pseudo_dir = "/path/to/pseudos"
```

## Notes

- The `default.json` file is a **copy** of the original configuration
- Modifying the original configuration does NOT update the default
- To update default: Set it as default again (overwrites `default.json`)
- The default configuration cannot be deleted directly - must clear it first
- This prevents accidental removal of the default setting

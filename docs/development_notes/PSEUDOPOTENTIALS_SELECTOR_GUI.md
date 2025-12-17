# Pseudopotentials Selector - GUI Integration

## Overview
The Calculation Setup and Workflow Builder pages now include an improved pseudopotentials selector that allows users to:

1. **Select from configured libraries** - Choose from pre-configured pseudopotential sets
2. **Filter by library and functional** - Filter options by library name (SSSP, PSLibrary, etc.) and functional type (PBE, LDA, etc.)
3. **Auto-populate elements** - Automatically assigns pseudopotentials for all elements in the structure
4. **Persistent selection** - The selected configuration persists across the interface using session state

## New Interface Components

### 1. Pseudopotentials Selector Section

Located in both:
- **Calculation Setup page** (key_prefix: "calc")
- **Workflow Builder page** (key_prefix: "workflow")

### 2. Filter Options

```
┌─────────────────────────────────────────────────────┐
│ 🧪 Pseudopotentials                                 │
├─────────────────────────────────────────────────────┤
│ Filter by Library:    [All ▼]  [SSSP, PSLibrary]   │
│ Filter by Functional: [All ▼]  [PBE, LDA, PBEsol]  │
└─────────────────────────────────────────────────────┘
```

### 3. Configuration Selector

```
┌─────────────────────────────────────────────────────────────────┐
│ Select Pseudopotential Configuration:                           │
│ [SSSP_1.1.2_PBE_efficiency (SSSP 1.1.2, PBE, 118 elements) ▼]  │
└─────────────────────────────────────────────────────────────────┘
```

### 4. Configuration Details (Expandable)

```
┌─────────────────────────────────────────────────────┐
│ 📋 Configuration Details                    [−]     │
├─────────────────────────────────────────────────────┤
│ Name: SSSP_1.1.2_PBE_efficiency                     │
│ Library: SSSP                                       │
│ Version: 1.1.2                                      │
│ Functional: PBE                                     │
│ Base Path: /home/user/pseudos/SSSP_1.1.2_PBE...    │
│ Location: Local                                     │
│ Elements: 118                                       │
│ Description: SSSP 1.1.2 PBE efficiency set...      │
└─────────────────────────────────────────────────────┘
```

### 5. Selected Pseudopotentials Table

```
┌─────────────────────────────────────────────────────┐
│ ✅ All required elements available                  │
│                                                      │
│ Selected Pseudopotentials:                          │
│ ┌─────────┬───────────────────────────────────────┐ │
│ │ Element │ File                                  │ │
│ ├─────────┼───────────────────────────────────────┤ │
│ │ Fe      │ Fe.pbe-spn-kjpaw_psl.0.2.1.UPF       │ │
│ │ O       │ O.pbe-n-kjpaw_psl.0.1.UPF            │ │
│ │ Si      │ Si.pbe-n-rrkjus_psl.1.0.0.UPF        │ │
│ └─────────┴───────────────────────────────────────┘ │
└─────────────────────────────────────────────────────┘
```

## Features

### 1. Filtering
- **By Library**: Filter configurations by pseudopotential library (SSSP, PSLibrary, Pseudo Dojo, etc.)
- **By Functional**: Filter by exchange-correlation functional (PBE, LDA, PBEsol, etc.)

### 2. Auto-Population
- Automatically detects elements in the loaded structure
- Assigns appropriate pseudopotentials from the selected configuration
- Shows which pseudopotentials will be used for each element

### 3. Validation
- Checks if all required elements are available in the selected configuration
- Shows error message if any elements are missing
- Allows manual entry for missing elements

### 4. Persistence
- Selected configuration is stored in session state
- Persists across page navigation
- Separate selections for Calculation Setup and Workflow Builder

### 5. Fallback Mode
If no configurations are available, the selector falls back to manual entry mode:

```
┌─────────────────────────────────────────────────────┐
│ ⚠️ No pseudopotential configurations found.         │
│                                                      │
│ Please configure pseudopotentials first:            │
│ 1. Enable "Show Configuration" in the sidebar       │
│ 2. Go to "🧪 Pseudopotentials Configuration"        │
│ 3. Auto-detect and save your pseudopotential library│
│                                                      │
│ Manual Entry:                                        │
│ Pseudopotential for Fe: [Fe.UPF        ]           │
│ Pseudopotential for O:  [O.UPF         ]           │
│ Pseudopotential for Si: [Si.UPF        ]           │
└─────────────────────────────────────────────────────┘
```

## Usage Flow

1. **Configure Pseudopotentials** (one-time setup):
   - Go to "🧪 Pseudopotentials Configuration" page
   - Auto-detect pseudopotentials from your library directory
   - Save the configuration with name, library, version, and functional

2. **Select in Calculation/Workflow**:
   - Load a structure
   - Navigate to Calculation Setup or Workflow Builder
   - Use filters to narrow down options (optional)
   - Select your pseudopotential configuration
   - Pseudopotentials are automatically assigned

3. **Verify Selection**:
   - Expand "Configuration Details" to see metadata
   - Review "Selected Pseudopotentials" table
   - Proceed with calculation setup

## Session State Variables

### Calculation Setup Page
- `calc_selected_pseudo_config`: Name of selected configuration
- `calc_filter_library`: Current library filter
- `calc_filter_functional`: Current functional filter
- `workflow_config['pseudopotentials']`: Dict of element -> filename
- `workflow_config['pseudo_dir']`: Base path for pseudopotentials
- `workflow_config['_selected_pseudo_config_name']`: Selected config name (for persistence)

### Workflow Builder Page  
- `workflow_selected_pseudo_config`: Name of selected configuration
- `workflow_filter_library`: Current library filter
- `workflow_filter_functional`: Current functional filter
- `workflow_config['pseudopotentials']`: Dict of element -> filename
- `workflow_config['pseudo_dir']`: Base path for pseudopotentials
- `workflow_config['_selected_pseudo_config_name']`: Selected config name (for persistence)

## Code Structure

```
qtgui/utils/pseudopotentials_selector.py
└── render_pseudopotentials_selector(elements, config_dict, key_prefix)
    ├── Load available configurations
    ├── Render filter controls (library, functional)
    ├── Apply filters to configuration list
    ├── Render configuration selector
    ├── Load selected configuration
    ├── Display configuration details (expandable)
    ├── Validate element availability
    ├── Populate config_dict with pseudopotentials
    └── Store pseudo_dir and selection for persistence
```

## Benefits

1. **User-Friendly**: No need to manually specify pseudopotentials for each element
2. **Consistent**: Uses the same pseudopotentials across calculations
3. **Flexible**: Supports filtering and selection from multiple libraries
4. **Validated**: Checks that all required elements are available
5. **Persistent**: Selection persists across the interface
6. **Backward Compatible**: Falls back to manual entry if no configurations available

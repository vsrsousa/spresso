# GUI Modular Architecture

This document describes the modular architecture of the xespresso GUI, following the principle that **calculation and workflow modules create objects, while job submission executes them**.

## Overview

The GUI has been refactored to follow a modular design pattern where:

1. **Calculation modules** (`gui/calculations/`) create and prepare Espresso and atoms objects
2. **Workflow modules** (`gui/workflows/`) orchestrate multiple calculation steps
3. **Job submission** (`gui/pages/job_submission.py`) receives prepared objects and executes them
4. All logic uses **xespresso's existing definitions and patterns**

## Architecture

```
xespresso/gui/
├── calculations/          # Calculation preparation modules
│   ├── __init__.py
│   ├── base.py           # Base calculation preparation class
│   └── preparation.py    # GUI config → Espresso calculator
│
├── workflows/            # Workflow orchestration modules
│   ├── __init__.py
│   └── base.py          # Multi-step workflow coordination
│
└── pages/               # GUI pages (Streamlit)
    ├── calculation_setup.py    # Configures and prepares calculations
    ├── workflow_builder.py     # Builds multi-step workflows
    └── job_submission.py       # Executes prepared calculations
```

## Key Modules

### 1. Calculation Modules (`gui/calculations/`)

**Purpose**: Create Espresso and atoms objects from GUI configuration

**Key Classes/Functions**:

- `BaseCalculationPreparation`: Base class for calculation preparation
- `CalculationPreparation`: Prepares calculations from GUI config
- `prepare_calculation_from_gui(atoms, config, label)`: Creates Espresso calculator from GUI config
- `dry_run_calculation(atoms, config, label)`: Prepares calculator and writes input files

**Example Usage**:
```python
from xespresso.gui.calculations import prepare_calculation_from_gui

# In GUI calculation_setup page:
config = st.session_state.workflow_config  # From GUI inputs
atoms = st.session_state.current_structure  # From structure viewer

# Calculation module creates objects
prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label='scf/fe')

# Store for job submission
st.session_state.espresso_calculator = calc
st.session_state.prepared_atoms = prepared_atoms
```

### 2. Workflow Modules (`gui/workflows/`)

**Purpose**: Orchestrate multiple calculation steps

**Key Classes**:

- `GUIWorkflow`: Coordinates multi-step calculations
  - Uses calculation modules to prepare each step
  - Manages workflow state and results
  - Executes steps in sequence

**Example Usage**:
```python
from xespresso.gui.workflows import GUIWorkflow

# Create workflow
workflow = GUIWorkflow(atoms, config, base_label='workflow')

# Add calculation steps (each uses calculation module internally)
workflow.add_calculation('scf', scf_config)
workflow.add_calculation('relax', relax_config)

# Execute steps
workflow.run_step('scf')
workflow.run_step('relax')
```

### 3. GUI Pages

#### Calculation Setup (`calculation_setup.py`)

**Responsibility**: Configure parameters and prepare calculation objects

**Flow**:
1. User configures calculation parameters (ecutwfc, kpts, pseudopotentials, etc.)
2. Click "Prepare Calculation" button
3. Calls `prepare_calculation_from_gui()` to create objects
4. Stores prepared atoms and calculator in session state
5. Objects are now ready for job submission

#### Workflow Builder (`workflow_builder.py`)

**Responsibility**: Build multi-step workflows

**Flow**:
1. User selects workflow type (SCF, SCF+Relax, etc.)
2. Configures common parameters
3. Click "Build Workflow" button
4. Creates `GUIWorkflow` instance
5. Adds calculation steps (each prepared by calculation modules)
6. Stores workflow in session state

#### Job Submission (`job_submission.py`)

**Responsibility**: Execute prepared calculations

**Flow**:
1. Check if calculator exists in session state (prepared by calculation_setup or workflow_builder)
2. If yes: Use prepared calculator
3. If no: Call `prepare_calculation_from_gui()` to create one
4. Execute the calculation
5. Display results

**Key Point**: Job submission does NOT create Espresso objects directly anymore. It either uses pre-prepared objects or delegates to calculation modules.

## Design Principles

### 1. Separation of Concerns

- **Calculation modules**: Object creation and preparation
- **Workflow modules**: Multi-step coordination
- **Job submission**: Execution only

### 2. Use xespresso's Patterns

All modules use xespresso's existing classes and patterns:
- `Espresso` class from xespresso
- `input_data` dictionary structure
- Parameter naming conventions (ecutwfc, kpts, etc.)
- Queue configuration for job submission

### 3. Reusability

Calculation and workflow modules can be used:
- In GUI pages
- In command-line scripts
- In automated workflows
- For testing

### 4. Testability

Each module can be tested independently:
```python
# Test calculation preparation
from xespresso.gui.calculations import prepare_calculation_from_gui

config = {'pseudopotentials': {'Fe': 'Fe.UPF'}, 'ecutwfc': 50, ...}
atoms, calc = prepare_calculation_from_gui(test_atoms, config, 'test')
assert calc is not None
assert calc.label == 'test'
```

## Migration Guide

### Before (Monolithic)

```python
# In job_submission.py (OLD - DO NOT USE)
calc_params = {
    'pseudopotentials': config['pseudopotentials'],
    'label': label,
}
input_data = {}
input_data['ecutwfc'] = config['ecutwfc']
# ... many lines of parameter building ...
calc = Espresso(**calc_params)  # Created in job submission!
```

### After (Modular)

```python
# In calculation_setup.py (NEW)
from xespresso.gui.calculations import prepare_calculation_from_gui
atoms, calc = prepare_calculation_from_gui(atoms, config, label)
st.session_state.espresso_calculator = calc

# In job_submission.py (NEW)
calc = st.session_state.espresso_calculator  # Just use prepared object
prepared_atoms.calc = calc
energy = prepared_atoms.get_potential_energy()
```

## Benefits

1. **Cleaner Code**: Each module has a single, clear responsibility
2. **Easier Maintenance**: Changes to calculation preparation only affect calculation module
3. **Better Testing**: Can test calculation preparation independently from GUI
4. **Consistency**: All calculations use the same preparation logic
5. **Follows xespresso Patterns**: Uses xespresso's definitions consistently

## Example

See `examples/gui_modular_example.py` for a complete example demonstrating:
- Direct use of calculation modules
- Dry run functionality
- Multi-step workflows
- How modules work together

## Future Extensions

The modular structure makes it easy to add:
- New calculation types (NSCF, bands, phonon, etc.)
- New workflow templates
- Calculation validation and checking
- Parameter optimization
- Remote execution support (already integrated via queue parameter)

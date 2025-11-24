# GUI Structures Module

## Overview

The `gui/structures/` module handles structure loading and exporting operations for the xespresso GUI, following the modular design pattern where dedicated modules handle specific functionality.

## Architecture

```
xespresso/gui/
├── structures/              # Structure handling modules
│   ├── __init__.py
│   ├── base.py             # Base structure handler class
│   ├── loader.py           # Structure loading operations
│   └── exporter.py         # Structure export operations
│
└── pages/
    └── structure_viewer.py  # Coordinates user interaction using structures module
```

## Modules

### Base (`base.py`)

**Purpose**: Provides base class for structure handling

**Key Class**:
- `BaseStructureHandler`: Base class with common structure operations

### Loader (`loader.py`)

**Purpose**: Handles loading structures from various sources

**Key Classes/Functions**:
- `StructureLoader`: Main class for loading operations
  - `load_from_file(filepath)`: Load from file path
  - `load_from_upload(file_bytes, filename)`: Load from uploaded bytes
  - `find_structure_files(directory)`: Find structure files in directory
- `load_structure_from_file(filepath)`: Convenience function
- `load_structure_from_upload(file_bytes, filename)`: Convenience function

**Features**:
- Path validation and security checks
- Support for multiple file formats (CIF, XYZ, PDB, VASP, etc.)
- Safe directory traversal with depth limits
- Temporary file handling for uploads

### Exporter (`exporter.py`)

**Purpose**: Handles exporting structures to various formats

**Key Classes/Functions**:
- `StructureExporter`: Main class for export operations
  - `export_to_file(filepath, format)`: Export to file
  - `export_to_bytes(format)`: Export to bytes for downloads
- `export_structure(atoms, format)`: Convenience function

**Features**:
- Support for multiple export formats
- Temporary file handling
- Format auto-detection from extensions

## Usage in GUI

### Structure Viewer Page

The `structure_viewer.py` page has been refactored to use the structures module:

**Before (Monolithic)**:
```python
# All loading logic inline in the page
with tempfile.NamedTemporaryFile(...) as tmp:
    tmp.write(uploaded_file.getvalue())
    atoms = ase_io.read(tmp.name)
    os.unlink(tmp.name)
```

**After (Modular)**:
```python
# Use structures module
from xespresso.gui.structures import load_structure_from_upload

atoms, loader = load_structure_from_upload(
    uploaded_file.getvalue(),
    uploaded_file.name
)
```

### Example: File Upload

```python
from xespresso.gui.structures import load_structure_from_upload

# In Streamlit page
uploaded_file = st.file_uploader("Upload structure")
if uploaded_file:
    try:
        atoms, loader = load_structure_from_upload(
            uploaded_file.getvalue(),
            uploaded_file.name
        )
        
        # Store in session state
        st.session_state.current_structure = atoms
        st.session_state.structure_info = loader.get_info()
        
        st.success(f"Loaded: {uploaded_file.name}")
    except Exception as e:
        st.error(f"Error: {e}")
```

### Example: Directory Browser

```python
from xespresso.gui.structures import StructureLoader, load_structure_from_file

# Find structure files
structure_files = StructureLoader.find_structure_files(
    workdir,
    max_depth=3,
    validate_safety=True
)

# Load selected file
if selected_file:
    atoms, loader = load_structure_from_file(selected_file)
    st.session_state.current_structure = atoms
```

### Example: Export

```python
from xespresso.gui.structures import export_structure

# Export to bytes for download
file_data = export_structure(atoms, format='cif')

st.download_button(
    label="Download CIF",
    data=file_data,
    file_name="structure.cif"
)
```

## Design Principles

### 1. Separation of Concerns

- **Structures module**: Handles loading/export operations
- **GUI page**: Coordinates user interaction

### 2. Security

All loading operations include:
- Path validation to prevent traversal attacks
- Safe directory restrictions (home directory, /tmp)
- Symlink attack prevention
- Depth limits for directory traversal

### 3. Error Handling

- Comprehensive error handling with logging
- Clean temporary file cleanup
- Informative error messages

### 4. Reusability

The structures module can be used:
- In GUI pages
- In command-line scripts
- In automated workflows
- For testing

## Benefits

1. **Cleaner Code**: Loading logic separated from GUI display logic
2. **Easier Testing**: Can test loading/export independently
3. **Better Security**: Centralized security validation
4. **Consistency**: All structure operations use the same patterns
5. **Maintainability**: Changes to loading logic only affect structures module

## Integration with Other Modules

The structures module works alongside other GUI modules:

```
User Action → Structure Viewer Page
    ↓
Structures Module (load structure)
    ↓
Session State (store atoms)
    ↓
Calculation Setup Page
    ↓
Calculations Module (create Espresso + atoms)
    ↓
Job Submission (execute)
```

## Future Extensions

The modular structure makes it easy to add:
- New file format support
- Structure validation and checking
- Structure manipulation operations
- Conversion between formats
- Structure database integration

# Pseudopotentials Module

The pseudopotentials module provides a comprehensive system for managing pseudopotential configurations in xespresso, following the same architecture pattern as the machines and codes modules.

## Overview

Pseudopotentials are **primarily stored on the local user computer** after downloading them from sources like SSSP, PSLibrary, or Pseudo Dojo. During calculations, xespresso automatically copies the needed pseudopotential files to remote machines.

## Key Features

- **Auto-detection** of .UPF files from directories
- **Local-first approach** - pseudopotentials stored on your computer
- **Version tracking** for libraries (e.g., SSSP 1.1.2, PSLibrary 1.0.0)
- **Metadata extraction** from UPF files (element, functional, type, z_valence)
- **JSON storage** in `~/.xespresso/pseudopotentials/`
- **GUI integration** for easy configuration
- **Optional remote storage** if needed

## Quick Start

### Using Python API

```python
from xespresso.pseudopotentials import create_pseudopotentials_config

# Configure SSSP pseudopotentials (stored locally)
config = create_pseudopotentials_config(
    name="SSSP_efficiency",
    base_path="/home/user/pseudopotentials/SSSP_1.1.2_PBE_efficiency",
    description="SSSP 1.1.2 PBE efficiency pseudopotentials",
    functional="PBE",
    library="SSSP",
    version="1.1.2",
    save=True
)
```

### Using the GUI

1. Start the xespresso GUI: `streamlit run -m qtgui`
2. Check **"Show Configuration"** in the sidebar
3. Navigate to **"🧪 Pseudopotentials Configuration"**
4. Fill in the configuration form:
   - **Configuration Name**: e.g., "SSSP_efficiency"
   - **Pseudopotentials Directory**: path to .UPF files
   - **Functional**: e.g., "PBE"
   - **Library Name**: e.g., "SSSP"
   - **Library Version**: e.g., "1.1.2"
5. Click **"🔍 Auto-Detect Pseudopotentials"**
6. Review and click **"💾 Save"**

## Module Structure

```
xespresso/pseudopotentials/
├── __init__.py          # Module initialization and public API
├── config.py            # Data structures (Pseudopotential, PseudopotentialsConfig)
├── manager.py           # Management utilities (create, load, save, delete)
└── detector.py          # Auto-detection logic for .UPF files
```

## Data Structures

### Pseudopotential Class

Represents a single pseudopotential file:

```python
@dataclass
class Pseudopotential:
    element: str              # Chemical symbol (e.g., 'Fe', 'Si')
    filename: str             # Filename (e.g., 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF')
    path: str                 # Full path to file
    functional: Optional[str] # Exchange-correlation functional (e.g., 'PBE')
    type: Optional[str]       # Type (e.g., 'PAW', 'Ultrasoft', 'Norm-conserving')
    z_valence: Optional[float] # Number of valence electrons
```

### PseudopotentialsConfig Class

Configuration for a set of pseudopotentials:

```python
@dataclass
class PseudopotentialsConfig:
    name: str                          # Configuration name
    base_path: str                     # Directory containing .UPF files
    pseudopotentials: Dict[str, Pseudopotential]
    machine_name: Optional[str] = None # Optional remote machine name
    description: Optional[str] = None  # Description
    functional: Optional[str] = None   # Primary functional (e.g., 'PBE')
    library: Optional[str] = None      # Library name (e.g., 'SSSP')
    version: Optional[str] = None      # Library version (e.g., '1.1.2')
```

## Usage Examples

### Load and Use in Calculations

```python
from xespresso.pseudopotentials import load_pseudopotentials_config
from xespresso import Espresso
from ase.build import bulk

# Load configuration
config = load_pseudopotentials_config("SSSP_efficiency")

# Get pseudopotentials dictionary
pseudopotentials = config.get_pseudopotentials_dict()
# Returns: {'Fe': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF', ...}

# Use in calculation
atoms = bulk('Fe', 'bcc', a=2.87)
calc = Espresso(
    pseudopotentials=pseudopotentials,
    pseudo_dir=config.base_path,
    # ... other parameters
)
atoms.calc = calc
```

### List Available Configurations

```python
from xespresso.pseudopotentials import PseudopotentialsManager

# List all saved configurations
configs = PseudopotentialsManager.list_configs()
print(f"Available: {configs}")

# Load specific configuration
config = PseudopotentialsManager.load_config("SSSP_efficiency")

# Get available elements
elements = config.list_elements()
print(f"Elements: {elements}")
```

### Configure Remote Pseudopotentials (Optional)

If pseudopotentials are stored on a remote machine:

```python
config = create_pseudopotentials_config(
    name="remote_sssp",
    base_path="/opt/pseudopotentials/SSSP",
    machine_name="my_cluster",  # Reference to configured machine
    library="SSSP",
    version="1.1.2",
    functional="PBE"
)
```

## Best Practices

1. **Store locally**: Keep pseudopotentials on your local computer for easier management
2. **Use descriptive names**: Include library, version, and type (e.g., "SSSP_1.1.2_PBE_efficiency")
3. **Track versions**: Always specify library version for reproducibility
4. **Organize by functional**: Keep different functionals in separate configurations
5. **Document choices**: Use the description field to explain your selection
6. **Test first**: Verify with a small calculation before production runs

## Supported Pseudopotential Libraries

The module works with any .UPF format pseudopotential files, including:

- **SSSP** (Standard Solid State Pseudopotentials)
  - Download: https://www.materialscloud.org/sssp
  - Efficiency and precision sets available
  
- **PSLibrary** (Plane-wave Self-consistent Field Library)
  - Download: https://dalcorso.github.io/pslibrary/
  - Multiple functionals available

- **Pseudo Dojo**
  - Download: http://www.pseudo-dojo.org/
  - Optimized for precision

- **Custom pseudopotentials**
  - Any .UPF format files

## File Organization

Configurations are stored as JSON files in `~/.xespresso/pseudopotentials/`:

```
~/.xespresso/pseudopotentials/
├── SSSP_efficiency.json
├── SSSP_precision.json
├── PSLibrary_pbe.json
└── custom_pseudos.json
```

Each JSON file contains:
- Configuration metadata (name, library, version, functional)
- Base path to pseudopotential files
- Detected pseudopotentials with element-specific information
- Optional machine name for remote storage

## API Reference

### Main Functions

- `create_pseudopotentials_config()` - Auto-detect and create configuration
- `load_pseudopotentials_config()` - Load saved configuration
- `PseudopotentialsManager.list_configs()` - List all configurations
- `PseudopotentialsManager.delete_config()` - Delete a configuration

### Configuration Methods

- `config.get_pseudopotential(element)` - Get pseudopotential for element
- `config.get_pseudopotentials_dict()` - Get element->filename dictionary
- `config.list_elements()` - Get list of available elements
- `config.to_json(filepath)` - Save to JSON file
- `PseudopotentialsConfig.from_json(filepath)` - Load from JSON file

## See Also

- `examples/pseudopotentials_example.py` - Comprehensive usage examples
- Machine configuration: `xespresso/machines/`
- Codes configuration: `xespresso/codes/`

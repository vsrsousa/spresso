# spresso GUI - Streamlit Interface

A graphical user interface for configuring and running Quantum ESPRESSO calculations with spresso (xespresso fork).

## Features

### 🖥️ Machine Configuration
- Configure local and remote execution environments
- Support for multiple schedulers (SLURM, PBS, SGE, direct)
- SSH connection setup for remote clusters
- Environment modules management

### ⚙️ Codes Configuration
- Auto-detect Quantum ESPRESSO executables
- Support for multiple QE versions
- Both local and remote code detection
- Save and load code configurations

### 🔬 Structure Viewer
- Upload structure files (CIF, XYZ, POSCAR, PDB, etc.)
- Build simple structures (bulk crystals, molecules)
- Interactive 3D visualization with Plotly
- Display detailed structure information
- Export structures in various formats

### 📊 Calculation Setup
- Select calculation type (SCF, relax, bands, DOS, etc.)
- Configure pseudopotentials
- Set basic DFT parameters (cutoffs, convergence)
- K-point configuration (Monkhorst-Pack or k-spacing)
- Spin polarization settings

### 🔄 Workflow Builder
- Quality presets (fast, moderate, accurate)
- Workflow templates (quick SCF, quick relax)
- Integrated parameter management

### 🚀 Job Submission
- Configuration summary and validation
- Dry-run option
- Job submission to configured machines

## Installation

### Install with GUI support:

```bash
# Install spresso with GUI dependencies
pip install spresso[gui]

# Or install dependencies manually
pip install spresso streamlit plotly
```

**Note:** The PyPI package is named `spresso`, but the Python module is `xespresso` (for backwards compatibility).
Install with `pip install spresso`, but import as `from xespresso import ...`

### From source:

```bash
git clone https://github.com/vsrsousa/spresso.git
cd spresso
pip install -e ".[gui]"
```

## Usage

### Launch the GUI:

```bash
# Using the installed command
spresso-gui
# or for backwards compatibility
xespresso-gui

# Or using Python module
python -m xespresso.gui

# Or with streamlit directly
streamlit run xespresso/gui/streamlit_app.py
```

The GUI will open in your default web browser at `http://localhost:8501`

## Workflow

### Typical workflow:

1. **Configure Machine** (🖥️)
   - Set up local or remote execution environment
   - Configure scheduler if needed
   - Save machine configuration

2. **Configure Codes** (⚙️)
   - Auto-detect or manually specify QE executables
   - Save codes configuration for the selected machine

3. **Load Structure** (🔬)
   - Upload a structure file (CIF, XYZ, etc.)
   - Or build a simple structure
   - Visualize in 3D

4. **Setup Calculation** (📊)
   - Choose calculation type
   - Configure pseudopotentials
   - Set DFT parameters
   - Configure k-points

5. **Build Workflow** (🔄)
   - Select quality preset
   - Choose workflow type
   - Set calculation label/directory

6. **Submit Job** (🚀)
   - Review configuration
   - Submit to configured machine
   - Monitor job status

## Configuration Files

The GUI uses the standard configuration directories:

- **Machines**: `~/.xespresso/machines/`
- **Codes**: `~/.xespresso/codes/`
- **Pseudopotentials**: `~/.xespresso/pseudo_configs/`

**Note:** Configuration directories still use the `xespresso` name for backwards compatibility.

## Screenshots

### Machine Configuration
Configure computational environments with support for local and remote execution, schedulers, and environment modules.

### Structure Viewer
Upload and visualize molecular and crystal structures in 3D with detailed information display.

### Calculation Setup
Easy configuration of DFT parameters, pseudopotentials, and k-points.

### Workflow Builder
Use quality presets and workflow templates for quick calculation setup.

## Requirements

- Python >= 3.6
- spresso (provides xespresso module)
- streamlit >= 1.28.0
- plotly >= 5.17.0
- ASE >= 3.22.0
- numpy, scipy, matplotlib

## Tips

1. **First-time setup**: Start with Machine Configuration, then Codes Configuration
2. **Testing**: Use quality='fast' preset for quick testing
3. **Remote clusters**: Make sure SSH keys are properly configured
4. **Structure formats**: Most ASE-compatible formats work (CIF, XYZ, POSCAR, etc.)

## Troubleshooting

### GUI doesn't start
- Make sure streamlit is installed: `pip install streamlit`
- Check that spresso is installed: `pip install spresso`
- The Python module should be accessible as `xespresso`

### Module loading issues
- Verify that modules are available: `module avail`
- Check env_setup configuration for non-interactive shells

### Remote connection issues
- Test SSH connection manually first
- Verify SSH key permissions (should be 600)
- Check firewall settings

## Development

To contribute to the GUI:

```bash
# Clone repository
git clone https://github.com/vsrsousa/spresso.git
cd spresso

# Install in development mode
pip install -e ".[gui]"

# Run GUI
streamlit run xespresso/gui/streamlit_app.py
```

## License

GPL v3 - Same as spresso/xespresso

## Support

- [GitHub Issues](https://github.com/vsrsousa/spresso/issues)
- [Documentation](https://github.com/vsrsousa/spresso)

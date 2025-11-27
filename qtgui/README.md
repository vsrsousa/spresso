# xespresso Qt GUI

PyQt5-based graphical user interface for xespresso - Quantum ESPRESSO calculations.

## Overview

This is an alternative GUI for xespresso using PyQt5 instead of Streamlit. It provides the same functionality as the Streamlit GUI but with a native desktop application feel.

## Features

- **Machine Configuration**: Create and edit computational machine configurations
- **Codes Configuration**: Auto-detect and configure Quantum ESPRESSO executables
- **Pseudopotentials Configuration**: Browse and configure pseudopotential files
- **Structure Viewer**: Load and visualize atomic structures with 3D visualization
- **Calculation Setup**: Configure calculation parameters
- **Workflow Builder**: Build multi-step calculation workflows
- **Job Submission**: Generate input files and run calculations
- **Results & Post-Processing**: View results and perform post-processing analysis

## Installation

Make sure PyQt5 is installed:

```bash
pip install PyQt5>=5.15.0
```

Or install all requirements:

```bash
pip install -r requirements.txt
```

## Usage

### Running the GUI

From the repository root:

```bash
python -m qtgui
```

Or directly:

```bash
python qtgui/main_app.py
```

### From Python

```python
from qtgui.main_app import main
main()
```

## Architecture

The Qt GUI follows a similar modular architecture as the Streamlit GUI:

```
qtgui/
├── __init__.py          # Package initialization
├── __main__.py          # Entry point for python -m qtgui
├── main_app.py          # Main application window
├── pages/               # Page modules
│   ├── __init__.py
│   ├── machine_config.py
│   ├── codes_config.py
│   ├── pseudopotentials_config.py
│   ├── structure_viewer.py
│   ├── calculation_setup.py
│   ├── workflow_builder.py
│   ├── job_submission.py
│   └── results_postprocessing.py
└── utils/               # Utility modules
    ├── __init__.py
    └── validation.py
```

### Session State

The Qt GUI uses a `SessionState` class to manage application state across pages, similar to `st.session_state` in Streamlit:

```python
from qtgui.main_app import session_state

# Get/set state
session_state['current_structure'] = atoms
atoms = session_state.get('current_structure')
```

## Comparison with Streamlit GUI

| Feature | Streamlit GUI | Qt GUI |
|---------|--------------|--------|
| Web-based | Yes | No |
| Native desktop | No | Yes |
| Installation | streamlit>=1.28.0 | PyQt5>=5.15.0 |
| State management | st.session_state | SessionState class |
| Visualization | Plotly, py3Dmol | Matplotlib |
| File dialogs | Browser-based | Native OS dialogs |

## Requirements

- Python 3.8+
- PyQt5 >= 5.15.0
- ASE >= 3.22.0 (for structure handling)
- matplotlib >= 3.4.0 (for 3D visualization)
- xespresso modules (for machine/codes configuration)

## Development

### Adding new pages

1. Create a new page module in `qtgui/pages/`
2. Inherit from `QWidget`
3. Accept `session_state` in constructor
4. Implement `refresh()` method for page updates
5. Add import and export to `qtgui/pages/__init__.py`
6. Add to page list in `main_app.py`

### Example page structure

```python
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel

class MyPage(QWidget):
    def __init__(self, session_state):
        super().__init__()
        self.session_state = session_state
        self._setup_ui()
    
    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("My Page"))
    
    def refresh(self):
        # Called when page needs to be refreshed
        pass
```

## License

Same as the main xespresso project.

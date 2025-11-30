# xespresso Qt GUI

PyQt5-based graphical user interface for xespresso - Quantum ESPRESSO calculations.

## Overview

This is an alternative GUI for xespresso using PyQt5 instead of Streamlit. It provides the same functionality as the Streamlit GUI but with a native desktop application feel.

## Features

### Main Window
- **Workflow Navigation**: Quick access to Structure Viewer, Calculation Setup, Workflow Builder, Job Submission, and Results pages
- **Session Management**: Create, save, load, rename, and switch between multiple sessions
- **Working Directory**: Browse and set working directory for calculations
- **Toolbar**: Quick access to Configuration and Session management

### Configuration Dialog (Non-Blocking)
- **Machine Configuration**: Create and edit computational machine configurations
- **Codes Configuration**: Auto-detect and configure Quantum ESPRESSO executables
- **Pseudopotentials Configuration**: Browse and configure pseudopotential files
- Opens as a separate window that doesn't block the main application
- Access via the "⚙️ Open Configuration..." button or Edit menu (Ctrl+,)

### Session Management
- **Multiple Sessions**: Work with multiple calculation sessions simultaneously
- **Session Persistence**: Sessions are automatically saved to disk (~/.xespresso/sessions/)
- **Session Switching**: Switch between sessions without losing data
- **Session Names**: Give your sessions meaningful names for easy identification

### Workflow Pages
- **Structure Viewer**: Load and visualize atomic structures with 3D visualization
- **Calculation Setup**: Configure calculation parameters (cutoffs, k-points, pseudopotentials)
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

The Qt GUI follows a modular architecture:

```
qtgui/
├── __init__.py          # Package initialization
├── __main__.py          # Entry point for python -m qtgui
├── main_app.py          # Main application window and SessionState
├── dialogs/             # Non-blocking dialogs
│   ├── __init__.py
│   └── configuration_dialog.py  # Configuration dialog with tabs
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

### SessionState

The Qt GUI uses an enhanced `SessionState` class to manage application state across pages:

```python
from qtgui.main_app import session_state

# Get/set state
session_state['current_structure'] = atoms
atoms = session_state.get('current_structure')

# Session management
session_state.create_session("My Session")
session_state.save_session()
session_state.switch_session(session_id)
session_state.list_sessions()
session_state.get_session_name()
session_state.get_current_session_id()

# State change listeners
def on_change():
    print("State changed!")
session_state.add_listener(on_change)
```

### Configuration Dialog

The configuration dialog opens as a non-blocking window:

```python
from qtgui.dialogs import ConfigurationDialog
from qtgui.main_app import session_state

# Create and show dialog
dialog = ConfigurationDialog(session_state, parent=main_window)
dialog.show()  # Non-blocking

# Show specific tab
dialog.show_machine_tab()
dialog.show_codes_tab()
dialog.show_pseudopotentials_tab()
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
| Configuration | In main pages | Non-blocking dialog |
| Multiple sessions | Yes | Yes (with persistence) |

## Requirements

- Python 3.8+
- PyQt5 >= 5.15.0
- ASE >= 3.22.0 (for structure handling)
- matplotlib >= 3.4.0 (for 3D visualization)
- xespresso modules (for machine/codes configuration)

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| Ctrl+N | New Session |
| Ctrl+S | Save Session |
| Ctrl+, | Open Configuration Dialog |
| Ctrl+Q | Exit |
| Ctrl+1 | Structure Viewer |
| Ctrl+2 | Calculation Setup |
| Ctrl+3 | Workflow Builder |
| Ctrl+4 | Job Submission |
| Ctrl+5 | Results |

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

### Adding state listeners

```python
def on_state_change():
    # React to state changes
    self._update_display()

session_state.add_listener(on_state_change)
```

## Changelog

### Version 1.1.0
- Added non-blocking Configuration Dialog
- Improved SessionState with multiple session support
- Added session persistence (save/load to disk)
- Added session switching without data loss
- Added toolbar for quick access
- Reorganized sidebar with cleaner workflow navigation

### Version 1.0.0
- Initial Qt GUI implementation
- Basic page structure
- Session state management

## License

Same as the main xespresso project.

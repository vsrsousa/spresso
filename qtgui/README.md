# xespresso Qt GUI

PySide6-based graphical user interface for xespresso - Quantum ESPRESSO calculations.

## Overview

This is an alternative GUI for xespresso using PySide6 (Qt 6) instead of Streamlit. It provides the same functionality as the Streamlit GUI but with a native desktop application feel.

**Version 2.0.0**: Added simplified stable version alongside the full-featured version.

## Two Versions Available

### 1. Simplified Version (Recommended for Stability)

A streamlined, stable interface with:
- Simple tab-based navigation
- No complex state listeners (avoids signal recursion issues)
- Direct state management
- Lighter resource usage

```bash
python -m qtgui --simple
```

### 2. Full Version (More Features)

The full-featured version with:
- Advanced session management
- Configuration dialogs
- Multi-session support
- Full workflow pages

```bash
python -m qtgui
```

## Features

### Simplified Version
- **Structure Tab**: Load and view atomic structures
- **Calculation Tab**: Configure basic calculation parameters
- **Files Tab**: Browse calculation files
- **Info Tab**: Quick start guide

### Full Version
- **Workflow Navigation**: Quick access to Structure Viewer, Calculation Setup, Workflow Builder, Job Submission, and Results pages
- **Session Management**: Create, save, load, rename, and switch between multiple sessions
- **Working Directory**: Browse and set working directory for calculations
- **Toolbar**: Quick access to Configuration and Session management
- **Configuration Dialog**: Non-blocking dialog for Machine, Codes, and Pseudopotentials configuration

## Installation

Make sure PySide6 is installed:

```bash
pip install PySide6>=6.5.0
```

Or install all requirements:

```bash
pip install -r requirements.txt
```

## Usage

### Running the GUI

Simplified (stable) version:
```bash
python -m qtgui --simple
```

Full version:
```bash
python -m qtgui
```

### From Python

```python
# Simplified version
from qtgui.main_app_simple import main
main()

# Full version
from qtgui.main_app import main
main()
```

## Architecture

The Qt GUI has two main modules:

```
qtgui/
├── __init__.py           # Package initialization
├── __main__.py           # Entry point (--simple flag for simplified version)
├── main_app.py           # Full application with pages and session management
├── main_app_simple.py    # Simplified application with tabs
├── dialogs/              # Non-blocking dialogs (full version)
├── pages/                # Page modules (full version)
└── utils/                # Utility modules
```

## Troubleshooting

### GUI Crashes or Hangs

If you experience crashes or hangs with the full version, try the simplified version:

```bash
python -m qtgui --simple
```

The simplified version avoids the complex state management that can cause signal recursion issues.

### Common Issues

1. **Import errors**: Make sure PySide6 is installed: `pip install PySide6`
2. **Display errors**: On headless systems, set `QT_QPA_PLATFORM=offscreen`
3. **File browser slow**: The simplified version uses a flat directory listing instead of recursive scanning

## Requirements

- Python 3.8+
- PySide6 >= 6.5.0
- ASE >= 3.22.0 (for structure handling)
- matplotlib >= 3.4.0 (for 3D visualization, optional)
- xespresso modules (for machine/codes configuration, optional)

## Keyboard Shortcuts (Full Version)

| Shortcut | Action |
|----------|--------|
| Ctrl+N | New Session |
| Ctrl+S | Save Session |
| Ctrl+, | Open Configuration Dialog |
| Ctrl+Q | Exit |
| Ctrl+1-5 | Navigate to pages |

## Changelog

### Version 2.0.0
- **Added Simplified Version**: New simplified, stable GUI using `--simple` flag
- **Fixed File Browser**: Improved file browser with limits to prevent blocking
- **Fixed Signal Recursion**: Added guards to prevent infinite loops in session handling
- **Bug Fixes**: Fixed various stability issues in the full version

### Version 1.2.0
- **Migration to PySide6**: Migrated from PyQt5 to PySide6 (Qt 6) for improved performance
- **Lazy imports**: Page modules are now loaded lazily for faster startup time

### Version 1.1.0
- Added non-blocking Configuration Dialog
- Improved SessionState with multiple session support
- Added session persistence (save/load to disk)

### Version 1.0.0
- Initial Qt GUI implementation

## License

Same as the main xespresso project.

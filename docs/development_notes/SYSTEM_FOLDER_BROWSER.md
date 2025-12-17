# System Folder Browser Enhancement

## Overview

The working directory browser in xespresso GUI has been enhanced with a native system folder dialog option using tkinter. This allows users to browse anywhere on their file system using the familiar operating system file picker.

## Features

### New: "Browse System Folders" Button

When tkinter is available on the system, a new "Browse System Folders" button appears at the top of the sidebar navigation:

- **Location**: Sidebar, above the Quick Access buttons
- **Icon**: 📂 Browse System Folders
- **Function**: Opens native OS folder picker dialog
- **Fallback**: Gracefully disabled if tkinter is not available

### Native Dialog Benefits

1. **Full System Access**: Browse anywhere on the file system, not just subdirectories
2. **Familiar Interface**: Uses the operating system's native folder picker
3. **Quick Navigation**: Jump directly to any location without multiple clicks
4. **Keyboard Navigation**: Use OS shortcuts to navigate quickly

### Existing Features (Still Available)

The original directory browser features remain unchanged:

- Quick access buttons for common directories (Home, Documents, Desktop, etc.)
- Parent directory navigation (⬆️ button)
- Subfolder navigation with selectbox
- Custom path input with create directory option

## Usage

### For Users

Simply click the "📂 Browse System Folders" button in the sidebar to open the native folder picker. Select your desired folder and it will be set as the working directory.

If the button doesn't appear, it means tkinter is not available on your system (see Installation below).

### For Developers

The functionality is implemented in `qtgui/utils/directory_browser.py`:

```python
from qtgui.utils.directory_browser import render_directory_browser

# Render the enhanced directory browser
selected_path = render_directory_browser(
    key="my_browser",
    initial_path="/home/user",
    help_text="Choose your working directory"
)
```

## Installation

### Linux

Install python3-tk package:

```bash
# Ubuntu/Debian
sudo apt-get install python3-tk

# Fedora/RHEL
sudo dnf install python3-tkinter

# Arch
sudo pacman -S tk
```

### macOS

tkinter comes pre-installed with Python on macOS.

### Windows

tkinter comes pre-installed with Python on Windows.

## Implementation Details

### File: `qtgui/utils/directory_browser.py`

#### New Function: `open_folder_dialog(initial_path)`

Opens a native folder selection dialog:

```python
def open_folder_dialog(initial_path: str) -> Optional[str]:
    """
    Open a native tkinter folder selection dialog.
    
    Args:
        initial_path: Initial directory to show in the dialog
        
    Returns:
        Selected folder path or None if cancelled/unavailable
    """
```

**Features**:
- Creates a hidden root window (no visible window frame)
- Sets dialog to topmost (appears above Streamlit window)
- Returns None if cancelled or if tkinter is unavailable
- Handles exceptions gracefully

#### Modified Function: `render_directory_browser()`

Enhanced to include the "Browse System Folders" button when tkinter is available:

```python
# Add "Browse System Folders" button if tkinter is available
if TKINTER_AVAILABLE:
    st.sidebar.markdown("**System Folder Browser:**")
    if st.sidebar.button("📂 Browse System Folders", ...):
        selected_folder = open_folder_dialog(current_path)
        if selected_folder:
            st.session_state[f'{key}_current_path'] = selected_folder
            st.rerun()
```

### Graceful Fallback

The implementation includes proper fallback handling:

```python
# Try to import tkinter for native file dialog
try:
    import tkinter as tk
    from tkinter import filedialog
    TKINTER_AVAILABLE = True
except ImportError:
    TKINTER_AVAILABLE = False
```

When tkinter is not available:
- The "Browse System Folders" button is not displayed
- All other directory browser features work normally
- No errors or warnings are shown to users

## Testing

### Manual Testing

1. Start the xespresso GUI: `xespresso-gui` or `python -m qtgui`
2. Look for the "📂 Browse System Folders" button in the sidebar
3. Click the button to open the native folder picker
4. Select a folder and verify it becomes the working directory

### Automated Testing

The directory browser can be imported and tested:

```python
from qtgui.utils.directory_browser import (
    render_directory_browser,
    open_folder_dialog,
    TKINTER_AVAILABLE
)

# Check if tkinter is available
print(f"tkinter available: {TKINTER_AVAILABLE}")

# Test folder dialog (returns None in headless mode)
result = open_folder_dialog("/tmp")
```

## Screenshots

### With tkinter Available

```
┌─ Working Directory ─────────────────────┐
│ System Folder Browser:                  │
│ [📂 Browse System Folders]              │
│ ─────────────────────────────────────── │
│ Quick Access:                           │
│ [🏠] [📂] [📊] [📄] [🖥️]              │
│ ─────────────────────────────────────── │
│ Current Path:                           │
│ 📁 /home/user/calculations              │
│ ...                                      │
└─────────────────────────────────────────┘
```

### Without tkinter (Fallback)

```
┌─ Working Directory ─────────────────────┐
│ Quick Access:                           │
│ [🏠] [📂] [📊] [📄] [🖥️]              │
│ ─────────────────────────────────────── │
│ Current Path:                           │
│ 📁 /home/user/calculations              │
│ ...                                      │
└─────────────────────────────────────────┘
```

## Advantages

1. **Better UX**: Users can navigate to any folder on the system quickly
2. **Familiar**: Uses OS-native dialogs that users already know
3. **Safe**: Proper error handling and graceful fallback
4. **No Breaking Changes**: Existing functionality remains intact
5. **Optional**: Works without tkinter, just with reduced features

## Known Limitations

1. **Display Requirement**: Requires a display/GUI environment (won't work in pure SSH sessions without X11 forwarding)
2. **System Dependency**: Requires python3-tk to be installed on Linux
3. **Modal Dialog**: Blocks Streamlit while dialog is open (expected behavior for native dialogs)

## Future Enhancements

Possible future improvements:

1. Add bookmarks/favorites for frequently used directories
2. Add recent directories list
3. Add directory creation from the native dialog
4. Add support for other native file pickers (e.g., zenity on Linux)

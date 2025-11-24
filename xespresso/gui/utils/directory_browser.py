"""
Directory browsing utilities for xespresso GUI.

Provides functionality to browse and select directories with subfolder navigation.
"""

import streamlit as st
import os
from pathlib import Path
from typing import Optional, List
import threading

# Try to import tkinter for native file dialog
try:
    import tkinter as tk
    from tkinter import filedialog

    TKINTER_AVAILABLE = True
except ImportError:
    TKINTER_AVAILABLE = False


def get_subdirectories(path: str) -> List[str]:
    """
    Get list of subdirectories in the given path.

    Args:
        path: Directory path to scan

    Returns:
        List of subdirectory paths
    """
    if not os.path.exists(path) or not os.path.isdir(path):
        return []

    try:
        subdirs = []
        for item in os.listdir(path):
            item_path = os.path.join(path, item)
            if os.path.isdir(item_path) and not item.startswith("."):
                subdirs.append(item_path)
        return sorted(subdirs)
    except (PermissionError, OSError):
        return []


def open_folder_dialog(initial_path: str) -> Optional[str]:
    """
    Open a native tkinter folder selection dialog.

    Args:
        initial_path: Initial directory to show in the dialog

    Returns:
        Selected folder path or None if cancelled/unavailable
    """
    if not TKINTER_AVAILABLE:
        return None

    try:
        # On macOS, we need to handle tkinter differently
        import sys
        import platform

        # Create a hidden root window
        root = tk.Tk()
        root.withdraw()

        # On macOS, we need to bring the dialog to front
        if platform.system() == "Darwin":
            # macOS specific: lift the window and make it topmost
            root.lift()
            root.attributes("-topmost", True)
            root.focus_force()

            # Call update to process the window management commands
            root.update()
        else:
            root.attributes("-topmost", True)

        # Open folder dialog
        folder_path = filedialog.askdirectory(
            initialdir=initial_path, title="Select Working Directory", parent=root
        )

        # Clean up
        root.destroy()

        # Return the selected path (empty string if cancelled)
        return folder_path if folder_path else None

    except Exception as e:
        # If tkinter fails for any reason, return None
        # On some systems, especially macOS with certain Python installations,
        # tkinter may not be properly configured
        return None


def render_directory_browser(
    key: str = "directory_browser",
    initial_path: Optional[str] = None,
    help_text: str = "Choose the base directory where calculation folders will be created",
) -> str:
    """
    Render an enhanced directory browser with subfolder navigation.

    Args:
        key: Unique key for the component
        initial_path: Initial directory path
        help_text: Help text to display

    Returns:
        Selected directory path
    """
    if initial_path is None:
        initial_path = os.path.expanduser("~")

    # Initialize browsing state
    if f"{key}_current_path" not in st.session_state:
        st.session_state[f"{key}_current_path"] = initial_path

    current_path = st.session_state[f"{key}_current_path"]

    # Add "Browse System Folders" button if tkinter is available
    if TKINTER_AVAILABLE:
        st.sidebar.markdown("**System Folder Browser:**")
        if st.sidebar.button(
            "📂 Browse System Folders",
            key=f"{key}_tkinter_browse",
            help="Open native file dialog to browse anywhere on your system",
            use_container_width=True,
        ):
            selected_folder = open_folder_dialog(current_path)
            if selected_folder:
                st.session_state[f"{key}_current_path"] = selected_folder
                st.rerun()
            else:
                # Check if we're on macOS and provide helpful message
                import platform

                if platform.system() == "Darwin":
                    st.sidebar.warning(
                        "⚠️ Could not open file dialog. On macOS, this feature may not work "
                        "with certain Python installations. Please use the alternative navigation "
                        "methods below (Quick Access, Custom Path, or subfolder navigation)."
                    )
                else:
                    st.sidebar.info("No folder selected")
        st.sidebar.markdown("---")
    else:
        # If tkinter is not available, show a note about alternatives
        st.sidebar.info(
            "💡 **Note**: System file browser not available. "
            "Use Quick Access buttons or Custom Path entry below."
        )

    # Common directories as quick shortcuts
    st.sidebar.markdown("**Quick Access:**")
    common_dirs = {
        "🏠 Home": os.path.expanduser("~"),
        "📂 Current Directory": os.getcwd(),
        "📊 Calculations": os.path.join(os.path.expanduser("~"), "calculations"),
        "📄 Documents": os.path.join(os.path.expanduser("~"), "Documents"),
        "🖥️ Desktop": os.path.join(os.path.expanduser("~"), "Desktop"),
    }

    cols = st.sidebar.columns(len(common_dirs))
    for idx, (label, path) in enumerate(common_dirs.items()):
        with cols[idx]:
            if st.button(
                label.split()[0],
                key=f"{key}_quick_{idx}",
                help=label,
                use_container_width=True,
            ):
                if os.path.exists(path):
                    st.session_state[f"{key}_current_path"] = path
                    st.rerun()

    st.sidebar.markdown("---")

    # Current path display with parent navigation
    st.sidebar.markdown("**Current Path:**")

    # Show path with ability to go up
    parent_path = os.path.dirname(current_path)
    if parent_path and parent_path != current_path:
        col1, col2 = st.sidebar.columns([1, 4])
        with col1:
            if st.button("⬆️", key=f"{key}_up", help="Go to parent directory"):
                st.session_state[f"{key}_current_path"] = parent_path
                st.rerun()
        with col2:
            st.caption(f"📁 {current_path}")
    else:
        st.sidebar.caption(f"📁 {current_path}")

    # Browse subdirectories
    subdirs = get_subdirectories(current_path)

    if subdirs:
        st.sidebar.markdown("**Subfolders:**")

        # Show subdirectories in a selectbox for better UX
        subdir_names = [os.path.basename(d) for d in subdirs]

        # Add option to stay in current directory
        display_options = ["(Stay here)"] + subdir_names

        selected_subdir = st.sidebar.selectbox(
            "Select subfolder:",
            options=display_options,
            key=f"{key}_subdir_selector",
            help="Choose a subfolder or stay in current directory",
        )

        if selected_subdir != "(Stay here)":
            idx = subdir_names.index(selected_subdir)
            col1, col2 = st.sidebar.columns([3, 1])
            with col2:
                if st.button("→", key=f"{key}_enter_subdir", help="Enter this folder"):
                    st.session_state[f"{key}_current_path"] = subdirs[idx]
                    st.rerun()
    else:
        st.sidebar.info("📂 No subfolders in this directory")

    st.sidebar.markdown("---")

    # Custom path input
    with st.sidebar.expander("✏️ Enter Custom Path", expanded=False):
        custom_path = st.text_input(
            "Enter path:",
            value=current_path,
            key=f"{key}_custom_path",
            help="Enter a custom directory path",
        )

        col1, col2 = st.sidebar.columns(2)
        with col1:
            if st.button("Go", key=f"{key}_go_custom", use_container_width=True):
                # Validate path
                expanded_path = os.path.expanduser(custom_path)
                expanded_path = os.path.abspath(expanded_path)

                if os.path.exists(expanded_path) and os.path.isdir(expanded_path):
                    st.session_state[f"{key}_current_path"] = expanded_path
                    st.rerun()
                else:
                    st.error("❌ Invalid directory path")

        with col2:
            if st.button("Create", key=f"{key}_create_dir", use_container_width=True):
                expanded_path = os.path.expanduser(custom_path)
                expanded_path = os.path.abspath(expanded_path)

                try:
                    os.makedirs(expanded_path, exist_ok=True)
                    st.session_state[f"{key}_current_path"] = expanded_path
                    st.success(f"✅ Created")
                    st.rerun()
                except Exception as e:
                    st.error(f"❌ Error: {e}")

    st.sidebar.markdown("---")
    st.sidebar.info(help_text)

    return current_path

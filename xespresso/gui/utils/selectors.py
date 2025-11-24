"""
Selector utilities for machine and code selection in the GUI.

These functions provide reusable components for selecting already-configured
machines and codes, separate from the configuration/setup pages.
"""

import streamlit as st
import os

try:
    from xespresso.machines.config.loader import (
        load_machine,
        list_machines,
        DEFAULT_CONFIG_PATH,
        DEFAULT_MACHINES_DIR,
    )
    from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR

    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False


def render_machine_selector(key="machine_selector", help_text=None):
    """
    Render a machine selector dropdown.

    This is separate from machine configuration - it only allows selecting
    from already-configured machines.

    Args:
        key: Unique key for the selector widget
        help_text: Optional help text for the selector

    Returns:
        tuple: (machine_name, machine_object) or (None, None) if no selection
    """
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available.")
        return None, None

    try:
        machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)

        if not machines_list:
            st.warning(
                "⚠️ No machines configured. Please configure a machine first in the Machine Configuration page."
            )
            return None, None

        # Use session state for persistent selection
        default_idx = 0
        if (
            st.session_state.get("current_machine_name")
            and st.session_state.current_machine_name in machines_list
        ):
            default_idx = machines_list.index(st.session_state.current_machine_name)

        selected_machine = st.selectbox(
            "Select Machine:",
            machines_list,
            index=default_idx,
            key=key,
            help=help_text or "Choose a configured machine for your calculations",
        )

        if selected_machine:
            try:
                machine = load_machine(
                    DEFAULT_CONFIG_PATH,
                    selected_machine,
                    DEFAULT_MACHINES_DIR,
                    return_object=True,
                )
                st.session_state.current_machine_name = selected_machine
                st.session_state.current_machine = machine

                # Show machine info
                with st.expander("📋 Machine Details", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"**Name:** {machine.name}")
                        st.write(f"**Execution:** {machine.execution}")
                        st.write(f"**Scheduler:** {machine.scheduler}")
                    with col2:
                        st.write(f"**Working Dir:** {machine.workdir}")
                        st.write(f"**Processors:** {machine.nprocs}")
                        if machine.is_remote:
                            st.write(f"**Host:** {machine.host}")

                return selected_machine, machine
            except Exception as e:
                st.error(f"Error loading machine: {e}")
                return None, None
    except Exception as e:
        st.warning(f"Could not load machines list: {e}")
        return None, None


def render_codes_selector(machine_name, key="codes_selector", help_text=None):
    """
    Render a codes/version selector for a specific machine.

    This is separate from codes configuration - it only allows selecting
    from already-configured codes and versions.

    Args:
        machine_name: Name of the machine to load codes for
        key: Unique key for the selector widget
        help_text: Optional help text for the selector

    Returns:
        CodesConfig object or None if no selection
    """
    if not XESPRESSO_AVAILABLE or not machine_name:
        return None

    try:
        codes_config = load_codes_config(machine_name, DEFAULT_CODES_DIR)

        if not codes_config:
            st.warning(
                f"⚠️ No codes configured for machine '{machine_name}'. Please configure codes first in the Codes Configuration page."
            )
            return None

        # Check if multiple versions are available
        if codes_config.versions:
            available_versions = codes_config.list_versions()

            # Show info about multiple versions if available
            if len(available_versions) > 1:
                st.info(
                    f"📦 Multiple QE versions available: {', '.join(available_versions)}"
                )

            # Version selector
            default_idx = 0
            if (
                st.session_state.get("selected_code_version")
                and st.session_state.selected_code_version in available_versions
            ):
                default_idx = available_versions.index(
                    st.session_state.selected_code_version
                )

            selected_version = st.selectbox(
                "Select QE Version:",
                available_versions,
                index=default_idx,
                key=key,
                help=help_text or "Choose which Quantum ESPRESSO version to use",
            )

            if selected_version:
                # Load codes for the selected version
                version_config = load_codes_config(
                    machine_name, DEFAULT_CODES_DIR, version=selected_version
                )
                st.session_state.selected_code_version = selected_version
                st.session_state.current_codes = version_config

                # Show version info
                with st.expander("⚙️ Version Details", expanded=False):
                    version_codes = version_config.get_all_codes(
                        version=selected_version
                    )
                    st.write(f"**Version:** {selected_version}")
                    if (
                        version_config.versions
                        and selected_version in version_config.versions
                    ):
                        label = version_config.versions[selected_version].get("label")
                        if label:
                            st.write(f"**Label:** {label}")
                        
                        # ALWAYS show modules if they exist in the codes JSON
                        # regardless of how many versions there are
                        modules = version_config.versions[selected_version].get('modules')
                        if modules:
                            st.write(f"**Modules:** {', '.join(modules)}")
                    st.write(f"**Codes:** {len(version_codes)} executables configured")

                    # Show code list
                    code_names = list(version_codes.keys())
                    st.write(f"**Available codes:** {', '.join(code_names)}")

                return version_config
        else:
            # Single version or no version structure
            st.session_state.current_codes = codes_config

            with st.expander("⚙️ Codes Details", expanded=False):
                all_codes = codes_config.get_all_codes()
                st.write(f"**Codes:** {len(all_codes)} executables configured")
                code_names = list(all_codes.keys())
                st.write(f"**Available codes:** {', '.join(code_names)}")

            return codes_config

    except Exception as e:
        st.error(f"Error loading codes: {e}")
        return None


def render_workdir_browser_with_button(current_dir=None, key="workdir_browser", label="Output Directory"):
    """
    Render a folder browser that looks similar to st.file_uploader.
    
    This provides a button-based interface for folder selection that mimics
    the appearance of the "Browse Files" button in file upload widgets.
    
    Args:
        current_dir: Current working directory path
        key: Unique key for the widget
        label: Label to display for the directory selector
        
    Returns:
        str: Selected directory path
    """
    import streamlit as st
    import os
    
    # Use per-key session state to avoid conflicts between multiple instances
    workdir_key = f"{key}_workdir"
    browse_mode_key = f"{key}_browse_mode"
    
    if current_dir is None:
        current_dir = st.session_state.get(workdir_key, os.path.join(os.getcwd(), "calculations"))
    
    # Initialize session state if not set
    if workdir_key not in st.session_state:
        st.session_state[workdir_key] = current_dir
    if browse_mode_key not in st.session_state:
        st.session_state[browse_mode_key] = False
    
    # Main container with similar styling to file_uploader
    st.markdown(f"**{label}**")
    
    # Create a button that looks like file uploader's browse button
    col1, col2 = st.columns([3, 1])
    
    with col1:
        # Display current selection
        current_path = st.session_state[workdir_key]
        st.text_input(
            "Selected folder:",
            value=current_path,
            key=f"{key}_display",
            disabled=True,
            label_visibility="collapsed"
        )
    
    with col2:
        if st.button("📁 Browse Folders", key=f"{key}_browse_btn", use_container_width=True):
            st.session_state[browse_mode_key] = not st.session_state[browse_mode_key]
            st.rerun()
    
    # Show folder browser when button is clicked (similar to file dialog)
    if st.session_state[browse_mode_key]:
        st.markdown("---")
        st.markdown("**📂 Select a folder:**")
        
        # Quick navigation buttons
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("🏠 Home", key=f"{key}_home", use_container_width=True):
                st.session_state[workdir_key] = os.path.expanduser("~")
                st.rerun()
        
        with col2:
            if st.button("📂 Current Dir", key=f"{key}_current", use_container_width=True):
                st.session_state[workdir_key] = os.getcwd()
                st.rerun()
        
        with col3:
            if st.button("⬆️ Parent", key=f"{key}_parent", use_container_width=True):
                current = st.session_state[workdir_key]
                st.session_state[workdir_key] = os.path.dirname(current)
                st.rerun()
        
        # Current path display
        workdir = st.session_state[workdir_key]
        st.info(f"📍 Current: `{workdir}`")
        
        # Manual path entry - use on_change callback to prevent infinite loop
        # Initialize the manual input key if not set
        manual_input_key = f"{key}_manual_input"
        if manual_input_key not in st.session_state:
            st.session_state[manual_input_key] = workdir
        
        new_path = st.text_input(
            "Or enter path manually:",
            value=st.session_state[manual_input_key],
            key=manual_input_key,
            help="Type or paste the full path to your desired folder"
        )
        
        # Only update workdir if the user actually changed the path
        # and it's different from current workdir (prevents infinite loop)
        if new_path and new_path != st.session_state[manual_input_key]:
            st.session_state[manual_input_key] = new_path
        
        if new_path:
            try:
                normalized_new_path = os.path.abspath(os.path.expanduser(new_path))
                normalized_workdir = os.path.abspath(os.path.expanduser(workdir))
                
                if normalized_new_path != normalized_workdir:
                    if os.path.exists(normalized_new_path) and os.path.isdir(normalized_new_path):
                        st.session_state[workdir_key] = normalized_new_path
                        st.session_state[manual_input_key] = normalized_new_path
                        st.rerun()
            except Exception:
                pass
        
        # Folder list for navigation
        try:
            # Security validation
            if os.path.isabs(workdir):
                real_workdir = os.path.realpath(workdir)
                
                if os.path.exists(real_workdir) and os.path.isdir(real_workdir):
                    # Get subdirectories
                    contents = os.listdir(real_workdir)
                    subdirs = []
                    for d in contents:
                        if d.startswith("."):
                            continue
                        subdir_path = os.path.join(real_workdir, d)
                        if os.path.isdir(subdir_path):
                            try:
                                # Validate path doesn't escape parent
                                if os.path.commonpath([real_workdir, os.path.realpath(subdir_path)]) == real_workdir:
                                    subdirs.append(d)
                            except (ValueError, OSError):
                                continue
                    subdirs.sort()
                    
                    if subdirs:
                        st.markdown("**Folders in current directory:**")
                        
                        # Use dropdown for folder selection
                        selected_subdir = st.selectbox(
                            "Select a folder:",
                            options=[""] + subdirs,
                            format_func=lambda x: "-- Select a folder --" if x == "" else f"📁 {x}",
                            key=f"{key}_folder_dropdown",
                            help="Choose a folder to navigate into"
                        )
                        
                        # Navigate button (only enabled when a folder is selected)
                        if selected_subdir:
                            if st.button("→ Open Folder", key=f"{key}_open_selected", 
                                        type="primary",
                                        use_container_width=True):
                                # Validate to prevent path traversal
                                if ".." not in selected_subdir and "/" not in selected_subdir and "\\" not in selected_subdir:
                                    new_workdir = os.path.realpath(os.path.join(real_workdir, selected_subdir))
                                    # Ensure new path is within parent
                                    try:
                                        if os.path.commonpath([real_workdir, new_workdir]) == real_workdir:
                                            st.session_state[workdir_key] = new_workdir
                                            st.rerun()
                                    except ValueError:
                                        pass
                    else:
                        st.info("ℹ️ No subfolders in current directory")
        except PermissionError:
            st.warning("⚠️ Permission denied to list directory contents")
        except Exception as e:
            st.warning(f"⚠️ Could not list folders: {e}")
        
        # Done button to close browser
        if st.button("✅ Use This Folder", key=f"{key}_done", type="primary", use_container_width=True):
            st.session_state[browse_mode_key] = False
            st.rerun()
        
        st.markdown("---")
    
    return st.session_state[workdir_key]


def render_workdir_browser(current_dir=None, key="workdir_browser"):
    """
    Render a clean and simple working directory browser/selector.

    This is a streamlined version that provides a cleaner interface
    similar to system file browsers, with direct folder navigation.

    Features:
    - Text input for direct path entry  
    - Quick access buttons (Current, Home, Parent)
    - Simple folder list for direct selection (like system file browser)

    Args:
        current_dir: Current working directory path
        key: Unique key for the widget

    Returns:
        str: Selected directory path
    """
    # Use per-key session state to avoid conflicts between multiple instances
    workdir_key = f"{key}_workdir"
    
    if current_dir is None:
        current_dir = st.session_state.get(workdir_key, st.session_state.get("local_workdir", os.getcwd()))

    # Initialize session state if not set
    if workdir_key not in st.session_state:
        st.session_state[workdir_key] = current_dir

    st.subheader("📁 Working Directory")

    # Path input and quick access buttons in a clean layout
    col1, col2, col3, col4 = st.columns([4, 1, 1, 1])

    with col1:
        # Don't use a separate key for the text input to avoid state conflicts
        # Just use value parameter and read from the return value
        workdir_input = st.text_input(
            "Directory Path:",
            value=st.session_state[workdir_key],
            help="Enter the path to your working directory",
        )
    
    # Update session state only if value actually changed
    if workdir_input:
        normalized_workdir = os.path.abspath(os.path.expanduser(workdir_input))
        # Only update and rerun if the normalized path is different
        if st.session_state[workdir_key] != normalized_workdir:
            st.session_state[workdir_key] = normalized_workdir
            # Don't call st.rerun() here - let Streamlit handle it naturally

    with col2:
        if st.button(
            "📂 Current", key=f"{key}_current", help="Go to current working directory"
        ):
            st.session_state[workdir_key] = os.getcwd()
            st.rerun()

    with col3:
        if st.button("🏠 Home", key=f"{key}_home", help="Go to home directory"):
            st.session_state[workdir_key] = os.path.expanduser("~")
            st.rerun()

    with col4:
        if st.button("⬆️ Up", key=f"{key}_parent", help="Go to parent directory"):
            current_dir = st.session_state[workdir_key]
            st.session_state[workdir_key] = os.path.dirname(current_dir)
            st.rerun()

    # Validate and normalize directory (use the value from session state)
    workdir = st.session_state[workdir_key]
    if workdir:
        try:
            if os.path.exists(workdir) and os.path.isdir(workdir):
                st.success(f"✅ `{workdir}`")

                # Clean folder browser - similar to system file browser
                st.markdown("---")
                st.subheader("📂 Browse Folders")

                try:
                    # Security validation
                    if not os.path.isabs(workdir):
                        st.error("❌ Invalid path: must be absolute")
                        return current_dir

                    # Resolve symlinks
                    real_workdir = os.path.realpath(workdir)

                    # Get subdirectories
                    contents = os.listdir(real_workdir)
                    subdirs = []
                    for d in contents:
                        # Skip hidden directories
                        if d.startswith("."):
                            continue
                        subdir_path = os.path.join(real_workdir, d)
                        # Validate path doesn't escape parent
                        if (
                            os.path.isdir(subdir_path)
                            and os.path.commonpath(
                                [real_workdir, os.path.realpath(subdir_path)]
                            )
                            == real_workdir
                        ):
                            subdirs.append(d)
                    subdirs.sort()

                    if subdirs:
                        st.info("💡 Select a folder from the dropdown to navigate into it")
                        
                        # Use dropdown for folder selection
                        selected_subdir = st.selectbox(
                            "Available folders:",
                            options=[""] + subdirs,  # Empty option to allow no selection
                            format_func=lambda x: "-- Select a folder --" if x == "" else f"📁 {x}",
                            key=f"{key}_folder_dropdown",
                            help="Choose a folder to navigate into"
                        )
                        
                        # Navigate button (only enabled when a folder is selected)
                        col1, col2 = st.columns([3, 1])
                        with col2:
                            if st.button("→ Open", key=f"{key}_open_folder", 
                                        disabled=(selected_subdir == ""),
                                        type="primary" if selected_subdir else "secondary",
                                        use_container_width=True):
                                # Validate to prevent path traversal
                                if selected_subdir and ".." not in selected_subdir and "/" not in selected_subdir and "\\" not in selected_subdir:
                                    new_workdir = os.path.realpath(
                                        os.path.join(real_workdir, selected_subdir)
                                    )
                                    # Ensure new path is within parent
                                    if os.path.commonpath([real_workdir, new_workdir]) == real_workdir:
                                        st.session_state[workdir_key] = new_workdir
                                        st.rerun()
                                else:
                                    st.error("❌ Invalid folder name")
                        
                        with col1:
                            if selected_subdir:
                                st.caption(f"Selected: {selected_subdir}")
                    else:
                        st.info("ℹ️ No subfolders in current directory")

                except PermissionError:
                    st.warning("⚠️ Permission denied to list directory contents")
                except Exception as e:
                    st.warning(f"⚠️ Could not list subfolders: {e}")

                return workdir
            else:
                st.error(f"❌ Directory does not exist: {workdir}")
                st.info("💡 Use the Home or Current button to navigate to a valid directory")
                return current_dir
        except Exception as e:
            st.error(f"❌ Invalid path: {e}")
            return current_dir

    return current_dir

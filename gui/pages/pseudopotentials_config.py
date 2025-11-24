"""
Pseudopotentials Configuration Page for xespresso GUI.

This module handles the pseudopotentials configuration interface,
allowing users to:
- Auto-detect pseudopotential files on machines
- Configure multiple pseudopotential libraries
- Save and load pseudopotential configurations
"""

import streamlit as st
import traceback
import os

try:
    from xespresso.machines.config.loader import (
        list_machines,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    from xespresso.pseudopotentials.manager import (
        create_pseudopotentials_config,
        load_pseudopotentials_config,
        PseudopotentialsManager,
        DEFAULT_PSEUDOPOTENTIALS_DIR
    )
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False


def render_pseudopotentials_config_page():
    """Render the pseudopotentials configuration page."""
    st.header("⚙️ Pseudopotentials Configuration")
    st.markdown("""
    **Configure and auto-detect** pseudopotential files for Quantum ESPRESSO calculations.
    
    Pseudopotentials are primarily stored on your **local computer** after downloading them.
    During calculations, xespresso automatically copies the needed pseudopotentials to remote machines.
    
    This page is for **configuration only** - once pseudopotentials are saved, you can use them
    in the Calculation Setup or Workflow Builder pages.
    """)
    
    st.info("""
    💡 **How Pseudopotentials Work:**
    - **Download** pseudopotentials to your local computer (e.g., from SSSP, PSLibrary, Pseudo Dojo)
    - **Configure** them here by pointing to the directory containing .UPF files
    - **Machine selection is optional** - only needed if pseudopotentials are stored on a remote machine
    - During calculations, xespresso automatically copies pseudopotentials to the remote machine
    - Configurations are saved to `~/.xespresso/pseudopotentials/`
    """)
    
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available. Cannot configure pseudopotentials.")
        return
    
    st.subheader("Pseudopotentials Configuration")
    
    # Option to specify if pseudopotentials are on a remote machine
    use_remote = st.checkbox(
        "Pseudopotentials are stored on a remote machine",
        value=False,
        help="Check this only if pseudopotentials are on a remote machine, not your local computer"
    )
    
    # Machine selection - only if pseudopotentials are on remote machine
    selected_machine = None
    if use_remote:
        try:
            machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            if machines_list:
                selected_machine = st.selectbox(
                    "Select Remote Machine:",
                    machines_list,
                    help="Choose the remote machine where pseudopotentials are stored"
                )
            else:
                st.warning("⚠️ No machines configured. Please configure a machine first, or uncheck the remote option to use local pseudopotentials.")
                selected_machine = None
        except Exception as e:
            st.warning(f"Could not load machines: {e}")
            selected_machine = None
    
    # Only proceed if we don't need a machine or if a machine is selected
    if not use_remote or selected_machine:
        st.subheader(f"Pseudopotentials Configuration{' for: ' + selected_machine if selected_machine else ''}")
        
        # Auto-detection section
        st.subheader("Auto-Detect Pseudopotentials")
        
        with st.form("detect_pseudos_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                config_name = st.text_input(
                    "Configuration Name *",
                    placeholder="e.g., SSSP_efficiency, pbe_standard, my_pseudos",
                    help="Unique name for this pseudopotential configuration"
                )
                base_path = st.text_input(
                    "Pseudopotentials Directory *",
                    placeholder="e.g., /home/user/pseudopotentials/SSSP, ~/pseudo/pbe" if not use_remote else "e.g., /opt/pseudopotentials/SSSP",
                    help="Directory containing .UPF pseudopotential files (local path or remote path)"
                )
                functional = st.text_input(
                    "Functional (optional)",
                    placeholder="e.g., PBE, LDA, PBEsol",
                    help="Exchange-correlation functional"
                )
            
            with col2:
                library = st.text_input(
                    "Library Name (optional)",
                    placeholder="e.g., SSSP, PSLibrary, Pseudo Dojo",
                    help="Name of the pseudopotential library"
                )
                library_version = st.text_input(
                    "Library Version (optional)",
                    placeholder="e.g., 1.1.2, 1.0.0, 0.4",
                    help="Version of the pseudopotential library"
                )
                description = st.text_area(
                    "Description (optional)",
                    placeholder="e.g., SSSP 1.1.2 PBE efficiency set for production calculations",
                    help="Brief description of this pseudopotential set"
                )
                recursive = st.checkbox(
                    "Search subdirectories recursively",
                    value=True,
                    help="Search for .UPF files in subdirectories"
                )
            
            st.info("""
            **💡 Tips:**
            - The directory should contain .UPF or .upf pseudopotential files
            - Files are typically named like `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF`
            - Auto-detection will extract element symbols from filenames
            - **Local pseudopotentials** are recommended (stored on your computer)
            - For remote pseudopotentials, ensure the path exists on the remote system
            """)
            
            detect_button = st.form_submit_button("🔍 Auto-Detect Pseudopotentials")
        
        if detect_button:
            if not config_name or not base_path:
                st.error("❌ Configuration name and base path are required!")
            else:
                with st.spinner("Detecting pseudopotential files..."):
                    try:
                        # Determine if we need SSH connection
                        ssh_conn = None
                        if use_remote and selected_machine:
                            # Load machine config to get SSH details
                            try:
                                from xespresso.machines.config.loader import load_machine
                                from xespresso.machines.machine import Machine
                                
                                machine_obj = load_machine(
                                    DEFAULT_CONFIG_PATH,
                                    selected_machine,
                                    DEFAULT_MACHINES_DIR,
                                    return_object=True
                                )
                                
                                if machine_obj and isinstance(machine_obj, Machine) and machine_obj.is_remote:
                                    ssh_conn = {
                                        'host': machine_obj.host,
                                        'username': machine_obj.username,
                                        'port': machine_obj.port if hasattr(machine_obj, 'port') else 22
                                    }
                            except Exception as e:
                                st.warning(f"Could not load machine SSH details: {e}")
                        
                        pseudo_config = create_pseudopotentials_config(
                            name=config_name.strip(),
                            base_path=base_path.strip(),
                            machine_name=selected_machine if use_remote else None,
                            ssh_connection=ssh_conn,
                            recursive=recursive,
                            description=description.strip() if description and description.strip() else None,
                            functional=functional.strip() if functional and functional.strip() else None,
                            library=library.strip() if library and library.strip() else None,
                            version=library_version.strip() if library_version and library_version.strip() else None,
                            save=False  # We'll save after user confirmation
                        )
                        
                        if pseudo_config and pseudo_config.pseudopotentials:
                            st.success(f"✅ Detected {len(pseudo_config.pseudopotentials)} pseudopotentials!")
                            st.session_state.detected_pseudos = pseudo_config
                        else:
                            st.warning("⚠️ No pseudopotentials detected. Check the path and try again.")
                            st.session_state.detected_pseudos = None
                    except Exception as e:
                        st.error(f"❌ Error detecting pseudopotentials: {e}")
                        st.code(traceback.format_exc())
                        st.session_state.detected_pseudos = None
        
        # Display detected pseudopotentials
        if hasattr(st.session_state, 'detected_pseudos') and st.session_state.detected_pseudos:
            pseudo_config = st.session_state.detected_pseudos
            
            st.subheader("Detected Pseudopotentials")
            
            # Display configuration details
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Configuration:** {pseudo_config.name}")
                if pseudo_config.machine_name:
                    st.markdown(f"**Remote Machine:** {pseudo_config.machine_name}")
                else:
                    st.markdown(f"**Location:** Local (your computer)")
                st.markdown(f"**Base Path:** {pseudo_config.base_path}")
            with col2:
                if pseudo_config.functional:
                    st.markdown(f"**Functional:** {pseudo_config.functional}")
                if pseudo_config.library:
                    lib_str = f"**Library:** {pseudo_config.library}"
                    if pseudo_config.version:
                        lib_str += f" v{pseudo_config.version}"
                    st.markdown(lib_str)
                if pseudo_config.description:
                    st.markdown(f"**Description:** {pseudo_config.description}")
            
            # Display pseudopotentials table
            st.markdown("### Available Elements")
            pseudo_data = []
            for element in sorted(pseudo_config.pseudopotentials.keys()):
                pseudo = pseudo_config.get_pseudopotential(element)
                pseudo_data.append({
                    "Element": element,
                    "Filename": pseudo.filename,
                    "Type": pseudo.type or "Unknown",
                    "Functional": pseudo.functional or "Unknown",
                })
            
            st.dataframe(pseudo_data, use_container_width=True)
            
            # Save option
            st.info("""
            **💾 Saving Pseudopotentials:**
            - Configuration will be saved to `~/.xespresso/pseudopotentials/{name}.json`
            - Pseudopotential files remain in their original location (local or remote)
            - During calculations, xespresso automatically copies needed files to remote machines
            - Multiple configurations can be saved
            """)
            
            if st.button("💾 Save Pseudopotentials Configuration"):
                try:
                    filepath = PseudopotentialsManager.save_config(
                        pseudo_config,
                        output_dir=DEFAULT_PSEUDOPOTENTIALS_DIR,
                        overwrite=True  # Allow overwrite since user explicitly clicked save
                    )
                    location = f"on remote machine '{pseudo_config.machine_name}'" if pseudo_config.machine_name else "locally"
                    st.success(f"✅ Pseudopotentials saved to: {filepath}")
                    st.info(f"📦 Saved configuration: **{pseudo_config.name}** ({location}) with {len(pseudo_config.pseudopotentials)} elements")
                    
                    # Clear the detected pseudos after successful save
                    st.session_state.detected_pseudos = None
                    st.session_state.current_pseudos = pseudo_config
                except Exception as e:
                    st.error(f"Error saving pseudopotentials: {e}")
                    st.code(traceback.format_exc())
        
        # Load existing configurations
        st.subheader("Existing Pseudopotential Configurations")
        
        try:
            existing_configs = PseudopotentialsManager.list_configs(DEFAULT_PSEUDOPOTENTIALS_DIR)
            
            if existing_configs:
                st.success(f"✅ Found {len(existing_configs)} saved configuration(s)")
                
                # Configuration selector
                selected_config = st.selectbox(
                    "Select Configuration to View:",
                    existing_configs,
                    help="Choose a pseudopotential configuration to view or use"
                )
                
                if st.button("Load Configuration"):
                    with st.spinner(f"Loading configuration '{selected_config}'..."):
                        try:
                            loaded_config = load_pseudopotentials_config(
                                selected_config,
                                DEFAULT_PSEUDOPOTENTIALS_DIR
                            )
                            
                            if loaded_config:
                                st.success(f"✅ Loaded configuration: {selected_config}")
                                st.session_state.current_pseudos = loaded_config
                                
                                # Display configuration details
                                st.markdown("### Configuration Details")
                                col1, col2 = st.columns(2)
                                with col1:
                                    st.markdown(f"**Name:** {loaded_config.name}")
                                    if loaded_config.machine_name:
                                        st.markdown(f"**Remote Machine:** {loaded_config.machine_name}")
                                    else:
                                        st.markdown(f"**Location:** Local (your computer)")
                                    st.markdown(f"**Base Path:** {loaded_config.base_path}")
                                with col2:
                                    if loaded_config.functional:
                                        st.markdown(f"**Functional:** {loaded_config.functional}")
                                    if loaded_config.library:
                                        lib_str = f"**Library:** {loaded_config.library}"
                                        if loaded_config.version:
                                            lib_str += f" v{loaded_config.version}"
                                        st.markdown(lib_str)
                                    st.markdown(f"**Elements:** {len(loaded_config.pseudopotentials)}")
                                
                                if loaded_config.description:
                                    st.markdown(f"**Description:** {loaded_config.description}")
                                
                                # Display pseudopotentials
                                st.markdown("### Available Pseudopotentials")
                                elements = loaded_config.list_elements()
                                st.markdown(f"**{len(elements)} elements:** {', '.join(elements)}")
                                
                                # Show detailed table in expander
                                with st.expander("📋 View Detailed Pseudopotentials Table"):
                                    pseudo_table_data = []
                                    for element in elements:
                                        pseudo = loaded_config.get_pseudopotential(element)
                                        pseudo_table_data.append({
                                            "Element": element,
                                            "Filename": pseudo.filename,
                                            "Type": pseudo.type or "-",
                                            "Functional": pseudo.functional or "-",
                                            "Path": pseudo.path,
                                        })
                                    st.dataframe(pseudo_table_data, use_container_width=True)
                                
                                # Option to set as default or delete
                                st.markdown("---")
                                
                                # Check if this is already the default
                                is_default = selected_config == "default"
                                has_default = PseudopotentialsManager.has_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
                                
                                col1, col2 = st.columns(2)
                                
                                with col1:
                                    # Set as default button (only for non-default configs)
                                    if not is_default:
                                        if st.button("⭐ Set as Default", help="Set this configuration as the default for calculations"):
                                            try:
                                                PseudopotentialsManager.set_default_config(
                                                    selected_config,
                                                    DEFAULT_PSEUDOPOTENTIALS_DIR
                                                )
                                                st.success(f"✅ Set '{selected_config}' as default configuration")
                                                st.info("💡 This configuration will be automatically selected in Calculation Setup and Workflow Builder")
                                                st.rerun()
                                            except Exception as e:
                                                st.error(f"Error setting default: {e}")
                                    else:
                                        st.info("ℹ️ This is the default configuration")
                                
                                with col2:
                                    # Clear default button (only if a default exists and viewing a non-default)
                                    if has_default and not is_default:
                                        if st.button("🚫 Clear Default", help="Remove the current default configuration"):
                                            try:
                                                PseudopotentialsManager.clear_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
                                                st.success("✅ Default configuration cleared")
                                                st.rerun()
                                            except Exception as e:
                                                st.error(f"Error clearing default: {e}")
                                
                                # Delete button (cannot delete 'default' directly, must clear it first)
                                st.markdown("---")
                                if not is_default:
                                    if st.button("🗑️ Delete This Configuration", type="secondary"):
                                        try:
                                            PseudopotentialsManager.delete_config(
                                                selected_config,
                                                DEFAULT_PSEUDOPOTENTIALS_DIR
                                            )
                                            st.success(f"✅ Deleted configuration: {selected_config}")
                                            if hasattr(st.session_state, 'current_pseudos'):
                                                del st.session_state.current_pseudos
                                            st.rerun()
                                        except Exception as e:
                                            st.error(f"Error deleting configuration: {e}")
                                else:
                                    st.warning("⚠️ Cannot delete the default configuration. Clear it first using the button above.")
                            else:
                                st.error(f"Failed to load configuration: {selected_config}")
                        except Exception as e:
                            st.error(f"Error loading configuration: {e}")
                            st.code(traceback.format_exc())
            else:
                st.info("No pseudopotential configurations saved yet. Create one above!")
        
        except Exception as e:
            st.warning(f"Could not load existing configurations: {e}")

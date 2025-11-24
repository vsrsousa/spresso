"""
Codes Configuration Page for xespresso GUI.

This module handles the Quantum ESPRESSO codes configuration interface,
allowing users to:
- Auto-detect QE executables on machines
- Configure multiple versions
- Save and load code configurations
"""

import streamlit as st
import traceback
import os

try:
    from xespresso.machines.config.loader import (
        list_machines,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    from xespresso.codes.manager import (
        detect_qe_codes, load_codes_config, CodesManager,
        DEFAULT_CODES_DIR
    )
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False


def render_codes_config_page():
    """Render the codes configuration page."""
    st.header("⚙️ Quantum ESPRESSO Codes Configuration")
    st.markdown("""
    **Configure and auto-detect** Quantum ESPRESSO executable paths.
    
    This page is for **configuration only** - once codes are saved, you can select versions
    in the Calculation Setup or Workflow Builder pages.
    
    Auto-detection is supported for both local and remote systems.
    """)
    
    st.info("""
    💡 **Configuration vs. Selection:**
    - **Configure** codes here (auto-detect and save)
    - **Select** code versions in Calculation Setup or Workflow Builder
    - Multiple QE versions can coexist for the same machine
    - Configurations are saved to `~/.xespresso/codes/`
    """)
    
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available. Cannot configure codes.")
        return
    
    # Machine selection
    try:
        machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
        if machines_list:
            selected_machine = st.selectbox(
                "Select Machine to Configure:",
                machines_list,
                help="Choose the machine to configure QE codes for"
            )
        else:
            st.warning("⚠️ No machines configured. Please configure a machine first.")
            selected_machine = None
    except Exception as e:
        st.warning(f"Could not load machines: {e}")
        selected_machine = None
    
    if selected_machine:
        st.subheader(f"Codes Configuration for: {selected_machine}")
        
        # Auto-detection section
        st.subheader("Auto-Detect Codes")
        
        with st.form("detect_codes_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                qe_prefix = st.text_input(
                    "QE Installation Prefix (optional)",
                    help="e.g., /opt/qe-7.2/bin"
                )
                # Feature 1: Explicit QE Version Specification
                qe_version = st.text_input(
                    "QE Version (optional but recommended)",
                    placeholder="e.g., 7.2, 7.1, 6.8",
                    help="Explicitly specify Quantum ESPRESSO version to avoid confusion with compiler versions"
                )
                label = st.text_input(
                    "Label (optional)",
                    help="Custom label for this version (e.g., 'production', 'dev', 'test')"
                )
                modules_str = st.text_area(
                    "Modules to Load (optional, one per line)",
                    help="Version-specific modules (e.g., 'qe/7.2' or 'quantum_espresso-7.4.1')"
                )
            
            with col2:
                search_paths_str = st.text_area(
                    "Additional Search Paths (optional, one per line)",
                    help="Additional directories to search for executables"
                )
                st.info("""
                **💡 Tip: Explicit Version**
                
                Auto-detection may pick up compiler versions.
                It's recommended to specify the QE version explicitly!
                """)
            
            detect_button = st.form_submit_button("🔍 Auto-Detect Codes")
        
        if detect_button:
            with st.spinner("Detecting Quantum ESPRESSO codes..."):
                try:
                    modules = [m.strip() for m in modules_str.split("\n") if m.strip()] if modules_str else None
                    search_paths = [p.strip() for p in search_paths_str.split("\n") if p.strip()] if search_paths_str else None
                    
                    codes_config = detect_qe_codes(
                        machine_name=selected_machine,
                        qe_prefix=qe_prefix if qe_prefix else None,
                        search_paths=search_paths,
                        modules=modules,
                        auto_load_machine=True,
                        qe_version=qe_version.strip() if qe_version and qe_version.strip() else None,
                        label=label.strip() if label and label.strip() else None
                    )
                    
                    # Check if any codes were detected (in main codes dict or versions structure)
                    if codes_config and codes_config.has_any_codes():
                        # Get all detected codes to count them
                        all_codes = codes_config.get_all_codes()
                        st.success(f"✅ Detected {len(all_codes)} codes!")
                        # Store in session state with a flag to show detected codes were just found
                        st.session_state.detected_codes = codes_config
                        st.session_state.detected_machine = selected_machine
                    else:
                        st.warning("⚠️ No codes detected. Check paths and modules.")
                        st.session_state.detected_codes = None
                except Exception as e:
                    st.error(f"❌ Error detecting codes: {e}")
                    st.code(traceback.format_exc())
                    st.session_state.detected_codes = None
        
        # Display detected codes (outside the detect_button block so it persists across reruns)
        if hasattr(st.session_state, 'detected_codes') and st.session_state.detected_codes:
            codes_config = st.session_state.detected_codes
            
            st.subheader("Detected Codes")
            
            # Get all codes from the appropriate location (main codes or version-specific)
            all_codes = codes_config.get_all_codes()
            
            codes_data = []
            for name, code in all_codes.items():
                codes_data.append({
                    "Code": name,
                    "Path": code.path,
                    "Version": code.version or "Unknown",
                    "Label": codes_config.label or "default"
                })
            st.table(codes_data)
            
            # Save option with clear explanation
            st.info("""
            **💾 Saving Codes:**
            - Configuration will be saved to `{machine}.json`
            - If **version/label** is specified, codes are stored in the `versions` structure within the file
            - With merge enabled, new versions are added to existing configurations without overwriting
            - Multiple QE versions with different labels can coexist in the same file
            - This prevents accidentally losing other version configurations
            """)
            
            if st.button("💾 Save Codes Configuration"):
                try:
                    filepath = CodesManager.save_config(
                        codes_config,
                        output_dir=DEFAULT_CODES_DIR,
                        overwrite=False,
                        merge=True,
                        interactive=False
                    )
                    st.success(f"✅ Codes saved to: {filepath}")
                    
                    # Show what was saved
                    if codes_config.qe_version:
                        st.info(f"📦 Added/updated QE version: **{codes_config.qe_version}**" + 
                               (f" with label **{codes_config.label}**" if codes_config.label else ""))
                    else:
                        st.info("💡 Configuration saved to default (no version specified)")
                    
                    # Clear the detected codes after successful save
                    st.session_state.detected_codes = None
                    st.session_state.current_codes = codes_config
                except Exception as e:
                    st.error(f"Error saving codes: {e}")
                    st.code(traceback.format_exc())
        
        # Load existing configuration
        st.subheader("Existing Codes Configuration")
        
        try:
            existing_codes = load_codes_config(selected_machine, DEFAULT_CODES_DIR)
            if existing_codes:
                st.success(f"✅ Loaded existing configuration for '{selected_machine}'")
                
                # Feature 2: Version Selection - Show available versions
                if existing_codes.versions:
                    available_versions = existing_codes.list_versions()
                    st.info(f"📦 Available QE versions in this configuration: {', '.join(available_versions)}")
                    
                    # Show labels for each version if available
                    version_labels = {}
                    for version in available_versions:
                        if version in existing_codes.versions and 'label' in existing_codes.versions[version]:
                            version_labels[version] = existing_codes.versions[version]['label']
                    
                    if version_labels:
                        st.markdown("**Version Labels:**")
                        for version, label in version_labels.items():
                            st.markdown(f"- Version {version}: `{label}`")
                    
                    # Version selector
                    st.subheader("Select QE Version for Calculations")
                    selected_version = st.selectbox(
                        "Choose QE Version:",
                        available_versions,
                        help="Select which Quantum ESPRESSO version to use for your calculations"
                    )
                    
                    # Load codes for the selected version
                    if st.button(f"Load QE {selected_version} Configuration"):
                        with st.spinner(f"Loading QE {selected_version}..."):
                            try:
                                version_config = load_codes_config(
                                    selected_machine, 
                                    DEFAULT_CODES_DIR, 
                                    version=selected_version
                                )
                                
                                if version_config:
                                    st.success(f"✅ Loaded QE {selected_version} configuration!")
                                    st.session_state.current_codes = version_config
                                    st.session_state.selected_qe_version = selected_version
                                    
                                    # Display codes for this version
                                    st.markdown(f"**Codes for QE {selected_version}:**")
                                    
                                    # Get codes from the version-specific structure
                                    version_codes = version_config.get_all_codes(version=selected_version)
                                    version_codes_data = []
                                    for name, code in version_codes.items():
                                        version_codes_data.append({
                                            "Code": name,
                                            "Path": code.path,
                                            "Version": code.version or selected_version,
                                        })
                                    st.table(version_codes_data)
                                    
                                    # Show version-specific modules if available
                                    # Check both top-level (backward compat) and version structure
                                    modules = None
                                    if version_config.versions and selected_version in version_config.versions:
                                        modules = version_config.versions[selected_version].get('modules')
                                    elif hasattr(version_config, 'modules') and version_config.modules:
                                        modules = version_config.modules
                                    
                                    if modules:
                                        st.markdown(f"**Modules for QE {selected_version}:**")
                                        for module in modules:
                                            st.markdown(f"- `{module}`")
                                else:
                                    st.warning(f"⚠️ Could not load QE {selected_version} configuration.")
                            except Exception as e:
                                st.error(f"❌ Error loading version: {e}")
                                st.code(traceback.format_exc())
                else:
                    # No version-specific configuration, show all codes
                    codes_data = []
                    for name, code in existing_codes.codes.items():
                        codes_data.append({
                            "Code": name,
                            "Path": code.path,
                            "Version": code.version or "Unknown"
                        })
                    st.table(codes_data)
                    
                    st.session_state.current_codes = existing_codes
                    
                    # Code/version selection for calculations
                    st.subheader("Select Code Version for Calculations")
                    if existing_codes.codes:
                        code_options = list(existing_codes.codes.keys())
                        selected_code = st.selectbox(
                            "Select QE code to use:",
                            code_options,
                            help="Choose which QE code to use for your calculations"
                        )
                        st.session_state.selected_code_version = selected_code
                        
                        selected_code_obj = existing_codes.codes[selected_code]
                        st.info(f"""
                        **Selected Code Details:**
                        - Path: `{selected_code_obj.path}`
                        - Version: {selected_code_obj.version or 'Unknown'}
                        """)
            else:
                st.info("ℹ️ No codes configuration found for this machine.")
        except Exception as e:
            st.warning(f"Could not load codes configuration: {e}")

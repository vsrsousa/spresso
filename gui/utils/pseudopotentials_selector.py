"""
Pseudopotentials selector utility for GUI pages.

This module provides a reusable component for selecting pseudopotentials
from configured libraries.
"""

import streamlit as st
from typing import Optional, Dict


def render_pseudopotentials_selector(
    elements: set,
    config_dict: dict,
    key_prefix: str = "calc"
) -> None:
    """
    Render a pseudopotentials selector component.
    
    This component allows users to select a pseudopotential configuration
    and automatically populates pseudopotentials for the required elements.
    The selection persists across the interface via session state.
    
    Args:
        elements: Set of element symbols needed for the calculation
        config_dict: Dictionary to store the selected pseudopotentials
        key_prefix: Prefix for widget keys to avoid conflicts
    """
    st.subheader("🧪 Pseudopotentials")
    
    try:
        from xespresso.pseudopotentials import (
            PseudopotentialsManager,
            DEFAULT_PSEUDOPOTENTIALS_DIR
        )
        
        # List available configurations
        available_configs = PseudopotentialsManager.list_configs(DEFAULT_PSEUDOPOTENTIALS_DIR)
        
        if not available_configs:
            st.warning("""
            ⚠️ No pseudopotential configurations found.
            
            Please configure pseudopotentials first:
            1. Enable "Show Configuration" in the sidebar
            2. Go to "🧪 Pseudopotentials Configuration"
            3. Auto-detect and save your pseudopotential library
            """)
            
            # Fallback to manual entry
            st.markdown("**Manual Entry:**")
            if "pseudopotentials" not in config_dict:
                config_dict["pseudopotentials"] = {}
            
            for elem in sorted(elements):
                pseudo = st.text_input(
                    f"Pseudopotential for {elem}:",
                    value=config_dict["pseudopotentials"].get(elem, f"{elem}.UPF"),
                    key=f"{key_prefix}_pseudo_manual_{elem}",
                )
                config_dict["pseudopotentials"][elem] = pseudo
            return
        
        # Initialize session state for selected config
        if f"{key_prefix}_selected_pseudo_config" not in st.session_state:
            # Check if there's a default configuration
            if PseudopotentialsManager.has_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR):
                default_config = PseudopotentialsManager.get_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
                if default_config:
                    st.session_state[f"{key_prefix}_selected_pseudo_config"] = "default"
                    st.info("ℹ️ Using default pseudopotentials configuration")
                else:
                    st.session_state[f"{key_prefix}_selected_pseudo_config"] = None
            else:
                st.session_state[f"{key_prefix}_selected_pseudo_config"] = None
        
        # Load all configurations to show metadata (excluding default from list)
        configs_info = []
        for config_name in available_configs:
            # Skip 'default' in the main list - it will be handled separately
            if config_name == 'default':
                continue
                
            try:
                pseudo_config = PseudopotentialsManager.load_config(
                    config_name,
                    DEFAULT_PSEUDOPOTENTIALS_DIR
                )
                if pseudo_config:
                    configs_info.append({
                        'name': config_name,
                        'config': pseudo_config,
                        'library': pseudo_config.library or 'Unknown',
                        'version': pseudo_config.version or '',
                        'functional': pseudo_config.functional or 'Unknown',
                        'elements': len(pseudo_config.pseudopotentials),
                        'location': 'Local' if not pseudo_config.machine_name else pseudo_config.machine_name
                    })
            except Exception as e:
                st.warning(f"Could not load config '{config_name}': {e}")
        
        if not configs_info:
            st.warning("⚠️ No valid pseudopotential configurations found.")
            return
        
        # Filter options
        col1, col2 = st.columns(2)
        
        with col1:
            # Filter by library
            libraries = sorted(set(c['library'] for c in configs_info))
            selected_library = st.selectbox(
                "Filter by Library:",
                ["All"] + libraries,
                key=f"{key_prefix}_filter_library",
                help="Filter configurations by pseudopotential library"
            )
        
        with col2:
            # Filter by functional
            functionals = sorted(set(c['functional'] for c in configs_info))
            selected_functional = st.selectbox(
                "Filter by Functional:",
                ["All"] + functionals,
                key=f"{key_prefix}_filter_functional",
                help="Filter configurations by exchange-correlation functional"
            )
        
        # Apply filters
        filtered_configs = configs_info
        if selected_library != "All":
            filtered_configs = [c for c in filtered_configs if c['library'] == selected_library]
        if selected_functional != "All":
            filtered_configs = [c for c in filtered_configs if c['functional'] == selected_functional]
        
        # Check if default configuration exists and add it as first option
        has_default = PseudopotentialsManager.has_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
        default_config_info = None
        
        if has_default:
            try:
                default_config = PseudopotentialsManager.get_default_config(DEFAULT_PSEUDOPOTENTIALS_DIR)
                if default_config:
                    default_config_info = {
                        'name': 'default',
                        'config': default_config,
                        'library': default_config.library or 'Unknown',
                        'version': default_config.version or '',
                        'functional': default_config.functional or 'Unknown',
                        'elements': len(default_config.pseudopotentials),
                        'location': 'Local' if not default_config.machine_name else default_config.machine_name
                    }
            except Exception as e:
                st.warning(f"Could not load default config: {e}")
        
        if not filtered_configs and not default_config_info:
            st.warning("⚠️ No configurations match the selected filters.")
            return
        
        # Build configuration options list
        # Add default as first option if it exists
        all_configs_for_selection = []
        config_options = []
        
        if default_config_info:
            all_configs_for_selection.append(default_config_info)
            config_options.append(
                f"⭐ Default ({default_config_info['library']} {default_config_info['version']}, "
                f"{default_config_info['functional']}, {default_config_info['elements']} elements)"
            )
        
        # Add filtered configs
        all_configs_for_selection.extend(filtered_configs)
        config_options.extend([
            f"{c['name']} ({c['library']} {c['version']}, {c['functional']}, {c['elements']} elements)"
            for c in filtered_configs
        ])
        
        # Find current selection index
        current_selection = st.session_state.get(f"{key_prefix}_selected_pseudo_config")
        default_idx = 0
        if current_selection:
            try:
                default_idx = [c['name'] for c in all_configs_for_selection].index(current_selection)
            except ValueError:
                # If previously selected config not in filtered list, default to first
                default_idx = 0
        
        selected_idx = st.selectbox(
            "Select Pseudopotential Configuration:",
            range(len(all_configs_for_selection)),
            format_func=lambda i: config_options[i],
            index=default_idx,
            key=f"{key_prefix}_pseudo_selector",
            help="Choose a pseudopotential configuration for this calculation"
        )
        
        selected_config_info = all_configs_for_selection[selected_idx]
        selected_config_name = selected_config_info['name']
        
        # Update session state
        st.session_state[f"{key_prefix}_selected_pseudo_config"] = selected_config_name
        
        # Load the selected configuration
        pseudo_config = selected_config_info['config']
        
        # Display configuration details
        with st.expander("📋 Configuration Details", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Name:** {pseudo_config.name}")
                st.markdown(f"**Library:** {pseudo_config.library or 'Not specified'}")
                if pseudo_config.version:
                    st.markdown(f"**Version:** {pseudo_config.version}")
                st.markdown(f"**Functional:** {pseudo_config.functional or 'Not specified'}")
            with col2:
                st.markdown(f"**Base Path:** {pseudo_config.base_path}")
                st.markdown(f"**Location:** {selected_config_info['location']}")
                st.markdown(f"**Elements:** {len(pseudo_config.pseudopotentials)}")
                if pseudo_config.description:
                    st.markdown(f"**Description:** {pseudo_config.description}")
        
        # Get pseudopotentials dictionary
        pseudo_dict = pseudo_config.get_pseudopotentials_dict()
        
        # Check if all required elements are available
        missing_elements = [elem for elem in elements if elem not in pseudo_dict]
        
        if missing_elements:
            st.error(f"""
            ❌ Missing pseudopotentials for: {', '.join(missing_elements)}
            
            The selected configuration does not contain pseudopotentials for all elements in your structure.
            Please select a different configuration or add the missing pseudopotentials manually.
            """)
            
            # Allow manual override for missing elements
            st.markdown("**Add Missing Pseudopotentials:**")
            if "pseudopotentials" not in config_dict:
                config_dict["pseudopotentials"] = {}
            
            for elem in missing_elements:
                pseudo = st.text_input(
                    f"Pseudopotential for {elem}:",
                    value=config_dict["pseudopotentials"].get(elem, f"{elem}.UPF"),
                    key=f"{key_prefix}_pseudo_missing_{elem}",
                )
                config_dict["pseudopotentials"][elem] = pseudo
            
            # Copy available elements from config
            for elem in elements:
                if elem in pseudo_dict:
                    config_dict["pseudopotentials"][elem] = pseudo_dict[elem]
        else:
            # All elements available - populate config
            st.success(f"✅ All required elements available in this configuration")
            
            if "pseudopotentials" not in config_dict:
                config_dict["pseudopotentials"] = {}
            
            # Show which pseudopotentials will be used
            st.markdown("**Selected Pseudopotentials:**")
            pseudo_table = []
            for elem in sorted(elements):
                pseudo_table.append({
                    "Element": elem,
                    "File": pseudo_dict[elem]
                })
                config_dict["pseudopotentials"][elem] = pseudo_dict[elem]
            
            st.table(pseudo_table)
        
        # Set ESPRESSO_PSEUDO environment variable to the config's base_path
        # This allows xespresso to find the pseudopotential files without
        # explicitly setting pseudo_dir in the calculator parameters
        # This respects xespresso's original implementation for local execution
        # and works correctly with remote execution (remote_mixin handles the transfer)
        import os
        os.environ["ESPRESSO_PSEUDO"] = pseudo_config.base_path
        
        # Also store in session state for tracking and persistence
        st.session_state["_espresso_pseudo_path"] = pseudo_config.base_path
        
        # Show current ESPRESSO_PSEUDO setting
        with st.expander("ℹ️ Environment Variable Info", expanded=False):
            st.code(f"ESPRESSO_PSEUDO = {os.environ.get('ESPRESSO_PSEUDO', 'Not set')}", language="bash")
            st.caption("This environment variable is automatically updated when you change pseudopotential configurations.")
        
        # Store the selected config name for persistence
        config_dict["_selected_pseudo_config_name"] = selected_config_name
        
    except ImportError as e:
        st.error(f"❌ Could not import pseudopotentials module: {e}")
        st.info("Falling back to manual pseudopotential entry")
        
        # Fallback to manual entry
        if "pseudopotentials" not in config_dict:
            config_dict["pseudopotentials"] = {}
        
        for elem in sorted(elements):
            pseudo = st.text_input(
                f"Pseudopotential for {elem}:",
                value=config_dict["pseudopotentials"].get(elem, f"{elem}.UPF"),
                key=f"{key_prefix}_pseudo_fallback_{elem}",
            )
            config_dict["pseudopotentials"][elem] = pseudo
    
    except Exception as e:
        st.error(f"❌ Error loading pseudopotentials: {e}")
        import traceback
        st.code(traceback.format_exc())

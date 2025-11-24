"""
Magnetic configuration selector utility for GUI pages.

This module provides a reusable component for configuring magnetic properties
using xespresso's setup_magnetic_config functionality.
"""

import streamlit as st
from typing import Dict, Set, Optional


# Predefined magnetic configurations for common elements
PREDEFINED_MAGNETIC_MOMENTS = {
    # 3d transition metals
    'Fe': {'ferromagnetic': 2.2, 'antiferromagnetic': [2.2, -2.2]},
    'Co': {'ferromagnetic': 1.7, 'antiferromagnetic': [1.7, -1.7]},
    'Ni': {'ferromagnetic': 0.6, 'antiferromagnetic': [0.6, -0.6]},
    'Mn': {'ferromagnetic': 5.0, 'antiferromagnetic': [5.0, -5.0]},
    'Cr': {'ferromagnetic': 3.0, 'antiferromagnetic': [3.0, -3.0]},
    'V': {'ferromagnetic': 2.0, 'antiferromagnetic': [2.0, -2.0]},
    'Ti': {'ferromagnetic': 1.0, 'antiferromagnetic': [1.0, -1.0]},
    # 4f rare earths
    'Gd': {'ferromagnetic': 7.0, 'antiferromagnetic': [7.0, -7.0]},
    'Nd': {'ferromagnetic': 3.0, 'antiferromagnetic': [3.0, -3.0]},
    'Sm': {'ferromagnetic': 5.0, 'antiferromagnetic': [5.0, -5.0]},
    'Eu': {'ferromagnetic': 7.0, 'antiferromagnetic': [7.0, -7.0]},
    'Tb': {'ferromagnetic': 9.0, 'antiferromagnetic': [9.0, -9.0]},
    'Dy': {'ferromagnetic': 10.0, 'antiferromagnetic': [10.0, -10.0]},
    'Ho': {'ferromagnetic': 10.0, 'antiferromagnetic': [10.0, -10.0]},
    'Er': {'ferromagnetic': 9.0, 'antiferromagnetic': [9.0, -9.0]},
}


def render_magnetic_selector(
    elements: Set[str], config_dict: dict, key_prefix: str = "calc"
) -> None:
    """
    Render a magnetic configuration selector component.

    This component allows users to configure magnetic properties per element
    following xespresso's setup_magnetic_config pattern.

    Args:
        elements: Set of element symbols in the structure
        config_dict: Dictionary to store the magnetic configuration
        key_prefix: Prefix for widget keys to avoid conflicts
    """
    st.subheader("🧲 Magnetic Configuration")

    # Initialize magnetic config in config_dict if not present
    if "enable_magnetism" not in config_dict:
        config_dict["enable_magnetism"] = False
    if "magnetic_config" not in config_dict:
        config_dict["magnetic_config"] = {}

    # Checkbox to enable magnetic configuration
    enable_magnetism = st.checkbox(
        "Enable Magnetic Configuration",
        value=config_dict.get("enable_magnetism", False),
        help="Enable to configure magnetic moments per element. Required for magnetic systems.",
        key=f"{key_prefix}_enable_magnetism",
    )
    config_dict["enable_magnetism"] = enable_magnetism

    if enable_magnetism:
        st.info(
            """
        💡 **How it works:**
        - Specify magnetic moment(s) for each element
        - Multiple values create non-equivalent atoms (different species)
        - Example: Fe = [2.2, -2.2] creates antiferromagnetic configuration
        - Example: Fe = [2.2] makes all Fe atoms equivalent with moment 2.2
        """
        )

        # Show predefined configuration selector
        st.markdown("**Quick Setup - Predefined Configurations:**")
        col1, col2 = st.columns(2)
        with col1:
            preset_type = st.selectbox(
                "Configuration Type:",
                ["Custom", "Ferromagnetic", "Antiferromagnetic"],
                key=f"{key_prefix}_preset_type",
                help="Select a predefined magnetic configuration"
            )
        with col2:
            if preset_type != "Custom" and st.button("Apply Preset", key=f"{key_prefix}_apply_preset"):
                for element in sorted(elements):
                    if element in PREDEFINED_MAGNETIC_MOMENTS:
                        presets = PREDEFINED_MAGNETIC_MOMENTS[element]
                        if preset_type == "Ferromagnetic" and 'ferromagnetic' in presets:
                            config_dict["magnetic_config"][element] = [presets['ferromagnetic']]
                        elif preset_type == "Antiferromagnetic" and 'antiferromagnetic' in presets:
                            config_dict["magnetic_config"][element] = presets['antiferromagnetic']
                    else:
                        # Default non-magnetic for elements without presets
                        config_dict["magnetic_config"][element] = [0]
                st.success(f"Applied {preset_type} preset!")
                st.rerun()

        st.markdown("---")
        st.markdown("**Configure Magnetic Moments per Element:**")
        
        # Create table-like header
        header_cols = st.columns([2, 3, 2])
        with header_cols[0]:
            st.markdown("**Species**")
        with header_cols[1]:
            st.markdown("**Magnetic Moments**")
        with header_cols[2]:
            st.markdown("**Result**")
        
        st.markdown("---")

        # Configure each element in table format (always visible)
        for element in sorted(elements):
            # Get current config for this element
            current_config = config_dict["magnetic_config"].get(element, [0])

            # Extract magnetic moments (handle both simple list and old dict format)
            if isinstance(current_config, dict) and "mag" in current_config:
                mag_moments = current_config.get("mag", [0])
            elif isinstance(current_config, list):
                mag_moments = current_config
            else:
                mag_moments = [0]

            # Create row for this element
            cols = st.columns([2, 3, 2])
            
            with cols[0]:
                st.markdown(f"**{element}**")
                # Show suggested value if available
                if element in PREDEFINED_MAGNETIC_MOMENTS:
                    fm_val = PREDEFINED_MAGNETIC_MOMENTS[element].get('ferromagnetic', 0)
                    st.caption(f"Typical: {fm_val}")

            with cols[1]:
                # Magnetic moments input
                moments_str = st.text_input(
                    f"Moments for {element}",
                    value=", ".join(map(str, mag_moments)),
                    help=f"Example: '2.2' for FM, '2.2, -2.2' for AFM, '0' for non-magnetic",
                    key=f"{key_prefix}_mag_moments_{element}",
                    label_visibility="collapsed"
                )

                # Parse moments
                try:
                    moments = [
                        float(m.strip()) for m in moments_str.split(",") if m.strip()
                    ]
                    if not moments:
                        moments = [0]
                except ValueError:
                    st.error(f"Invalid input for {element}")
                    moments = [0]

                # Store as simple list
                config_dict["magnetic_config"][element] = moments

            with cols[2]:
                # Show what will be created
                num_species = len(moments)
                if num_species == 1:
                    if moments[0] == 0:
                        st.caption("Non-magnetic")
                    else:
                        st.caption(f"1 species (m={moments[0]})")
                else:
                    species_names = [element] + [f"{element}{i}" for i in range(1, num_species)]
                    st.caption(f"{num_species} species")

        # Additional options
        st.markdown("---")
        st.markdown("**Additional Options:**")

        col1, col2 = st.columns(2)
        with col1:
            expand_cell = st.checkbox(
                "Auto-expand cell",
                value=config_dict.get("expand_cell", False),
                help="Automatically expand cell if needed to accommodate magnetic configuration",
                key=f"{key_prefix}_expand_cell",
            )
            config_dict["expand_cell"] = expand_cell

        with col2:
            # QE version for format selection
            current_qe_version = config_dict.get("qe_version", "auto")

            if current_qe_version is None or current_qe_version == "auto":
                display_value = "auto"
            elif isinstance(current_qe_version, str):
                try:
                    major_version = int(current_qe_version.split(".")[0])
                    if major_version >= 7:
                        display_value = "7.x"
                    else:
                        display_value = "6.x"
                except (ValueError, IndexError):
                    display_value = "auto"
            else:
                display_value = "auto"

            qe_version_display = st.selectbox(
                "QE Version Format:",
                ["auto", "6.x", "7.x"],
                index=["auto", "6.x", "7.x"].index(display_value),
                help="Select QE version format for compatibility",
                key=f"{key_prefix}_qe_version",
            )

            if qe_version_display == "auto":
                config_dict["qe_version"] = None
            elif current_qe_version and current_qe_version not in ["auto", "6.x", "7.x"]:
                config_dict["qe_version"] = current_qe_version
            else:
                config_dict["qe_version"] = qe_version_display

        # Show summary
        with st.expander("📋 Configuration Summary", expanded=False):
            st.json(config_dict["magnetic_config"])
    else:
        # Clear magnetic config if disabled
        config_dict["magnetic_config"] = {}

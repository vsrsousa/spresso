"""
Magnetic configuration selector utility for GUI pages.

This module provides a reusable component for configuring magnetic properties
using xespresso's setup_magnetic_config functionality.
"""

import streamlit as st
from typing import Dict, Set, Optional


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
        st.markdown("**Configure Magnetic Moments:**")
        st.info(
            """
        💡 **How it works:**
        - Specify magnetic moment(s) for each element
        - Multiple values create non-equivalent atoms (different species)
        - Example: Fe = [1, -1] creates antiferromagnetic configuration
        - Example: Fe = [1] makes all Fe atoms equivalent with moment 1
        """
        )

        # Configure each element
        for element in sorted(elements):
            # Auto-expand if magnetism is already configured for this element
            # A configured element has non-zero moments or multiple values
            current_mag = config_dict["magnetic_config"].get(element, [0])
            is_configured = bool(
                current_mag 
                and current_mag != [0] 
                and (len(current_mag) > 1 or (len(current_mag) == 1 and current_mag[0] != 0))
            )
            
            with st.expander(f"**{element}** Magnetization", expanded=is_configured):
                # Get current config for this element
                current_config = config_dict["magnetic_config"].get(element, [0])

                # Extract magnetic moments (handle both simple list and old dict format)
                if isinstance(current_config, dict) and "mag" in current_config:
                    # Old format with Hubbard - extract just magnetic moments
                    mag_moments = current_config.get("mag", [0])
                elif isinstance(current_config, list):
                    # Simple list format
                    mag_moments = current_config
                else:
                    # Default
                    mag_moments = [0]

                # Magnetic moments input
                st.markdown("**Magnetic Moments:**")
                moments_str = st.text_input(
                    f"Moments for {element} (comma-separated):",
                    value=", ".join(map(str, mag_moments)),
                    help="Example: '1, -1' for AFM, '1' for all equivalent, '0' for non-magnetic",
                    key=f"{key_prefix}_mag_moments_{element}",
                )

                # Parse moments
                try:
                    moments = [
                        float(m.strip()) for m in moments_str.split(",") if m.strip()
                    ]
                    if not moments:
                        moments = [0]
                except ValueError:
                    st.error(
                        f"Invalid input for {element}. Use numbers separated by commas."
                    )
                    moments = [0]

                # Store as simple list (magnetic configuration only)
                config_dict["magnetic_config"][element] = moments

                # Show what will be created
                num_species = len(moments)
                if num_species == 1:
                    st.caption(
                        f"✓ Will create 1 species: {element} with moment {moments[0]}"
                    )
                else:
                    species_names = [element] + [
                        f"{element}{i}" for i in range(1, num_species)
                    ]
                    st.caption(
                        f"✓ Will create {num_species} species: {', '.join(species_names)}"
                    )

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
            # Get the current qe_version and determine display value
            current_qe_version = config_dict.get("qe_version", "auto")

            # Determine which category it falls into for display
            if current_qe_version is None or current_qe_version == "auto":
                display_value = "auto"
            elif isinstance(current_qe_version, str):
                # Try to parse version to determine major version
                try:
                    major_version = int(current_qe_version.split(".")[0])
                    if major_version >= 7:
                        display_value = "7.x"
                    else:
                        display_value = "6.x"
                except (ValueError, IndexError):
                    # If parsing fails, default to auto
                    display_value = "auto"
            else:
                display_value = "auto"

            qe_version_display = st.selectbox(
                "QE Version Format:",
                ["auto", "6.x", "7.x"],
                index=["auto", "6.x", "7.x"].index(display_value),
                help="Select QE version format for compatibility (affects input file format)",
                key=f"{key_prefix}_qe_version",
            )

            # Store the version appropriately
            # If "auto" is selected, set to None to let xespresso auto-detect
            # Otherwise, preserve the actual version from config or use the category
            if qe_version_display == "auto":
                config_dict["qe_version"] = None
            elif current_qe_version and current_qe_version not in [
                "auto",
                "6.x",
                "7.x",
            ]:
                # Preserve actual version string if it exists (e.g., "7.4")
                config_dict["qe_version"] = current_qe_version
            else:
                # Use the category selection
                config_dict["qe_version"] = qe_version_display

        # Show summary
        with st.expander("📋 Configuration Summary", expanded=False):
            st.json(config_dict["magnetic_config"])
    else:
        # Clear magnetic config if disabled
        config_dict["magnetic_config"] = {}

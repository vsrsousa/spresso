"""
Standalone Hubbard configuration selector utility for GUI pages.

This module provides a component for configuring Hubbard parameters
independently of magnetic configuration.
"""

import streamlit as st
from typing import Dict, Set, Optional


# Element-to-orbital mapping for common elements requiring Hubbard U
ELEMENT_ORBITAL_MAP = {
    # 3d transition metals (first row)
    "Sc": ["3d"],
    "Ti": ["3d"],
    "V": ["3d"],
    "Cr": ["3d"],
    "Mn": ["3d"],
    "Fe": ["3d"],
    "Co": ["3d"],
    "Ni": ["3d"],
    "Cu": ["3d"],
    "Zn": ["3d"],
    # 4d transition metals (second row)
    "Y": ["4d"],
    "Zr": ["4d"],
    "Nb": ["4d"],
    "Mo": ["4d"],
    "Tc": ["4d"],
    "Ru": ["4d"],
    "Rh": ["4d"],
    "Pd": ["4d"],
    "Ag": ["4d"],
    "Cd": ["4d"],
    # 5d transition metals (third row)
    "La": ["5d"],
    "Hf": ["5d"],
    "Ta": ["5d"],
    "W": ["5d"],
    "Re": ["5d"],
    "Os": ["5d"],
    "Ir": ["5d"],
    "Pt": ["5d"],
    "Au": ["5d"],
    "Hg": ["5d"],
    # Lanthanides (4f)
    "Ce": ["4f"],
    "Pr": ["4f"],
    "Nd": ["4f"],
    "Pm": ["4f"],
    "Sm": ["4f"],
    "Eu": ["4f"],
    "Gd": ["4f"],
    "Tb": ["4f"],
    "Dy": ["4f"],
    "Ho": ["4f"],
    "Er": ["4f"],
    "Tm": ["4f"],
    "Yb": ["4f"],
    "Lu": ["4f"],
    # Actinides (5f)
    "Th": ["5f"],
    "Pa": ["5f"],
    "U": ["5f"],
    "Np": ["5f"],
    "Pu": ["5f"],
    "Am": ["5f"],
    "Cm": ["5f"],
    "Bk": ["5f"],
    "Cf": ["5f"],
    # p-block elements (sometimes need U)
    "O": ["2p"],
    "N": ["2p"],
    "C": ["2p"],
    "S": ["3p"],
    "P": ["3p"],
    "Si": ["3p"],
}


def get_suggested_orbitals(element: str) -> list:
    """
    Get suggested orbital types for a given element.

    Args:
        element: Element symbol (e.g., 'Fe', 'O')

    Returns:
        List of suggested orbital strings
    """
    return ELEMENT_ORBITAL_MAP.get(element, ["3d"])  # Default to 3d if unknown


def render_hubbard_selector(
    elements: Set[str], config_dict: dict, key_prefix: str = "calc"
) -> None:
    """
    Render a Hubbard configuration selector component.

    This component allows users to configure Hubbard U and V parameters
    for DFT+U calculations.

    Args:
        elements: Set of element symbols in the structure
        config_dict: Dictionary to store the Hubbard configuration
        key_prefix: Prefix for widget keys to avoid conflicts
    """
    st.subheader("⚛️ Hubbard (DFT+U) Configuration")

    # Initialize Hubbard config in config_dict if not present
    if "enable_hubbard" not in config_dict:
        config_dict["enable_hubbard"] = False
    if "hubbard_u" not in config_dict:
        config_dict["hubbard_u"] = {}

    # Checkbox to enable Hubbard configuration
    enable_hubbard = st.checkbox(
        "Enable Hubbard (DFT+U) Corrections",
        value=config_dict.get("enable_hubbard", False),
        help="Enable DFT+U for strongly correlated systems (transition metals, rare earths)",
        key=f"{key_prefix}_enable_hubbard",
    )
    config_dict["enable_hubbard"] = enable_hubbard

    if enable_hubbard:
        st.markdown("**Configure Hubbard U Parameters:**")
        st.info(
            """
        💡 **DFT+U Corrections:**
        - Used for strongly correlated systems (3d/4f electrons)
        - Common for transition metals (Fe, Mn, Co, Ni, etc.)
        - U values typically range from 2-8 eV
        - Check literature for appropriate U values for your system
        """
        )

        # Format selection - default to NEW format (QE >= 7.0)
        format_type = st.radio(
            "Hubbard Format:",
            ["New Format (QE >= 7.0) - Recommended", "Old Format (QE < 7.0)"],
            index=0 if config_dict.get("hubbard_format", "new") == "new" else 1,
            help="New format (QE 7.0+) is recommended and allows orbital-specific U parameters with different correction schemes",
            key=f"{key_prefix}_hubbard_format",
        )
        config_dict["hubbard_format"] = "new" if "New" in format_type else "old"

        # Hubbard Flavor Selection (U, U+J, U+V, etc.)
        st.markdown("---")
        st.markdown("**Hubbard Correction Flavor:**")
        st.info(
            """
        💡 **Different Hubbard correction schemes:**
        - **U only**: Simple on-site correction (most common)
        - **U+J**: On-site U with Hund's exchange J (for more accurate treatment)
        - **U+V**: On-site U with inter-site V interactions (for extended systems)
        - **U+J+V**: Complete correction with all parameters
        """
        )
        
        hubbard_flavor = st.selectbox(
            "Select correction type:",
            ["U only", "U+J", "U+V", "U+J+V"],
            index=["U only", "U+J", "U+V", "U+J+V"].index(
                config_dict.get("hubbard_flavor", "U only")
            ),
            help="Choose which Hubbard parameters to use. U only is most common, U+J for better accuracy, U+V for inter-site interactions.",
            key=f"{key_prefix}_hubbard_flavor",
        )
        config_dict["hubbard_flavor"] = hubbard_flavor
        
        st.markdown("---")

        # Configure U for each element
        for element in sorted(elements):
            # Auto-expand if Hubbard U is already configured for this element
            is_configured = element in config_dict["hubbard_u"] and config_dict["hubbard_u"][element] != 0.0
            
            with st.expander(f"**{element}** Hubbard U", expanded=is_configured):
                # Current U value
                current_u = config_dict["hubbard_u"].get(element, 0.0)

                if config_dict["hubbard_format"] == "new":
                    # New format with orbital specification
                    st.markdown("**Orbital-Specific U:**")

                    # Get suggested orbitals for this element
                    suggested_orbitals = get_suggested_orbitals(element)
                    default_orbital = suggested_orbitals[0]

                    # Show info about suggested orbitals
                    if len(suggested_orbitals) > 1:
                        st.info(
                            f"💡 Common orbitals for {element}: {', '.join(suggested_orbitals)}"
                        )
                    else:
                        st.info(f"💡 Common orbital for {element}: {default_orbital}")

                    # Get orbital type (use selectbox for common orbitals, allow custom)
                    current_orbital = config_dict.get(
                        f"hubbard_orbital_{element}", default_orbital
                    )

                    # Provide selectbox with suggested orbitals plus "Custom" option
                    orbital_options = suggested_orbitals + ["Custom..."]

                    # Determine which option to select
                    if current_orbital in suggested_orbitals:
                        orbital_select_idx = suggested_orbitals.index(current_orbital)
                    else:
                        orbital_select_idx = len(
                            suggested_orbitals
                        )  # "Custom..." option

                    orbital_selection = st.selectbox(
                        f"Select orbital for {element}:",
                        options=orbital_options,
                        index=orbital_select_idx,
                        help=f"Choose orbital for Hubbard U. Common choices for {element}: {', '.join(suggested_orbitals)}",
                        key=f"{key_prefix}_hubbard_orbital_select_{element}",
                    )

                    # If "Custom..." is selected, show text input
                    if orbital_selection == "Custom...":
                        orbital = st.text_input(
                            f"Custom orbital for {element}:",
                            value=(
                                current_orbital
                                if current_orbital not in suggested_orbitals
                                else ""
                            ),
                            help="Examples: '3d', '4f', '2p'",
                            key=f"{key_prefix}_hubbard_orbital_custom_{element}",
                        )
                    else:
                        orbital = orbital_selection

                    config_dict[f"hubbard_orbital_{element}"] = orbital

                    u_value = st.number_input(
                        f"U value (eV) for {element}-{orbital}:",
                        value=float(current_u),
                        min_value=0.0,
                        max_value=20.0,
                        step=0.1,
                        help="Hubbard U parameter in eV",
                        key=f"{key_prefix}_hubbard_u_{element}",
                    )
                else:
                    # Old format - just U value
                    u_value = st.number_input(
                        f"U value (eV) for {element}:",
                        value=float(current_u),
                        min_value=0.0,
                        max_value=20.0,
                        step=0.1,
                        help="Hubbard U parameter in eV",
                        key=f"{key_prefix}_hubbard_u_{element}",
                    )

                # Store U value
                if u_value > 0:
                    config_dict["hubbard_u"][element] = u_value
                elif element in config_dict["hubbard_u"]:
                    del config_dict["hubbard_u"][element]
                
                # Add J parameter if flavor includes J
                if "J" in config_dict.get("hubbard_flavor", "U only"):
                    st.markdown("---")
                    st.markdown("**J Parameter (Hund's Exchange):**")
                    
                    # Initialize J dict if not present
                    if "hubbard_j" not in config_dict:
                        config_dict["hubbard_j"] = {}
                    
                    current_j = config_dict["hubbard_j"].get(element, 0.0)
                    
                    j_value = st.number_input(
                        f"J value (eV) for {element}:",
                        value=float(current_j),
                        min_value=0.0,
                        max_value=5.0,
                        step=0.1,
                        help="Hund's exchange J parameter in eV (typically 0.5-1.0 eV)",
                        key=f"{key_prefix}_hubbard_j_{element}",
                    )
                    
                    # Store J value
                    if j_value > 0:
                        config_dict["hubbard_j"][element] = j_value
                    elif element in config_dict.get("hubbard_j", {}):
                        del config_dict["hubbard_j"][element]
        
        # Add V parameters if flavor includes V (inter-site interactions)
        if "V" in config_dict.get("hubbard_flavor", "U only"):
            st.markdown("---")
            st.markdown("**V Parameters (Inter-site Interactions):**")
            st.info(
                """
            💡 **Inter-site Hubbard V:**
            - V describes interactions between orbitals on different sites
            - Useful for extended systems or when nearest-neighbor interactions are important
            - Requires specifying pairs of species and orbitals
            """
            )
            
            # Initialize V list if not present
            if "hubbard_v" not in config_dict:
                config_dict["hubbard_v"] = []
            
            # Show existing V parameters
            if config_dict["hubbard_v"]:
                st.markdown("**Current V Parameters:**")
                for idx, v_param in enumerate(config_dict["hubbard_v"]):
                    col1, col2 = st.columns([4, 1])
                    with col1:
                        st.caption(
                            f"{idx+1}. {v_param.get('species1', '?')}-{v_param.get('orbital1', '?')} ↔ "
                            f"{v_param.get('species2', '?')}-{v_param.get('orbital2', '?')}: {v_param.get('value', 0)} eV"
                        )
                    with col2:
                        if st.button("🗑️", key=f"{key_prefix}_remove_v_{idx}", help="Remove this V parameter"):
                            config_dict["hubbard_v"].pop(idx)
                            st.rerun()
            
            # Add new V parameter
            with st.expander("➕ Add V Parameter", expanded=False):
                elements_list = sorted(list(elements))
                
                col1, col2 = st.columns(2)
                with col1:
                    v_species1 = st.selectbox(
                        "Species 1:",
                        elements_list,
                        key=f"{key_prefix}_v_species1",
                    )
                    v_orbital1 = st.text_input(
                        "Orbital 1:",
                        value="3d",
                        key=f"{key_prefix}_v_orbital1",
                        help="e.g., 3d, 4f, 2p",
                    )
                
                with col2:
                    v_species2 = st.selectbox(
                        "Species 2:",
                        elements_list,
                        key=f"{key_prefix}_v_species2",
                    )
                    v_orbital2 = st.text_input(
                        "Orbital 2:",
                        value="2p",
                        key=f"{key_prefix}_v_orbital2",
                        help="e.g., 3d, 4f, 2p",
                    )
                
                v_value = st.number_input(
                    "V value (eV):",
                    value=1.0,
                    min_value=0.0,
                    max_value=10.0,
                    step=0.1,
                    key=f"{key_prefix}_v_value",
                    help="Inter-site V parameter in eV",
                )
                
                if st.button("Add V Parameter", key=f"{key_prefix}_add_v"):
                    config_dict["hubbard_v"].append({
                        "species1": v_species1,
                        "orbital1": v_orbital1,
                        "species2": v_species2,
                        "orbital2": v_orbital2,
                        "value": v_value,
                        "i": 1,
                        "j": 1,
                    })
                    st.success(f"Added V: {v_species1}-{v_orbital1} ↔ {v_species2}-{v_orbital2}")
                    st.rerun()

        # Additional Hubbard options - Different flavors/projectors for new format
        if config_dict["hubbard_format"] == "new":
            st.markdown("---")
            st.markdown("**Projector Type (Flavor):**")
            st.info(
                "💡 **Different projector flavors** affect how Hubbard corrections are applied:\n"
                "- **ortho-atomic**: Orthogonalized atomic orbitals (recommended by QE)\n"
                "- **atomic**: Raw atomic orbitals\n"
                "- **norm-atomic**: Normalized atomic orbitals\n"
                "- **wf**: Wannier functions\n"
                "- **pseudo**: Pseudopotential projectors"
            )
            projector = st.selectbox(
                "Projector:",
                ["ortho-atomic", "atomic", "norm-atomic", "wf", "pseudo"],
                index=["ortho-atomic", "atomic", "norm-atomic", "wf", "pseudo"].index(
                    config_dict.get("hubbard_projector", "ortho-atomic")
                ),
                help="Choose the projector flavor for Hubbard calculations. Each flavor has different convergence and accuracy characteristics.",
                key=f"{key_prefix}_hubbard_projector",
            )
            config_dict["hubbard_projector"] = projector

        # Show summary
        with st.expander("📋 Hubbard Configuration Summary", expanded=False):
            if config_dict["hubbard_u"]:
                st.json(
                    {
                        "format": config_dict["hubbard_format"],
                        "U_parameters": config_dict["hubbard_u"],
                        "projector": (
                            config_dict.get("hubbard_projector")
                            if config_dict["hubbard_format"] == "new"
                            else "N/A"
                        ),
                    }
                )
            else:
                st.warning("No Hubbard U parameters configured")
    else:
        # Clear Hubbard config if disabled
        config_dict["hubbard_u"] = {}

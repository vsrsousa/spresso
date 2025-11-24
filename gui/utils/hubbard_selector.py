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
        st.info(
            """
        💡 **DFT+U Corrections:**
        - Used for strongly correlated systems (3d/4f electrons)
        - Common for transition metals (Fe, Mn, Co, Ni, etc.)
        - U values typically range from 2-8 eV
        """
        )

        # Format and Flavor selection in columns
        col1, col2 = st.columns(2)
        
        with col1:
            format_type = st.selectbox(
                "Hubbard Format:",
                ["New (QE >= 7.0)", "Old (QE < 7.0)"],
                index=0 if config_dict.get("hubbard_format", "new") == "new" else 1,
                help="New format (QE 7.0+) is recommended",
                key=f"{key_prefix}_hubbard_format",
            )
            config_dict["hubbard_format"] = "new" if "New" in format_type else "old"

        with col2:
            hubbard_flavor = st.selectbox(
                "Correction Type:",
                ["U only", "U+J", "U+V", "U+J+V"],
                index=["U only", "U+J", "U+V", "U+J+V"].index(
                    config_dict.get("hubbard_flavor", "U only")
                ),
                help="U only is most common, U+J for better accuracy",
                key=f"{key_prefix}_hubbard_flavor",
            )
            config_dict["hubbard_flavor"] = hubbard_flavor

        st.markdown("---")
        st.markdown("**Configure U Parameters per Element:**")
        
        # Create table-like header
        if config_dict["hubbard_format"] == "new":
            header_cols = st.columns([2, 2, 2, 2])
            with header_cols[0]:
                st.markdown("**Species**")
            with header_cols[1]:
                st.markdown("**Orbital**")
            with header_cols[2]:
                st.markdown("**U (eV)**")
            with header_cols[3]:
                st.markdown("**Typical**")
        else:
            header_cols = st.columns([2, 3, 2])
            with header_cols[0]:
                st.markdown("**Species**")
            with header_cols[1]:
                st.markdown("**U (eV)**")
            with header_cols[2]:
                st.markdown("**Typical**")
        
        st.markdown("---")

        # Configure U for each element in table format (always visible)
        for element in sorted(elements):
            current_u = config_dict["hubbard_u"].get(element, 0.0)
            suggested_orbitals = get_suggested_orbitals(element)
            default_orbital = suggested_orbitals[0]
            
            # Typical U values for common elements
            typical_u = {
                'Fe': 4.0, 'Co': 3.5, 'Ni': 3.0, 'Mn': 4.0, 'Cr': 3.5,
                'V': 3.0, 'Ti': 2.5, 'Cu': 4.0, 'Zn': 4.0,
                'Gd': 6.0, 'Nd': 5.0, 'Ce': 5.0, 'U': 4.0,
                'O': 8.0,  # For oxygen p orbitals
            }.get(element, 0.0)

            if config_dict["hubbard_format"] == "new":
                # New format with orbital specification
                cols = st.columns([2, 2, 2, 2])
                
                with cols[0]:
                    st.markdown(f"**{element}**")
                
                with cols[1]:
                    current_orbital = config_dict.get(
                        f"hubbard_orbital_{element}", default_orbital
                    )
                    orbital = st.selectbox(
                        f"Orbital for {element}",
                        options=suggested_orbitals + ["Custom"],
                        index=(
                            suggested_orbitals.index(current_orbital)
                            if current_orbital in suggested_orbitals
                            else 0
                        ),
                        key=f"{key_prefix}_hubbard_orbital_{element}",
                        label_visibility="collapsed"
                    )
                    if orbital == "Custom":
                        orbital = st.text_input(
                            f"Custom orbital for {element}",
                            value=current_orbital if current_orbital not in suggested_orbitals else "3d",
                            key=f"{key_prefix}_hubbard_orbital_custom_{element}",
                            label_visibility="collapsed"
                        )
                    config_dict[f"hubbard_orbital_{element}"] = orbital
                
                with cols[2]:
                    u_value = st.number_input(
                        f"U for {element}",
                        value=float(current_u),
                        min_value=0.0,
                        max_value=20.0,
                        step=0.5,
                        key=f"{key_prefix}_hubbard_u_{element}",
                        label_visibility="collapsed"
                    )
                
                with cols[3]:
                    if typical_u > 0:
                        st.caption(f"~{typical_u} eV")
                    else:
                        st.caption("N/A")
            else:
                # Old format - just U value
                cols = st.columns([2, 3, 2])
                
                with cols[0]:
                    st.markdown(f"**{element}**")
                
                with cols[1]:
                    u_value = st.number_input(
                        f"U for {element}",
                        value=float(current_u),
                        min_value=0.0,
                        max_value=20.0,
                        step=0.5,
                        key=f"{key_prefix}_hubbard_u_{element}",
                        label_visibility="collapsed"
                    )
                
                with cols[2]:
                    if typical_u > 0:
                        st.caption(f"~{typical_u} eV")
                    else:
                        st.caption("N/A")

            # Store U value
            if u_value > 0:
                config_dict["hubbard_u"][element] = u_value
            elif element in config_dict["hubbard_u"]:
                del config_dict["hubbard_u"][element]

        # Add J parameter section if flavor includes J
        if "J" in config_dict.get("hubbard_flavor", "U only"):
            st.markdown("---")
            st.markdown("**J Parameters (Hund's Exchange):**")
            
            if "hubbard_j" not in config_dict:
                config_dict["hubbard_j"] = {}
            
            # J parameter table header
            j_header = st.columns([2, 3, 2])
            with j_header[0]:
                st.markdown("**Species**")
            with j_header[1]:
                st.markdown("**J (eV)**")
            with j_header[2]:
                st.markdown("**Typical**")
            
            st.markdown("---")
            
            for element in sorted(elements):
                if element in config_dict["hubbard_u"] and config_dict["hubbard_u"][element] > 0:
                    current_j = config_dict["hubbard_j"].get(element, 0.0)
                    
                    cols = st.columns([2, 3, 2])
                    with cols[0]:
                        st.markdown(f"**{element}**")
                    with cols[1]:
                        j_value = st.number_input(
                            f"J for {element}",
                            value=float(current_j),
                            min_value=0.0,
                            max_value=5.0,
                            step=0.1,
                            key=f"{key_prefix}_hubbard_j_{element}",
                            label_visibility="collapsed"
                        )
                    with cols[2]:
                        st.caption("0.5-1.0 eV")
                    
                    if j_value > 0:
                        config_dict["hubbard_j"][element] = j_value
                    elif element in config_dict.get("hubbard_j", {}):
                        del config_dict["hubbard_j"][element]
        
        # Add V parameters if flavor includes V
        if "V" in config_dict.get("hubbard_flavor", "U only"):
            st.markdown("---")
            st.markdown("**V Parameters (Inter-site Interactions):**")
            
            if "hubbard_v" not in config_dict:
                config_dict["hubbard_v"] = []
            
            # Show existing V parameters
            if config_dict["hubbard_v"]:
                for idx, v_param in enumerate(config_dict["hubbard_v"]):
                    col1, col2 = st.columns([4, 1])
                    with col1:
                        st.caption(
                            f"{v_param.get('species1', '?')}-{v_param.get('orbital1', '?')} ↔ "
                            f"{v_param.get('species2', '?')}-{v_param.get('orbital2', '?')}: {v_param.get('value', 0)} eV"
                        )
                    with col2:
                        if st.button("🗑️", key=f"{key_prefix}_remove_v_{idx}"):
                            config_dict["hubbard_v"].pop(idx)
                            st.rerun()
            
            # Add new V parameter
            with st.expander("➕ Add V Parameter", expanded=False):
                elements_list = sorted(list(elements))
                
                col1, col2 = st.columns(2)
                with col1:
                    v_species1 = st.selectbox("Species 1:", elements_list, key=f"{key_prefix}_v_species1")
                    v_orbital1 = st.text_input("Orbital 1:", value="3d", key=f"{key_prefix}_v_orbital1")
                
                with col2:
                    v_species2 = st.selectbox("Species 2:", elements_list, key=f"{key_prefix}_v_species2")
                    v_orbital2 = st.text_input("Orbital 2:", value="2p", key=f"{key_prefix}_v_orbital2")
                
                v_value = st.number_input("V value (eV):", value=1.0, min_value=0.0, max_value=10.0, step=0.1, key=f"{key_prefix}_v_value")
                
                if st.button("Add V Parameter", key=f"{key_prefix}_add_v"):
                    config_dict["hubbard_v"].append({
                        "species1": v_species1, "orbital1": v_orbital1,
                        "species2": v_species2, "orbital2": v_orbital2,
                        "value": v_value, "i": 1, "j": 1,
                    })
                    st.rerun()

        # Projector type for new format
        if config_dict["hubbard_format"] == "new":
            st.markdown("---")
            projector = st.selectbox(
                "Projector Type:",
                ["ortho-atomic", "atomic", "norm-atomic", "wf", "pseudo"],
                index=["ortho-atomic", "atomic", "norm-atomic", "wf", "pseudo"].index(
                    config_dict.get("hubbard_projector", "ortho-atomic")
                ),
                help="ortho-atomic is recommended by QE",
                key=f"{key_prefix}_hubbard_projector",
            )
            config_dict["hubbard_projector"] = projector

        # Show summary
        with st.expander("📋 Hubbard Configuration Summary", expanded=False):
            if config_dict["hubbard_u"]:
                st.json({
                    "format": config_dict["hubbard_format"],
                    "U_parameters": config_dict["hubbard_u"],
                    "projector": config_dict.get("hubbard_projector") if config_dict["hubbard_format"] == "new" else "N/A",
                })
            else:
                st.warning("No Hubbard U parameters configured")
    else:
        # Clear Hubbard config if disabled
        config_dict["hubbard_u"] = {}

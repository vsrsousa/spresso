"""
Structure visualization component.
"""

import streamlit as st
from typing import Optional
from ase import Atoms


def show_structure_info(atoms: Optional[Atoms]):
    """Display structure information."""
    if atoms is None:
        st.error("No structure data available")
        return
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Number of Atoms", len(atoms))
    
    with col2:
        formula = atoms.get_chemical_formula()
        st.metric("Formula", formula)
    
    with col3:
        volume = atoms.get_volume()
        st.metric("Volume (Ų)", f"{volume:.2f}")
    
    # Atomic positions
    with st.expander("Atomic Positions"):
        positions_data = {
            'Atom': atoms.get_chemical_symbols(),
            'X (Å)': atoms.positions[:, 0],
            'Y (Å)': atoms.positions[:, 1],
            'Z (Å)': atoms.positions[:, 2],
        }
        st.dataframe(positions_data, width='stretch')
    
    # Cell parameters
    with st.expander("Cell Parameters"):
        cell = atoms.get_cell()
        st.write(f"**Cell vectors (Å):**")
        for i, vec in enumerate(cell):
            st.write(f"a{i+1}: [{vec[0]:.4f}, {vec[1]:.4f}, {vec[2]:.4f}]")


def structure_comparison(atoms1: Optional[Atoms], atoms2: Optional[Atoms],
                        label1: str = 'Structure 1', 
                        label2: str = 'Structure 2'):
    """Compare two structures."""
    if atoms1 is None or atoms2 is None:
        st.error("Both structures required for comparison")
        return
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader(label1)
        show_structure_info(atoms1)
    
    with col2:
        st.subheader(label2)
        show_structure_info(atoms2)

"""
Filter and search widgets.
"""

import streamlit as st
from ..config import COLORS
from ..utils.constants import CALCULATION_METHODS


def structure_filters():
    """Create structure filter widgets."""
    with st.expander("🔍 Filters", expanded=False):
        col1, col2 = st.columns(2)
        
        with col1:
            formula = st.text_input(
                "Chemical Formula",
                placeholder="e.g., Si, MnO, Al2O3",
                help="Filter by chemical formula (partial match)"
            )
        
        with col2:
            natoms = st.slider(
                "Number of Atoms",
                min_value=1,
                max_value=100,
                value=(1, 100),
                help="Filter by number of atoms"
            )
        
        return {
            'formula': formula if formula else None,
            'natoms_min': natoms[0],
            'natoms_max': natoms[1],
        }


def calculation_filters():
    """Create calculation filter widgets."""
    with st.expander("🔍 Filters", expanded=False):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            method = st.multiselect(
                "Calculation Method",
                options=list(CALCULATION_METHODS.keys()),
                help="Filter by calculation type"
            )
        
        with col2:
            status = st.multiselect(
                "Status",
                options=['Success', 'Failed'],
                default=['Success'],
                help="Filter by calculation status"
            )
        
        with col3:
            machine = st.text_input(
                "Machine",
                placeholder="e.g., medusa, local",
                help="Filter by execution machine"
            )
        
        return {
            'method': method if method else None,
            'status': status,
            'machine': machine if machine else None,
        }


def search_filters():
    """Create advanced search filter widgets."""
    st.subheader("Advanced Search")
    
    search_type = st.radio(
        "Search Type",
        ['Structures', 'Calculations', 'Both'],
        horizontal=True
    )
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        formula = st.text_input("Formula")
    
    with col2:
        energy_min = st.number_input("Energy Min (eV)", value=None)
        energy_max = st.number_input("Energy Max (eV)", value=None)
    
    with col3:
        calc_method = st.multiselect(
            "Calculation Method",
            options=list(CALCULATION_METHODS.keys())
        )
    
    return {
        'search_type': search_type,
        'formula': formula if formula else None,
        'energy_min': energy_min,
        'energy_max': energy_max,
        'method': calc_method if calc_method else None,
    }

"""
Advanced Search page.
"""

import streamlit as st
from ..database_interface import (
    search_structures, get_all_calculations, get_structure_details
)
from ..components.filters import search_filters
from ..utils.formatters import format_energy


def main():
    """Advanced Search main page."""
    st.set_page_config(page_title="Search", page_icon="🔎", layout="wide")
    
    st.title("🔎 Advanced Search")
    st.write("Search across structures and calculations with multiple criteria")
    
    # Initialize session state
    if 'search_performed' not in st.session_state:
        st.session_state.search_performed = False
    if 'struct_results' not in st.session_state:
        st.session_state.struct_results = None
    if 'calc_results' not in st.session_state:
        st.session_state.calc_results = None
    if 'search_type' not in st.session_state:
        st.session_state.search_type = None
    
    # Create search interface
    filters = search_filters()
    
    if st.button("🔍 Search", width='stretch', key="search_button"):
        st.session_state.search_performed = True
        st.session_state.search_type = filters['search_type']
        
        # Search structures
        if filters['search_type'] in ['Structures', 'Both']:
            struct_results = search_structures(
                formula=filters['formula'],
                energy_min=filters['energy_min'],
                energy_max=filters['energy_max'],
            )
            st.session_state.struct_results = struct_results
        
        # Search calculations
        if filters['search_type'] in ['Calculations', 'Both']:
            all_calcs = get_all_calculations()
            
            if not all_calcs.empty:
                calcs_results = all_calcs.copy()
                if filters['method']:
                    calcs_results = calcs_results[calcs_results['method'].isin(filters['method'])]
                st.session_state.calc_results = calcs_results
    
    # Display results if search was performed
    if st.session_state.search_performed:
        st.divider()
        
        search_type = st.session_state.search_type
        
        # Search structures
        if search_type in ['Structures', 'Both']:
            st.subheader("Structures Search")
            
            results = st.session_state.struct_results
            
            if results is None or results.empty:
                st.info("No structures found matching your criteria")
            else:
                st.success(f"Found {len(results)} structures")
                st.dataframe(results, width='stretch', hide_index=True)
                
                # Allow selection for detailed view
                if len(results) > 0:
                    selected_id = st.selectbox(
                        "View structure details:",
                        options=results['id'].tolist(),
                        format_func=lambda x: f"ID {x} - {results[results['id']==x]['formula'].values[0]}"
                    )
                    
                    if selected_id:
                        atoms = get_structure_details(selected_id)
                        if atoms:
                            st.write(f"**Formula:** {atoms.get_chemical_formula()}")
                            st.write(f"**Number of atoms:** {len(atoms)}")
                            st.write(f"**Volume:** {atoms.get_volume():.2f} Ų")
        
        # Search calculations
        if search_type in ['Calculations', 'Both']:
            st.subheader("Calculations Search")
            
            results = st.session_state.calc_results
            
            if results is None or results.empty:
                st.info("No calculations found matching your criteria")
            else:
                st.success(f"Found {len(results)} calculations")
                st.dataframe(
                    results[['hash', 'method', 'energy', 'steps', 'success']],
                    width='stretch',
                    hide_index=True
                )
    else:
        st.info("Configure your search criteria and click the Search button")


if __name__ == "__main__":
    main()

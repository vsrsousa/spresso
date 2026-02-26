"""
Derivation Lineage page.
"""

import streamlit as st
from ..database_interface import get_derivations, get_all_structures, get_structure_details
from ..components.derivation_graph import (
    show_derivation_tree, show_derivation_history, plot_energy_evolution
)
from ..components.structure_viewer import structure_comparison
from ..components.metrics import info_box


def main():
    """Derivation Lineage main page."""
    st.set_page_config(page_title="Lineage", page_icon="🌳", layout="wide")
    
    st.title("🌳 Derivation Lineage")
    st.write("Explore derivation relationships between structures")
    
    # Initialize session state
    if 'selected_source_id' not in st.session_state:
        st.session_state.selected_source_id = None
    if 'current_derivations' not in st.session_state:
        st.session_state.current_derivations = None
    
    # Get all structures for reference
    all_structs = get_all_structures()
    
    if all_structs.empty:
        st.warning("No structures in database")
        return
    
    # Select source structure
    st.subheader("Select Source Structure")
    selected_id = st.selectbox(
        "Choose a structure to see its derivations:",
        options=all_structs['id'].tolist(),
        format_func=lambda x: f"ID {x} - {all_structs[all_structs['id']==x]['formula'].values[0]}"
    )
    
    # Store in session and fetch derivations
    if selected_id != st.session_state.selected_source_id:
        st.session_state.selected_source_id = selected_id
        st.session_state.current_derivations = get_derivations(selected_id)
    
    if st.session_state.selected_source_id:
        derivations = st.session_state.current_derivations
        
        if derivations is None or derivations.empty:
            info_box(f"No derivations found from Structure {selected_id}")
        else:
            st.divider()
            
            # Display statistics
            st.subheader("Derivation Statistics")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric("Total Derivations", len(derivations))
            
            with col2:
                methods = derivations['method'].unique()
                st.metric("Calculation Types", len(methods))
            
            with col3:
                energy_change_total = derivations['energy_change'].sum()
                st.metric("Total Energy Change (eV)", f"{energy_change_total:.4f}")
            
            st.divider()
            
            # Derivation tree
            st.subheader("Derivation Tree")
            show_derivation_tree(derivations)
            
            st.divider()
            
            # Derivation history
            st.subheader("Derivation History")
            show_derivation_history(derivations)
            
            st.divider()
            
            # Energy evolution
            st.subheader("Energy Evolution")
            plot_energy_evolution(derivations)
            
            st.divider()
            
            # Detailed comparison
            st.subheader("Structure Comparison")
            
            col1, col2 = st.columns(2)
            
            with col1:
                compare_with = st.selectbox(
                    "Compare with derived structure:",
                    options=derivations['derived'].tolist(),
                    key="derived_select",
                    format_func=lambda x: f"Structure {x}"
                )
            
            with col2:
                if st.button("Compare", width='stretch'):
                    source_atoms = get_structure_details(selected_id)
                    derived_atoms = get_structure_details(compare_with)
                    
                    if source_atoms and derived_atoms:
                        structure_comparison(
                            source_atoms, derived_atoms,
                            label1=f"Source (ID {selected_id})",
                            label2=f"Derived (ID {compare_with})"
                        )


if __name__ == "__main__":
    main()

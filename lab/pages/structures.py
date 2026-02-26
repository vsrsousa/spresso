"""
Structures Explorer page.
"""

import streamlit as st
from ..database_interface import (
    get_all_structures, get_structure_details, search_structures, get_structure_creator,
    get_execution_history
)
from ..components.structure_viewer import show_structure_info
from ..components.filters import structure_filters
from ..components.metrics import stats_row
from ..utils.formatters import format_energy, format_timestamp


def main():
    """Structures Explorer main page."""
    st.set_page_config(page_title="Structures", page_icon="🏗️", layout="wide")
    
    st.title("🏗️ Structures Explorer")
    st.write("Browse all structures in the database")
    
    # Initialize session state
    if 'search_performed' not in st.session_state:
        st.session_state.search_performed = True  # Load all structures by default
    if 'search_results' not in st.session_state:
        st.session_state.search_results = None
    
    # Get filters from user
    filters = structure_filters()
    
    # Load structures on first run
    if st.session_state.search_results is None:
        st.session_state.search_results = get_all_structures()
    
    # Search button to apply filters
    if st.button("Search / Filter", width='stretch'):
        with st.spinner("Searching..."):
            if filters['formula']:
                results = search_structures(
                    formula=filters['formula'],
                )
            else:
                results = get_all_structures()
            
            st.session_state.search_results = results
            st.session_state.search_performed = True
    
    # Display results
    if st.session_state.search_performed and st.session_state.search_results is not None:
        results = st.session_state.search_results
        
        if results.empty:
            st.warning("No structures found matching the criteria")
        else:
            # Display statistics
            stats_row({
                "Total Structures": len(results),
                "Avg Energy": f"{results['energy'].mean():.4f} eV" if results['energy'].notna().any() else "N/A",
                "Min Energy": f"{results['energy'].min():.4f} eV" if results['energy'].notna().any() else "N/A",
            })
            
            st.divider()
            
            # Display table
            st.subheader("Structures")
            
            display_cols = ['id', 'formula', 'natoms', 'energy']
            st.dataframe(
                results[display_cols],
                width='stretch',
                hide_index=True,
            )
            
            # Select structure for detailed view
            selected_id = st.selectbox(
                "Select structure for details:",
                options=results['id'].tolist(),
                format_func=lambda x: f"ID {x} - {results[results['id']==x]['formula'].values[0]}"
            )
            
            if selected_id:
                st.divider()
                st.subheader(f"Structure Details - ID {selected_id}")
                
                # Check if this structure was derived
                creator = get_structure_creator(selected_id)
                
                if creator:
                    st.info("📊 This structure was derived from another structure")
                    
                    with st.expander("View derivation origin", expanded=True):
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.write("**Input Structure**")
                            input_atoms = get_structure_details(creator['source_structure_id'])
                            if input_atoms:
                                st.write(f"ID: {creator['source_structure_id']}")
                                st.write(f"Formula: {input_atoms.get_chemical_formula()}")
                            else:
                                st.write(f"ID: {creator['source_structure_id']}")
                        
                        with col2:
                            st.write("**Calculation**")
                            st.write(f"Method: {creator['method']}")
                            st.write(f"Steps: {creator['steps']}")
                            st.write(f"Status: {'✅ Success' if creator['success'] else '❌ Failed'}")
                        
                        with col3:
                            st.write("**Output (This structure)**")
                            st.write(f"ID: {selected_id}")
                            st.write(f"Energy Change: {creator['energy_change']:.4f} eV")
                            st.write(f"Final Energy: {creator['energy']:.6f} eV")
                        
                        # Execution history
                        if st.button("Show execution history", key=f"creator_exec_{selected_id}"):
                            exec_hist = get_execution_history(creator['calc_id'])
                            
                            if not exec_hist.empty:
                                st.write("**Execution History**")
                                for _, row in exec_hist.iterrows():
                                    st.write(f"""
                                    - **Machine**: {row['machine']}
                                    - **Wall Time**: {format_wall_time(row['wall_time'])}
                                    - **QE Version**: {row['qe_version'] or 'N/A'}
                                    - **Status**: {'✅ Success' if row['success'] else '❌ Failed'}
                                    - **Timestamp**: {row['timestamp']}
                                    """)
                            else:
                                st.info("No execution history available")
                    
                    st.divider()
                
                atoms = get_structure_details(selected_id)
                if atoms:
                    show_structure_info(atoms)


if __name__ == "__main__":
    main()


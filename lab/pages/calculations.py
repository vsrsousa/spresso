"""
Calculations Timeline page.
"""

import streamlit as st
import plotly.express as px
from ..database_interface import (
    get_all_calculations, get_execution_history, get_all_structures, get_structure_details,
    get_structure_creator
)
from ..components.filters import calculation_filters
from ..components.metrics import stats_row, success_box, error_box, info_box
from ..utils.formatters import format_energy, format_timestamp, format_wall_time


def get_calculations_for_structure(structure_id: int):
    """Get all calculations for a given input structure."""
    try:
        from lab.database_interface import get_databases
        wf = get_databases()
        if not wf:
            return None
        
        cursor = wf.provenance.conn.execute('''
            SELECT 
                id, calculation_hash, calculation_method, energy, 
                convergence_steps, success, timestamp, output_structure_id,
                convergence_params
            FROM calculations
            WHERE input_structure_id = ?
            ORDER BY timestamp DESC
        ''', (structure_id,))
        
        calcs = []
        for row in cursor:
            calcs.append({
                'id': row[0],
                'hash': row[1][:8] + '...',
                'method': row[2],
                'energy': row[3],
                'steps': row[4],
                'success': row[5],
                'timestamp': row[6],
                'output_structure_id': row[7],
                'convergence_params': row[8],
            })
        
        import pandas as pd
        return pd.DataFrame(calcs) if calcs else pd.DataFrame()
    except Exception as e:
        return None


def main():
    """Calculations Timeline main page."""
    st.set_page_config(page_title="Calculations", page_icon="⚙️", layout="wide")
    
    st.title("⚙️ Calculations Timeline")
    st.write("View calculations organized by input structure")
    
    # Initialize session state
    if 'selected_struct_for_calcs' not in st.session_state:
        st.session_state.selected_struct_for_calcs = None
    if 'struct_calcs_displayed' not in st.session_state:
        st.session_state.struct_calcs_displayed = None
    
    # Get all structures
    all_structs = get_all_structures()
    
    if all_structs.empty:
        st.warning("No structures in database")
        return
    
    # Select structure
    st.subheader("Select Input Structure")
    selected_struct = st.selectbox(
        "Choose a structure to see its calculations:",
        options=all_structs['id'].tolist(),
        format_func=lambda x: f"ID {x} - {all_structs[all_structs['id']==x]['formula'].values[0]}"
    )
    
    # Get input structure formula
    input_struct_row = all_structs[all_structs['id'] == selected_struct]
    input_formula = input_struct_row['formula'].values[0] if not input_struct_row.empty else 'N/A'
    
    # Fetch calculations for selected structure
    if selected_struct != st.session_state.selected_struct_for_calcs:
        st.session_state.selected_struct_for_calcs = selected_struct
        st.session_state.struct_calcs_displayed = get_calculations_for_structure(selected_struct)
    
    if st.session_state.struct_calcs_displayed is not None and not st.session_state.struct_calcs_displayed.empty:
        calcs = st.session_state.struct_calcs_displayed
        
        st.divider()
        
        # Check if this structure is derived (has a creator calculation)
        creator = get_structure_creator(selected_struct)
        
        # DERIVATION ORIGIN (if applicable)
        if creator:
            st.subheader("📊 Derivation Origin - How this structure was created")
            
            with st.expander(f"✅ {creator['method'].upper()} | Energy change: {creator['energy_change']:.4f} eV", expanded=True):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.write("**Input Structure**")
                    input_struct = creator['source_structure_id']
                    input_atoms = get_structure_details(input_struct)
                    if input_atoms:
                        st.write(f"ID: {input_struct}")
                        st.write(f"Formula: {input_atoms.get_chemical_formula()}")
                    else:
                        st.write(f"ID: {input_struct}")
                
                with col2:
                    st.write("**Calculation**")
                    st.write(f"Method: {creator['method']}")
                    st.write(f"Steps: {creator['steps']}")
                    st.write(f"Final Energy: {creator['energy']:.6f} eV")
                    st.write(f"Status: {'✅ Success' if creator['success'] else '❌ Failed'}")
                
                with col3:
                    st.write("**Output Structure (This one)**")
                    st.write(f"ID: {selected_struct}")
                    st.write(f"Formula: {input_formula}")
                    st.write(f"Energy Change: {creator['energy_change']:.4f} eV")
            
            st.divider()
            st.subheader("📈 Subsequent Calculations - Using this structure as input")
        else:
            st.subheader("⚙️ Calculations - Using this structure as input")
        
        # SUBSEQUENT CALCULATIONS
        if not calcs.empty:
            # Statistics
            success_count = calcs[calcs['success'] == True].shape[0]
            failed_count = calcs[calcs['success'] == False].shape[0]
            
            stats_row({
                "Total": len(calcs),
                "Success": success_count,
                "Failed": failed_count,
            })
            
            st.divider()
            
            # Timeline chart
            st.subheader("Calculation Timeline")
            
            if not calcs.empty and 'timestamp' in calcs.columns:
                fig = px.bar(
                    calcs,
                    x='timestamp',
                    y='steps',
                    color='success',
                    title=f'Calculations from {input_formula} (ID {selected_struct})',
                    labels={'steps': 'SCF Iterations', 'success': 'Status'},
                    color_discrete_map={True: '#5FB878', False: '#D0021B'},
                    hover_data=['method', 'energy']
                )
                st.plotly_chart(fig, width='stretch')
            
            st.divider()
            
            # Calculations with expandable details
            st.subheader("Calculations Details")
            
            for idx, (_, calc) in enumerate(calcs.iterrows()):
                with st.expander(
                    f"{'✅' if calc['success'] else '❌'} {calc['method'].upper()} | "
                    f"Energy: {calc['energy']:.4f} eV | Output: Struct {int(calc['output_structure_id']) if not (isinstance(calc['output_structure_id'], float) and calc['output_structure_id'] != calc['output_structure_id']) else 'N/A'} | "
                    f"{calc['timestamp']}"
                ):
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.write("**Input Structure**")
                        st.write(f"ID: {selected_struct}")
                        st.write(f"Formula: {input_formula}")
                        
                        st.write("\n**Calculation Data**")
                        st.write(f"Method: {calc['method']}")
                        st.write(f"Hash: {calc['hash']}")
                        st.write(f"SCF Iterations: {calc['steps']}")
                        
                    with col2:
                        st.write("**Output Structure**")
                        output_struct_id = calc['output_structure_id']
                        
                        if output_struct_id is not None and not (isinstance(output_struct_id, float) and output_struct_id != output_struct_id):  # Check for NaN
                            output_struct_id_int = int(output_struct_id)
                            output_atoms = get_structure_details(output_struct_id_int)
                            if output_atoms:
                                st.write(f"ID: {output_struct_id_int}")
                                st.write(f"Formula: {output_atoms.get_chemical_formula()}")
                                st.write(f"Atoms: {len(output_atoms)}")
                                st.write(f"Volume: {output_atoms.get_volume():.2f} Ų")
                            else:
                                st.write(f"ID: {output_struct_id_int}")
                                st.write("Could not load structure data")
                        else:
                            st.write("No output structure (calculation may have failed)")
                        
                        st.write("\n**Energy**")
                        st.write(f"Final: {calc['energy']:.6f} eV")
                    
                    # Execution history
                    if st.button(f"Show execution history", key=f"exec_{calc['id']}"):
                        exec_hist = get_execution_history(calc['id'])
                        
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
        else:
            info_box(f"No subsequent calculations for Structure {selected_struct}")
    else:
        if st.session_state.selected_struct_for_calcs is not None:
            info_box(f"No calculations found for Structure {selected_struct} ({input_formula})")


if __name__ == "__main__":
    main()

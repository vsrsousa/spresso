"""
Derivation graph visualization.
"""

import streamlit as st
import plotly.graph_objects as go
import pandas as pd
from ..config import COLORS


def show_derivation_tree(derivations_df: pd.DataFrame):
    """Display derivation tree as a hierarchical visualization."""
    if derivations_df.empty:
        st.info("No derivations found")
        return
    
    # Create a simple text-based tree
    st.write("**Derivation Tree:**")
    
    for _, row in derivations_df.iterrows():
        source = row['source']
        derived = row['derived']
        method = row['method']
        energy_change = row.get('energy_change', 'N/A')
        
        st.write(f"""
        ```
        Structure {source}
            ↓ [{method}]
        Structure {derived}
        Energy change: {energy_change}
        ```
        """)


def show_derivation_history(derivations_df: pd.DataFrame):
    """Show derivation history as timeline."""
    if derivations_df.empty:
        st.info("No derivations found")
        return
    
    # Sort by timestamp if available
    if 'timestamp' in derivations_df.columns:
        derivations_df = derivations_df.sort_values('timestamp')
    
    # Display as table
    st.dataframe(
        derivations_df[['source', 'derived', 'method', 'energy_change']],
        width='stretch',
        hide_index=True
    )


def plot_energy_evolution(derivations_df: pd.DataFrame):
    """Plot energy evolution across derivations."""
    if derivations_df.empty or 'energy_change' not in derivations_df.columns:
        st.info("No energy data available")
        return
    
    # Create cumulative energy
    energy_data = derivations_df[['derived', 'energy_change']].copy()
    energy_data['cumulative'] = energy_data['energy_change'].cumsum()
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=energy_data['derived'].astype(str),
        y=energy_data['energy_change'],
        name='Energy Change',
        marker=dict(color=COLORS['primary']),
    ))
    
    fig.add_trace(go.Scatter(
        x=energy_data['derived'].astype(str),
        y=energy_data['cumulative'],
        name='Cumulative Energy',
        mode='lines+markers',
        marker=dict(color=COLORS['warning']),
        line=dict(width=2),
    ))
    
    fig.update_layout(
        title="Energy Evolution Across Derivations",
        xaxis_title="Structure ID",
        yaxis_title="Energy (eV)",
        hovermode='x unified',
        height=400,
    )
    
    st.plotly_chart(fig, width='stretch')

"""
XESPRESSO Lab - Main Application

Entry point for the Streamlit interface.
Run with: streamlit run lab/app.py
"""

import streamlit as st
from pathlib import Path
import sys

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lab.config import STREAMLIT_CONFIG, DATABASES_DIR
from lab.database_interface import get_statistics
from lab.components.metrics import stats_row
from lab.pages import (
    structures as page_structures,
    calculations as page_calculations,
    lineage as page_lineage,
    search as page_search,
)


def configure_page():
    """Configure Streamlit page settings."""
    st.set_page_config(
        page_title="XESPRESSO Lab",
        page_icon="🔬",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def sidebar_navigation():
    """Create sidebar navigation."""
    with st.sidebar:
        st.title("🔬 XESPRESSO Lab")
        st.divider()
        
        page = st.radio(
            "Navigate to:",
            ["Home", "Structures", "Calculations", "Lineage", "Search"],
        )
        
        st.divider()
        
        # Database info
        with st.expander("📊 Database Info"):
            st.write(f"**Location:** `{DATABASES_DIR}`")
            
            stats = get_statistics()
            if stats:
                st.metric("Structures", stats.get('total_structures', 0))
                st.metric("Calculations", stats.get('total_calculations', 0))
                st.metric("Derivations", stats.get('total_derivations', 0))
                st.metric("Executions", stats.get('total_executions', 0))
            else:
                st.warning("Could not load database statistics")
        
        st.divider()
        st.caption("v0.1.0 - Provenance Explorer")
    
    return page


def home_page():
    """Home page."""
    st.title("🔬 XESPRESSO Lab")
    st.subheader("Provenance Explorer & Analysis Interface")
    
    st.write("""
    Welcome to XESPRESSO Lab! This interface allows you to explore:
    
    - **Structures**: Browse and search all atomic structures in your database
    - **Calculations**: View calculation history and execution details
    - **Lineage**: Explore derivation relationships between structures
    - **Search**: Advanced search with multiple criteria
    """)
    
    st.divider()
    
    # Quick statistics
    st.subheader("Quick Statistics")
    
    stats = get_statistics()
    if stats:
        stats_row({
            "Structures": stats.get('total_structures', 0),
            "Calculations": stats.get('total_calculations', 0),
            "Derivations": stats.get('total_derivations', 0),
        })
    else:
        st.warning("Could not load database. Please ensure databases exist at:")
        st.code(f"{DATABASES_DIR}")
    
    st.divider()
    
    # Getting started
    st.subheader("Getting Started")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Structures**")
        st.write("Browse all structures and view their properties")
    
    with col2:
        st.write("**Lineage**")
        st.write("Explore how structures evolved through calculations")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Calculations**")
        st.write("Check calculation history and execution logs")
    
    with col2:
        st.write("**Search**")
        st.write("Find structures and calculations by criteria")


def main():
    """Main application."""
    configure_page()
    
    page = sidebar_navigation()
    
    if page == "Home":
        home_page()
    elif page == "Structures":
        page_structures.main()
    elif page == "Calculations":
        page_calculations.main()
    elif page == "Lineage":
        page_lineage.main()
    elif page == "Search":
        page_search.main()


if __name__ == "__main__":
    main()

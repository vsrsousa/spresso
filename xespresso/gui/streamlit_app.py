"""
Streamlit GUI for xespresso - Quantum ESPRESSO Configuration Interface

This application provides a user-friendly interface for:
- Configuring machines (local/remote execution environments)
- Setting up Quantum ESPRESSO codes
- Viewing and selecting molecular structures
- Configuring calculations and workflows
- Submitting computational jobs

This is the main entry point for the modular GUI.
"""

import streamlit as st
import os
import json
import tempfile
from pathlib import Path
import traceback

# Configure page
st.set_page_config(
    page_title="xespresso GUI",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Import session manager functions early to get current session
try:
    from xespresso.gui.utils.session_manager import get_current_session_id, get_active_sessions
    SESSION_MANAGER_AVAILABLE = True
except ImportError:
    SESSION_MANAGER_AVAILABLE = False

# Title and description with session info
if SESSION_MANAGER_AVAILABLE:
    current_id = get_current_session_id()
    active_sessions = get_active_sessions()
    session_name = active_sessions.get(current_id, {}).get('name', 'Session 1')
    st.title(f"⚛️ xespresso - {session_name}")
else:
    st.title("⚛️ xespresso Configuration GUI")

st.markdown("""
Welcome to the xespresso graphical interface for Quantum ESPRESSO calculations.
Configure your computational environment, select structures, and submit jobs easily.
""")

# Import modular page components
try:
    from xespresso.gui.pages import (
        render_machine_config_page,
        render_codes_config_page,
        render_pseudopotentials_config_page,
        render_structure_viewer_page,
        render_calculation_setup_page,
        render_workflow_builder_page,
        render_job_submission_page,
        render_results_postprocessing_page
    )
    PAGES_AVAILABLE = True
except ImportError as e:
    st.error(f"⚠️ Error importing page modules: {e}")
    PAGES_AVAILABLE = False

# Import xespresso modules for pages that still need them inline
try:
    from xespresso.machines.machine import Machine
    from xespresso.machines.config.loader import (
        load_machine, save_machine, list_machines,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    from xespresso.codes.manager import (
        detect_qe_codes, load_codes_config,
        DEFAULT_CODES_DIR
    )
    from xespresso.workflow import (
        CalculationWorkflow, quick_scf, quick_relax, PRESETS
    )
    from xespresso import Espresso
    XESPRESSO_AVAILABLE = True
except ImportError as e:
    st.error(f"⚠️ Error importing xespresso modules: {e}")
    st.info("Make sure xespresso is properly installed.")
    XESPRESSO_AVAILABLE = False

# Import visualization modules
try:
    from ase import io
    from ase.build import bulk, molecule
    import numpy as np
    ASE_AVAILABLE = True
except ImportError:
    st.warning("⚠️ ASE not available. Structure visualization will be limited.")
    ASE_AVAILABLE = False

try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    st.warning("⚠️ Plotly not available. 3D visualization will be limited.")
    PLOTLY_AVAILABLE = False

# Import utility functions
try:
    from xespresso.gui.utils import validate_path, create_3d_structure_plot, display_structure_info
    from xespresso.gui.utils.session_manager import render_session_manager, get_current_session_id, get_active_sessions
    from xespresso.gui.utils.directory_browser import render_directory_browser
    UTILS_AVAILABLE = True
except ImportError as e:
    st.warning(f"⚠️ GUI utilities not fully available: {e}")
    UTILS_AVAILABLE = False
    # Fallback implementations
    def validate_path(path, allow_creation=False):
        """Fallback path validation."""
        return True, path, None
    
    def create_3d_structure_plot(atoms):
        """Fallback plot function."""
        return None
    
    def display_structure_info(atoms):
        """Fallback structure info display."""
        st.write(f"Structure: {atoms.get_chemical_formula()}")
    
    def render_session_manager(key="session_manager"):
        """Fallback session manager."""
        pass
    
    def render_directory_browser(key="directory_browser", initial_path=None, help_text=""):
        """Fallback directory browser."""
        return initial_path or os.path.expanduser("~")

# Initialize session state
if 'current_structure' not in st.session_state:
    st.session_state.current_structure = None
if 'current_machine' not in st.session_state:
    st.session_state.current_machine = None
if 'current_machine_name' not in st.session_state:
    st.session_state.current_machine_name = None
if 'current_codes' not in st.session_state:
    st.session_state.current_codes = None
if 'selected_code_version' not in st.session_state:
    st.session_state.selected_code_version = None
if 'workflow_config' not in st.session_state:
    st.session_state.workflow_config = {}

# Working directory is now session-specific
# Initialize only if not already set (will be restored from session state when switching)
if 'working_directory' not in st.session_state:
    st.session_state.working_directory = os.path.expanduser("~")

# Add session manager to sidebar first
if UTILS_AVAILABLE:
    render_session_manager()

st.sidebar.markdown("---")

# Working Directory selector (for calculation pages only)
st.sidebar.subheader("📁 Working Directory")

# Use enhanced directory browser
if UTILS_AVAILABLE:
    selected_workdir = render_directory_browser(
        key="workdir_browser",
        initial_path=st.session_state.get('working_directory', os.path.expanduser("~")),
        help_text="Choose the base directory where calculation folders will be created"
    )
    
    # Update session state when directory changes
    if selected_workdir != st.session_state.get('working_directory'):
        st.session_state.working_directory = selected_workdir
else:
    # Fallback to simple selectbox if utils not available
    common_dirs = [
        os.path.expanduser("~"),
        os.getcwd(),
        os.path.join(os.path.expanduser("~"), "calculations"),
        os.path.join(os.path.expanduser("~"), "Documents"),
        os.path.join(os.path.expanduser("~"), "Desktop"),
    ]
    
    if st.session_state.working_directory not in common_dirs:
        common_dirs.insert(0, st.session_state.working_directory)
    
    def format_dir(path):
        if path == os.path.expanduser("~"):
            return "🏠 Home"
        elif path == os.getcwd():
            return "📂 Current Directory"
        elif path.endswith("calculations"):
            return "📊 Calculations"
        elif path.endswith("Documents"):
            return "📄 Documents"
        elif path.endswith("Desktop"):
            return "🖥️ Desktop"
        else:
            return f"📁 {os.path.basename(path)}"
    
    selected_workdir = st.sidebar.selectbox(
        "Select working directory:",
        options=common_dirs,
        format_func=format_dir,
        key="workdir_selector",
        help="Choose the base directory where calculation folders will be created"
    )
    
    if selected_workdir != st.session_state.working_directory:
        st.session_state.working_directory = selected_workdir
    
    st.sidebar.caption(f"📍 {st.session_state.working_directory}")
    st.sidebar.info("💡 Calculation folders will be created here based on calc/label")

# Sidebar navigation
st.sidebar.markdown("---")
st.sidebar.title("Navigation")

# Add toggle for configuration pages
show_config = st.sidebar.checkbox(
    "⚙️ Show Configuration",
    value=False,
    help="Show machine and codes configuration (needed only for initial setup)"
)

st.sidebar.markdown("---")

# Build page list based on whether config is shown
if show_config:
    page = st.sidebar.radio(
        "Select Page:",
        [
            "🖥️ Machine Configuration",
            "⚙️ Codes Configuration",
            "🧪 Pseudopotentials Configuration",
            "🔬 Structure Viewer",
            "📊 Calculation Setup",
            "🔄 Workflow Builder",
            "🚀 Job Submission & Files",
            "📈 Results & Post-Processing"
        ]
    )
else:
    page = st.sidebar.radio(
        "Select Page:",
        [
            "🔬 Structure Viewer",
            "📊 Calculation Setup",
            "🔄 Workflow Builder",
            "🚀 Job Submission & Files",
            "📈 Results & Post-Processing"
        ]
    )

# Determine if this is a calculation page (for any page-specific logic later)
is_calculation_page = page not in ["🖥️ Machine Configuration", "⚙️ Codes Configuration", "🧪 Pseudopotentials Configuration"]

# Page routing
if page == "🖥️ Machine Configuration":
    if PAGES_AVAILABLE:
        render_machine_config_page()
    else:
        st.error("Page modules not available. Please check installation.")

elif page == "⚙️ Codes Configuration":
    if PAGES_AVAILABLE:
        render_codes_config_page()
    else:
        st.error("Page modules not available. Please check installation.")

elif page == "🧪 Pseudopotentials Configuration":
    if PAGES_AVAILABLE:
        render_pseudopotentials_config_page()
    else:
        st.error("Page modules not available. Please check installation.")

elif page == "🔬 Structure Viewer":
    if PAGES_AVAILABLE:
        render_structure_viewer_page()
    else:
        st.error("Page modules not available. Please check installation.")

elif page == "📊 Calculation Setup":
    if PAGES_AVAILABLE:
        render_calculation_setup_page()
    else:
        st.error("Page modules not available. Please check installation.")

elif page == "🔄 Workflow Builder":
    if PAGES_AVAILABLE:
        render_workflow_builder_page()
    else:
        st.error("Page modules not available. Please check installation.")

# Page 6: Job Submission & File Management
elif page == "🚀 Job Submission & Files":
    if PAGES_AVAILABLE:
        render_job_submission_page()
    else:
        st.error("Page modules not available. Please check installation.")

# Page 7: Results & Post-Processing
elif page == "📈 Results & Post-Processing":
    if PAGES_AVAILABLE:
        render_results_postprocessing_page()
    else:
        st.error("Page modules not available. Please check installation.")

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("""
### About
**xespresso GUI** - Streamlit interface for Quantum ESPRESSO calculations

Version: 1.0.0

[Documentation](https://github.com/superstar54/xespresso) | 
[Report Issue](https://github.com/superstar54/xespresso/issues)
""")

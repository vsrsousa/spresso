"""
Streamlit GUI for xespresso - Quantum ESPRESSO Configuration Interface

This application provides a user-friendly interface for:
- Configuring machines (local/remote execution environments)
- Setting up Quantum ESPRESSO codes
- Viewing and selecting molecular structures
- Configuring calculations and workflows
- Submitting computational jobs
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

# Title and description
st.title("⚛️ xespresso Configuration GUI")
st.markdown("""
Welcome to the xespresso graphical interface for Quantum ESPRESSO calculations.
Configure your computational environment, select structures, and submit jobs easily.
""")

# Import xespresso modules
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
if 'local_workdir' not in st.session_state:
    st.session_state.local_workdir = os.getcwd()

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio(
    "Select Configuration Step:",
    [
        "🖥️ Machine Configuration",
        "⚙️ Codes Configuration", 
        "🔬 Structure Viewer",
        "📊 Calculation Setup",
        "🔄 Workflow Builder",
        "🚀 Job Submission",
        "📈 Results & Post-Processing"
    ]
)

# Helper function for path validation
def validate_path(path, allow_creation=False):
    """
    Validate and sanitize file paths to prevent path injection.
    
    Args:
        path: Path to validate
        allow_creation: If True, allow non-existent paths (for file creation)
    
    Returns:
        Tuple of (is_valid, normalized_path, error_message)
    """
    if not path:
        return False, None, "Path cannot be empty"
    
    try:
        # Normalize and resolve the path
        normalized = os.path.abspath(os.path.expanduser(path))
        
        # Check for path traversal attempts
        if '..' in os.path.relpath(normalized, os.path.expanduser('~')):
            # Allow if it's an absolute path or in allowed directories
            allowed_dirs = ['/tmp', '/home', '/Users', os.path.expanduser('~')]
            if not any(normalized.startswith(d) for d in allowed_dirs):
                return False, None, "Path traversal not allowed"
        
        # Check if path exists (if required)
        if not allow_creation and not os.path.exists(normalized):
            return False, normalized, f"Path does not exist: {normalized}"
        
        return True, normalized, None
        
    except Exception as e:
        return False, None, f"Invalid path: {str(e)}"

# Helper functions for structure visualization
def create_3d_structure_plot(atoms):
    """Create a 3D plotly visualization of atomic structure."""
    if not PLOTLY_AVAILABLE:
        return None
    
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Color map for common elements
    color_map = {
        'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red',
        'F': 'green', 'P': 'orange', 'S': 'yellow',
        'Cl': 'green', 'Fe': 'brown', 'Cu': 'brown',
        'Al': 'silver', 'Si': 'pink', 'Pt': 'silver'
    }
    
    colors = [color_map.get(s, 'purple') for s in symbols]
    
    # Create scatter plot
    fig = go.Figure(data=[go.Scatter3d(
        x=positions[:, 0],
        y=positions[:, 1],
        z=positions[:, 2],
        mode='markers+text',
        marker=dict(
            size=12,
            color=colors,
            line=dict(color='black', width=1)
        ),
        text=symbols,
        textposition="top center",
        hovertemplate='<b>%{text}</b><br>x: %{x:.2f}<br>y: %{y:.2f}<br>z: %{z:.2f}<extra></extra>'
    )])
    
    # Add cell visualization if present
    if atoms.cell is not None and atoms.pbc.any():
        cell = atoms.cell.array
        # Draw cell edges
        edges = [
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],  # bottom
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 1],  # top
            [1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0], [0, 1, 0], [0, 1, 1]
        ]
        
        edge_points = np.array([np.dot(edge, cell) for edge in edges])
        
        fig.add_trace(go.Scatter3d(
            x=edge_points[:, 0],
            y=edge_points[:, 1],
            z=edge_points[:, 2],
            mode='lines',
            line=dict(color='black', width=2),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    fig.update_layout(
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            aspectmode='data'
        ),
        margin=dict(l=0, r=0, t=0, b=0),
        height=500
    )
    
    return fig

def display_structure_info(atoms):
    """Display information about the atomic structure."""
    st.subheader("Structure Information")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Number of Atoms", len(atoms))
        st.metric("Chemical Formula", atoms.get_chemical_formula())
    
    with col2:
        symbols = atoms.get_chemical_symbols()
        unique_elements = list(set(symbols))
        st.metric("Unique Elements", len(unique_elements))
        st.write("**Elements:**", ", ".join(sorted(unique_elements)))
    
    with col3:
        if atoms.cell is not None:
            st.metric("Cell Volume", f"{atoms.get_volume():.2f} Å³")
            pbc_str = "".join(["T" if p else "F" for p in atoms.pbc])
            st.metric("PBC", pbc_str)
    
    # Display cell parameters
    if atoms.cell is not None and atoms.pbc.any():
        st.subheader("Cell Parameters")
        cell_params = atoms.cell.cellpar()
        col1, col2, col3 = st.columns(3)
        with col1:
            st.write(f"**a:** {cell_params[0]:.3f} Å")
            st.write(f"**b:** {cell_params[1]:.3f} Å")
            st.write(f"**c:** {cell_params[2]:.3f} Å")
        with col2:
            st.write(f"**α:** {cell_params[3]:.2f}°")
            st.write(f"**β:** {cell_params[4]:.2f}°")
            st.write(f"**γ:** {cell_params[5]:.2f}°")

# Page 1: Machine Configuration
if page == "🖥️ Machine Configuration":
    st.header("Machine Configuration")
    st.markdown("""
    Configure the computational machine/cluster where calculations will run.
    Supports both local and remote (SSH) execution environments.
    
    **Note:** Saving a machine configuration creates/updates machine-specific JSON files in `~/.xespresso/machines/`.
    Pre-configured machines in `machines.json` are not modified.
    """)
    
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available. Cannot configure machines.")
    else:
        # List existing machines
        st.subheader("Existing Machines")
        try:
            machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            if machines_list:
                # Use session state for persistent selection
                default_idx = 0
                if st.session_state.current_machine_name and st.session_state.current_machine_name in machines_list:
                    default_idx = machines_list.index(st.session_state.current_machine_name) + 1
                
                selected_machine = st.selectbox(
                    "Select a machine to edit or view:",
                    ["[Create New]"] + machines_list,
                    index=default_idx,
                    key="machine_selector"
                )
                
                # Update session state
                if selected_machine != "[Create New]":
                    st.session_state.current_machine_name = selected_machine
            else:
                st.info("No machines configured yet. Create your first machine below.")
                selected_machine = "[Create New]"
        except Exception as e:
            st.warning(f"Could not load machines list: {e}")
            selected_machine = "[Create New]"
        
        # Create or edit machine
        st.subheader("Machine Configuration")
        
        # Load existing machine if selected
        if selected_machine != "[Create New]":
            try:
                machine = load_machine(DEFAULT_CONFIG_PATH, selected_machine, DEFAULT_MACHINES_DIR, return_object=True)
                st.success(f"✅ Loaded machine: {selected_machine}")
                st.session_state.current_machine = machine
                
                # Display current configuration
                with st.expander("📋 View Current Configuration", expanded=False):
                    st.json({
                        "name": machine.name,
                        "execution": machine.execution,
                        "scheduler": machine.scheduler,
                        "workdir": machine.workdir,
                        "nprocs": machine.nprocs,
                        "launcher": machine.launcher,
                        "use_modules": machine.use_modules if hasattr(machine, 'use_modules') else False,
                        "modules": machine.modules if hasattr(machine, 'modules') else [],
                    })
            except Exception as e:
                st.error(f"Error loading machine: {e}")
                machine = None
        else:
            machine = None
        
        # Configuration form
        with st.form("machine_config_form"):
            col1, col2 = st.columns(2)
            
            with col1:
                machine_name = st.text_input(
                    "Machine Name",
                    value=machine.name if machine else "",
                    help="Unique identifier for this machine"
                )
                
                execution = st.selectbox(
                    "Execution Mode",
                    ["local", "remote"],
                    index=0 if not machine or machine.execution == "local" else 1
                )
                
                scheduler = st.selectbox(
                    "Scheduler Type",
                    ["direct", "slurm", "pbs", "sge"],
                    index=["direct", "slurm", "pbs", "sge"].index(machine.scheduler) if machine else 0,
                    help="Job scheduler system"
                )
            
            with col2:
                workdir = st.text_input(
                    "Working Directory",
                    value=machine.workdir if machine else "./calculations",
                    help="Directory for calculation files"
                )
                
                nprocs = st.number_input(
                    "Number of Processors",
                    min_value=1,
                    value=machine.nprocs if machine else 1,
                    help="Default number of processors"
                )
                
                launcher = st.text_input(
                    "MPI Launcher",
                    value=machine.launcher if machine else "mpirun -np {nprocs}",
                    help="MPI launch command template"
                )
            
            # Remote configuration
            if execution == "remote":
                st.subheader("Remote Connection Settings")
                col1, col2 = st.columns(2)
                
                with col1:
                    host = st.text_input(
                        "Host",
                        value=machine.host if machine and machine.is_remote else "",
                        help="Remote hostname or IP"
                    )
                    username = st.text_input(
                        "Username",
                        value=machine.username if machine and machine.is_remote else "",
                        help="SSH username"
                    )
                
                with col2:
                    port = st.number_input(
                        "SSH Port",
                        min_value=1,
                        max_value=65535,
                        value=machine.port if machine and machine.is_remote else 22
                    )
                    ssh_key = st.text_input(
                        "SSH Key Path",
                        value=machine.auth.get("ssh_key", "~/.ssh/id_rsa") if machine and machine.is_remote else "~/.ssh/id_rsa",
                        help="Path to SSH private key"
                    )
            
            # Module configuration
            st.subheader("Environment Modules")
            use_modules = st.checkbox(
                "Use Environment Modules",
                value=machine.use_modules if machine else False
            )
            
            if use_modules:
                modules_str = st.text_area(
                    "Modules to Load (one per line)",
                    value="\n".join(machine.modules) if machine and machine.modules else "",
                    help="Environment modules to load before execution"
                )
            
            # Advanced settings
            with st.expander("Advanced Settings"):
                prepend = st.text_area(
                    "Prepend Commands",
                    value="\n".join(machine.prepend) if machine and isinstance(machine.prepend, list) else (machine.prepend if machine else ""),
                    help="Commands to run before calculation"
                )
                postpend = st.text_area(
                    "Postpend Commands", 
                    value="\n".join(machine.postpend) if machine and isinstance(machine.postpend, list) else (machine.postpend if machine else ""),
                    help="Commands to run after calculation"
                )
                env_setup = st.text_input(
                    "Environment Setup",
                    value=machine.env_setup if machine and hasattr(machine, 'env_setup') else "",
                    help="Shell commands to setup environment (e.g., 'source /etc/profile')"
                )
            
            # Scheduler resources
            if scheduler != "direct":
                st.subheader("Scheduler Resources")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    nodes = st.number_input("Nodes", min_value=1, value=1)
                with col2:
                    ntasks = st.number_input("Tasks per Node", min_value=1, value=20)
                with col3:
                    time = st.text_input("Wall Time", value="24:00:00")
                
                partition = st.text_input("Partition/Queue", value="")
            
            # Submit buttons
            col1, col2 = st.columns(2)
            with col1:
                submit = st.form_submit_button("💾 Save Machine Configuration")
            with col2:
                test = st.form_submit_button("🔍 Test Connection")
        
        # Handle form submission
        if submit:
            try:
                # Build machine config
                machine_config = {
                    "name": machine_name,
                    "execution": execution,
                    "scheduler": scheduler,
                    "workdir": workdir,
                    "nprocs": nprocs,
                    "launcher": launcher,
                    "use_modules": use_modules,
                }
                
                if use_modules:
                    machine_config["modules"] = [m.strip() for m in modules_str.split("\n") if m.strip()]
                
                if prepend:
                    machine_config["prepend"] = [p.strip() for p in prepend.split("\n") if p.strip()]
                if postpend:
                    machine_config["postpend"] = [p.strip() for p in postpend.split("\n") if p.strip()]
                if env_setup:
                    machine_config["env_setup"] = env_setup
                
                if execution == "remote":
                    machine_config["host"] = host
                    machine_config["username"] = username
                    machine_config["port"] = port
                    machine_config["auth"] = {
                        "method": "key",
                        "ssh_key": ssh_key
                    }
                
                if scheduler != "direct":
                    machine_config["resources"] = {
                        "nodes": nodes,
                        "ntasks-per-node": ntasks,
                        "time": time,
                    }
                    if partition:
                        machine_config["resources"]["partition"] = partition
                
                # Create Machine object
                new_machine = Machine(**machine_config)
                
                # Save machine
                save_machine(new_machine, DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
                
                st.success(f"✅ Machine '{machine_name}' saved successfully!")
                st.info(f"💾 Configuration saved to: ~/.xespresso/machines/{machine_name}.json")
                st.session_state.current_machine = new_machine
                st.session_state.current_machine_name = machine_name
                
            except Exception as e:
                st.error(f"❌ Error saving machine: {e}")
                st.code(traceback.format_exc())
        
        # Handle test connection
        if test:
            st.subheader("Connection Test Results")
            try:
                # Build minimal machine config for testing
                test_config = {
                    "name": machine_name,
                    "execution": execution,
                    "scheduler": "direct",  # Use direct for testing
                    "workdir": workdir,
                    "nprocs": 1,
                }
                
                if execution == "remote":
                    test_config["host"] = host
                    test_config["username"] = username
                    test_config["port"] = port
                    test_config["auth"] = {
                        "method": "key",
                        "ssh_key": ssh_key
                    }
                
                test_machine = Machine(**test_config)
                
                # Test connection
                if execution == "local":
                    st.success("✅ Local machine - connection OK")
                    st.info(f"Working directory: {workdir}")
                    st.info(f"Current user: {os.environ.get('USER', 'unknown')}")
                else:
                    # Test remote connection
                    with st.spinner("Testing SSH connection..."):
                        try:
                            # Try to establish SSH connection
                            import paramiko
                            ssh = paramiko.SSHClient()
                            # Load known hosts for security
                            try:
                                ssh.load_system_host_keys()
                            except Exception:
                                pass  # Known hosts file may not exist
                            
                            # Use WarningPolicy - warns but allows connection for testing
                            # Note: For production, use RejectPolicy and pre-configure host keys
                            ssh.set_missing_host_key_policy(paramiko.WarningPolicy())
                            
                            # Expand ssh_key path and validate
                            key_path = os.path.expanduser(ssh_key)
                            
                            if not os.path.isfile(key_path):
                                st.error(f"❌ SSH key not found: {key_path}")
                                st.info("💡 Check the SSH key path")
                            else:
                                ssh.connect(
                                    hostname=host,
                                    username=username,
                                    port=port,
                                    key_filename=key_path,
                                    timeout=10
                                )
                            
                            # Test command execution
                            stdin, stdout, stderr = ssh.exec_command('echo "Connection test successful"')
                            output = stdout.read().decode().strip()
                            
                            ssh.close()
                            
                            st.success(f"✅ SSH connection successful!")
                            st.info(f"Connected to: {username}@{host}:{port}")
                            st.info(f"Test output: {output}")
                        except paramiko.AuthenticationException:
                            st.error("❌ Authentication failed. Check username and SSH key.")
                        except paramiko.SSHException as e:
                            st.error(f"❌ SSH error: {e}")
                        except FileNotFoundError:
                            st.error(f"❌ SSH key not found: {key_path}")
                        except Exception as e:
                            st.error(f"❌ Connection failed: {e}")
                            
            except Exception as e:
                st.error(f"❌ Test failed: {e}")
                st.code(traceback.format_exc())

# Page 2: Codes Configuration
elif page == "⚙️ Codes Configuration":
    st.header("Quantum ESPRESSO Codes Configuration")
    st.markdown("""
    Configure Quantum ESPRESSO executable paths for different machines.
    Auto-detection is supported for both local and remote systems.
    """)
    
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available. Cannot configure codes.")
    else:
        # Machine selection
        try:
            machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
            if machines_list:
                selected_machine = st.selectbox(
                    "Select Machine:",
                    machines_list,
                    help="Choose the machine to configure codes for"
                )
            else:
                st.warning("⚠️ No machines configured. Please configure a machine first.")
                selected_machine = None
        except Exception as e:
            st.warning(f"Could not load machines: {e}")
            selected_machine = None
        
        if selected_machine:
            st.subheader(f"Codes Configuration for: {selected_machine}")
            
            # Auto-detection section
            st.subheader("Auto-Detect Codes")
            
            with st.form("detect_codes_form"):
                col1, col2 = st.columns(2)
                
                with col1:
                    qe_prefix = st.text_input(
                        "QE Installation Prefix (optional)",
                        help="e.g., /opt/qe-7.2/bin"
                    )
                    version_label = st.text_input(
                        "Version Label (optional)",
                        help="Custom label for this version (e.g., 'qe-7.2', 'qe-dev')"
                    )
                    modules_str = st.text_area(
                        "Modules to Load (optional, one per line)",
                        help="Version-specific modules (e.g., 'qe/7.2' or 'quantum_espresso-7.4.1')"
                    )
                
                with col2:
                    search_paths_str = st.text_area(
                        "Additional Search Paths (optional, one per line)",
                        help="Additional directories to search for executables"
                    )
                
                detect_button = st.form_submit_button("🔍 Auto-Detect Codes")
            
            if detect_button:
                with st.spinner("Detecting Quantum ESPRESSO codes..."):
                    try:
                        modules = [m.strip() for m in modules_str.split("\n") if m.strip()] if modules_str else None
                        search_paths = [p.strip() for p in search_paths_str.split("\n") if p.strip()] if search_paths_str else None
                        
                        codes_config = detect_qe_codes(
                            machine_name=selected_machine,
                            qe_prefix=qe_prefix if qe_prefix else None,
                            search_paths=search_paths,
                            modules=modules,
                            auto_load_machine=True
                        )
                        
                        # Check if any codes were detected (in main codes dict or versions structure)
                        if codes_config and codes_config.has_any_codes():
                            # Get all detected codes to count them
                            all_codes = codes_config.get_all_codes()
                            st.success(f"✅ Detected {len(all_codes)} codes!")
                            
                            # Add version label if provided
                            if version_label:
                                codes_config.version_label = version_label
                            
                            st.session_state.current_codes = codes_config
                            
                            # Display detected codes
                            st.subheader("Detected Codes")
                            codes_data = []
                            for name, code in all_codes.items():
                                codes_data.append({
                                    "Code": name,
                                    "Path": code.path,
                                    "Version": code.version or "Unknown",
                                    "Label": version_label or "default"
                                })
                            st.table(codes_data)
                            
                            # Save option with clear explanation
                            st.info("""
                            **💾 Saving Codes:**
                            - Detected codes will be **merged** with existing configurations
                            - Multiple versions on the same machine are supported
                            - Existing codes with different paths/versions will be kept
                            """)
                            
                            if st.button("💾 Save Codes Configuration"):
                                try:
                                    from xespresso.codes.manager import CodesManager
                                    filepath = CodesManager.save_config(
                                        codes_config,
                                        output_dir=DEFAULT_CODES_DIR,
                                        overwrite=False,
                                        merge=True
                                    )
                                    st.success(f"✅ Codes saved to: {filepath}")
                                    st.info("Multiple versions are preserved. Reload the page to see all versions.")
                                except Exception as e:
                                    st.error(f"Error saving codes: {e}")
                                    st.code(traceback.format_exc())
                        else:
                            st.warning("⚠️ No codes detected. Check paths and modules.")
                    except Exception as e:
                        st.error(f"❌ Error detecting codes: {e}")
                        st.code(traceback.format_exc())
            
            # Load existing configuration
            st.subheader("Existing Codes Configuration")
            try:
                existing_codes = load_codes_config(selected_machine, DEFAULT_CODES_DIR)
                if existing_codes:
                    st.success(f"✅ Loaded existing configuration")
                    
                    codes_data = []
                    for name, code in existing_codes.codes.items():
                        codes_data.append({
                            "Code": name,
                            "Path": code.path,
                            "Version": code.version or "Unknown",
                            "Modules": ", ".join(code.modules) if hasattr(code, 'modules') and code.modules else "None"
                        })
                    st.table(codes_data)
                    
                    st.session_state.current_codes = existing_codes
                    
                    # Code/version selection for calculations
                    st.subheader("Select Code Version for Calculations")
                    if existing_codes.codes:
                        code_options = list(existing_codes.codes.keys())
                        selected_code = st.selectbox(
                            "Select QE version to use:",
                            code_options,
                            help="Choose which version of QE to use for your calculations"
                        )
                        st.session_state.selected_code_version = selected_code
                        
                        selected_code_obj = existing_codes.codes[selected_code]
                        st.info(f"""
                        **Selected Code Details:**
                        - Path: `{selected_code_obj.path}`
                        - Version: {selected_code_obj.version or 'Unknown'}
                        """)
                else:
                    st.info("ℹ️ No codes configuration found for this machine.")
            except Exception as e:
                st.warning(f"Could not load codes configuration: {e}")

# Page 3: Structure Viewer
elif page == "🔬 Structure Viewer":
    st.header("Structure Viewer")
    st.markdown("""
    Upload or select molecular/crystal structures for your calculations.
    Supports CIF, POSCAR, XYZ, and other ASE-compatible formats.
    """)
    
    # Show currently loaded structure if available
    if st.session_state.current_structure is not None:
        st.success(f"✅ Current structure: {st.session_state.current_structure.get_chemical_formula()} ({len(st.session_state.current_structure)} atoms)")
    
    if not ASE_AVAILABLE:
        st.error("ASE not available. Structure viewing is disabled.")
    else:
        # Structure source selection
        structure_source = st.radio(
            "Structure Source:",
            ["Upload File", "Build Structure", "Load from File", "ASE Database"]
        )
        
        atoms = None
        
        if structure_source == "Upload File":
            uploaded_file = st.file_uploader(
                "Upload Structure File",
                type=['cif', 'xyz', 'pdb', 'poscar', 'vasp', 'traj'],
                help="Supported formats: CIF, XYZ, PDB, POSCAR, etc."
            )
            
            if uploaded_file is not None:
                try:
                    # Save to temporary file
                    with tempfile.NamedTemporaryFile(delete=False, suffix=Path(uploaded_file.name).suffix) as tmp_file:
                        tmp_file.write(uploaded_file.getvalue())
                        tmp_path = tmp_file.name
                    
                    # Read structure
                    atoms = io.read(tmp_path)
                    st.success(f"✅ Loaded structure from {uploaded_file.name}")
                    os.unlink(tmp_path)
                    
                except Exception as e:
                    st.error(f"❌ Error reading structure: {e}")
        
        elif structure_source == "Build Structure":
            st.subheader("Build Simple Structure")
            
            build_type = st.selectbox(
                "Structure Type:",
                ["Bulk Crystal", "Molecule"]
            )
            
            if build_type == "Bulk Crystal":
                col1, col2 = st.columns(2)
                with col1:
                    element = st.text_input("Element", value="Fe")
                    crystal_structure = st.selectbox(
                        "Crystal Structure",
                        ["fcc", "bcc", "hcp", "diamond", "sc"]
                    )
                with col2:
                    a_param = st.number_input("Lattice Parameter (Å)", value=3.6, step=0.1)
                    cubic = st.checkbox("Cubic Cell", value=True)
                
                if st.button("Build Crystal"):
                    try:
                        atoms = bulk(
                            element,
                            crystal_structure,
                            a=a_param,
                            cubic=cubic
                        )
                        st.success(f"✅ Built {element} {crystal_structure} structure")
                    except Exception as e:
                        st.error(f"❌ Error building structure: {e}")
            
            else:  # Molecule
                molecule_name = st.text_input(
                    "Molecule Name",
                    value="H2O",
                    help="Common molecules: H2O, CO2, CH4, etc."
                )
                
                if st.button("Build Molecule"):
                    try:
                        atoms = molecule(molecule_name)
                        st.success(f"✅ Built {molecule_name} molecule")
                    except Exception as e:
                        st.error(f"❌ Error building molecule: {e}")
        
        elif structure_source == "Load from File":
            file_path = st.text_input(
                "File Path",
                help="Enter full path to structure file"
            )
            
            if st.button("Load File") and file_path:
                try:
                    atoms = io.read(file_path)
                    st.success(f"✅ Loaded structure from {file_path}")
                except Exception as e:
                    st.error(f"❌ Error loading file: {e}")
        
        else:  # ASE Database
            st.subheader("ASE Database")
            
            db_path = st.text_input(
                "Database Path",
                value=st.session_state.get('ase_db_path', os.path.expanduser("~/.xespresso/structures.db")),
                help="Path to ASE database file"
            )
            
            # Validate database path
            is_valid, normalized_db_path, error_msg = validate_path(db_path, allow_creation=True)
            if not is_valid:
                st.error(f"❌ Invalid database path: {error_msg}")
                normalized_db_path = None
            else:
                st.session_state['ase_db_path'] = normalized_db_path
            
            # Database operations
            db_operation = st.radio(
                "Operation:",
                ["Load from Database", "Save to Database"]
            )
            
            if normalized_db_path and is_valid:
                if db_operation == "Load from Database":
                    if os.path.exists(normalized_db_path):
                        try:
                            from ase.db import connect
                            db = connect(normalized_db_path)
                            
                            # List structures in database
                            rows = list(db.select())
                            if rows:
                                st.write(f"Found {len(rows)} structures in database")
                            
                            # Create selection table
                            structures_info = []
                            for row in rows:
                                structures_info.append({
                                    "ID": row.id,
                                    "Formula": row.formula,
                                    "Atoms": row.natoms,
                                    "Tags": ", ".join(row.key_value_pairs.keys()) if row.key_value_pairs else ""
                                })
                            
                            st.table(structures_info)
                            
                            selected_id = st.number_input(
                                "Select structure ID to load:",
                                min_value=1,
                                max_value=len(rows),
                                value=1
                            )
                            
                            if st.button("Load Selected Structure"):
                                try:
                                    row = db.get(id=selected_id)
                                    atoms = row.toatoms()
                                    st.success(f"✅ Loaded structure ID {selected_id}: {atoms.get_chemical_formula()}")
                                except Exception as e:
                                    st.error(f"❌ Error loading structure: {e}")
                            else:
                                st.info("Database is empty. Save structures to start building your library.")
                        except Exception as e:
                            st.error(f"❌ Error reading database: {e}")
                    else:
                        st.info(f"Database does not exist yet. It will be created when you save your first structure.")
                
                else:  # Save to Database
                    if st.session_state.current_structure is not None:
                        current_atoms = st.session_state.current_structure
                        st.info(f"Ready to save: {current_atoms.get_chemical_formula()} ({len(current_atoms)} atoms)")
                        
                        # Add metadata
                        save_tags = st.text_input(
                            "Tags (comma-separated)",
                            help="Add tags to help identify this structure later"
                        )
                        
                        save_description = st.text_area(
                            "Description (optional)",
                            help="Add notes about this structure"
                        )
                        
                        if st.button("💾 Save to Database"):
                            try:
                                from ase.db import connect
                                db = connect(normalized_db_path)
                                
                                # Parse tags
                                key_value_pairs = {}
                                if save_tags:
                                    for tag in save_tags.split(','):
                                        tag = tag.strip()
                                        if tag:
                                            key_value_pairs[tag] = True
                                
                                if save_description:
                                    key_value_pairs['description'] = save_description
                                
                                # Save to database
                                db.write(current_atoms, **key_value_pairs)
                                st.success(f"✅ Structure saved to database: {normalized_db_path}")
                            except Exception as e:
                                st.error(f"❌ Error saving to database: {e}")
                    else:
                        st.warning("⚠️ No structure loaded. Load a structure first before saving to database.")
        
        # Display structure if loaded (or show current structure)
        display_atoms = atoms if atoms is not None else st.session_state.current_structure
        
        if display_atoms is not None:
            # Update session state if new structure was loaded
            if atoms is not None:
                st.session_state.current_structure = atoms
            
            # Display structure info
            display_structure_info(display_atoms)
            
            # 3D Visualization
            st.subheader("3D Visualization")
            if PLOTLY_AVAILABLE:
                fig = create_3d_structure_plot(display_atoms)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("⚠️ Plotly not available. Install plotly for 3D visualization.")
            
            # Export structure
            st.subheader("Export Structure")
            col1, col2 = st.columns(2)
            
            with col1:
                export_format = st.selectbox(
                    "Export Format",
                    ["cif", "xyz", "poscar", "pdb"]
                )
            
            with col2:
                export_filename = st.text_input(
                    "Filename",
                    value=f"structure.{export_format}"
                )
            
            if st.button("💾 Export Structure"):
                try:
                    export_path = os.path.join(tempfile.gettempdir(), export_filename)
                    io.write(export_path, display_atoms, format=export_format)
                    
                    with open(export_path, 'r') as f:
                        file_content = f.read()
                    
                    st.download_button(
                        "📥 Download File",
                        file_content,
                        file_name=export_filename,
                        mime="text/plain"
                    )
                    
                    st.success(f"✅ Structure exported as {export_format}")
                except Exception as e:
                    st.error(f"❌ Error exporting: {e}")

# Page 4: Calculation Setup
elif page == "📊 Calculation Setup":
    st.header("Calculation Setup")
    st.markdown("""
    Configure the type of calculation and basic parameters.
    """)
    
    if st.session_state.current_structure is None:
        st.warning("⚠️ No structure loaded. Please load a structure first in the Structure Viewer.")
    else:
        st.info(f"✅ Working with: {st.session_state.current_structure.get_chemical_formula()}")
        
        # Calculation type
        calc_type = st.selectbox(
            "Calculation Type",
            [
                "SCF (Self-Consistent Field)",
                "Relaxation (Geometry Optimization)",
                "VC-Relax (Cell + Geometry Optimization)",
                "Bands (Band Structure)",
                "DOS (Density of States)",
                "NSCF (Non-Self-Consistent)",
                "Phonon",
                "NEB (Nudged Elastic Band)"
            ]
        )
        
        st.session_state.workflow_config['calc_type'] = calc_type.split()[0].lower()
        
        # Pseudopotentials
        st.subheader("Pseudopotentials")
        
        atoms = st.session_state.current_structure
        unique_elements = list(set(atoms.get_chemical_symbols()))
        
        st.write(f"**Elements in structure:** {', '.join(unique_elements)}")
        
        pseudo_method = st.radio(
            "Pseudopotential Selection:",
            ["Manual Entry", "Load Configuration"]
        )
        
        pseudopotentials = {}
        
        if pseudo_method == "Manual Entry":
            st.write("Enter pseudopotential file for each element:")
            for element in unique_elements:
                pseudo = st.text_input(
                    f"{element}",
                    value=f"{element}.pbe-n-kjpaw_psl.1.0.0.UPF",
                    key=f"pseudo_{element}"
                )
                pseudopotentials[element] = pseudo
        else:
            try:
                from xespresso.utils import list_pseudo_configs, load_pseudo_config
                configs = list_pseudo_configs()
                
                if configs:
                    selected_config = st.selectbox(
                        "Select Configuration:",
                        configs
                    )
                    
                    if st.button("Load Configuration"):
                        config = load_pseudo_config(selected_config)
                        pseudopotentials = config.get('pseudopotentials', {})
                        st.success(f"✅ Loaded pseudopotentials from {selected_config}")
                        st.json(pseudopotentials)
                else:
                    st.warning("No saved pseudopotential configurations found.")
            except Exception as e:
                st.error(f"Error loading configurations: {e}")
        
        st.session_state.workflow_config['pseudopotentials'] = pseudopotentials
        
        # Basic parameters
        st.subheader("Calculation Parameters")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Use session state for persistence
            default_ecutwfc = st.session_state.workflow_config.get('ecutwfc', 50.0)
            ecutwfc = st.number_input(
                "Kinetic Energy Cutoff (ecutwfc) [Ry]",
                min_value=10.0,
                max_value=200.0,
                value=float(default_ecutwfc),
                step=5.0,
                help="Plane-wave cutoff energy"
            )
            
            # Add dual parameter
            default_dual = st.session_state.workflow_config.get('dual', 4.0)
            dual = st.number_input(
                "Dual Parameter (ecutrho/ecutwfc ratio)",
                min_value=1.0,
                max_value=12.0,
                value=float(default_dual),
                step=0.5,
                help="Ratio between charge density and wavefunction cutoffs (typically 4-8)"
            )
            
            if calc_type.startswith("SCF") or calc_type.startswith("Relaxation"):
                default_conv_thr = st.session_state.workflow_config.get('conv_thr', 1e-6)
                conv_thr = st.number_input(
                    "Convergence Threshold",
                    min_value=1e-10,
                    max_value=1e-4,
                    value=float(default_conv_thr),
                    format="%.1e",
                    help="SCF convergence threshold"
                )
        
        with col2:
            # Calculate ecutrho from dual
            ecutrho = ecutwfc * dual
            st.number_input(
                "Charge Density Cutoff (ecutrho) [Ry]",
                min_value=10.0,
                max_value=800.0,
                value=ecutrho,
                step=10.0,
                help="Charge density cutoff = dual × ecutwfc",
                disabled=True
            )
        
        st.session_state.workflow_config.update({
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
            'dual': dual,
        })
        
        if calc_type.startswith("SCF") or calc_type.startswith("Relaxation"):
            st.session_state.workflow_config['conv_thr'] = conv_thr
        
        # Smearing options
        st.subheader("Electronic Occupations")
        
        col1, col2 = st.columns(2)
        with col1:
            occupations = st.selectbox(
                "Occupation Type",
                ["smearing", "fixed", "tetrahedra"],
                index=["smearing", "fixed", "tetrahedra"].index(
                    st.session_state.workflow_config.get('occupations', 'smearing')
                ),
                help="Method for determining electronic occupations"
            )
            st.session_state.workflow_config['occupations'] = occupations
        
        with col2:
            if occupations == "smearing":
                smearing_type = st.selectbox(
                    "Smearing Type",
                    ["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
                    index=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"].index(
                        st.session_state.workflow_config.get('smearing', 'gaussian')
                    ),
                    help="Type of smearing function"
                )
                st.session_state.workflow_config['smearing'] = smearing_type
                
                degauss = st.number_input(
                    "Smearing Width (degauss) [Ry]",
                    min_value=0.001,
                    max_value=0.1,
                    value=st.session_state.workflow_config.get('degauss', 0.02),
                    step=0.001,
                    format="%.3f",
                    help="Width of smearing (typically 0.01-0.03 Ry)"
                )
                st.session_state.workflow_config['degauss'] = degauss
        
        # K-points
        st.subheader("K-point Sampling")
        
        kpt_method = st.radio(
            "K-point Method:",
            ["K-spacing", "Monkhorst-Pack Grid"]
        )
        
        if kpt_method == "K-spacing":
            kspacing = st.slider(
                "K-spacing (Å⁻¹)",
                min_value=0.1,
                max_value=1.0,
                value=0.3,
                step=0.05,
                help="Smaller values = denser k-point mesh"
            )
            st.session_state.workflow_config['kspacing'] = kspacing
            
            # Show equivalent grid
            try:
                from xespresso import kpts_from_spacing
                kpts = kpts_from_spacing(atoms, kspacing)
                st.info(f"Equivalent Monkhorst-Pack grid: {kpts[0]} × {kpts[1]} × {kpts[2]}")
            except:
                pass
        else:
            col1, col2, col3 = st.columns(3)
            with col1:
                k1 = st.number_input("k₁", min_value=1, value=4)
            with col2:
                k2 = st.number_input("k₂", min_value=1, value=4)
            with col3:
                k3 = st.number_input("k₃", min_value=1, value=4)
            
            st.session_state.workflow_config['kpts'] = (k1, k2, k3)
        
        # Spin polarization
        st.subheader("Spin Polarization")
        default_nspin = st.session_state.workflow_config.get('nspin', 1)
        nspin = st.selectbox(
            "Spin Treatment",
            [1, 2, 4],
            index=[1, 2, 4].index(default_nspin),
            format_func=lambda x: {
                1: "Non-spin-polarized",
                2: "Spin-polarized (collinear)",
                4: "Non-collinear + spin-orbit"
            }[x],
            help="Spin treatment for magnetic systems"
        )
        st.session_state.workflow_config['nspin'] = nspin
        
        # DFT+U section
        st.subheader("DFT+U Configuration")
        use_dft_u = st.checkbox(
            "Enable DFT+U",
            value=st.session_state.workflow_config.get('use_dft_u', False),
            help="Add Hubbard U correction for strongly correlated systems"
        )
        st.session_state.workflow_config['use_dft_u'] = use_dft_u
        
        if use_dft_u:
            st.info("Configure Hubbard U parameters for each element")
            
            atoms = st.session_state.current_structure
            unique_elements = list(set(atoms.get_chemical_symbols()))
            
            hubbard_u = st.session_state.workflow_config.get('hubbard_u', {})
            
            for element in unique_elements:
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.write(f"**{element}**")
                with col2:
                    u_val = st.number_input(
                        f"U value (eV)",
                        min_value=0.0,
                        max_value=10.0,
                        value=hubbard_u.get(element, {}).get('U', 0.0),
                        step=0.5,
                        key=f"hubbard_u_{element}"
                    )
                with col3:
                    orbital = st.selectbox(
                        f"Orbital",
                        ["2p", "3d", "4f"],
                        index=["2p", "3d", "4f"].index(hubbard_u.get(element, {}).get('orbital', '3d')),
                        key=f"hubbard_orbital_{element}"
                    )
                
                if u_val > 0:
                    hubbard_u[element] = {'U': u_val, 'orbital': orbital}
            
            st.session_state.workflow_config['hubbard_u'] = hubbard_u
        
        # Calculation-specific options
        if calc_type.startswith("Bands"):
            st.subheader("Band Structure Settings")
            
            band_path_method = st.radio(
                "K-path Selection:",
                ["Automatic (seekpath)", "Custom Path"]
            )
            
            if band_path_method == "Automatic (seekpath)":
                st.info("Will use automatic k-path detection based on crystal symmetry")
                st.session_state.workflow_config['band_path'] = 'auto'
            else:
                st.write("Define custom k-path (e.g., 'GXMGRX' for cubic systems)")
                custom_path = st.text_input(
                    "K-path",
                    value=st.session_state.workflow_config.get('custom_band_path', 'GXMGRX'),
                    help="Specify high-symmetry points"
                )
                st.session_state.workflow_config['band_path'] = 'custom'
                st.session_state.workflow_config['custom_band_path'] = custom_path
                
            npoints = st.number_input(
                "Number of k-points along path",
                min_value=10,
                max_value=500,
                value=st.session_state.workflow_config.get('band_npoints', 100),
                help="Total number of k-points along the band path"
            )
            st.session_state.workflow_config['band_npoints'] = npoints
        
        st.success("✅ Calculation parameters configured!")

# Page 5: Workflow Builder
elif page == "🔄 Workflow Builder":
    st.header("Workflow Builder")
    st.markdown("""
    Build complete calculation workflows using quality presets.
    """)
    
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available.")
    elif st.session_state.current_structure is None:
        st.warning("⚠️ No structure loaded. Please load a structure first.")
    else:
        st.info(f"✅ Working with: {st.session_state.current_structure.get_chemical_formula()}")
        
        # Quality presets
        st.subheader("Quality Presets")
        
        quality = st.select_slider(
            "Quality Level",
            options=["fast", "moderate", "accurate"],
            value="moderate",
            help="Predefined parameter sets for different accuracy/speed tradeoffs"
        )
        
        # Display preset info
        if XESPRESSO_AVAILABLE:
            preset_info = PRESETS.get(quality, {})
            st.info(f"""
            **{quality.upper()} preset:**
            - ecutwfc: {preset_info.get('ecutwfc', 'N/A')} Ry
            - ecutrho: {preset_info.get('ecutrho', 'N/A')} Ry
            - conv_thr: {preset_info.get('conv_thr', 'N/A')}
            - Default k-spacing: {preset_info.get('kspacing', 'N/A')} Å⁻¹
            """)
        
        st.session_state.workflow_config['quality'] = quality
        
        # Workflow type
        st.subheader("Workflow Configuration")
        
        workflow_type = st.selectbox(
            "Workflow Type",
            ["Quick SCF", "Quick Relax", "Custom Workflow"]
        )
        
        if workflow_type == "Quick Relax":
            relax_type = st.selectbox(
                "Relaxation Type",
                ["relax", "vc-relax"],
                format_func=lambda x: {
                    "relax": "Relax atoms only",
                    "vc-relax": "Relax atoms + cell"
                }[x]
            )
            st.session_state.workflow_config['relax_type'] = relax_type
        
        # Label/directory
        st.subheader("Calculation Settings")
        
        label = st.text_input(
            "Calculation Label",
            value="calc/structure",
            help="Directory path for calculation files"
        )
        st.session_state.workflow_config['label'] = label
        
        # Summary
        st.subheader("Workflow Summary")
        
        st.json(st.session_state.workflow_config)
        
        if st.button("✅ Create Workflow"):
            st.success("✅ Workflow configured! Proceed to Job Submission to run the calculation.")

# Page 6: Job Submission
elif page == "🚀 Job Submission":
    st.header("Job Submission")
    st.markdown("""
    Submit your configured calculation to the selected machine.
    """)
    
    # Check prerequisites
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available.")
    elif st.session_state.current_structure is None:
        st.warning("⚠️ No structure loaded.")
    elif not st.session_state.workflow_config:
        st.warning("⚠️ No workflow configured.")
    else:
        # Display configuration summary
        st.subheader("Configuration Summary")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Structure:**")
            atoms = st.session_state.current_structure
            st.write(f"- Formula: {atoms.get_chemical_formula()}")
            st.write(f"- Atoms: {len(atoms)}")
            
            st.write("**Machine:**")
            if st.session_state.current_machine:
                machine = st.session_state.current_machine
                st.write(f"- Name: {machine.name}")
                st.write(f"- Type: {machine.execution}")
            else:
                st.write("- Not configured")
        
        with col2:
            st.write("**Workflow:**")
            config = st.session_state.workflow_config
            st.write(f"- Quality: {config.get('quality', 'N/A')}")
            st.write(f"- Type: {config.get('calc_type', 'N/A')}")
            st.write(f"- Label: {config.get('label', 'N/A')}")
            
            st.write("**Codes:**")
            if st.session_state.current_codes:
                codes = st.session_state.current_codes
                st.write(f"- Configured: {len(codes.codes)} codes")
            else:
                st.write("- Not configured")
        
        # Submission options
        st.subheader("Submission Options")
        
        # Local working directory selection
        st.write("**Local Working Directory:**")
        col1, col2 = st.columns([3, 1])
        with col1:
            local_workdir = st.text_input(
                "Working Directory",
                value=st.session_state.local_workdir,
                help="Directory where calculation files will be saved locally"
            )
            st.session_state.local_workdir = local_workdir
        with col2:
            if st.button("📁 Use Current"):
                st.session_state.local_workdir = os.getcwd()
                st.rerun()
        
        # Show where files will be saved  
        # Validate local workdir path
        is_valid_workdir, normalized_workdir, error_msg = validate_path(local_workdir, allow_creation=True)
        if not is_valid_workdir:
            st.error(f"❌ Invalid working directory: {error_msg}")
        else:
            st.info(f"📂 Files will be saved to: `{normalized_workdir}`")
            
            # Create directory if it doesn't exist
            if not os.path.exists(normalized_workdir):
                if st.checkbox("Create directory if it doesn't exist", value=True):
                    try:
                        os.makedirs(normalized_workdir, exist_ok=True)
                        st.success(f"✅ Directory created: {normalized_workdir}")
                    except Exception as e:
                        st.error(f"❌ Could not create directory: {e}")
        
        dry_run = st.checkbox(
            "Dry Run (don't actually submit)",
            value=True,
            help="Generate input files without running"
        )
        
        # Submit button
        if st.button("🚀 Submit Job", type="primary"):
            with st.spinner("Submitting job..."):
                try:
                    st.info("📝 Job submission functionality will be implemented")
                    st.info("This will create input files and submit to the configured machine")
                    
                    # Show what would be done
                    st.subheader("Submission Details")
                    st.write("**Steps that would be performed:**")
                    st.write("1. ✓ Create calculation directory")
                    st.write("2. ✓ Write structure file")
                    st.write("3. ✓ Generate Quantum ESPRESSO input")
                    st.write("4. ✓ Submit to scheduler (if configured)")
                    st.write("5. ✓ Monitor job status")
                    
                    # Show file locations (use normalized path if valid)
                    display_workdir = normalized_workdir if is_valid_workdir else local_workdir
                    st.subheader("📂 File Locations")
                    st.info(f"""
                    **Input files:** `{display_workdir}/`
                    - Structure file: `{display_workdir}/structure.cif`
                    - QE input: `{display_workdir}/espresso.pwi`
                    - Job script: `{display_workdir}/run.sh`
                    """)
                    
                    if dry_run:
                        st.success("✅ Dry run completed - no job submitted")
                    else:
                        st.success("✅ Job submitted successfully!")
                        st.balloons()
                    
                except Exception as e:
                    st.error(f"❌ Error submitting job: {e}")
                    st.code(traceback.format_exc())

# Page 7: Results & Post-Processing
elif page == "📈 Results & Post-Processing":
    st.header("Results & Post-Processing")
    st.markdown("""
    View calculation results, analyze output files, and perform post-processing.
    """)
    
    # Working directory selection
    st.subheader("Select Calculation Directory")
    
    results_dir = st.text_input(
        "Results Directory",
        value=st.session_state.local_workdir,
        help="Path to directory containing calculation results"
    )
    
    # Validate results directory path
    is_valid_results, normalized_results_dir, error_msg = validate_path(results_dir, allow_creation=False)
    
    if not is_valid_results:
        st.error(f"❌ Invalid results directory: {error_msg}")
    elif os.path.exists(normalized_results_dir):
        st.success(f"✅ Directory found: {normalized_results_dir}")
        
        # List output files
        st.subheader("Output Files")
        
        try:
            files = os.listdir(normalized_results_dir)
            output_files = [f for f in files if f.endswith(('.out', '.pwo', '.xml', '.log'))]
            
            if output_files:
                selected_file = st.selectbox(
                    "Select output file to view:",
                    output_files
                )
                
                # Validate selected filename (no path traversal)
                if '..' in selected_file or '/' in selected_file or '\\' in selected_file:
                    st.error("❌ Invalid filename")
                else:
                    file_path = os.path.join(normalized_results_dir, selected_file)
                
                    # Display file info
                    file_size = os.path.getsize(file_path)
                    st.info(f"File: {selected_file} | Size: {file_size / 1024:.2f} KB")
                    
                    # View file content
                    if st.button("📄 View File Content"):
                        try:
                            with open(file_path, 'r') as f:
                                content = f.read()
                            
                            # Show in expandable text area
                            with st.expander("File Content", expanded=True):
                                st.text_area(
                                    "Output",
                                    value=content,
                                    height=400,
                                    key="file_content"
                                )
                            
                            # Parse for key information
                            st.subheader("Extracted Information")
                            
                            # Simple parsing for common outputs
                            if "Final energy" in content or "!" in content:
                                st.write("**Energy Information:**")
                                for line in content.split('\n'):
                                    if "Final energy" in line or (line.strip().startswith("!") and "total energy" in line.lower()):
                                        st.code(line.strip())
                            
                            if "convergence has been achieved" in content.lower():
                                st.success("✅ Calculation converged successfully")
                            elif "convergence NOT achieved" in content.lower():
                                st.warning("⚠️ Calculation did not converge")
                            
                        except Exception as e:
                            st.error(f"❌ Error reading file: {e}")
                    
                    # Download button
                    try:
                        with open(file_path, 'r') as f:
                            file_content = f.read()
                        
                        st.download_button(
                            "📥 Download Output File",
                            file_content,
                            file_name=selected_file,
                            mime="text/plain"
                        )
                    except Exception as e:
                        st.error(f"❌ Error preparing download: {e}")
            else:
                st.warning("⚠️ No output files found in this directory.")
        
        except Exception as e:
            st.error(f"❌ Error listing files: {e}")
        
        # Structure visualization from results
        st.subheader("Structure Visualization")
        
        try:
            structure_files = [f for f in files if f.endswith(('.cif', '.xyz', '.pdb', '.poscar', 'CONTCAR'))]
            
            if structure_files and ASE_AVAILABLE:
                selected_structure = st.selectbox(
                    "Select structure file:",
                    structure_files
                )
                
                # Validate structure filename
                if '..' in selected_structure or '/' in selected_structure or '\\' in selected_structure:
                    st.error("❌ Invalid structure filename")
                elif st.button("🔬 Visualize Structure"):
                    try:
                        struct_path = os.path.join(normalized_results_dir, selected_structure)
                        atoms = io.read(struct_path)
                        
                        st.success(f"✅ Loaded: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
                        
                        # Display structure info
                        display_structure_info(atoms)
                        
                        # 3D Visualization
                        if PLOTLY_AVAILABLE:
                            fig = create_3d_structure_plot(atoms)
                            if fig:
                                st.plotly_chart(fig, use_container_width=True)
                        
                    except Exception as e:
                        st.error(f"❌ Error visualizing structure: {e}")
            else:
                st.info("No structure files found for visualization.")
        
        except Exception as e:
            st.error(f"❌ Error searching for structure files: {e}")
        
        # Post-processing tools
        st.subheader("Post-Processing Tools")
        
        post_tool = st.selectbox(
            "Select Tool:",
            [
                "Energy Analysis",
                "DOS Plotting",
                "Band Structure Plotting",
                "Structure Comparison"
            ]
        )
        
        if post_tool == "Energy Analysis":
            st.info("📊 Energy analysis tools will extract and plot total energy convergence.")
            st.write("*Feature coming soon*")
        
        elif post_tool == "DOS Plotting":
            st.info("📈 DOS plotting tools will visualize density of states from dos.x output.")
            st.write("*Feature coming soon*")
        
        elif post_tool == "Band Structure Plotting":
            st.info("📉 Band structure plotting from bands.x output.")
            st.write("*Feature coming soon*")
        
        elif post_tool == "Structure Comparison":
            st.info("🔄 Compare initial and final structures from relaxation calculations.")
            st.write("*Feature coming soon*")
    
    else:
        st.error(f"❌ Directory not found: {results_dir}")
        st.info("Please check the path or complete a calculation first.")

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("""
### About
**xespresso GUI** - Streamlit interface for Quantum ESPRESSO calculations

Version: 1.0.0

[Documentation](https://github.com/superstar54/xespresso) | 
[Report Issue](https://github.com/superstar54/xespresso/issues)
""")

"""
Machine Configuration Page for xespresso GUI.

This module handles the machine configuration interface, allowing users to:
- Create and edit machine configurations
- Test connections to remote machines
- Save configurations for later use
"""

import streamlit as st
import os
import traceback

try:
    from xespresso.machines.machine import Machine
    from xespresso.machines.config.loader import (
        load_machine, save_machine, list_machines,
        DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
    )
    XESPRESSO_AVAILABLE = True
except ImportError:
    XESPRESSO_AVAILABLE = False

# SSH connection testing is optional (requires paramiko)
try:
    from xespresso.utils.auth import test_ssh_connection
    SSH_TEST_AVAILABLE = True
except ImportError:
    SSH_TEST_AVAILABLE = False
    def test_ssh_connection(username, host, key_path=None, port=22):
        """Fallback SSH connection test when paramiko is not available."""
        import subprocess
        import os
        
        key_path = os.path.expanduser(key_path) if key_path else None
        cmd = ["ssh", "-p", str(port), "-o", "PasswordAuthentication=no", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5"]
        if key_path:
            cmd += ["-i", key_path]
        cmd += [f"{username}@{host}", "echo 'Connection successful'"]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False


def render_machine_config_page():
    """Render the machine configuration page."""
    st.header("⚙️ Machine Configuration")
    st.markdown("""
    **Create and configure** computational machines/clusters for your calculations.
    
    This page is for **configuration only** - once machines are saved, you can select them 
    in the Calculation Setup or Workflow Builder pages.
    
    Supports both local and remote (SSH) execution environments.
    """)
    
    st.info("""
    💡 **Configuration vs. Selection:**
    - **Configure** machines here (one-time setup or updates)
    - **Select** configured machines in Calculation Setup or Workflow Builder
    - Configurations are saved to `~/.xespresso/machines/`
    """)
    
    if not XESPRESSO_AVAILABLE:
        st.error("xespresso modules not available. Cannot configure machines.")
        return
    
    # List existing machines
    st.subheader("Existing Machines")
    try:
        machines_list = list_machines(DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR)
        if machines_list:
            st.success(f"✅ {len(machines_list)} machine(s) configured: {', '.join(machines_list)}")
            
            # Use session state for persistent selection
            default_idx = 0
            if st.session_state.current_machine_name and st.session_state.current_machine_name in machines_list:
                default_idx = machines_list.index(st.session_state.current_machine_name) + 1
            
            selected_machine = st.selectbox(
                "Select a machine to edit:",
                ["[Create New Machine]"] + machines_list,
                index=default_idx,
                key="machine_selector",
                help="Choose an existing machine to edit or create a new one"
            )
            
            # Update session state
            if selected_machine != "[Create New Machine]":
                st.session_state.current_machine_name = selected_machine
        else:
            st.info("No machines configured yet. Create your first machine below.")
            selected_machine = "[Create New Machine]"
    except Exception as e:
        st.warning(f"Could not load machines list: {e}")
        selected_machine = "[Create New Machine]"
    
    # Create or edit machine
    st.subheader("Machine Configuration Form")
    
    # Load existing machine if selected
    if selected_machine != "[Create New Machine]":
        try:
            machine = load_machine(DEFAULT_CONFIG_PATH, selected_machine, DEFAULT_MACHINES_DIR, return_object=True)
            st.success(f"✅ Editing machine: {selected_machine}")
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
        host = ""
        username = ""
        port = 22
        ssh_key = "~/.ssh/id_rsa"
        
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
        
        modules_str = ""
        if use_modules:
            modules_str = st.text_area(
                "Modules to Load (one per line)",
                value="\n".join(machine.modules) if machine and machine.modules else "",
                help="Environment modules to load before execution"
            )
        
        # Advanced settings
        prepend = ""
        postpend = ""
        env_setup = ""
        
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
        nodes = 1
        ntasks = 20
        time = "24:00:00"
        partition = ""
        
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
    
    # Handle test connection - now using xespresso's test_ssh_connection
    if test:
        st.subheader("Connection Test Results")
        try:
            if execution == "local":
                st.success("✅ Local machine - connection OK")
                st.info(f"Working directory: {workdir}")
                st.info(f"Current user: {os.environ.get('USER', 'unknown')}")
            else:
                # Test remote connection using xespresso's built-in function
                with st.spinner("Testing SSH connection..."):
                    key_path = os.path.expanduser(ssh_key)
                    
                    if not os.path.isfile(key_path):
                        st.error(f"❌ SSH key not found: {key_path}")
                        st.info("💡 Check the SSH key path")
                    else:
                        # Use xespresso's test_ssh_connection function
                        success = test_ssh_connection(username, host, key_path, port)
                        
                        if success:
                            st.success(f"✅ SSH connection successful!")
                            st.info(f"Connected to: {username}@{host}:{port}")
                        else:
                            st.error("❌ SSH connection failed.")
                            st.info("💡 Check your credentials and SSH key configuration")
                        
        except Exception as e:
            st.error(f"❌ Test failed: {e}")
            st.code(traceback.format_exc())

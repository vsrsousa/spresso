"""
Workflow Builder Page for xespresso GUI.

This page creates multi-step workflows using the GUI workflow modules,
which coordinate multiple calculations following xespresso's design patterns.
"""

import streamlit as st


def render_workflow_builder_page():
    """
    Render the workflow builder page.

    This page uses the GUI workflow modules to orchestrate multiple
    calculations, with each calculation prepared by calculation modules.
    """
    st.header("🔄 Workflow Builder")
    st.markdown(
        """
    Build multi-step workflows using **workflow modules** that orchestrate calculations.
    
    **Modular Design:**
    - Workflow modules (`gui/workflows/`) coordinate multiple calculations
    - Each calculation step uses calculation modules to prepare objects
    - Job submission executes the prepared workflow steps
    """
    )

    # Check if structure is loaded
    if (
        "current_structure" not in st.session_state
        or st.session_state.current_structure is None
    ):
        st.warning(
            "⚠️ No structure loaded. Please load a structure first in the Structure Viewer page."
        )
        return

    atoms = st.session_state.current_structure
    st.success(
        f"✅ Structure loaded: {atoms.get_chemical_formula()} ({len(atoms)} atoms)"
    )

    st.markdown("---")

    # Initialize workflow_config if needed
    if "workflow_config" not in st.session_state:
        st.session_state.workflow_config = {}

    config = st.session_state.workflow_config

    # ===== MOVED TO TOP: Machine and Code Selection =====
    st.subheader("🖥️ Execution Environment")
    st.info(
        """
    **First, select the machine and code version** for this workflow.
    These settings will be used for all steps in the workflow.
    """
    )

    col1, col2 = st.columns(2)

    with col1:
        # Machine Selection
        st.write("**Machine:**")
        try:
            from xespresso.machines.config.loader import (
                list_machines,
                load_machine,
                DEFAULT_CONFIG_PATH,
                DEFAULT_MACHINES_DIR,
            )

            available_machines = list_machines()

            if available_machines:
                # Determine default index for machine selector
                default_machine_idx = 0
                if (
                    st.session_state.get("selected_machine")
                    and st.session_state.selected_machine in available_machines
                ):
                    default_machine_idx = available_machines.index(
                        st.session_state.selected_machine
                    )

                selected_machine_name = st.selectbox(
                    "Select Machine:",
                    options=available_machines,
                    index=default_machine_idx,
                    help="Machine where the workflow will run",
                    key="workflow_machine_selector",
                )
                # Store selection in shared state variable for cross-page persistence
                st.session_state.selected_machine = selected_machine_name
                config["machine_name"] = selected_machine_name

                # Load the machine object
                try:
                    machine = load_machine(
                        DEFAULT_CONFIG_PATH,
                        selected_machine_name,
                        DEFAULT_MACHINES_DIR,
                        return_object=True,
                    )
                    st.session_state.workflow_machine = machine

                    # Show machine info
                    st.caption(f"Type: {machine.execution}")
                    if machine.scheduler:
                        st.caption(f"Scheduler: {machine.scheduler}")
                except Exception as e:
                    st.warning(f"Could not load machine: {e}")
            else:
                st.warning(
                    "⚠️ No machines configured. Please configure a machine first in the Machine Configuration page."
                )
                st.session_state.selected_machine = None
                config["machine_name"] = None
        except ImportError:
            st.error("❌ Machine configuration modules not available")
            st.session_state.selected_machine = None
            config["machine_name"] = None

    with col2:
        # Code/Version Selection using the proper selector
        st.write("**Code Version:**")
        if st.session_state.get("selected_machine"):
            try:
                from xespresso.codes.manager import load_codes_config, DEFAULT_CODES_DIR

                codes = load_codes_config(
                    st.session_state.selected_machine, DEFAULT_CODES_DIR
                )

                if codes and codes.has_any_codes():
                    # Check if versions are available
                    available_versions = codes.list_versions()

                    if codes.versions and available_versions:
                        # Show version selector (whether single or multiple versions)
                        if len(available_versions) > 1:
                            st.info(
                                f"📦 Multiple QE versions available: {', '.join(available_versions)}"
                            )

                        # Version selector
                        default_idx = 0
                        if (
                            st.session_state.get("selected_version")
                            and st.session_state.selected_version
                            in available_versions
                        ):
                            default_idx = available_versions.index(
                                st.session_state.selected_version
                            )

                        selected_version = st.selectbox(
                            "Select QE Version:",
                            available_versions,
                            index=default_idx,
                            key="workflow_version_selector",
                            help="Choose which Quantum ESPRESSO version to use for this workflow",
                        )

                        # Store selected version in shared state for cross-page persistence
                        st.session_state.selected_version = selected_version
                        config["qe_version"] = selected_version

                        # Get codes for selected version
                        version_codes = codes.get_all_codes(version=selected_version)

                        # Show version details and ALWAYS retrieve modules if defined
                        with st.expander("⚙️ Version Details", expanded=False):
                            st.write(f"**Version:** {selected_version}")
                            if codes.versions and selected_version in codes.versions:
                                version_config = codes.versions[selected_version]
                                if version_config.get("label"):
                                    st.write(f"**Label:** {version_config['label']}")
                                if version_config.get("qe_prefix"):
                                    st.write(
                                        f"**Prefix:** {version_config['qe_prefix']}"
                                    )
                                # ALWAYS show and store modules if they exist in codes JSON
                                if version_config.get("modules"):
                                    modules = version_config["modules"]
                                    st.write(f"**Modules:** {', '.join(modules)}")
                                    config["modules"] = modules
                            st.write(
                                f"**Available codes:** {', '.join(version_codes.keys())}"
                            )

                    else:
                        # Single version or no version structure
                        version_codes = codes.get_all_codes()
                        selected_version = codes.qe_version
                        if selected_version:
                            st.caption(f"QE Version: {selected_version}")
                            config["qe_version"] = selected_version

                    # Show available codes and allow selection
                    code_names = list(version_codes.keys())
                    if code_names:
                        st.caption(
                            f"✓ {len(code_names)} codes configured: {', '.join(code_names[:3])}{' ...' if len(code_names) > 3 else ''}"
                        )

                        # Individual code selection
                        st.markdown("**Select Code:**")
                        default_code_idx = 0
                        # Try to select 'pw' by default if available
                        if "pw" in code_names:
                            default_code_idx = code_names.index("pw")
                        elif (
                            st.session_state.get("selected_code")
                            and st.session_state.selected_code in code_names
                        ):
                            default_code_idx = code_names.index(
                                st.session_state.selected_code
                            )

                        selected_code = st.selectbox(
                            "Choose code executable:",
                            code_names,
                            index=default_code_idx,
                            key="workflow_code_selector",
                            help="Select which Quantum ESPRESSO executable to use (typically pw for workflows)",
                        )

                        # Store selected code in shared state for cross-page persistence
                        st.session_state.selected_code = selected_code
                        config["selected_code"] = selected_code

                        # Show code details
                        selected_code_obj = version_codes[selected_code]
                        st.caption(f"📍 Path: {selected_code_obj.path}")
                        if (
                            hasattr(selected_code_obj, "version")
                            and selected_code_obj.version
                        ):
                            st.caption(f"📦 Version: {selected_code_obj.version}")

                else:
                    st.warning(
                        f"⚠️ No codes configured for machine '{st.session_state.selected_machine}'. Please configure codes in the Codes Configuration page."
                    )
                    config["qe_version"] = None
                    config["selected_code"] = None
            except Exception as e:
                st.warning(f"Could not load codes: {e}")
                config["qe_version"] = None
        else:
            st.info("Select a machine first")
            config["qe_version"] = None
            config["selected_code"] = None

    st.markdown("---")

    # Workflow Type Selection
    st.subheader("🔧 Workflow Type")
    workflow_type = st.selectbox(
        "Select Workflow:",
        [
            "Single SCF",
            "SCF + Relaxation",
            "SCF + Relax + SCF (on relaxed structure)",
            "Custom Multi-Step",
        ],
        help="Choose a predefined workflow or build a custom one",
    )

    # Initialize workflow_config if needed
    if "workflow_config" not in st.session_state:
        st.session_state.workflow_config = {}

    config = st.session_state.workflow_config

    # Common Parameters
    st.subheader("⚙️ Common Parameters")

    col1, col2 = st.columns(2)
    with col1:
        ecutwfc = st.number_input(
            "Energy Cutoff (Ry):", value=50.0, min_value=10.0, max_value=200.0, step=5.0
        )
        config["ecutwfc"] = ecutwfc

    with col2:
        ecutrho = st.number_input(
            "Charge Density Cutoff (Ry):", value=ecutwfc * 8, min_value=40.0, step=20.0
        )
        config["ecutrho"] = ecutrho

    # K-points
    st.subheader("🔷 K-points")
    kpts_mode = st.radio(
        "K-points Mode:",
        ["K-spacing", "Explicit Grid"],
        horizontal=True,
        key="workflow_kpts_mode",
    )

    if kpts_mode == "K-spacing":
        # K-spacing mode: use slider to adjust density, then convert to kpts
        kspacing_value = st.slider(
            "K-spacing (Å⁻¹):",
            min_value=0.1,
            max_value=1.0,
            value=float(config.get("kspacing_ui", 0.3)),
            step=0.05,
            help="K-point density in reciprocal space. This will be converted to k-point grid.",
            key="workflow_kspacing_slider",
        )

        # Store the UI value for persistence (not in actual config)
        config["kspacing_ui"] = kspacing_value

        # Convert kspacing to kpts using the structure
        from xespresso import kpts_from_spacing

        computed_kpts = kpts_from_spacing(atoms, kspacing_value)

        # Store the computed kpts in config (not kspacing)
        config["kpts"] = computed_kpts

        # Display the computed k-points to the user
        st.info(
            f"ℹ️ Computed k-point grid: {computed_kpts[0]} × {computed_kpts[1]} × {computed_kpts[2]}"
        )

        # Remove kspacing from config as it should not be passed as a parameter
        if "kspacing" in config:
            del config["kspacing"]
    else:
        # Explicit Grid mode
        col1, col2, col3 = st.columns(3)
        with col1:
            k1 = st.number_input(
                "k₁:",
                value=config.get("kpts", (4, 4, 4))[0],
                min_value=1,
                max_value=20,
                key="workflow_k1",
            )
        with col2:
            k2 = st.number_input(
                "k₂:",
                value=config.get("kpts", (4, 4, 4))[1],
                min_value=1,
                max_value=20,
                key="workflow_k2",
            )
        with col3:
            k3 = st.number_input(
                "k₃:",
                value=config.get("kpts", (4, 4, 4))[2],
                min_value=1,
                max_value=20,
                key="workflow_k3",
            )
        config["kpts"] = (int(k1), int(k2), int(k3))
        # Remove kspacing from config as we're using explicit grid
        if "kspacing" in config:
            del config["kspacing"]

    # Occupations
    occupations = st.selectbox(
        "Occupations:", ["smearing", "fixed", "tetrahedra"], key="workflow_occupations"
    )
    config["occupations"] = occupations

    # Pseudopotentials - using reusable selector component
    from xespresso.gui.utils.pseudopotentials_selector import render_pseudopotentials_selector
    
    elements = set(atoms.get_chemical_symbols())
    render_pseudopotentials_selector(elements, config, key_prefix="workflow")

    st.markdown("---")

    # Magnetic Configuration - optional expandable section
    from xespresso.gui.utils.magnetic_selector import render_magnetic_selector
    render_magnetic_selector(elements, config, key_prefix="workflow")

    st.markdown("---")

    # Hubbard Configuration - optional expandable section
    from xespresso.gui.utils.hubbard_selector import render_hubbard_selector
    render_hubbard_selector(elements, config, key_prefix="workflow")

    st.markdown("---")

    # Workflow-specific configuration
    if workflow_type == "Single SCF":
        st.subheader("📊 Single SCF Calculation")
        st.info("This workflow performs a single self-consistent field calculation.")
        config["calc_type"] = "scf"

    elif workflow_type == "SCF + Relaxation":
        st.subheader("🔄 SCF + Relaxation Workflow")
        st.info(
            """
        This workflow:
        1. Performs an initial SCF calculation
        2. Relaxes the structure geometry
        """
        )

        relax_type = st.radio(
            "Relaxation Type:", ["relax", "vc-relax"], horizontal=True
        )
        config["relax_type"] = relax_type

        forc_conv_thr = st.number_input(
            "Force Convergence (Ry/bohr):", value=1.0e-3, format="%.1e"
        )
        config["forc_conv_thr"] = forc_conv_thr

    elif workflow_type == "SCF + Relax + SCF (on relaxed structure)":
        st.subheader("🔄 Complete Relaxation Workflow")
        st.info(
            """
        This workflow:
        1. Performs initial SCF calculation
        2. Relaxes the structure
        3. Performs final SCF on relaxed structure for accurate energy
        """
        )

        relax_type = st.radio(
            "Relaxation Type:",
            ["relax", "vc-relax"],
            horizontal=True,
            key="wf_relax_type2",
        )
        config["relax_type"] = relax_type

    else:  # Custom Multi-Step
        st.subheader("🔧 Custom Multi-Step Workflow")
        st.info("Build a custom workflow with multiple calculation steps.")

        num_steps = st.number_input(
            "Number of Steps:", min_value=1, max_value=10, value=2
        )

        if "custom_steps" not in st.session_state:
            st.session_state.custom_steps = []

        for i in range(num_steps):
            with st.expander(f"Step {i+1}"):
                step_type = st.selectbox(
                    "Calculation Type:",
                    ["scf", "relax", "vc-relax", "nscf", "bands"],
                    key=f"step_{i}_type",
                )
                st.session_state.custom_steps.append({"type": step_type, "index": i})

    st.markdown("---")

    st.markdown("---")

    # Resources Configuration
    st.subheader("⚙️ Resources Configuration")
    st.info(
        """
    Configure computational resources for this workflow.
    By default, resources are taken from the machine configuration.
    Enable "Adjust Resources" to customize values for this workflow.
    """
    )

    # Initialize resources in config if not present
    if "resources" not in config:
        config["resources"] = {}

    # Checkbox to enable custom resources
    adjust_resources = st.checkbox(
        "Adjust Resources",
        value=config.get("adjust_resources", False),
        help="Enable to customize resource values for this workflow. Otherwise, defaults from machine configuration are used.",
        key="workflow_adjust_resources",
    )
    config["adjust_resources"] = adjust_resources

    if adjust_resources:
        st.markdown("**Custom Resources:**")

        # Get scheduler type and default resources from machine if available
        scheduler_type = "direct"
        default_resources = {}
        default_nprocs = 1
        default_launcher = "mpirun -np {nprocs}"

        if st.session_state.get("workflow_machine"):
            machine = st.session_state.workflow_machine
            scheduler_type = getattr(machine, "scheduler", "direct")
            default_nprocs = getattr(machine, "nprocs", 1)
            default_launcher = getattr(machine, "launcher", "mpirun -np {nprocs}")
            if hasattr(machine, "resources"):
                default_resources = machine.resources or {}

        # For direct execution, only show nprocs
        if scheduler_type == "direct":
            st.info(
                "ℹ️ **Direct Execution Mode**: Only processor count is configurable. Scheduler resources (nodes, memory, time, etc.) are not applicable for direct execution."
            )

            # Number of processors for direct execution
            nprocs = st.number_input(
                "Number of Processors (nprocs):",
                value=int(config.get("nprocs", default_nprocs)),
                min_value=1,
                max_value=256,
                help="Number of processor cores to use for the calculation",
                key="workflow_nprocs",
            )
            config["nprocs"] = nprocs

            # Update launcher to use the adjusted nprocs value
            # Handles both template placeholders {nprocs} and hardcoded values
            import re

            if "{nprocs}" in default_launcher:
                # Template placeholder - replace it
                resolved_launcher = default_launcher.replace("{nprocs}", str(nprocs))
                config[
                    "launcher"
                ] = default_launcher  # Store template for future adjustments
                st.caption(f"💡 Launcher will be: `{resolved_launcher}`")
            else:
                # Check for hardcoded nprocs values and replace them
                # Pattern matches: -np <number>, -n <number>, --np <number>
                patterns = [
                    (r"(-np\s+)\d+", r"\g<1>" + str(nprocs)),  # -np 16 -> -np 8
                    (r"(-n\s+)\d+", r"\g<1>" + str(nprocs)),  # -n 16 -> -n 8
                    (r"(--np\s+)\d+", r"\g<1>" + str(nprocs)),  # --np 16 -> --np 8
                ]

                resolved_launcher = default_launcher
                for pattern, replacement in patterns:
                    resolved_launcher = re.sub(pattern, replacement, resolved_launcher)

                # Store the updated launcher
                config["launcher"] = resolved_launcher

                if resolved_launcher != default_launcher:
                    st.caption(
                        f"💡 Launcher will be: `{resolved_launcher}` (updated from machine default)"
                    )
                else:
                    st.caption(f"💡 Launcher: `{resolved_launcher}`")
        else:
            # For schedulers (slurm, pbs, sge), show full resource configuration
            st.info(
                f"ℹ️ **Scheduler Mode ({scheduler_type.upper()})**: Configure resources for job scheduler submission."
            )

            # Resource inputs in columns
            col1, col2 = st.columns(2)

            with col1:
                nodes = st.number_input(
                    "Nodes:",
                    value=int(
                        config["resources"].get(
                            "nodes", default_resources.get("nodes", 1)
                        )
                    ),
                    min_value=1,
                    max_value=1000,
                    help="Number of compute nodes to use",
                    key="workflow_nodes",
                )
                config["resources"]["nodes"] = nodes

                ntasks_per_node = st.number_input(
                    "Tasks per Node:",
                    value=int(
                        config["resources"].get(
                            "ntasks-per-node",
                            default_resources.get("ntasks-per-node", 16),
                        )
                    ),
                    min_value=1,
                    max_value=256,
                    help="Number of MPI tasks per node",
                    key="workflow_ntasks",
                )
                config["resources"]["ntasks-per-node"] = ntasks_per_node

                mem = st.text_input(
                    "Memory:",
                    value=config["resources"].get(
                        "mem", default_resources.get("mem", "32G")
                    ),
                    help="Memory per node (e.g., 32G, 64GB)",
                    key="workflow_mem",
                )
                config["resources"]["mem"] = mem

            with col2:
                time = st.text_input(
                    "Time Limit:",
                    value=config["resources"].get(
                        "time", default_resources.get("time", "02:00:00")
                    ),
                    help="Wall time limit (format: HH:MM:SS)",
                    key="workflow_time",
                )
                config["resources"]["time"] = time

                partition = st.text_input(
                    "Partition/Queue:",
                    value=config["resources"].get(
                        "partition", default_resources.get("partition", "compute")
                    ),
                    help="Scheduler partition or queue name",
                    key="workflow_partition",
                )
                config["resources"]["partition"] = partition

                # Additional resource options
                account = st.text_input(
                    "Account (optional):",
                    value=config["resources"].get(
                        "account", default_resources.get("account", "")
                    ),
                    help="Account or project code for billing",
                    key="workflow_account",
                )
                if account:
                    config["resources"]["account"] = account

            st.caption(
                "💡 These custom resources will override the machine defaults for this workflow."
            )
    else:
        # Show default resources from machine
        if st.session_state.get("workflow_machine"):
            machine = st.session_state.workflow_machine
            scheduler_type = getattr(machine, "scheduler", "direct")

            if scheduler_type == "direct":
                # For direct execution, only show nprocs
                st.info("**Using default configuration from machine:**")
                nprocs = getattr(machine, "nprocs", 1)
                launcher = getattr(machine, "launcher", "mpirun -np {nprocs}")
                st.caption(f"Number of Processors: {nprocs}")
                st.caption(
                    f"Launcher: {launcher.replace('{nprocs}', str(nprocs)) if '{nprocs}' in launcher else launcher}"
                )
            elif hasattr(machine, "resources") and machine.resources:
                # For schedulers, show full resource configuration
                st.info("**Using default resources from machine configuration:**")
                col1, col2 = st.columns(2)
                with col1:
                    if "nodes" in machine.resources:
                        st.caption(f"Nodes: {machine.resources['nodes']}")
                    if "ntasks-per-node" in machine.resources:
                        st.caption(
                            f"Tasks per node: {machine.resources['ntasks-per-node']}"
                        )
                    if "mem" in machine.resources:
                        st.caption(f"Memory: {machine.resources['mem']}")
                with col2:
                    if "time" in machine.resources:
                        st.caption(f"Time limit: {machine.resources['time']}")
                    if "partition" in machine.resources:
                        st.caption(f"Partition: {machine.resources['partition']}")
        else:
            st.caption("Select a machine to see default resources")

    st.markdown("---")

    # Build Workflow Button
    st.subheader("✨ Build Workflow")
    st.info(
        """
    Click below to build the workflow using **workflow modules**.
    The workflow will use calculation modules to prepare each step.
    """
    )

    if st.button("🔄 Build Workflow", type="primary"):
        try:
            from xespresso.gui.workflows import GUIWorkflow

            # Validate configuration
            if not config.get("pseudopotentials"):
                st.error("❌ Please specify pseudopotentials for all elements")
                return

            # Validate machine selection
            if not st.session_state.get("workflow_machine"):
                st.error("❌ Please select a machine")
                return

            # Add machine to config as queue parameter (for backwards compatibility)
            config["queue"] = st.session_state.workflow_machine

            # Apply custom resources if enabled
            if config.get("adjust_resources"):
                if config.get("resources"):
                    config["queue"]["resources"] = config["resources"]
                    st.info(f"   Using custom resources: {config['resources']}")
                
                # Apply custom nprocs if adjusted (for direct execution)
                if "nprocs" in config:
                    config["queue"]["nprocs"] = config["nprocs"]
                    st.info(f"   Using custom nprocs: {config['nprocs']}")
                
                # Apply custom launcher if adjusted
                if "launcher" in config:
                    # If launcher still has {nprocs} placeholder, substitute it
                    launcher = config["launcher"]
                    if "{nprocs}" in launcher and "nprocs" in config:
                        launcher = launcher.replace("{nprocs}", str(config["nprocs"]))
                    config["queue"]["launcher"] = launcher
                    st.info(f"   Using custom launcher: {launcher}")

            # Create workflow using workflow module
            st.info("📦 Creating workflow using workflow module...")
            st.info(f"   Machine: {st.session_state.selected_machine}")
            if st.session_state.get("selected_code"):
                st.info(f"   Code: {st.session_state.selected_code}")

            base_label = "workflow"
            workflow = GUIWorkflow(atoms, config, base_label=base_label)

            # Add calculation steps based on workflow type
            if workflow_type == "Single SCF":
                scf_config = config.copy()
                scf_config["calc_type"] = "scf"
                workflow.add_calculation("scf", scf_config)

            elif workflow_type == "SCF + Relaxation":
                # Step 1: SCF
                scf_config = config.copy()
                scf_config["calc_type"] = "scf"
                workflow.add_calculation("scf", scf_config)

                # Step 2: Relax
                relax_config = config.copy()
                relax_config["calc_type"] = config.get("relax_type", "relax")
                workflow.add_calculation("relax", relax_config)

            elif workflow_type == "SCF + Relax + SCF (on relaxed structure)":
                # Step 1: Initial SCF
                scf_config = config.copy()
                scf_config["calc_type"] = "scf"
                workflow.add_calculation("scf_initial", scf_config)

                # Step 2: Relax
                relax_config = config.copy()
                relax_config["calc_type"] = config.get("relax_type", "relax")
                workflow.add_calculation("relax", relax_config)

                # Step 3: Final SCF
                scf_final_config = config.copy()
                scf_final_config["calc_type"] = "scf"
                workflow.add_calculation("scf_final", scf_final_config)

            # Store workflow in session state
            st.session_state.gui_workflow = workflow

            st.success(
                f"""
            ✅ **Workflow built successfully!**
            
            The workflow module has created a {workflow_type} workflow with:
            - {len(workflow.calculations)} calculation step(s)
            - Each step prepared using calculation modules
            
            The workflow is stored in session state and ready for execution.
            """
            )

            # Show workflow summary
            with st.expander("📋 Workflow Summary"):
                st.write(f"**Workflow Type:** {workflow_type}")
                st.write(f"**Base Label:** {base_label}")
                st.write(f"**Number of Steps:** {len(workflow.calculations)}")
                st.write("**Steps:**")
                for i, (name, calc_info) in enumerate(workflow.calculations.items(), 1):
                    st.write(
                        f"{i}. **{name}** - {calc_info['config'].get('calc_type', 'unknown')}"
                    )

        except Exception as e:
            st.error(f"❌ Error building workflow: {e}")
            import traceback

            with st.expander("Error Details"):
                st.code(traceback.format_exc())

    # Show current workflow if exists
    if "gui_workflow" in st.session_state and st.session_state.gui_workflow:
        st.markdown("---")
        st.subheader("📊 Current Workflow")
        workflow = st.session_state.gui_workflow

        st.write(f"**Steps:** {len(workflow.calculations)}")
        for name, calc_info in workflow.calculations.items():
            st.write(f"- {name}: {calc_info['config'].get('calc_type', 'unknown')}")

    st.markdown("---")
    st.info(
        """
    **Next Steps:**
    1. Configure your workflow parameters above
    2. Click "Build Workflow" to create the workflow using workflow modules
    3. Go to "Job Submission" page to execute the workflow steps
    """
    )

"""Job Submission Page for xespresso GUI."""
import streamlit as st
import os
from pathlib import Path


def render_job_submission_page():
    """Render the job submission page with enhanced job file viewer, editor, and submission."""
    st.header("Job Submission & File Management")
    st.markdown(
        """
    **Generate calculation files, browse directories, and run calculations using xespresso.**
    
    - **Generate Files (Dry Run)**: Creates Espresso calculator and generates input files with `calc.write_input(atoms)` for review
    - **File Browser**: Browse, view, and edit existing calculation files
    - **Run Calculation**: Creates Espresso calculator and runs with `calc.get_potential_energy()`
    """
    )

    # Create tabs for different functionalities
    tab1, tab2, tab3 = st.tabs(
        ["📂 File Browser", "🧪 Generate Files (Dry Run)", "🚀 Run Calculation"]
    )

    with tab1:
        render_file_browser_tab()

    with tab2:
        render_dry_run_tab()

    with tab3:
        render_job_submission_tab()


def render_dry_run_tab():
    """Render the dry run tab for generating input and job files without submission."""
    st.subheader("🧪 Generate Calculation Files (Dry Run)")
    st.markdown(
        """
    Generate input files and job scripts for your calculation **without submitting the job**.
    This is useful for:
    - Testing your configuration before submission
    - Manually reviewing and editing files before running
    - Creating files to transfer to another system
    
    **Note:** This uses the calculator prepared in Calculation Setup or creates a new one if needed.
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

    # Safety check: Ensure atoms is actually an Atoms object, not a string
    try:
        from ase import Atoms

        if not isinstance(atoms, Atoms):
            st.error(
                f"❌ Structure is not a valid ASE Atoms object (type: {type(atoms).__name__}). Please reload the structure in the Structure Viewer page."
            )
            if isinstance(atoms, str):
                st.info(f"Debug: Structure appears to be a string: {atoms[:100]}...")
            st.info(
                "💡 Tip: Go to Structure Viewer and reload your structure, then try again."
            )
            return
    except ImportError:
        st.error("❌ ASE not available. Cannot verify structure.")
        return

    st.success(
        f"✅ Structure loaded: {atoms.get_chemical_formula()} ({len(atoms)} atoms)"
    )

    st.markdown("---")

    # Check if calculation is configured
    if (
        "workflow_config" not in st.session_state
        or not st.session_state.workflow_config.get("pseudopotentials")
    ):
        st.warning(
            "⚠️ No calculation configured. Please configure your calculation first in the Calculation Setup page."
        )
        st.info(
            """
        **Required configuration:**
        - Pseudopotentials for all elements
        - Calculation parameters (ecutwfc, kpts, etc.)
        - Machine and codes selection
        
        Go to **📊 Calculation Setup** page to configure these.
        """
        )
        return

    # Show current configuration summary
    st.subheader("📋 Current Configuration")
    config = st.session_state.workflow_config

    col1, col2, col3 = st.columns(3)
    with col1:
        st.write("**Calculation:**")
        st.metric("Type", config.get("calc_type", "scf").upper())
        st.metric("Energy Cutoff", f"{config.get('ecutwfc', 50)} Ry")
        if "kspacing" in config:
            st.metric("K-spacing", f"{config.get('kspacing')} Å⁻¹")
        elif "kpts" in config:
            kpts = config.get("kpts")
            st.metric("K-points", f"{kpts[0]}×{kpts[1]}×{kpts[2]}")

    with col2:
        st.write("**Pseudopotentials:**")
        for elem, pseudo in config.get("pseudopotentials", {}).items():
            st.text(f"  {elem}: {pseudo}")

    with col3:
        st.write("**Machine & Resources:**")
        if "machine_name" in config:
            st.text(f"Machine: {config['machine_name']}")
        if "qe_version" in config:
            st.text(f"QE Version: {config['qe_version']}")
        if "selected_code" in config:
            st.text(f"Code: {config['selected_code']}")

        # Show resources if configured
        if config.get("adjust_resources") and "resources" in config:
            st.text("**Custom Resources:**")
            res = config["resources"]
            if "nodes" in res:
                st.text(f"  Nodes: {res['nodes']}")
            if "ntasks-per-node" in res:
                st.text(f"  Tasks/node: {res['ntasks-per-node']}")
            if "time" in res:
                st.text(f"  Time: {res['time']}")
        elif "machine_name" in config:
            # Show machine defaults
            machine = st.session_state.get("calc_machine") or st.session_state.get(
                "workflow_machine"
            )
            if machine:
                if hasattr(machine, "scheduler") and machine.scheduler:
                    st.text(f"Scheduler: {machine.scheduler}")
                if hasattr(machine, "nprocs"):
                    st.text(f"Processors: {machine.nprocs}")
                elif hasattr(machine, "resources") and machine.resources:
                    res = machine.resources
                    if "nodes" in res:
                        st.text(f"Nodes: {res['nodes']}")
                    if "ntasks-per-node" in res:
                        st.text(f"Tasks/node: {res['ntasks-per-node']}")

    st.markdown("---")

    # Use working directory from session state + calculation label
    st.subheader("📁 Output Directory")

    # Get base working directory from session state
    base_workdir = st.session_state.get("working_directory", os.path.expanduser("~"))
    st.info(f"📍 Base directory: `{base_workdir}`")

    # Label/subfolder for this calculation - inherited from Calculation Setup
    # If label is in workflow_config, use it (set in Calculation Setup)
    # Otherwise, generate a default
    if "label" in config and config["label"]:
        label = config["label"]
        st.success(f"📋 Using label from Calculation Setup: `{label}`")
        st.caption("💡 To change the label, go back to Calculation Setup page")
    else:
        # Fallback: generate default label
        default_label = (
            f"{config.get('calc_type', 'scf')}/{atoms.get_chemical_formula()}"
        )
        st.warning("⚠️ No label set in Calculation Setup. Using default.")
        label = st.text_input(
            "Calculation Label (subfolder):",
            value=default_label,
            help="Label for this calculation - will create subfolder under working directory",
            key="dry_run_label",
        )

    # Full path where files will be created
    full_path = os.path.join(base_workdir, label)

    # Validate the full path to prevent path traversal in the label
    try:
        full_path = os.path.realpath(full_path)
        base_workdir_real = os.path.realpath(base_workdir)

        # Ensure full_path is under base_workdir (prevent path traversal via label)
        if not full_path.startswith(base_workdir_real):
            st.error("❌ Invalid calculation label - path traversal detected")
            return
    except (OSError, ValueError) as e:
        st.error(f"❌ Invalid path: {e}")
        return

    st.success(f"✅ Files will be created in: `{full_path}`")

    st.markdown("---")

    # Generate button
    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        generate_button = st.button(
            "🧪 Generate Files (Dry Run)",
            type="primary",
            help="Generate input and job files without submitting",
            key="generate_files_button",
        )

    if generate_button:
        st.info("🧪 **Dry Run Mode** - Generating input files...")

        with st.spinner("Generating files..."):
            try:
                from gui.calculations import (
                    dry_run_calculation,
                    prepare_calculation_from_gui,
                )
                from ase import io as ase_io

                # Create output directory if it doesn't exist
                os.makedirs(full_path, exist_ok=True)

                # Save structure file
                structure_filename = f"{atoms.get_chemical_formula()}.cif"
                structure_path = os.path.join(full_path, structure_filename)
                ase_io.write(structure_path, atoms)
                st.info(f"💾 Saved structure: {structure_filename}")

                # Pass the label directly to Espresso calculator
                # Espresso expects just a label like "scf/Gd2" and will handle local/remote paths correctly
                if (
                    "espresso_calculator" in st.session_state
                    and st.session_state.espresso_calculator is not None
                ):
                    st.info(
                        "📦 Using pre-configured calculator from Calculation Setup..."
                    )
                    calc = st.session_state.espresso_calculator
                    prepared_atoms = st.session_state.get("prepared_atoms", atoms)

                    # Update the label - Espresso will handle the path resolution
                    calc.set_label(label, calc.prefix)

                    # Write input files using xespresso's method
                    calc.write_input(prepared_atoms)
                else:
                    # Use calculation module to prepare atoms and calculator, then generate files
                    st.info("🔧 Creating Espresso calculator and generating files...")
                    prepared_atoms, calc = dry_run_calculation(
                        atoms, config, label=label
                    )

                st.success("✅ Files generated successfully using xespresso!")

                # Display results
                st.subheader("📄 Generated Files")

                # List generated files
                generated_files = []
                if os.path.exists(full_path):
                    for f in os.listdir(full_path):
                        if f.endswith(
                            (".pwi", ".asei", ".cif", "job_file", ".sh", ".slurm")
                        ):
                            generated_files.append(f)

                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Files created:**")
                    for f in generated_files:
                        st.markdown(f"- ✅ {f}")

                with col2:
                    st.write("**Location:**")
                    st.code(full_path)

                st.markdown("---")

                # Preview input file
                input_file_path = os.path.join(full_path, "espresso.pwi")
                if os.path.exists(input_file_path):
                    st.subheader("👁️ Input File Preview")
                    try:
                        # Validate path is under full_path for security
                        input_path = os.path.realpath(input_file_path)
                        if not input_path.startswith(full_path):
                            st.error(
                                "❌ Security error: input file path is outside expected directory"
                            )
                        else:
                            with open(input_path, "r") as f:
                                input_content = f.read()

                            with st.expander("View Input File", expanded=True):
                                st.code(
                                    input_content, language="fortran", line_numbers=True
                                )

                                # Download button
                                st.download_button(
                                    label="⬇️ Download Input File",
                                    data=input_content,
                                    file_name="espresso.pwi",
                                    mime="text/plain",
                                )
                    except Exception as e:
                        st.error(f"Error reading input file: {e}")

                # Preview job file if it exists
                job_file_path = os.path.join(full_path, "job_file")
                if os.path.exists(job_file_path):
                    st.subheader("👁️ Job Script Preview")
                    try:
                        # Validate path is under full_path for security
                        job_path = os.path.realpath(job_file_path)
                        if not job_path.startswith(full_path):
                            st.error(
                                "❌ Security error: job file path is outside expected directory"
                            )
                        else:
                            with open(job_path, "r") as f:
                                job_content = f.read()

                            with st.expander("View Job Script", expanded=False):
                                st.code(job_content, language="bash", line_numbers=True)

                                # Download button
                                st.download_button(
                                    label="⬇️ Download Job Script",
                                    data=job_content,
                                    file_name="job_file",
                                    mime="text/plain",
                                )
                    except Exception as e:
                        st.error(f"Error reading job file: {e}")

                st.markdown("---")

                # Next steps
                st.subheader("✨ Next Steps")
                st.info(
                    """
                **Files have been generated!** You can now:
                
                1. **Review the files** using the File Browser tab above
                2. **Edit the files** if needed (use Edit mode in File Browser)
                3. **Run the calculation** using the Run Calculation tab
                4. **Transfer files** to another system if needed
                """
                )

            except Exception as e:
                st.error(f"❌ Error generating files: {e}")
                import traceback

                with st.expander("Error Details"):
                    st.code(traceback.format_exc())


def render_file_browser_tab():
    """Render the file browser tab."""
    # Working Directory Browser Section
    st.subheader("📁 Browse Calculation Folders")

    # Use base working directory from session state
    base_workdir = st.session_state.get("working_directory", os.path.expanduser("~"))
    st.info(f"📍 Browsing in: `{base_workdir}`")

    if not os.path.exists(base_workdir) or not os.path.isdir(base_workdir):
        st.error(f"❌ Invalid working directory: {base_workdir}")
        return

    st.markdown("---")

    # Calculation Folder Navigation
    st.subheader("📂 Calculation Folders")
    st.info(
        "💡 xespresso organizes calculations in label-based folders (e.g., 'scf/formula')"
    )

    # Find calculation folders (those with input files)
    calc_folders = []
    input_file_extensions = [".in", ".pwi", ".phi", ".ppi", ".bandi"]

    try:
        for root, dirs, files in os.walk(base_workdir, topdown=True):
            # Limit depth to avoid too much recursion
            depth = root[len(base_workdir) :].count(os.sep)
            if depth < 4:
                has_input_files = any(
                    f.endswith(tuple(input_file_extensions))
                    or f == "job_file"
                    or f.endswith(".sh")
                    or f.endswith(".slurm")
                    for f in files
                )
                if has_input_files:
                    calc_folders.append(root)

        if calc_folders:
            st.success(f"✅ Found {len(calc_folders)} calculation folder(s)")

            # Select calculation folder
            selected_folder = st.selectbox(
                "Select Calculation Folder:",
                calc_folders,
                format_func=lambda x: os.path.relpath(x, base_workdir)
                if x != base_workdir
                else ".",
            )

            if selected_folder:
                st.info(
                    f"📍 Selected: `{os.path.relpath(selected_folder, base_workdir)}`"
                )

                # List files in the selected folder
                try:
                    files_in_folder = [
                        f
                        for f in os.listdir(selected_folder)
                        if os.path.isfile(os.path.join(selected_folder, f))
                    ]

                    # Categorize files
                    input_files = [
                        f
                        for f in files_in_folder
                        if f.endswith(tuple(input_file_extensions))
                    ]
                    job_files = [
                        f
                        for f in files_in_folder
                        if f == "job_file" or f.endswith((".sh", ".slurm"))
                    ]
                    output_files = [
                        f
                        for f in files_in_folder
                        if f.endswith((".out", ".pwo", ".pho", ".ppo"))
                    ]
                    other_files = [
                        f
                        for f in files_in_folder
                        if f not in input_files + job_files + output_files
                    ]

                    # Display file categories
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("Input Files", len(input_files))
                    with col2:
                        st.metric("Job Scripts", len(job_files))
                    with col3:
                        st.metric("Output Files", len(output_files))
                    with col4:
                        st.metric("Other Files", len(other_files))

                    st.markdown("---")

                    # File Viewer/Editor Section
                    st.subheader("📄 File Viewer & Editor")

                    # File type selector
                    file_type = st.radio(
                        "File Category:",
                        ["Input Files", "Job Scripts", "Output Files", "All Files"],
                        horizontal=True,
                    )

                    # Select file to view based on category
                    if file_type == "Input Files":
                        available_files = input_files
                    elif file_type == "Job Scripts":
                        available_files = job_files
                    elif file_type == "Output Files":
                        available_files = output_files
                    else:
                        available_files = files_in_folder

                    if available_files:
                        selected_file = st.selectbox("Select File:", available_files)

                        if selected_file:
                            file_path = os.path.join(selected_folder, selected_file)

                            # Security check
                            if (
                                not os.path.commonpath([selected_folder, file_path])
                                == selected_folder
                            ):
                                st.error("❌ Invalid file path")
                                return

                            # View/Edit mode selector
                            mode = st.radio(
                                "Mode:",
                                ["View", "Edit"],
                                horizontal=True,
                                key="file_mode",
                            )

                            try:
                                with open(file_path, "r") as f:
                                    file_content = f.read()

                                # Show file info
                                file_stat = os.stat(file_path)
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("File Size", f"{file_stat.st_size} bytes")
                                with col2:
                                    st.metric("Lines", len(file_content.splitlines()))
                                with col3:
                                    from datetime import datetime

                                    mod_time = datetime.fromtimestamp(
                                        file_stat.st_mtime
                                    )
                                    st.metric(
                                        "Modified", mod_time.strftime("%Y-%m-%d %H:%M")
                                    )

                                if mode == "View":
                                    # View mode - display as code
                                    st.subheader("File Content")
                                    st.code(
                                        file_content,
                                        language="bash"
                                        if selected_file.endswith(
                                            (".sh", ".slurm", "job_file")
                                        )
                                        else "fortran",
                                        line_numbers=True,
                                    )

                                    # Download button
                                    st.download_button(
                                        label="⬇️ Download File",
                                        data=file_content,
                                        file_name=selected_file,
                                        mime="text/plain",
                                    )

                                    # Parse job file information
                                    if selected_file in job_files:
                                        st.subheader("Job File Summary")

                                        # Extract SLURM/PBS directives
                                        scheduler_directives = []
                                        commands = []
                                        for line in file_content.splitlines():
                                            line_stripped = line.strip()
                                            if line_stripped.startswith(
                                                "#SBATCH"
                                            ) or line_stripped.startswith("#PBS"):
                                                scheduler_directives.append(
                                                    line_stripped
                                                )
                                            elif (
                                                line_stripped
                                                and not line_stripped.startswith("#")
                                            ):
                                                commands.append(line_stripped)

                                        if scheduler_directives:
                                            st.markdown("**Scheduler Directives:**")
                                            for directive in scheduler_directives:
                                                st.markdown(f"- `{directive}`")

                                        if commands:
                                            st.markdown("**Key Commands:**")
                                            for cmd in commands[:10]:
                                                st.markdown(f"- `{cmd}`")
                                            if len(commands) > 10:
                                                st.markdown(
                                                    f"... and {len(commands) - 10} more commands"
                                                )

                                else:
                                    # Edit mode - text area for editing
                                    st.subheader("Edit File Content")
                                    st.warning(
                                        "⚠️ **Caution:** Editing files can affect your calculations. Make sure you know what you're doing!"
                                    )

                                    edited_content = st.text_area(
                                        "File Content:",
                                        value=file_content,
                                        height=400,
                                        key="file_editor",
                                    )

                                    col1, col2 = st.columns(2)
                                    with col1:
                                        if st.button("💾 Save Changes", type="primary"):
                                            try:
                                                with open(file_path, "w") as f:
                                                    f.write(edited_content)
                                                st.success(
                                                    f"✅ File saved successfully: {selected_file}"
                                                )
                                            except Exception as e:
                                                st.error(f"❌ Error saving file: {e}")

                                    with col2:
                                        if st.button("↩️ Revert Changes"):
                                            st.rerun()

                            except Exception as e:
                                st.error(f"❌ Error reading file: {e}")
                    else:
                        st.warning(f"⚠️ No files found in the '{file_type}' category")

                except Exception as e:
                    st.error(f"❌ Error listing files: {e}")
        else:
            st.warning("⚠️ No calculation folders found in the working directory")
            st.info(
                """
            **Tip:** Calculation folders typically contain:
            - Input files (`*.in`, `*.pwi`, etc.)
            - Job scripts (`job_file`, `*.sh`, `*.slurm`)
            - Output files (`*.out`, `*.pwo`, etc.)
            
            Make sure you've set up calculations first in the Calculation Setup page.
            """
            )

    except Exception as e:
        st.error(f"❌ Error scanning directory: {e}")

    # Additional info section
    st.markdown("---")
    st.subheader("💡 Tips")
    st.info(
        """
    **Working with xespresso calculations:**
    - xespresso organizes files in label-based folders (e.g., `calc/label/`)
    - Input files are named based on the calculation type (e.g., `*.pwi` for pw.x)
    - Job files are typically named `job_file` or have `.sh`/`.slurm` extensions
    - The results folder is the same as the calculation folder
    
    **File editing:**
    - Use the Edit mode to modify input parameters
    - Always backup important files before editing
    - Changes are saved immediately when you click "Save Changes"
    """
    )


def render_job_submission_tab():
    """Render the job submission tab for running calculations with xespresso."""
    st.subheader("🚀 Run Calculation")
    st.markdown(
        """
    Run a calculation using xespresso by calling `calc.get_potential_energy()`.
    
    This will:
    - Create an Espresso calculator from your configuration
    - Automatically generate input files if they don't exist
    - Run the calculation and return the energy
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

    # Safety check: Ensure atoms is actually an Atoms object, not a string
    try:
        from ase import Atoms

        if not isinstance(atoms, Atoms):
            st.error(
                f"❌ Structure is not a valid ASE Atoms object (type: {type(atoms).__name__}). Please reload the structure in the Structure Viewer page."
            )
            if isinstance(atoms, str):
                st.info(f"Debug: Structure appears to be a string: {atoms[:100]}...")
            st.info(
                "💡 Tip: Go to Structure Viewer and reload your structure, then try again."
            )
            return
    except ImportError:
        st.error("❌ ASE not available. Cannot verify structure.")
        return

    st.success(
        f"✅ Structure loaded: {atoms.get_chemical_formula()} ({len(atoms)} atoms)"
    )

    st.markdown("---")

    # Check if calculation is configured
    if (
        "workflow_config" not in st.session_state
        or not st.session_state.workflow_config.get("pseudopotentials")
    ):
        st.warning(
            "⚠️ No calculation configured. Please configure your calculation first in the Calculation Setup page."
        )
        st.info(
            """
        **Required configuration:**
        - Pseudopotentials for all elements
        - Calculation parameters (ecutwfc, kpts, etc.)
        - Machine and codes selection
        
        Go to **📊 Calculation Setup** page to configure these.
        """
        )
        return

    # Show current configuration summary
    st.subheader("📋 Current Configuration")
    config = st.session_state.workflow_config

    col1, col2, col3 = st.columns(3)
    with col1:
        st.write("**Calculation:**")
        st.metric("Type", config.get("calc_type", "scf").upper())
        st.metric("Energy Cutoff", f"{config.get('ecutwfc', 50)} Ry")
        if "kspacing" in config:
            st.metric("K-spacing", f"{config.get('kspacing')} Å⁻¹")
        elif "kpts" in config:
            kpts = config.get("kpts")
            st.metric("K-points", f"{kpts[0]}×{kpts[1]}×{kpts[2]}")

    with col2:
        st.write("**Pseudopotentials:**")
        for elem, pseudo in config.get("pseudopotentials", {}).items():
            st.text(f"  {elem}: {pseudo}")

    with col3:
        st.write("**Machine & Resources:**")
        if "machine_name" in config:
            st.text(f"Machine: {config['machine_name']}")
        if "qe_version" in config:
            st.text(f"QE Version: {config['qe_version']}")
        if "selected_code" in config:
            st.text(f"Code: {config['selected_code']}")

        # Show resources if configured
        if config.get("adjust_resources") and "resources" in config:
            st.text("**Custom Resources:**")
            res = config["resources"]
            if "nodes" in res:
                st.text(f"  Nodes: {res['nodes']}")
            if "ntasks-per-node" in res:
                st.text(f"  Tasks/node: {res['ntasks-per-node']}")
            if "time" in res:
                st.text(f"  Time: {res['time']}")
        elif "machine_name" in config:
            # Show machine defaults
            machine = st.session_state.get("calc_machine") or st.session_state.get(
                "workflow_machine"
            )
            if machine:
                if hasattr(machine, "scheduler") and machine.scheduler:
                    st.text(f"Scheduler: {machine.scheduler}")
                if hasattr(machine, "nprocs"):
                    st.text(f"Processors: {machine.nprocs}")
                elif hasattr(machine, "resources") and machine.resources:
                    res = machine.resources
                    if "nodes" in res:
                        st.text(f"Nodes: {res['nodes']}")
                    if "ntasks-per-node" in res:
                        st.text(f"Tasks/node: {res['ntasks-per-node']}")

    st.markdown("---")

    # Working directory and label
    st.subheader("📁 Output Location")

    # Use base working directory from session state
    base_workdir = st.session_state.get("working_directory", os.path.expanduser("~"))
    st.info(f"📍 Base directory: `{base_workdir}`")

    # Label for calculation - inherited from Calculation Setup
    # If label is in workflow_config, use it (set in Calculation Setup)
    # Otherwise, generate a default
    if "label" in config and config["label"]:
        label = config["label"]
        st.success(f"📋 Using label from Calculation Setup: `{label}`")
        st.caption("💡 To change the label, go back to Calculation Setup page")
    else:
        # Fallback: generate default label
        default_label = (
            f"{config.get('calc_type', 'scf')}/{atoms.get_chemical_formula()}"
        )
        st.warning("⚠️ No label set in Calculation Setup. Using default.")
        label = st.text_input(
            "Calculation Label (subfolder):",
            value=default_label,
            help="Label for this calculation - will create subfolder under working directory",
            key="run_calc_label",
        )

    # Full path where calculation will run
    full_path = os.path.join(base_workdir, label)

    # Validate the full path to prevent path traversal in the label
    try:
        full_path = os.path.realpath(full_path)
        base_workdir_real = os.path.realpath(base_workdir)

        # Ensure full_path is under base_workdir (prevent path traversal via label)
        if not full_path.startswith(base_workdir_real):
            st.error("❌ Invalid calculation label - path traversal detected")
            return
    except (OSError, ValueError) as e:
        st.error(f"❌ Invalid path: {e}")
        return

    st.success(f"✅ Calculation will run in: `{full_path}`")

    st.markdown("---")

    # Run calculation button
    st.subheader("▶️ Execute Calculation")

    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        run_button = st.button(
            "🚀 Run Calculation",
            type="primary",
            help="Run calculation using calc.get_potential_energy()",
            key="run_calculation_button",
        )

    if run_button:
        st.info(
            "🚀 **Running Calculation** - Using xespresso's calc.get_potential_energy()..."
        )

        with st.spinner("Running calculation..."):
            try:
                from gui.calculations import prepare_calculation_from_gui

                # Create output directory if it doesn't exist
                os.makedirs(full_path, exist_ok=True)

                # Pass the label directly to Espresso calculator
                # Espresso expects just a label like "scf/Gd2" and will handle local/remote paths correctly
                # Check if calculator already exists in session_state
                # (created by Calculation Setup or Workflow Builder using calculation modules)
                if (
                    "espresso_calculator" in st.session_state
                    and st.session_state.espresso_calculator is not None
                ):
                    st.info(
                        "📦 Using pre-configured calculator from Calculation Setup (prepared by calculation module)..."
                    )
                    calc = st.session_state.espresso_calculator
                    prepared_atoms = st.session_state.get("prepared_atoms", atoms)

                    # Update the label - Espresso will handle the path resolution
                    calc.set_label(label, calc.prefix)
                    
                    # Check if previous calculation failed and force recalculation if needed
                    # This allows users to rerun failed calculations from the GUI
                    if hasattr(calc, 'results') and 'convergence' in calc.results:
                        convergence_status = calc.results['convergence']
                        if convergence_status > 0:
                            st.warning(f"⚠️ Previous calculation failed (status: {convergence_status}). Forcing recalculation...")
                            # Clear cached results to force a fresh calculation
                            calc.reset()
                            st.info("✓ Calculator reset - ready for fresh calculation")
                else:
                    # Use calculation module to prepare atoms and Espresso calculator
                    # Following the principle: calculation modules create objects, job submission executes
                    st.info(
                        "🔧 Using calculation module to prepare atoms and Espresso calculator..."
                    )

                    prepared_atoms, calc = prepare_calculation_from_gui(
                        atoms, config, label=label
                    )

                    st.info("✅ Calculation module prepared objects from configuration!")
                    
                    # Check if this calculator loaded a failed previous calculation
                    # If so, reset it to allow rerun
                    if hasattr(calc, 'results') and 'convergence' in calc.results:
                        convergence_status = calc.results['convergence']
                        if convergence_status > 0:
                            st.warning(f"⚠️ Found previous failed calculation (status: {convergence_status}). Forcing recalculation...")
                            # Clear cached results to force a fresh calculation
                            calc.reset()
                            st.info("✓ Calculator reset - ready for fresh calculation")

                # Attach calculator to prepared atoms (relationship maintained by calculation module)
                st.info("🔗 Attaching calculator to atoms object...")
                prepared_atoms.calc = calc

                # Run calculation using get_potential_energy()
                st.info("⚡ Calling atoms.get_potential_energy()...")
                energy = prepared_atoms.get_potential_energy()

                # Display results
                st.success("✅ Calculation completed successfully!")
                st.markdown("---")
                st.subheader("📊 Results")

                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Total Energy", f"{energy:.6f} eV")
                with col2:
                    st.metric("Structure", prepared_atoms.get_chemical_formula())

                # Show calculation details
                with st.expander("📋 Calculation Details"):
                    st.write("**Input Parameters:**")
                    if hasattr(calc, "parameters") and "input_data" in calc.parameters:
                        st.json(calc.parameters["input_data"])
                    else:
                        st.json(
                            calc.input_data if hasattr(calc, "input_data") else config
                        )

                    st.write("**Pseudopotentials:**")
                    for species, pseudo in config["pseudopotentials"].items():
                        st.text(f"  {species}: {pseudo}")

                    if hasattr(calc, "kpts") and calc.kpts:
                        st.write("**K-points:**")
                        st.text(f"  {calc.kpts}")
                    elif "kspacing" in config:
                        st.write("**K-spacing:**")
                        st.text(f"  {config['kspacing']} Å⁻¹")

                # Show output location
                st.markdown("---")
                st.subheader("📁 Output Files")
                st.info(f"Calculation files saved in: `{full_path}`")

                # List generated files
                if os.path.exists(full_path):
                    generated_files = []
                    for f in os.listdir(full_path):
                        if f.endswith((".pwi", ".pwo", ".out", ".log")):
                            generated_files.append(f)

                    if generated_files:
                        st.write("**Generated files:**")
                        for f in generated_files:
                            st.markdown(f"- `{f}`")

            except Exception as e:
                st.error(f"❌ Calculation error: {e}")
                import traceback

                with st.expander("Error Details"):
                    st.code(traceback.format_exc())

"""
editor.py

Utility for interactively editing existing machine configurations for xespresso.

This module supports:
- Editing any field of a machine profile
- Reusing current values as defaults
- Preserving unchanged fields
- Compatible with machines created via creator.py
- Logging and warnings for traceability

Usage:
from xespresso.machines.config.editor import edit_machine
edit_machine("cluster_a")  # Edit an existing machine
"""

import os
import json
try:
    from xespresso.utils import warnings as warnings
    from xespresso.utils.logging import get_logger
    logger = get_logger()
    warnings.apply_custom_format()
except ImportError:
    # Fallback to standard Python warnings and logging
    import warnings
    import logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO)

DEFAULT_CONFIG_PATH = os.path.expanduser("~/.xespresso/machines.json")

def edit_machine(machine_name: str, path: str = DEFAULT_CONFIG_PATH):
    """
    Interactively edits an existing machine configuration.

    Parameters:
    - machine_name (str): Name of the machine to edit
    - path (str): Path to the config file
    """
    if not os.path.exists(path):
        print(f"❌ Config file not found at {path}")
        logger.error(f"Config file not found: {path}")
        return

    try:
        with open(path) as f:
            config = json.load(f)
        logger.info(f"Loaded config from {path}")
    except Exception as e:
        print("❌ Failed to load config file.")
        logger.error(f"Failed to parse config: {e}")
        return

    machines = config.get("machines", {})
    if machine_name not in machines:
        print(f"❌ Machine '{machine_name}' not found.")
        logger.warning(f"Machine '{machine_name}' not found.")
        return

    machine = machines[machine_name]
    print(f"✏️ Editing machine: {machine_name}")

    # Execution mode
    execution = input(f"Execution mode [local/remote] [{machine.get('execution', 'local')}]: ").strip() or machine.get("execution", "local")
    machine["execution"] = execution

    # Scheduler
    scheduler = input(f"Scheduler [direct/slurm] [{machine.get('scheduler', 'direct')}]: ").strip() or machine.get("scheduler", "direct")
    machine["scheduler"] = scheduler

    # Workdir
    workdir = input(f"Workdir path [{machine.get('workdir', './xespresso')}]: ").strip() or machine.get("workdir", "./xespresso")
    machine["workdir"] = workdir

    # Remote fields
    if execution == "remote":
        host = input(f"Remote host [{machine.get('host', '')}]: ").strip() or machine.get("host", "")
        machine["host"] = host

        port = input(f"SSH port [{machine.get('port', 22)}]: ").strip() or str(machine.get("port", 22))
        machine["port"] = int(port)

        username = input(f"SSH username [{machine.get('username', '')}]: ").strip() or machine.get("username", "")
        machine["username"] = username

        auth = machine.get("auth", {})
        ssh_key = input(f"Path to SSH key [{auth.get('ssh_key', '~/.ssh/id_rsa.pub')}]: ").strip() or auth.get("ssh_key", "~/.ssh/id_rsa.pub")
        machine["auth"] = {"method": "key", "ssh_key": ssh_key}

    # Resources
    if scheduler == "slurm":
        print("🧮 Edit job resources (press Enter to keep current):")
        res = machine.get("resources", {})
        nodes = input(f"Number of nodes [{res.get('nodes', '')}]: ").strip()
        ntasks = input(f"Tasks per node [{res.get('ntasks-per-node', '')}]: ").strip()
        time = input(f"Walltime [{res.get('time', '')}]: ").strip()
        partition = input(f"Partition [{res.get('partition', '')}]: ").strip()

        machine["resources"] = {
            "nodes": int(nodes) if nodes else res.get("nodes"),
            "ntasks-per-node": int(ntasks) if ntasks else res.get("ntasks-per-node"),
            "time": time if time else res.get("time"),
            "partition": partition if partition else res.get("partition")
        }

    # nprocs
    nprocs = input(f"Number of processes [{machine.get('nprocs', 1)}]: ").strip() or str(machine.get("nprocs", 1))
    machine["nprocs"] = int(nprocs)

    # Launcher
    print("🧭 Define the launcher command used to run Quantum ESPRESSO.")
    print("You may use the placeholder {nprocs}, which will be replaced at runtime.")
    print("Examples:")
    print(" - mpirun -np {nprocs}          (direct or manual MPI)")
    print(" - srun --mpi=pmi2              (Slurm with Intel MPI)")
    print("Note: For Slurm with Intel MPI, use 'srun --mpi=pmi2' without {nprocs}.")
    launcher = input(f"Launcher command [{machine.get('launcher', 'mpirun -np {nprocs}')}]: ").strip() or machine.get("launcher", "mpirun -np {nprocs}")
    machine["launcher"] = launcher

    # Modules
    use_modules = input(f"Use modules? [y/N] [{'y' if machine.get('use_modules', False) else 'N'}]: ").strip().lower()
    machine["use_modules"] = use_modules == "y"

    if machine["use_modules"]:
        modules_default = machine.get("modules", [])
        modules_input = input(f"Modules to load (comma-separated) [{', '.join(modules_default)}]: ").strip()
        machine["modules"] = [m.strip() for m in modules_input.split(",")] if modules_input else modules_default
    else:
        # Mantém a lista existente, mesmo que não seja usada no job_file
        logger.info("Modules disabled — keeping existing module list.")

    # Prepend / Postpend
    prepend_default = machine.get("prepend", [])
    prepend_input = input(f"Commands to prepend before job (comma-separated) [{', '.join(prepend_default)}]: ").strip()
    machine["prepend"] = [cmd.strip() for cmd in prepend_input.split(",")] if prepend_input else prepend_default

    postpend_default = machine.get("postpend", [])
    postpend_input = input(f"Commands to postpend after job (comma-separated) [{', '.join(postpend_default)}]: ").strip()
    machine["postpend"] = [cmd.strip() for cmd in postpend_input.split(",")] if postpend_input else postpend_default

    # Save
    config["machines"][machine_name] = machine
    try:
        with open(path, "w") as f:
            json.dump(config, f, indent=2)
        print(f"✅ Machine '{machine_name}' updated in {path}")
        logger.info(f"Machine '{machine_name}' updated successfully.")
    except Exception as e:
        print("❌ Failed to save updated machine.")
        logger.error(f"Failed to write config file: {e}")

"""
creator.py

Utility for interactively creating or editing machine configurations for xespresso workflows.

Supports:
- Creating new machine profiles
- Editing existing profiles via edit_machine()
- Local and remote execution modes
- Key-based SSH authentication only
- Optional job resources for Slurm
- Launcher command with {nprocs} placeholder
- Loading presets from templates/ or external paths
- Logging and warnings for traceability

Usage:
from xespresso.machines.config.creator import create_machine
create_machine()  # Launch interactive setup
create_machine(preset_path="/path/to/preset.json")  # Load preset automatically
"""

import os
import json
try:
    from xespresso.utils.auth import generate_ssh_key, install_ssh_key, test_ssh_connection
    AUTH_AVAILABLE = True
except ImportError:
    AUTH_AVAILABLE = False
from .editor import edit_machine
from .presets import list_presets, load_preset
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

def create_machine(path: str = DEFAULT_CONFIG_PATH, preset_path: str = None):
    logger.info("Starting interactive machine configuration.")
    os.makedirs(os.path.dirname(path), exist_ok=True)

    config = {"machines": {}}
    if os.path.exists(path):
        try:
            with open(path) as f:
                config = json.load(f)
            logger.info(f"Loaded existing config from {path}")
        except Exception as e:
            logger.error(f"Failed to load config: {e}")
            return

    machine_name = input("Machine name (e.g. local_desktop, cluster_a): ").strip()
    if not machine_name:
        logger.warning("Machine name was left empty.")
        print("❌ Machine name cannot be empty.")
        return

    if machine_name in config["machines"]:
        print(f"⚠️ Machine '{machine_name}' already exists.")
        print("Options:")
        print(" [1] Overwrite from scratch")
        print(" [2] Edit existing configuration")
        choice = input("Choose an option [1/2]: ").strip()
        if choice == "2":
            edit_machine(machine_name, path)
            return
        elif choice != "1":
            print("❌ Invalid choice. No changes made.")
            logger.warning(f"Invalid overwrite option selected: '{choice}'")
            return

    machine = {
        "execution": "local",
        "scheduler": "direct",
        "workdir": "./xespresso",
        "modules": [],
        "use_modules": False,
        "prepend": [],
        "postpend": [],
        "resources": {}
    }

    # Load preset from argument
    if preset_path:
        if os.path.isfile(preset_path):
            try:
                with open(preset_path) as f:
                    preset = json.load(f)
                machine.update(preset)
                logger.info(f"Preset loaded from argument: {preset_path}")
            except Exception as e:
                print("⚠️ Failed to load preset from argument.")
                logger.warning(f"Preset load failed: {e}")
        else:
            print(f"⚠️ Preset path '{preset_path}' not found.")
            logger.warning(f"Invalid preset path: {preset_path}")

    # Unified preset input
    available_presets = list_presets()
    if available_presets:
        print("🧭 Available machine presets:")
        for p in available_presets:
            print(f" - {p}")
    preset_input = input("Choose a preset name from above, or enter full path to a preset (.json). Press Enter to skip: ").strip()

    if preset_input:
        if preset_input in available_presets:
            try:
                preset = load_preset(preset_input)
                machine.update(preset)
                logger.info(f"Preset '{preset_input}' loaded successfully.")
            except Exception as e:
                print("⚠️ Failed to load preset.")
                logger.warning(f"Preset load failed: {e}")
        elif os.path.isfile(preset_input):
            try:
                with open(preset_input) as f:
                    external_preset = json.load(f)
                machine.update(external_preset)
                logger.info(f"Custom preset loaded from: {preset_input}")
            except Exception as e:
                print("⚠️ Failed to load custom preset.")
                logger.warning(f"Preset load failed: {e}")
        else:
            print(f"⚠️ Preset '{preset_input}' not found.")
            logger.warning(f"Invalid preset input: {preset_input}")

    # Interactive prompts
    machine["execution"] = input(f"Execution mode [local/remote] [{machine['execution']}]: ").strip().lower() or machine["execution"]
    machine["scheduler"] = input(f"Scheduler [direct/slurm] [{machine['scheduler']}]: ").strip() or machine["scheduler"]
    machine["workdir"] = input(f"Workdir path [{machine['workdir']}]: ").strip() or machine["workdir"]

    if machine["execution"] == "remote":
        machine["host"] = input(f"Remote host [{machine.get('host', '')}]: ").strip() or machine.get("host", "")
        machine["port"] = int(input(f"SSH port [22]: ").strip() or machine.get("port", 22))
        machine["username"] = input(f"SSH username [{machine.get('username', '')}]: ").strip() or machine.get("username", "")
        ssh_key = input(f"Path to SSH key [~/.ssh/id_rsa.pub]: ").strip() or "~/.ssh/id_rsa.pub"
        ssh_key = os.path.expanduser(ssh_key)

        if not os.path.isfile(ssh_key):
            logger.warning(f"SSH key not found at {ssh_key}")
            print(f"⚠️ SSH key '{ssh_key}' not found.")
            if AUTH_AVAILABLE:
                create_key = input("Generate new SSH key pair now? [y/N]: ").strip().lower()
                if create_key == "y":
                    try:
                        generate_ssh_key(ssh_key.replace(".pub", ""))
                        logger.info(f"SSH key generated at {ssh_key.replace('.pub', '')}")
                        print("✅ Key created.")
                        install = input("Install this key on the remote server now? [y/N]: ").strip().lower()
                        if install == "y":
                            install_ssh_key(machine["username"], machine["host"], ssh_key, machine["port"])
                    except Exception as e:
                        logger.error(f"Failed to generate SSH key: {e}")
                        print("❌ Failed to generate SSH key.")
            else:
                print("⚠️ SSH key management not available. Please generate manually.")
        else:
            machine["auth"] = {"method": "key", "ssh_key": ssh_key, "port": machine["port"]}
            logger.info(f"Using existing SSH key: {ssh_key}")
            if AUTH_AVAILABLE:
                test = input("Test SSH connectivity with this key? [y/N]: ").strip().lower()
                if test == "y":
                    success = test_ssh_connection(machine["username"], machine["host"], ssh_key.replace(".pub", ""), machine["port"])
                    if success:
                        print("✅ SSH key is already installed and working.")
                        logger.info("SSH key validated successfully.")
                    else:
                        print("⚠️ SSH connection failed. You may need to install the key.")
                        install = input("Install this key on the remote server now? [y/N]: ").strip().lower()
                        if install == "y":
                            install_ssh_key(machine["username"], machine["host"], ssh_key, machine["port"])
            else:
                print("⚠️ SSH testing not available. Please test manually.")

        machine["auth"] = {"method": "key", "ssh_key": ssh_key, "port": machine["port"]}

    if machine["scheduler"] == "slurm":
        print("🧮 Define job resources (press Enter to skip any):")
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

    if machine["execution"] == "remote":
        nodes = machine["resources"].get("nodes", 1)
        ntasks = machine["resources"].get("ntasks-per-node", 1)
        machine["nprocs"] = nodes * ntasks
    else:
        nprocs_input = input(f"Number of processes [{machine.get('nprocs', 1)}]: ").strip()
        try:
            machine["nprocs"] = int(nprocs_input) if nprocs_input else machine.get("nprocs", 1)
        except ValueError:
            print("⚠️ Invalid input. Using default nprocs = 1.")
            machine["nprocs"] = 1

    print("🧭 Define the launcher command used to run Quantum ESPRESSO.")
    print("You may use the placeholder {nprocs}, which will be replaced at runtime.")
    print("Examples:")
    print(" - mpirun -np {nprocs}          (direct or manual MPI)")
    print(" - srun --mpi=pmi2              (Slurm with Intel MPI)")
    print("Note: For Slurm with Intel MPI, use 'srun --mpi=pmi2' without {nprocs}.")
    default_launcher = machine.get("launcher", "mpirun -np {nprocs}")
    launcher = input(f"Launcher command [{default_launcher}]: ").strip() or default_launcher
    machine["launcher"] = launcher
    logger.info(f"Launcher set to: {launcher}")

    # Ask user if they want to save as individual file or in machines.json
    print("\n💾 Choose how to save the machine configuration:")
    print(" [1] Add to machines.json (traditional, all machines in one file)")
    print(" [2] Save as individual JSON file (recommended, one file per machine)")
    save_choice = input("Choose save format [1/2] [2]: ").strip() or "2"
    
    if save_choice == "2":
        # Save as individual file
        machines_dir = os.path.join(os.path.dirname(path), "machines")
        os.makedirs(machines_dir, exist_ok=True)
        individual_path = os.path.join(machines_dir, f"{machine_name}.json")
        
        try:
            with open(individual_path, "w") as f:
                json.dump(machine, f, indent=2)
            print(f"✅ Machine '{machine_name}' saved to {individual_path}")
            logger.info(f"Machine '{machine_name}' saved as individual file: {individual_path}")
        except Exception as e:
            print("❌ Failed to save machine configuration.")
            logger.error(f"Failed to write individual machine file: {e}")
    else:
        # Save to machines.json (traditional)
        config["machines"][machine_name] = machine
        
        try:
            with open(path, "w") as f:
                json.dump(config, f, indent=2)
            print(f"✅ Machine '{machine_name}' saved to {path}")
            logger.info(f"Machine '{machine_name}' saved successfully.")
        except Exception as e:
            print("❌ Failed to save machine configuration.")
            logger.error(f"Failed to write config file: {e}")

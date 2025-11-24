#!/usr/bin/env python3
"""
Example: Using env_setup for code detection on remote systems

This example demonstrates how to use the env_setup parameter to make
the module command available when detecting Quantum ESPRESSO codes
via SSH on remote systems.
"""

from xespresso.codes import detect_qe_codes, create_codes_config
from xespresso.machines.machine import Machine

def example_basic_env_setup():
    """Basic usage of env_setup parameter"""
    print("=" * 60)
    print("Example 1: Basic env_setup usage")
    print("=" * 60)
    
    # Detect codes with env_setup to make module command available
    config = detect_qe_codes(
        machine_name="my_cluster",
        modules=["quantum-espresso/7.2"],
        env_setup="source /etc/profile",
        ssh_connection={
            'host': 'cluster.example.edu',
            'username': 'myusername',
            'port': 22
        },
        auto_load_machine=False  # Don't try to load machine config
    )
    
    print(f"Detected {len(config.codes)} codes")


def example_machine_with_env_setup():
    """Using env_setup in Machine configuration"""
    print("\n" + "=" * 60)
    print("Example 2: Machine configuration with env_setup")
    print("=" * 60)
    
    # Create a machine with env_setup
    machine = Machine(
        name="remote_cluster",
        execution="remote",
        host="cluster.example.edu",
        username="myusername",
        port=22,
        auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
        scheduler="slurm",
        workdir="/scratch/myusername",
        use_modules=True,
        modules=["intel-compiler/2021", "quantum-espresso/7.2"],
        env_setup="source /etc/profile.d/modules.sh"
    )
    
    # Save machine configuration
    machine.to_file("~/.xespresso/machines/remote_cluster.json")
    print(f"Machine '{machine.name}' configuration saved")
    
    # Convert to queue for use with Espresso calculator
    queue = machine.to_queue()
    print(f"Queue contains env_setup: {'env_setup' in queue}")


def example_auto_loading():
    """Automatic loading of env_setup from machine config"""
    print("\n" + "=" * 60)
    print("Example 3: Automatic env_setup loading")
    print("=" * 60)
    
    # When you have a machine config file with env_setup defined,
    # detect_qe_codes will automatically use it
    
    # This assumes you have a machine config at:
    # ~/.xespresso/machines/my_cluster.json
    # with "env_setup": "source /etc/profile.d/modules.sh"
    
    try:
        config = detect_qe_codes(
            machine_name="my_cluster",
            # env_setup will be loaded from machine config automatically
            auto_load_machine=True  # This is the default
        )
        print(f"Detected {len(config.codes)} codes using auto-loaded env_setup")
    except Exception as e:
        print(f"Note: This example requires a machine config file: {e}")


def example_multiple_env_sources():
    """Using multiple environment sources"""
    print("\n" + "=" * 60)
    print("Example 4: Multiple environment sources")
    print("=" * 60)
    
    # You can chain multiple source commands
    config = detect_qe_codes(
        machine_name="advanced_cluster",
        modules=["quantum-espresso/7.2"],
        env_setup="source /etc/profile && source /opt/intel/bin/compilervars.sh intel64 && source ~/.bashrc",
        ssh_connection={
            'host': 'cluster.example.edu',
            'username': 'myusername',
        },
        auto_load_machine=False
    )
    
    print(f"Detected {len(config.codes)} codes with complex env_setup")


def example_local_with_modules():
    """Using env_setup for local detection (less common)"""
    print("\n" + "=" * 60)
    print("Example 5: Local detection with env_setup")
    print("=" * 60)
    
    # Even for local detection, env_setup can be useful if you need
    # to set up the environment before running commands
    config = detect_qe_codes(
        machine_name="local",
        modules=["quantum-espresso/7.2"],
        env_setup="source /usr/share/modules/init/bash",
        auto_load_machine=False
    )
    
    print(f"Detected {len(config.codes)} codes locally with env_setup")


def example_create_and_save_with_env_setup():
    """Create and save config with env_setup"""
    print("\n" + "=" * 60)
    print("Example 6: Create and save config with env_setup")
    print("=" * 60)
    
    # Detect and automatically save the configuration
    config = create_codes_config(
        machine_name="production_cluster",
        qe_prefix="/opt/qe-7.2/bin",
        modules=["quantum-espresso/7.2"],
        env_setup="source /etc/profile",
        ssh_connection={
            'host': 'prod.cluster.edu',
            'username': 'produser',
        },
        save=True,  # Automatically save to file
        merge=True,  # Merge with existing config if it exists
        auto_load_machine=False
    )
    
    print(f"Configuration saved with {len(config.codes)} codes")


def example_machine_config_file():
    """Example machine configuration file with env_setup"""
    print("\n" + "=" * 60)
    print("Example 7: Machine configuration file format")
    print("=" * 60)
    
    config_example = """
{
  "name": "my_cluster",
  "execution": "remote",
  "host": "cluster.example.edu",
  "username": "myusername",
  "port": 22,
  "auth": {
    "method": "key",
    "ssh_key": "~/.ssh/id_rsa"
  },
  "scheduler": "slurm",
  "workdir": "/scratch/myusername",
  "use_modules": true,
  "modules": ["quantum-espresso/7.2"],
  "env_setup": "source /etc/profile.d/modules.sh",
  "resources": {
    "nodes": 1,
    "ntasks": 16,
    "time": "01:00:00",
    "partition": "regular"
  }
}
"""
    
    print("Save this to ~/.xespresso/machines/my_cluster.json:")
    print(config_example)
    print("\nThen use with:")
    print("  config = detect_qe_codes(machine_name='my_cluster')")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Environment Setup (env_setup) Examples")
    print("=" * 60)
    print("\nThese examples show how to use env_setup to make the module")
    print("command available when detecting codes via SSH.")
    print("\nNote: These examples won't actually run without proper SSH")
    print("      credentials and remote access. They demonstrate the API.")
    print("=" * 60)
    
    # Run examples (they will fail without proper setup, but show the API)
    try:
        example_basic_env_setup()
    except Exception as e:
        print(f"Example 1 error (expected): {type(e).__name__}")
    
    try:
        example_machine_with_env_setup()
    except Exception as e:
        print(f"Example 2 error: {type(e).__name__}")
    
    try:
        example_auto_loading()
    except Exception as e:
        print(f"Example 3 error (expected): {type(e).__name__}")
    
    try:
        example_multiple_env_sources()
    except Exception as e:
        print(f"Example 4 error (expected): {type(e).__name__}")
    
    try:
        example_local_with_modules()
    except Exception as e:
        print(f"Example 5 error: {type(e).__name__}")
    
    try:
        example_create_and_save_with_env_setup()
    except Exception as e:
        print(f"Example 6 error (expected): {type(e).__name__}")
    
    # This one always works as it just prints
    example_machine_config_file()
    
    print("\n" + "=" * 60)
    print("For more information, see docs/ENV_SETUP.md")
    print("=" * 60)

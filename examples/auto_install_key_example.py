#!/usr/bin/env python3
"""
Example: Automatic SSH Key Installation on Connection Failure

This example demonstrates the new auto_install_key feature that automatically
installs SSH keys on remote servers when authentication fails.

The feature is useful when:
- Setting up a new remote machine for the first time
- The SSH key hasn't been copied to the remote server yet
- You want to automate the setup process

Note: This requires that:
1. Password authentication is enabled on the remote server (for ssh-copy-id)
2. The public key file exists (e.g., ~/.ssh/id_rsa.pub)
3. You have the password available for initial authentication
"""

from xespresso.utils.auth import RemoteAuth, generate_ssh_key, install_ssh_key
from xespresso.machines.machine import Machine


def example_basic_auto_install():
    """Basic usage of auto_install_key feature"""
    print("=" * 60)
    print("Example 1: Basic auto_install_key usage")
    print("=" * 60)
    
    # Configure authentication with auto_install_key enabled
    auth_config = {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa",
        "port": 22,
        "auto_install_key": True  # Enable automatic key installation
    }
    
    # Create RemoteAuth instance
    remote = RemoteAuth(
        username="myusername",
        host="cluster.example.edu",
        auth_config=auth_config
    )
    
    # When you call connect(), if authentication fails:
    # 1. It will automatically try to install the SSH key using ssh-copy-id
    # 2. Then retry the connection
    # 3. If successful, you're connected!
    try:
        remote.connect()
        print(f"✅ Connected to {remote.host}")
        remote.close()
    except Exception as e:
        print(f"❌ Connection failed: {e}")


def example_machine_with_auto_install():
    """Using auto_install_key in Machine configuration"""
    print("\n" + "=" * 60)
    print("Example 2: Machine configuration with auto_install_key")
    print("=" * 60)
    
    # Create a machine with auto_install_key enabled
    machine = Machine(
        name="new_cluster",
        execution="remote",
        host="newcluster.example.edu",
        username="myusername",
        port=22,
        auth={
            "method": "key",
            "ssh_key": "~/.ssh/id_rsa",
            "auto_install_key": True  # Enable automatic key installation
        },
        scheduler="slurm",
        workdir="/scratch/myusername"
    )
    
    # Save machine configuration
    machine.to_file("~/.xespresso/machines/new_cluster.json")
    print(f"Machine '{machine.name}' configuration saved with auto_install_key enabled")
    
    # When you use this machine, it will automatically install
    # the SSH key if the first connection fails


def example_without_auto_install():
    """Default behavior without auto_install_key (backward compatible)"""
    print("\n" + "=" * 60)
    print("Example 3: Default behavior (auto_install_key=False)")
    print("=" * 60)
    
    # By default, auto_install_key is False for backward compatibility
    auth_config = {
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa",
        "port": 22,
        # auto_install_key defaults to False
    }
    
    remote = RemoteAuth(
        username="myusername",
        host="cluster.example.edu",
        auth_config=auth_config
    )
    
    print(f"auto_install_key is: {remote.auto_install_key}")
    print("Connection will fail immediately if authentication fails")


def example_manual_key_setup():
    """Manual SSH key setup (alternative approach)"""
    print("\n" + "=" * 60)
    print("Example 4: Manual SSH key setup")
    print("=" * 60)
    
    # You can also manually install keys before connecting
    ssh_key = "~/.ssh/id_rsa"
    
    print(f"To generate SSH key pair manually:")
    print(f"  ssh-keygen -t rsa -b 4096 -f {ssh_key}")
    
    print(f"\nTo manually install the key on remote server:")
    print(f"  ssh-copy-id -p 22 -i {ssh_key}.pub myusername@cluster.example.edu")
    
    # Now connect normally
    auth_config = {
        "method": "key",
        "ssh_key": ssh_key,
        "port": 22,
        "auto_install_key": False  # Not needed since we installed manually
    }
    
    remote = RemoteAuth(
        username="myusername",
        host="cluster.example.edu",
        auth_config=auth_config
    )
    print("\n✅ Remote auth configured for manual key setup")


def example_machine_config_file():
    """Example machine configuration file with auto_install_key"""
    print("\n" + "=" * 60)
    print("Example 5: Machine configuration file format")
    print("=" * 60)
    
    config_example = """
{
  "name": "auto_setup_cluster",
  "execution": "remote",
  "host": "cluster.example.edu",
  "username": "myusername",
  "port": 22,
  "auth": {
    "method": "key",
    "ssh_key": "~/.ssh/id_rsa",
    "auto_install_key": true
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
    
    print("Save this to ~/.xespresso/machines/auto_setup_cluster.json:")
    print(config_example)
    print("\nWith this configuration, the first time you connect:")
    print("1. If authentication fails, it will run ssh-copy-id automatically")
    print("2. You'll be prompted for your password (once)")
    print("3. Future connections will work seamlessly with the key")


def example_workflow_with_auto_install():
    """Complete workflow example with auto_install_key"""
    print("\n" + "=" * 60)
    print("Example 6: Complete workflow with auto_install_key")
    print("=" * 60)
    
    workflow_example = """
from xespresso import Espresso
from ase.build import bulk
from xespresso.machines.machine import Machine

# Create machine with auto_install_key
machine = Machine(
    name="production_cluster",
    execution="remote",
    host="prod.cluster.edu",
    username="produser",
    port=22,
    auth={
        "method": "key",
        "ssh_key": "~/.ssh/id_rsa",
        "auto_install_key": True  # First-time setup made easy!
    },
    scheduler="slurm",
    workdir="/scratch/produser"
)

# Convert to queue
queue = machine.to_queue()

# Setup calculation
atoms = bulk('Al', 'fcc', a=4.05)
pseudopotentials = {'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF'}

calc = Espresso(
    label='al_bulk/scf',
    pseudopotentials=pseudopotentials,
    queue=queue,
    # ... other parameters ...
)

# Run calculation
# On first run, if SSH key isn't installed:
# 1. Connection attempt fails
# 2. auto_install_key triggers ssh-copy-id
# 3. You enter password once
# 4. Connection succeeds
# 5. Calculation runs
# Future runs work seamlessly!
atoms.calc = calc
energy = atoms.get_potential_energy()
"""
    
    print("Complete workflow code:")
    print(workflow_example)


def example_security_considerations():
    """Security considerations for auto_install_key"""
    print("\n" + "=" * 60)
    print("Example 7: Security Considerations")
    print("=" * 60)
    
    print("\n⚠️ Important Security Notes:")
    print("1. auto_install_key requires password authentication to be enabled")
    print("   on the remote server (for ssh-copy-id to work)")
    print("\n2. You will be prompted for your password when ssh-copy-id runs")
    print("   - This is normal and expected")
    print("   - The password is entered interactively and not stored")
    print("\n3. After key installation, password authentication can be disabled")
    print("   on the remote server for enhanced security")
    print("\n4. The feature is disabled by default (auto_install_key=False)")
    print("   - Maintains backward compatibility")
    print("   - Must be explicitly enabled")
    print("\n5. Public key must exist before auto-install can work")
    print("   - Generate with: ssh-keygen -t rsa -b 4096")
    print("   - Default location: ~/.ssh/id_rsa and ~/.ssh/id_rsa.pub")
    print("\n✅ Best Practice:")
    print("   - Enable auto_install_key for initial setup")
    print("   - Disable it after successful key installation")
    print("   - Use strong SSH keys (RSA 4096-bit or Ed25519)")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("Automatic SSH Key Installation Examples")
    print("=" * 60)
    print("\nThese examples show how to use the auto_install_key feature")
    print("to automatically install SSH keys when authentication fails.")
    print("\nNote: These examples won't actually run without proper SSH")
    print("      credentials and remote access. They demonstrate the API.")
    print("=" * 60)
    
    # Run examples
    try:
        example_basic_auto_install()
    except Exception as e:
        print(f"Example 1 error (expected): {type(e).__name__}")
    
    try:
        example_machine_with_auto_install()
    except Exception as e:
        print(f"Example 2 error: {type(e).__name__}")
    
    try:
        example_without_auto_install()
    except Exception as e:
        print(f"Example 3 error: {type(e).__name__}")
    
    try:
        example_manual_key_setup()
    except Exception as e:
        print(f"Example 4 error: {type(e).__name__}")
    
    # These always work as they just print
    example_machine_config_file()
    example_workflow_with_auto_install()
    example_security_considerations()
    
    print("\n" + "=" * 60)
    print("For more information, see:")
    print("- xespresso/utils/auth.py (RemoteAuth class)")
    print("- tests/test_auto_install_key.py (test examples)")
    print("=" * 60)

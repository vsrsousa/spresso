#!/usr/bin/env python
"""
Example: Using the new Machine class and modular configuration

This script demonstrates the new features:
1. Machine class for type-safe configuration
2. Modular configuration loading
3. Connection persistence verification
"""

import tempfile
import os
from xespresso.machines import Machine, load_machine
from xespresso.machines.config import list_machines

def example_1_machine_class():
    """Example 1: Creating and using Machine class"""
    print("\n" + "="*60)
    print("Example 1: Machine Class")
    print("="*60)
    
    # Create a local machine
    local_machine = Machine(
        name="my_local",
        execution="local",
        scheduler="direct",
        workdir="./calculations",
        nprocs=4
    )
    
    print(f"\nCreated machine: {local_machine}")
    print(f"Is remote? {local_machine.is_remote}")
    print(f"Is local? {local_machine.is_local}")
    
    # Convert to queue format (backward compatible)
    queue = local_machine.to_queue()
    print(f"\nQueue format: {queue['execution']}, {queue['scheduler']}")
    print(f"Local dir: {queue['local_dir']}")
    
    return local_machine


def example_2_modular_config():
    """Example 2: Modular configuration with individual files"""
    print("\n" + "="*60)
    print("Example 2: Modular Configuration")
    print("="*60)
    
    # Create temporary directory structure
    with tempfile.TemporaryDirectory() as tmpdir:
        machines_dir = os.path.join(tmpdir, "machines")
        os.makedirs(machines_dir)
        
        # Create individual machine files
        machine1 = Machine(
            name="cluster1",
            execution="remote",
            scheduler="slurm",
            host="cluster1.edu",
            username="user",
            workdir="/home/user/calc",
            auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            nprocs=16
        )
        machine1.to_file(os.path.join(machines_dir, "cluster1.json"))
        
        machine2 = Machine(
            name="cluster2",
            execution="remote",
            scheduler="slurm",
            host="cluster2.edu",
            username="user",
            workdir="/home/user/calc",
            auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            nprocs=32
        )
        machine2.to_file(os.path.join(machines_dir, "cluster2.json"))
        
        # Create default.json
        import json
        with open(os.path.join(machines_dir, "default.json"), 'w') as f:
            json.dump({"default": "cluster1"}, f)
        
        print(f"\nCreated machines directory: {machines_dir}")
        print("Files created:")
        for f in os.listdir(machines_dir):
            print(f"  - {f}")
        
        # List all machines
        machines = list_machines(machines_dir=machines_dir)
        print(f"\nAvailable machines: {machines}")
        
        # Load specific machine
        loaded = load_machine(
            machines_dir=machines_dir,
            machine_name="cluster1",
            return_object=True
        )
        print(f"\nLoaded machine: {loaded}")
        print(f"Host: {loaded.host}")
        print(f"Processors: {loaded.nprocs}")


def example_3_connection_persistence():
    """Example 3: Understanding connection persistence"""
    print("\n" + "="*60)
    print("Example 3: Connection Persistence")
    print("="*60)
    
    print("""
Connection persistence is AUTOMATIC and TRANSPARENT!

How it works:
1. First calculation creates connection to (host, user)
2. Connection is cached in RemoteExecutionMixin._remote_sessions
3. Subsequent calculations REUSE the cached connection
4. Only creates new connection for different (host, user)

Example flow:
    queue = load_machine("cluster1")
    
    # First job
    scheduler1 = get_scheduler(calc1, queue, command1)
    scheduler1.run()  # Creates connection to cluster1
    
    # Second job - REUSES connection!
    scheduler2 = get_scheduler(calc2, queue, command2)
    scheduler2.run()  # No new connection, just submits job
    
    # Third job - REUSES same connection!
    scheduler3 = get_scheduler(calc3, queue, command3)
    scheduler3.run()  # No new connection, just submits job

Benefits:
- ⚡ Faster job submission (no SSH handshake)
- 🎯 More reliable (established connection)
- 💪 Less server load (fewer connections)

See docs/REMOTE_CONNECTION_PERSISTENCE.md for details.
    """)


def example_4_migration():
    """Example 4: Migrating from old to new format"""
    print("\n" + "="*60)
    print("Example 4: Migration Guide")
    print("="*60)
    
    print("""
Old way (still works!):
    queue = load_machine("cluster1")
    # Returns dict with all configuration

New way (recommended):
    # Option 1: Get Machine object
    machine = load_machine("cluster1", return_object=True)
    queue = machine.to_queue()
    
    # Option 2: Create Machine directly
    machine = Machine(
        name="my_cluster",
        execution="remote",
        host="cluster.edu",
        username="user",
        ...
    )
    queue = machine.to_queue()

Benefits of new way:
- ✅ Type safety and validation
- ✅ Better IDE support
- ✅ Clearer error messages
- ✅ Easier to test
- ✅ Self-documenting code

No breaking changes! Choose when to migrate.
    """)


def main():
    """Run all examples"""
    print("\n" + "="*60)
    print("XEspresso Machine Configuration Examples")
    print("="*60)
    
    try:
        example_1_machine_class()
        example_2_modular_config()
        example_3_connection_persistence()
        example_4_migration()
        
        print("\n" + "="*60)
        print("✅ All examples completed successfully!")
        print("="*60)
        print("\nNext steps:")
        print("1. See docs/MACHINE_CONFIGURATION.md for complete guide")
        print("2. See docs/REMOTE_CONNECTION_PERSISTENCE.md for persistence details")
        print("3. Check examples/ directory for configuration templates")
        print()
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

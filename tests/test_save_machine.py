"""
Tests for the save_machine function in loader.py.

This test suite verifies that the save_machine function works correctly
with both individual files and machines.json format.
"""

import pytest
import tempfile
import os
import json
from pathlib import Path

from xespresso.machines.machine import Machine
from xespresso.machines.config.loader import (
    save_machine, load_machine, list_machines,
    DEFAULT_CONFIG_PATH, DEFAULT_MACHINES_DIR
)


class TestSaveMachine:
    """Test suite for save_machine function."""
    
    def test_save_machine_individual_file(self):
        """Test saving a machine as an individual file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Create test machine
            machine = Machine(
                name="test_cluster",
                execution="remote",
                scheduler="slurm",
                workdir="/scratch/work",
                host="cluster.example.com",
                username="user123",
                auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
                nprocs=48
            )
            
            # Save as individual file
            filepath = save_machine(machine, config_path, machines_dir, use_individual_file=True)
            
            # Verify file exists
            expected_path = os.path.join(machines_dir, "test_cluster.json")
            assert filepath == expected_path, f"Expected {expected_path}, got {filepath}"
            assert os.path.exists(filepath), f"File not created: {filepath}"
            
            # Verify content
            with open(filepath) as f:
                data = json.load(f)
            
            assert data['name'] == 'test_cluster'
            assert data['execution'] == 'remote'
            assert data['scheduler'] == 'slurm'
            assert data['host'] == 'cluster.example.com'
    
    def test_save_machine_to_machines_json(self):
        """Test saving a machine to machines.json."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Create test machine
            machine = Machine(
                name="local_workstation",
                execution="local",
                scheduler="direct",
                workdir="/home/user/calculations",
                nprocs=8
            )
            
            # Save to machines.json
            filepath = save_machine(machine, config_path, machines_dir, use_individual_file=False)
            
            # Verify file exists
            assert filepath == config_path
            assert os.path.exists(config_path)
            
            # Verify content
            with open(config_path) as f:
                config = json.load(f)
            
            assert 'machines' in config
            assert 'local_workstation' in config['machines']
            assert config['machines']['local_workstation']['execution'] == 'local'
    
    def test_save_and_load_machine(self):
        """Test the full save and load cycle."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Create and save a machine
            original = Machine(
                name="integration_test",
                execution="local",
                scheduler="slurm",
                workdir="/tmp/work",
                nprocs=16,
                launcher="srun --mpi=pmi2",
                use_modules=True,
                modules=["gcc/11.2.0", "openmpi/4.1.0"]
            )
            
            save_machine(original, config_path, machines_dir, use_individual_file=True)
            
            # Load it back
            loaded = load_machine(config_path, "integration_test", machines_dir, return_object=True)
            
            # Verify properties match
            assert loaded.name == original.name
            assert loaded.execution == original.execution
            assert loaded.scheduler == original.scheduler
            assert loaded.workdir == original.workdir
            assert loaded.nprocs == original.nprocs
            assert loaded.launcher == original.launcher
            assert loaded.use_modules == original.use_modules
            assert loaded.modules == original.modules
    
    def test_list_machines_after_save(self):
        """Test that list_machines includes saved machines."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Initially empty
            machines = list_machines(config_path, machines_dir)
            assert len(machines) == 0
            
            # Save a machine as individual file
            m1 = Machine(name="machine1", execution="local", scheduler="direct", workdir="/tmp")
            save_machine(m1, config_path, machines_dir, use_individual_file=True)
            
            # Save another to machines.json
            m2 = Machine(name="machine2", execution="local", scheduler="direct", workdir="/tmp")
            save_machine(m2, config_path, machines_dir, use_individual_file=False)
            
            # List should include both
            machines = list_machines(config_path, machines_dir)
            assert len(machines) == 2
            assert 'machine1' in machines
            assert 'machine2' in machines
    
    def test_save_machine_overwrites_existing(self):
        """Test that saving a machine with the same name overwrites the existing one."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Save initial machine
            m1 = Machine(
                name="test_machine",
                execution="local",
                scheduler="direct",
                workdir="/tmp/old",
                nprocs=4
            )
            save_machine(m1, config_path, machines_dir, use_individual_file=True)
            
            # Save updated machine with same name
            m2 = Machine(
                name="test_machine",
                execution="local",
                scheduler="slurm",
                workdir="/tmp/new",
                nprocs=8
            )
            save_machine(m2, config_path, machines_dir, use_individual_file=True)
            
            # Load and verify it's the updated version
            loaded = load_machine(config_path, "test_machine", machines_dir, return_object=True)
            assert loaded.workdir == "/tmp/new"
            assert loaded.nprocs == 8
            assert loaded.scheduler == "slurm"
    
    def test_save_machine_invalid_input(self):
        """Test that save_machine raises error for invalid input."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Try to save non-Machine object
            with pytest.raises(ValueError, match="Expected Machine object"):
                save_machine({"not": "a machine"}, config_path, machines_dir)
    
    def test_save_machine_with_complex_configuration(self):
        """Test saving a machine with complex configuration."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            machines_dir = os.path.join(tmpdir, "machines")
            
            # Create complex machine
            machine = Machine(
                name="complex_cluster",
                execution="remote",
                scheduler="slurm",
                workdir="/scratch/project/calculations",
                host="hpc.university.edu",
                username="researcher",
                port=2222,
                auth={"method": "key", "ssh_key": "~/.ssh/cluster_key"},
                nprocs=128,
                launcher="srun --mpi=pmi2",
                use_modules=True,
                modules=["intel/2023.1", "impi/2021.9", "qe/7.2"],
                prepend=["export OMP_NUM_THREADS=1", "ulimit -s unlimited"],
                postpend=["echo 'Job completed'"],
                resources={
                    "nodes": 4,
                    "ntasks-per-node": 32,
                    "time": "48:00:00",
                    "partition": "compute",
                    "account": "project123"
                },
                env_setup="source /etc/profile.d/modules.sh"
            )
            
            # Save and load
            save_machine(machine, config_path, machines_dir, use_individual_file=True)
            loaded = load_machine(config_path, "complex_cluster", machines_dir, return_object=True)
            
            # Verify all properties
            assert loaded.name == machine.name
            assert loaded.execution == machine.execution
            assert loaded.scheduler == machine.scheduler
            assert loaded.host == machine.host
            assert loaded.username == machine.username
            assert loaded.port == machine.port
            assert loaded.nprocs == machine.nprocs
            assert loaded.launcher == machine.launcher
            assert loaded.use_modules == machine.use_modules
            assert loaded.modules == machine.modules
            assert loaded.prepend == machine.prepend
            assert loaded.postpend == machine.postpend
            assert loaded.resources == machine.resources
            assert loaded.env_setup == machine.env_setup


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

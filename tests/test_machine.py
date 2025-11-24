"""
test_machine.py

Tests for the Machine class and configuration loading.
"""

import os
import json
import tempfile
import pytest
from xespresso.machines.machine import Machine
from xespresso.machines.config.loader import load_machine, list_machines, _get_default_machine_name


class TestMachine:
    """Test suite for Machine class."""
    
    def test_machine_creation_local(self):
        """Test creating a local machine."""
        machine = Machine(
            name="local_test",
            execution="local",
            scheduler="direct",
            workdir="/tmp/test"
        )
        
        assert machine.name == "local_test"
        assert machine.execution == "local"
        assert machine.is_local
        assert not machine.is_remote
        assert machine.scheduler == "direct"
        assert machine.workdir == "/tmp/test"
    
    def test_machine_creation_remote(self):
        """Test creating a remote machine."""
        machine = Machine(
            name="remote_test",
            execution="remote",
            scheduler="slurm",
            workdir="/home/user/calc",
            host="cluster.example.com",
            username="user",
            auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            port=22
        )
        
        assert machine.name == "remote_test"
        assert machine.execution == "remote"
        assert machine.is_remote
        assert not machine.is_local
        assert machine.host == "cluster.example.com"
        assert machine.username == "user"
    
    def test_machine_validation_remote_missing_host(self):
        """Test that remote machine requires host."""
        with pytest.raises(ValueError, match="requires 'host'"):
            Machine(
                name="invalid",
                execution="remote",
                username="user"
            )
    
    def test_machine_validation_remote_missing_username(self):
        """Test that remote machine requires username."""
        with pytest.raises(ValueError, match="requires 'username'"):
            Machine(
                name="invalid",
                execution="remote",
                host="example.com"
            )
    
    def test_machine_from_dict(self):
        """Test creating machine from dictionary."""
        config = {
            "execution": "local",
            "scheduler": "direct",
            "workdir": "/tmp/test",
            "nprocs": 4
        }
        
        machine = Machine.from_dict("test_machine", config)
        
        assert machine.name == "test_machine"
        assert machine.execution == "local"
        assert machine.nprocs == 4
    
    def test_machine_to_dict(self):
        """Test converting machine to dictionary."""
        machine = Machine(
            name="test",
            execution="local",
            scheduler="direct",
            workdir="/tmp",
            nprocs=8
        )
        
        config = machine.to_dict()
        
        assert config["execution"] == "local"
        assert config["scheduler"] == "direct"
        assert config["workdir"] == "/tmp"
        assert config["nprocs"] == 8
    
    def test_machine_to_queue_local(self):
        """Test converting local machine to queue dict."""
        machine = Machine(
            name="test",
            execution="local",
            scheduler="direct",
            workdir="/tmp/calc",
            prepend=["export PATH=/opt/qe/bin:$PATH"],
            postpend=["echo done"]
        )
        
        queue = machine.to_queue()
        
        assert queue["execution"] == "local"
        assert queue["scheduler"] == "direct"
        assert queue["local_dir"] == "/tmp/calc"
        assert "export PATH" in queue["prepend"]
        assert "echo done" in queue["postpend"]
    
    def test_machine_to_queue_remote(self):
        """Test converting remote machine to queue dict."""
        machine = Machine(
            name="test",
            execution="remote",
            scheduler="slurm",
            workdir="/home/user/calc",
            host="cluster.edu",
            username="user",
            auth={"method": "key", "ssh_key": "~/.ssh/id_rsa"}
        )
        
        queue = machine.to_queue()
        
        assert queue["execution"] == "remote"
        assert queue["remote_host"] == "cluster.edu"
        assert queue["remote_user"] == "user"
        assert queue["remote_dir"] == "/home/user/calc"
        assert queue["remote_auth"]["method"] == "key"
    
    def test_machine_to_queue_script_normalization(self):
        """Test that script blocks are normalized to strings."""
        machine = Machine(
            name="test",
            execution="local",
            prepend=["line1", "line2", "line3"],
            postpend=["cleanup1", "cleanup2"]
        )
        
        queue = machine.to_queue()
        
        assert isinstance(queue["prepend"], str)
        assert isinstance(queue["postpend"], str)
        assert "line1\nline2\nline3" == queue["prepend"]
        assert "cleanup1\ncleanup2" == queue["postpend"]
    
    def test_machine_file_operations(self):
        """Test saving and loading machine from file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_machine.json")
            
            # Create and save
            machine1 = Machine(
                name="test",
                execution="local",
                scheduler="direct",
                workdir="/tmp",
                nprocs=4
            )
            machine1.to_file(filepath)
            
            # Load and verify
            machine2 = Machine.from_file(filepath)
            
            assert machine2.name == "test"
            assert machine2.execution == "local"
            assert machine2.scheduler == "direct"
            assert machine2.nprocs == 4


class TestMachineLoader:
    """Test suite for machine loading functions."""
    
    def test_load_from_single_machines_json(self):
        """Test loading machine from traditional machines.json."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            
            config = {
                "machines": {
                    "test_machine": {
                        "execution": "local",
                        "scheduler": "direct",
                        "workdir": "/tmp/test",
                        "nprocs": 4
                    }
                }
            }
            
            with open(config_path, 'w') as f:
                json.dump(config, f)
            
            # Load as dict
            queue = load_machine(config_path, "test_machine")
            assert queue is not None
            assert queue["execution"] == "local"
            assert queue["local_dir"] == "/tmp/test"
            
            # Load as object
            machine = load_machine(config_path, "test_machine", return_object=True)
            assert machine is not None
            assert isinstance(machine, Machine)
            assert machine.name == "test_machine"
    
    def test_load_from_individual_file(self):
        """Test loading machine from individual JSON file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            machines_dir = os.path.join(tmpdir, "machines")
            os.makedirs(machines_dir)
            
            machine_file = os.path.join(machines_dir, "cluster1.json")
            config = {
                "execution": "remote",
                "scheduler": "slurm",
                "workdir": "/home/user/calc",
                "host": "cluster.edu",
                "username": "user",
                "auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"}
            }
            
            with open(machine_file, 'w') as f:
                json.dump(config, f)
            
            # Load machine
            machine = load_machine(
                config_path=os.path.join(tmpdir, "machines.json"),
                machine_name="cluster1",
                machines_dir=machines_dir,
                return_object=True
            )
            
            assert machine is not None
            assert isinstance(machine, Machine)
            assert machine.name == "cluster1"
            assert machine.is_remote
    
    def test_default_machine_from_machines_json(self):
        """Test loading default machine specified in machines.json."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "machines.json")
            
            config = {
                "default": "my_cluster",
                "machines": {
                    "my_cluster": {
                        "execution": "local",
                        "scheduler": "direct",
                        "workdir": "/tmp"
                    }
                }
            }
            
            with open(config_path, 'w') as f:
                json.dump(config, f)
            
            # Should use configured default
            default_name = _get_default_machine_name(config_path)
            assert default_name == "my_cluster"
    
    def test_list_machines_combined(self):
        """Test listing machines from both sources."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create machines.json
            config_path = os.path.join(tmpdir, "machines.json")
            config = {
                "machines": {
                    "machine1": {"execution": "local"},
                    "machine2": {"execution": "local"}
                }
            }
            with open(config_path, 'w') as f:
                json.dump(config, f)
            
            # Create individual machine files
            machines_dir = os.path.join(tmpdir, "machines")
            os.makedirs(machines_dir)
            
            for name in ["machine3", "machine4"]:
                machine_file = os.path.join(machines_dir, f"{name}.json")
                with open(machine_file, 'w') as f:
                    json.dump({"execution": "local"}, f)
            
            # List all machines
            machines = list_machines(config_path, machines_dir)
            
            assert len(machines) == 4
            assert "machine1" in machines
            assert "machine2" in machines
            assert "machine3" in machines
            assert "machine4" in machines


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

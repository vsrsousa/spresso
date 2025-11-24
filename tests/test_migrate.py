"""
test_migrate.py

Tests for machine configuration migration functionality.
"""

import os
import json
import tempfile
import pytest
from xespresso.machines.config.migrate import migrate_machines, rollback_migration


class TestMigrateMachines:
    """Test suite for migrate_machines function."""
    
    def test_migrate_all_machines(self):
        """Test migrating all machines from machines.json."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create machines.json
            machines_json = os.path.join(tmpdir, "machines.json")
            config = {
                "default": "local_test",
                "machines": {
                    "local_test": {
                        "execution": "local",
                        "scheduler": "direct",
                        "workdir": "/tmp/test",
                        "nprocs": 4
                    },
                    "cluster_test": {
                        "execution": "remote",
                        "scheduler": "slurm",
                        "workdir": "/home/user/calc",
                        "host": "cluster.edu",
                        "username": "user",
                        "auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"}
                    }
                }
            }
            
            with open(machines_json, 'w') as f:
                json.dump(config, f)
            
            # Migrate
            output_dir = os.path.join(tmpdir, "machines")
            results = migrate_machines(
                machines_json_path=machines_json,
                output_dir=output_dir
            )
            
            # Verify results
            assert results["success"]
            assert len(results["migrated"]) == 2
            assert "local_test" in results["migrated"]
            assert "cluster_test" in results["migrated"]
            assert len(results["failed"]) == 0
            assert results["default"] == "local_test"
            
            # Verify files created
            assert os.path.exists(os.path.join(output_dir, "local_test.json"))
            assert os.path.exists(os.path.join(output_dir, "cluster_test.json"))
            assert os.path.exists(os.path.join(output_dir, "default.json"))
            
            # Verify content
            with open(os.path.join(output_dir, "local_test.json")) as f:
                local_config = json.load(f)
            assert local_config["execution"] == "local"
            assert local_config["nprocs"] == 4
            
            with open(os.path.join(output_dir, "default.json")) as f:
                default_config = json.load(f)
            assert default_config["default"] == "local_test"
    
    def test_migrate_specific_machines(self):
        """Test migrating only specific machines."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create machines.json with 3 machines
            machines_json = os.path.join(tmpdir, "machines.json")
            config = {
                "machines": {
                    "machine1": {"execution": "local"},
                    "machine2": {"execution": "local"},
                    "machine3": {"execution": "local"}
                }
            }
            
            with open(machines_json, 'w') as f:
                json.dump(config, f)
            
            # Migrate only machine1 and machine3
            output_dir = os.path.join(tmpdir, "machines")
            results = migrate_machines(
                machines_json_path=machines_json,
                output_dir=output_dir,
                machine_names=["machine1", "machine3"]
            )
            
            # Verify results
            assert results["success"]
            assert len(results["migrated"]) == 2
            assert "machine1" in results["migrated"]
            assert "machine3" in results["migrated"]
            
            # Verify only specified files created
            assert os.path.exists(os.path.join(output_dir, "machine1.json"))
            assert os.path.exists(os.path.join(output_dir, "machine3.json"))
            assert not os.path.exists(os.path.join(output_dir, "machine2.json"))
    
    def test_migrate_skip_existing(self):
        """Test that existing files are skipped by default."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create machines.json
            machines_json = os.path.join(tmpdir, "machines.json")
            config = {
                "machines": {
                    "machine1": {"execution": "local", "nprocs": 4},
                    "machine2": {"execution": "local", "nprocs": 8}
                }
            }
            
            with open(machines_json, 'w') as f:
                json.dump(config, f)
            
            # Create existing file for machine1
            output_dir = os.path.join(tmpdir, "machines")
            os.makedirs(output_dir)
            existing_file = os.path.join(output_dir, "machine1.json")
            with open(existing_file, 'w') as f:
                json.dump({"execution": "local", "nprocs": 99}, f)
            
            # Migrate
            results = migrate_machines(
                machines_json_path=machines_json,
                output_dir=output_dir,
                overwrite=False
            )
            
            # Verify results
            assert results["success"]
            assert len(results["migrated"]) == 1
            assert "machine2" in results["migrated"]
            assert len(results["skipped"]) == 1
            assert "machine1" in results["skipped"]
            
            # Verify machine1 was not overwritten
            with open(existing_file) as f:
                content = json.load(f)
            assert content["nprocs"] == 99  # Original value preserved
    
    def test_migrate_overwrite_existing(self):
        """Test that existing files are overwritten when specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create machines.json
            machines_json = os.path.join(tmpdir, "machines.json")
            config = {
                "machines": {
                    "machine1": {"execution": "local", "nprocs": 4}
                }
            }
            
            with open(machines_json, 'w') as f:
                json.dump(config, f)
            
            # Create existing file
            output_dir = os.path.join(tmpdir, "machines")
            os.makedirs(output_dir)
            existing_file = os.path.join(output_dir, "machine1.json")
            with open(existing_file, 'w') as f:
                json.dump({"execution": "local", "nprocs": 99}, f)
            
            # Migrate with overwrite=True
            results = migrate_machines(
                machines_json_path=machines_json,
                output_dir=output_dir,
                overwrite=True
            )
            
            # Verify results
            assert results["success"]
            assert len(results["migrated"]) == 1
            assert "machine1" in results["migrated"]
            assert len(results["skipped"]) == 0
            
            # Verify machine1 was overwritten
            with open(existing_file) as f:
                content = json.load(f)
            assert content["nprocs"] == 4  # New value
    
    def test_migrate_missing_file(self):
        """Test migration with missing machines.json."""
        with tempfile.TemporaryDirectory() as tmpdir:
            nonexistent = os.path.join(tmpdir, "missing.json")
            output_dir = os.path.join(tmpdir, "machines")
            
            results = migrate_machines(
                machines_json_path=nonexistent,
                output_dir=output_dir
            )
            
            # Verify failure
            assert not results["success"]
            assert "file_not_found" in results["errors"]
    
    def test_migrate_no_machines(self):
        """Test migration with empty machines.json."""
        with tempfile.TemporaryDirectory() as tmpdir:
            machines_json = os.path.join(tmpdir, "machines.json")
            config = {"machines": {}}
            
            with open(machines_json, 'w') as f:
                json.dump(config, f)
            
            output_dir = os.path.join(tmpdir, "machines")
            results = migrate_machines(
                machines_json_path=machines_json,
                output_dir=output_dir
            )
            
            # Verify appropriate handling
            assert not results["success"]
            assert "no_machines" in results["errors"]


class TestRollbackMigration:
    """Test suite for rollback_migration function."""
    
    def test_rollback_all_machines(self):
        """Test rolling back all migrated machines."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "machines")
            os.makedirs(output_dir)
            
            # Create individual machine files
            for name in ["machine1", "machine2", "machine3"]:
                filepath = os.path.join(output_dir, f"{name}.json")
                with open(filepath, 'w') as f:
                    json.dump({"execution": "local"}, f)
            
            # Rollback
            results = rollback_migration(output_dir=output_dir)
            
            # Verify results
            assert results["success"]
            assert len(results["removed"]) == 3
            assert "machine1" in results["removed"]
            assert "machine2" in results["removed"]
            assert "machine3" in results["removed"]
            
            # Verify files removed
            assert not os.path.exists(os.path.join(output_dir, "machine1.json"))
            assert not os.path.exists(os.path.join(output_dir, "machine2.json"))
            assert not os.path.exists(os.path.join(output_dir, "machine3.json"))
    
    def test_rollback_specific_machines(self):
        """Test rolling back only specific machines."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "machines")
            os.makedirs(output_dir)
            
            # Create individual machine files
            for name in ["machine1", "machine2", "machine3"]:
                filepath = os.path.join(output_dir, f"{name}.json")
                with open(filepath, 'w') as f:
                    json.dump({"execution": "local"}, f)
            
            # Rollback only machine1 and machine3
            results = rollback_migration(
                output_dir=output_dir,
                machine_names=["machine1", "machine3"]
            )
            
            # Verify results
            assert results["success"]
            assert len(results["removed"]) == 2
            assert "machine1" in results["removed"]
            assert "machine3" in results["removed"]
            
            # Verify only specified files removed
            assert not os.path.exists(os.path.join(output_dir, "machine1.json"))
            assert os.path.exists(os.path.join(output_dir, "machine2.json"))
            assert not os.path.exists(os.path.join(output_dir, "machine3.json"))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

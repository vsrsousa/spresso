"""
Test that machine configurations properly sanitize embedded quotes in JSON values.

This test verifies that quotes embedded in JSON string values (e.g., '"export VAR=1"')
are properly removed when loading machine configurations, preventing quotes from
appearing in generated job scripts.
"""
import sys
import os
import json
import tempfile
import shutil

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_sanitize_string_value():
    """Test that string values with embedded quotes are sanitized."""
    from xespresso.machines.config.loader import sanitize_string_value
    
    # Test double quotes
    assert sanitize_string_value('"export OMP_NUM_THREADS=1"') == 'export OMP_NUM_THREADS=1'
    
    # Test single quotes
    assert sanitize_string_value("'export OMP_NUM_THREADS=1'") == 'export OMP_NUM_THREADS=1'
    
    # Test no quotes
    assert sanitize_string_value('export OMP_NUM_THREADS=1') == 'export OMP_NUM_THREADS=1'
    
    # Test with whitespace
    assert sanitize_string_value('  "export VAR=1"  ') == 'export VAR=1'
    
    # Test None
    assert sanitize_string_value(None) is None
    
    print("✓ String value sanitization works correctly")
    return True


def test_sanitize_list_values():
    """Test that list values with embedded quotes are sanitized."""
    from xespresso.machines.config.loader import sanitize_list_values
    
    # Test list with quoted strings
    input_list = ['"export OMP_NUM_THREADS=1"', '"module load gcc"']
    expected = ['export OMP_NUM_THREADS=1', 'module load gcc']
    assert sanitize_list_values(input_list) == expected
    
    # Test mixed list
    input_list = ['"quoted"', 'not quoted', '  "spaced"  ']
    expected = ['quoted', 'not quoted', 'spaced']
    assert sanitize_list_values(input_list) == expected
    
    # Test empty list
    assert sanitize_list_values([]) == []
    
    # Test non-list
    assert sanitize_list_values("not a list") == "not a list"
    
    print("✓ List value sanitization works correctly")
    return True


def test_machine_config_sanitization():
    """Test that machine configurations are properly sanitized when loaded."""
    from xespresso.machines.config.loader import load_machine
    from xespresso.schedulers.slurm import SlurmScheduler
    
    # Create a config file with embedded quotes
    config_dir = tempfile.mkdtemp()
    config_file = os.path.join(config_dir, "machines.json")
    
    try:
        config_with_quotes = {
            "machines": {
                "test_machine": {
                    "name": "test_machine",
                    "execution": "local",
                    "scheduler": "slurm",
                    "workdir": "./work",
                    "nprocs": 16,
                    "launcher": '"srun --mpi=pmi2"',  # Embedded quotes
                    "use_modules": True,
                    "modules": ['"quantum-espresso/7.4.1"'],  # Embedded quotes in list
                    "prepend": ['"export OMP_NUM_THREADS=1"'],  # Embedded quotes in list
                    "postpend": ['"echo Done"'],  # Embedded quotes in list
                    "resources": {
                        "nodes": 1,
                        "ntasks-per-node": 16,
                        "mem": "48G",
                        "time": "01:00:00",
                        "partition": "parallel"
                    }
                }
            }
        }
        
        # Save config
        with open(config_file, "w") as f:
            json.dump(config_with_quotes, f, indent=2)
        
        # Load machine
        machine = load_machine(config_file, "test_machine", return_object=True)
        
        # Check that quotes were removed
        assert machine.launcher == 'srun --mpi=pmi2', f"Expected 'srun --mpi=pmi2', got {repr(machine.launcher)}"
        assert machine.prepend == ['export OMP_NUM_THREADS=1'], f"Expected ['export OMP_NUM_THREADS=1'], got {repr(machine.prepend)}"
        assert machine.postpend == ['echo Done'], f"Expected ['echo Done'], got {repr(machine.postpend)}"
        assert machine.modules == ['quantum-espresso/7.4.1'], f"Expected ['quantum-espresso/7.4.1'], got {repr(machine.modules)}"
        
        print("✓ Machine configuration sanitization works correctly")
        return True
        
    finally:
        shutil.rmtree(config_dir, ignore_errors=True)


def test_job_file_no_quotes():
    """Test that generated job files don't contain unwanted quotes."""
    from xespresso.machines.config.loader import load_machine
    from xespresso.schedulers.slurm import SlurmScheduler
    
    # Create a config file with embedded quotes
    config_dir = tempfile.mkdtemp()
    config_file = os.path.join(config_dir, "machines.json")
    
    # Mock calculator
    class MockCalc:
        def __init__(self):
            self.prefix = "espresso"
            self.directory = tempfile.mkdtemp()
    
    try:
        config_with_quotes = {
            "machines": {
                "test_machine": {
                    "name": "test_machine",
                    "execution": "local",
                    "scheduler": "slurm",
                    "workdir": "./work",
                    "nprocs": 16,
                    "launcher": '"srun --mpi=pmi2"',
                    "use_modules": True,
                    "modules": ["quantum-espresso/7.4.1"],
                    "prepend": ['"export OMP_NUM_THREADS=1"'],
                    "resources": {
                        "nodes": 1,
                        "ntasks-per-node": 16,
                        "mem": "48G",
                        "time": "01:00:00",
                        "partition": "parallel"
                    }
                }
            }
        }
        
        # Save config
        with open(config_file, "w") as f:
            json.dump(config_with_quotes, f, indent=2)
        
        # Load machine and create queue
        machine = load_machine(config_file, "test_machine", return_object=True)
        queue = machine.to_queue()
        
        # Create calculator and scheduler
        calc = MockCalc()
        command = queue["launcher"] + " pw.x -in espresso.pwi > espresso.pwo"
        scheduler = SlurmScheduler(calc, queue, command)
        scheduler.write_script()
        
        # Read job file
        job_file_path = os.path.join(calc.directory, "job_file")
        with open(job_file_path, "r") as f:
            job_content = f.read()
        
        # Check for unwanted quotes
        assert '"export' not in job_content, "Job file contains quoted export command"
        assert '"srun' not in job_content, "Job file contains quoted srun command"
        assert 'export OMP_NUM_THREADS=1' in job_content, "Job file should contain unquoted export command"
        assert 'srun --mpi=pmi2' in job_content, "Job file should contain unquoted srun command"
        
        print("✓ Job file generation works correctly without quotes")
        return True
        
    finally:
        shutil.rmtree(config_dir, ignore_errors=True)
        if 'calc' in locals():
            shutil.rmtree(calc.directory, ignore_errors=True)


if __name__ == '__main__':
    print("=" * 60)
    print("Testing Machine Configuration Quotes Sanitization")
    print("=" * 60)
    
    tests = [
        test_sanitize_string_value,
        test_sanitize_list_values,
        test_machine_config_sanitization,
        test_job_file_no_quotes,
    ]
    
    results = []
    for test in tests:
        print(f"\nRunning: {test.__name__}")
        print("-" * 60)
        try:
            results.append(test())
        except Exception as e:
            print(f"✗ Test failed: {e}")
            import traceback
            traceback.print_exc()
            results.append(False)
    
    print("\n" + "=" * 60)
    print(f"Results: {sum(results)}/{len(results)} tests passed")
    print("=" * 60)
    
    sys.exit(0 if all(results) else 1)

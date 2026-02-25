"""
=============================================================================
MOCKING FEATURES FOR TESTING - Best Practices
=============================================================================

When components are not physically installed (like SLURM), you can mock them
for testing purposes. This guide shows different mocking strategies.
"""

# ============================================================================
# 1. MOCKING SYSTEM COMMANDS (like sbatch)
# ============================================================================

from unittest.mock import patch, MagicMock
import subprocess

def example_mock_system_command():
    """Mock a system command that may not be installed."""
    
    # Before: Would fail if sbatch is not installed
    # subprocess.run(['sbatch', 'job_file'], check=True)
    
    # After: Mock the subprocess call
    with patch('subprocess.run') as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        # Now subprocess.run works without needing real sbatch
        result = subprocess.run(['sbatch', 'job_file'], check=True)
        print(f"✅ Command mocked successfully")


# ============================================================================
# 2. MOCKING AVAILABILITY CHECKS (like check_slurm_available)
# ============================================================================

def example_mock_availability_check():
    """Mock availability check to simulate feature being installed."""
    
    from unittest.mock import patch
    
    # Scenario: SLURM not installed physically
    # Problem: Code tries to check if sbatch exists and raises error
    # Solution: Mock the availability check
    
    with patch('xespresso.scheduler.check_slurm_available') as mock_check:
        # Make mock return None (no exception = feature available)
        mock_check.return_value = None
        
        # Now code thinks SLURM is available
        # Create scheduler, generate job_file, etc.
        print(f"✅ Availability check mocked - feature simulated as available")


# ============================================================================
# 3. MOCKING MODULE IMPORTS (when library is not installed)
# ============================================================================

def example_mock_module_import():
    """Mock importing a module that may not be installed."""
    
    from unittest.mock import MagicMock, patch
    import sys
    
    # Scenario: Trying to import a quantum software that's not installed
    # Solution: Mock the import
    
    # Method 1: Mock sys.modules directly
    mock_qe = MagicMock()
    sys.modules['quantum_espresso'] = mock_qe
    
    # Now `import quantum_espresso` will work
    # from quantum_espresso import pw
    # pw will be a mock object
    print(f"✅ Module import mocked")
    
    # Method 2: Use patch for cleaner context management
    with patch.dict('sys.modules', {'quantum_espresso': MagicMock()}):
        # import quantum_espresso  # Would work here
        print(f"✅ Module import mocked in context")


# ============================================================================
# 4. MOCKING FILE OPERATIONS
# ============================================================================

def example_mock_file_operations():
    """Mock file operations for testing without real files."""
    
    from unittest.mock import patch, mock_open
    import os
    
    # Scenario: Code reads a config file from remote server via SSH
    # Solution: Mock the file reading
    
    file_content = """
    nodes=2
    ntasks-per-node=16
    time=04:00:00
    """
    
    with patch('builtins.open', mock_open(read_data=file_content)):
        with open('/remote/config.txt', 'r') as f:
            content = f.read()
            print(f"✅ File operation mocked, read: {len(content)} bytes")


# ============================================================================
# 5. MOCKING SSH/REMOTE CONNECTIONS
# ============================================================================

def example_mock_ssh_connection():
    """Mock SSH operations for testing without actual remote server."""
    
    from unittest.mock import patch, MagicMock
    
    # Scenario: Code connects via SSH to run remote job
    # Solution: Mock the SSH connection
    
    # Mock RemoteAuth or paramiko
    with patch('xespresso.schedulers.remote_connection.RemoteAuth') as mock_auth:
        # Setup mock return values
        mock_instance = MagicMock()
        mock_instance.run_command.return_value = ('stdout', 'stderr')
        mock_instance.retrieve_file.return_value = None
        mock_auth.return_value = mock_instance
        
        # Now code can "connect" and run commands without real SSH
        print(f"✅ SSH connection mocked")


# ============================================================================
# 6. MOCKING CALCULATOR EXECUTION
# ============================================================================

def example_mock_calculator_execution():
    """Mock DFT calculator execution."""
    
    from unittest.mock import patch, MagicMock
    from ase.build import bulk
    
    # Scenario: Full DFT calculation would take hours
    # Solution: Mock the results
    
    with patch('xespresso.Espresso.execute') as mock_execute:
        # Setup mock to return quickly
        mock_execute.return_value = None
        
        # Also mock read_results to return fake energies
        with patch('xespresso.Espresso.read_results') as mock_read:
            mock_read.return_value = None
            
            # Now calculation "runs" instantly
            print(f"✅ Calculator execution mocked")


# ============================================================================
# 7. COMPLETE TEST EXAMPLE: Mocking Multiple Components
# ============================================================================

def example_complete_workflow_mock():
    """Complete example: Mock multiple components for full workflow."""
    
    from unittest.mock import patch, MagicMock, mock_open
    from ase.build import bulk
    from xespresso import Espresso
    
    print("\n" + "=" * 70)
    print("COMPLETE WORKFLOW MOCK EXAMPLE")
    print("=" * 70)
    
    # Scenario: Test remote SLURM execution without:
    # - SLURM installed locally
    # - SSH access to real cluster
    # - Quantum ESPRESSO installed
    
    with patch('xespresso.scheduler.check_slurm_available'):
        with patch('xespresso.schedulers.remote_connection.RemoteAuth') as mock_auth:
            with patch('xespresso.Espresso.read_results') as mock_read:
                # Setup mocks
                mock_connection = MagicMock()
                mock_connection.run_command.return_value = ('Submitted batch job 12345', '')
                mock_auth.return_value = mock_connection
                
                # Now create workflow
                atoms = bulk("Si", cubic=True)
                calc = Espresso(
                    pseudopotentials={"Si": "Si.pbe.UPF"},
                    queue={
                        "execution": "remote",
                        "scheduler": "slurm",
                        "remote_host": "cluster.edu",
                        "remote_user": "user",
                    }
                )
                
                print("\n✅ Mocked workflow components:")
                print("   - check_slurm_available: Bypassed")
                print("   - RemoteAuth: Mocked SSH connection")
                print("   - read_results: Mocked calculator results")
                
                print("\n📝 With these mocks, you can:")
                print("   - Test all code paths without real infrastructure")
                print("   - Run tests in CI/CD pipelines")
                print("   - Test without waiting for long calculations")
                print("   - Test without specific HPC access")


# ============================================================================
# 8. PATCHING STRATEGIES
# ============================================================================

def example_patching_strategies():
    """Different ways to patch/mock code."""
    
    print("\n" + "=" * 70)
    print("PATCHING STRATEGIES")
    print("=" * 70)
    
    print("""
1. DECORATOR APPROACH:
   @patch('module.function')
   def test_something(mock_func):
       mock_func.return_value = 42
       assert function() == 42

2. CONTEXT MANAGER APPROACH:
   with patch('module.function') as mock_func:
       mock_func.return_value = 42
       assert function() == 42

3. DIRECT PATCHING:
   patcher = patch('module.function', return_value=42)
   mock_func = patcher.start()
   try:
       assert function() == 42
   finally:
       patcher.stop()

4. PATCH DICT (for environment/config):
   with patch.dict('os.environ', {'VAR': 'value'}):
       # VAR is now 'value'

5. PATCH MULTIPLE:
   with patch('mod.func1') as m1, patch('mod.func2') as m2:
       m1.return_value = 1
       m2.return_value = 2
       # Both patched
    """)


# ============================================================================
# 9. MOCK VERIFICATION
# ============================================================================

def example_mock_verification():
    """Verify that mocks were called correctly."""
    
    from unittest.mock import patch
    
    print("\n" + "=" * 70)
    print("MOCK VERIFICATION")
    print("=" * 70)
    
    with patch('xespresso.scheduler.check_slurm_available') as mock_check:
        # Simulate code that uses the mock
        try:
            mock_check()  # Simulate calling the mocked function
        except:
            pass
        
        # Verify it was called
        print(f"\n✅ Mock verification:")
        print(f"   - Called: {mock_check.called}")
        print(f"   - Call count: {mock_check.call_count}")
        print(f"   - Called with: {mock_check.call_args}")
        print(f"   - Call args list: {mock_check.call_args_list}")


# ============================================================================
# SUMMARY
# ============================================================================

"""
🎯 KEY TAKEAWAYS:

1. Use mock/patch for any external dependency:
   - System commands (sbatch, srun, etc.)
   - Network operations (SSH, HTTP)
   - File I/O (reading configs from remote servers)
   - Long-running operations (DFT calculations)

2. Common mock targets in xespresso:
   - check_slurm_available() - Mock when SLURM not installed
   - RemoteAuth - Mock SSH connections
   - Espresso.execute() - Mock job execution
   - read_results() - Mock parsing output files
   - subprocess.run() - Mock system commands

3. Testing benefits:
   ✅ Tests run in seconds instead of hours
   ✅ Tests work without target hardware/software
   ✅ Tests work in CI/CD pipelines (Docker containers, etc.)
   ✅ Tests are deterministic and reproducible
   ✅ Tests can simulate error conditions

4. When to use mocks:
   - External system not available (SLURM, remote cluster)
   - External system is expensive (DFT calculations)
   - Testing error handling from external systems
   - Unit testing (tests should be independent)

5. When NOT to use mocks:
   - Integration tests (should use real components)
   - End-to-end tests (should use real infrastructure)
   - Performance testing (mocks don't represent real performance)
"""

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("MOCKING FEATURES FOR TESTING")
    print("=" * 70)
    
    example_mock_availability_check()
    example_mock_file_operations()
    example_complete_workflow_mock()
    example_patching_strategies()
    example_mock_verification()
    
    print("\n" + "=" * 70)
    print("For more info, see: https://docs.python.org/3/library/unittest.mock.html")
    print("=" * 70)

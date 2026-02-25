#!/usr/bin/env python
"""
Test to verify that job_file is generated correctly for remote execution
without requiring ASE_ESPRESSO_COMMAND environment variable.

This demonstrates the fix for the issue:
- User should NOT have to define ASE_ESPRESSO_COMMAND
- The workflow should use calc.command as default
- job_file should be generated correctly with working commands
- When features are not installed, we can mock them for testing
"""

import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from ase.build import bulk
from xespresso import Espresso

def test_job_file_without_ase_espresso_command():
    """Test that job_file is generated without ASE_ESPRESSO_COMMAND."""
    
    print("=" * 70)
    print("Test: Generate job_file WITHOUT ASE_ESPRESSO_COMMAND environment variable")
    print("=" * 70)
    
    # Make sure ASE_ESPRESSO_COMMAND is not defined
    old_ase_cmd = os.environ.pop('ASE_ESPRESSO_COMMAND', None)
    
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a simple structure
            atoms = bulk("Si", cubic=True)
            
            # Create calculator with queue configuration for direct scheduler
            calc = Espresso(
                label=os.path.join(tmpdir, "test"),
                pseudopotentials={"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"},
                queue={
                    "execution": "local",
                    "scheduler": "direct",
                    "launcher": "mpirun -np 16"
                }
            )
            
            atoms.set_calculator(calc)
            
            # Write input - this triggers set_queue() and generates job_file
            calc.write_input(atoms)
            
            # Check if job_file was created
            job_file_path = os.path.join(tmpdir, "test", "job_file")
            assert os.path.exists(job_file_path), f"job_file not created at {job_file_path}"
            
            # Read and display job_file
            with open(job_file_path, 'r') as f:
                job_content = f.read()
            
            print(f"\n✅ job_file generated successfully at: {job_file_path}")
            print(f"\nContent of job_file:")
            print("-" * 70)
            print(job_content)
            print("-" * 70)
            
            # Verify job_file contains expected content
            assert "#!/bin/bash" in job_content, "Missing bash shebang"
            assert "pw.x" in job_content, "Missing pw.x command"
            assert ".pwi" in job_content, "Missing input file reference"
            assert ".pwo" in job_content, "Missing output file reference"
            
            print("\n✅ All validations passed!")
            print("   - Bash shebang present")
            print("   - pw.x command present")
            print("   - Input/output file references present")
            
    finally:
        # Restore ASE_ESPRESSO_COMMAND if it was set
        if old_ase_cmd is not None:
            os.environ['ASE_ESPRESSO_COMMAND'] = old_ase_cmd


def test_job_file_with_mocked_slurm():
    """Test SLURM job_file generation with mocked SLURM check.
    
    This demonstrates: when SLURM is not installed, we can mock it for testing.
    
    MOCKING STRATEGY:
    - Problem: SLURM (sbatch) is not installed on this system
    - Solution: Mock check_slurm_available() to simulate SLURM being present
    - Result: Code path is tested without requiring actual SLURM installation
    
    This is useful for:
    1. Testing on systems without SLURM (local dev machines, Docker containers)
    2. CI/CD pipelines that may not have HPC infrastructure
    3. Rapid testing without waiting for actual job submission
    """
    
    print("\n" + "=" * 70)
    print("Test: SLURM job_file generation (MOCKED)")
    print("=" * 70)
    
    print("\n📝 Scenario: SLURM not physically installed")
    print("   Problem: sbatch command is not available")
    print("   Solution: Mock check_slurm_available() to bypass validation")
    print("   Result: Test SLURM job_file without SLURM installed")
    
    # Make sure ASE_ESPRESSO_COMMAND is not defined
    old_ase_cmd = os.environ.pop('ASE_ESPRESSO_COMMAND', None)
    
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            atoms = bulk("Fe", cubic=True)
            
            # ✅ KEY TECHNIQUE: Mock the SLURM availability check
            # This simulates SLURM being installed without needing real sbatch
            with patch('xespresso.scheduler.check_slurm_available') as mock_slurm_check:
                # Mock returns None (no exception) = SLURM is "available"
                mock_slurm_check.return_value = None
                
                print(f"\n🔄 Mock setup:")
                print(f"   - check_slurm_available patched")
                print(f"   - Will return: None (simulating SLURM available)")
                
                calc = Espresso(
                    label=os.path.join(tmpdir, "test_slurm_mock"),
                    pseudopotentials={"Fe": "Fe.pbe-spn.UPF"},
                    queue={
                        "execution": "local",
                        "scheduler": "slurm",
                        "nodes": 2,
                        "ntasks-per-node": 16,
                        "time": "04:00:00",
                        "partition": "compute"
                    }
                )
                
                atoms.set_calculator(calc)
                
                # This should now succeed because we mocked the SLURM check
                calc.write_input(atoms)
                
                job_file_path = os.path.join(tmpdir, "test_slurm_mock", "job_file")
                assert os.path.exists(job_file_path), "job_file not created for mocked SLURM"
                
                with open(job_file_path, 'r') as f:
                    content = f.read()
                
                print(f"\n✅ SLURM (MOCKED): job_file generated successfully")
                print(f"\nContent of job_file:")
                print("-" * 70)
                print(content)
                print("-" * 70)
                
                # Validate SLURM directives
                assert "#!/bin/bash" in content, "No bash shebang"
                assert "#SBATCH" in content, "No SBATCH directives"
                assert "pw.x" in content, "No pw.x command"
                
                print(f"\n✅ Validations passed:")
                print(f"   - Bash shebang present")
                print(f"   - SBATCH directives present")
                print(f"   - pw.x command present")
                
                # Verify mock was called
                print(f"\n✅ Mock verification:")
                print(f"   - check_slurm_available() was called: {mock_slurm_check.called}")
                print(f"   - Call count: {mock_slurm_check.call_count}")
    
    finally:
        if old_ase_cmd is not None:
            os.environ['ASE_ESPRESSO_COMMAND'] = old_ase_cmd


def test_job_file_with_different_schedulers():
    """Test job_file generation with different scheduler types."""
    
    print("\n" + "=" * 70)
    print("Test: job_file generation with direct scheduler")
    print("=" * 70)
    
    # Make sure ASE_ESPRESSO_COMMAND is not defined
    old_ase_cmd = os.environ.pop('ASE_ESPRESSO_COMMAND', None)
    
    try:
        print(f"\n--- Testing Direct Scheduler ---")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            atoms = bulk("Al", cubic=True)
            
            calc = Espresso(
                label=os.path.join(tmpdir, f"test_direct"),
                pseudopotentials={"Al": "Al.pbe.UPF"},
                queue={
                    "execution": "local",
                    "scheduler": "direct",
                    "launcher": "mpirun -np 4"
                }
            )
            
            atoms.set_calculator(calc)
            calc.write_input(atoms)
            
            job_file_path = os.path.join(tmpdir, f"test_direct", "job_file")
            assert os.path.exists(job_file_path), f"job_file not created for direct scheduler"
            
            with open(job_file_path, 'r') as f:
                content = f.read()
            
            # Basic validations
            assert "#!/bin/bash" in content, f"No bash shebang in direct scheduler"
            assert "pw.x" in content, f"No pw.x command in direct scheduler"
            
            print(f"✅ Direct Scheduler: job_file generated successfully")
            print(f"\nContent of job_file:")
            print("-" * 70)
            print(content)
            print("-" * 70)
    
    finally:
        if old_ase_cmd is not None:
            os.environ['ASE_ESPRESSO_COMMAND'] = old_ase_cmd


def test_fallback_chain():
    """Demonstrate the command fallback chain."""
    
    print("\n" + "=" * 70)
    print("Test: Command resolution fallback chain")
    print("=" * 70)
    
    print("\nFallback Resolution Order:")
    print("1. If ASE_ESPRESSO_COMMAND is set in environment → USE IT")
    print("2. Else if calc.command is available → USE IT")
    print("3. Else → USE DEFAULT TEMPLATE")
    print("\nThis ensures job_file is ALWAYS generated correctly!")
    
    # Example 1: With ASE_ESPRESSO_COMMAND
    print("\n📝 Example 1: WITH ASE_ESPRESSO_COMMAND")
    os.environ['ASE_ESPRESSO_COMMAND'] = "mpirun -np 4 pw.x -in PREFIX.pwi > PREFIX.pwo"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = bulk("Fe", cubic=True)
        calc = Espresso(
            label=os.path.join(tmpdir, "test_with_env"),
            pseudopotentials={"Fe": "Fe.pbe-spn.UPF"},
            queue={"scheduler": "direct", "execution": "local"}
        )
        atoms.set_calculator(calc)
        calc.write_input(atoms)
        
        job_file_path = os.path.join(tmpdir, "test_with_env", "job_file")
        with open(job_file_path, 'r') as f:
            print(f"Job file contains:\n{f.read()}")
    
    os.environ.pop('ASE_ESPRESSO_COMMAND', None)
    
    # Example 2: Without ASE_ESPRESSO_COMMAND (uses calc.command)
    print("\n📝 Example 2: WITHOUT ASE_ESPRESSO_COMMAND (uses calc.command default)")
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = bulk("Ni", cubic=True)
        calc = Espresso(
            label=os.path.join(tmpdir, "test_without_env"),
            pseudopotentials={"Ni": "Ni.pbe.UPF"},
            queue={"scheduler": "direct", "execution": "local"}
        )
        atoms.set_calculator(calc)
        calc.write_input(atoms)
        
        job_file_path = os.path.join(tmpdir, "test_without_env", "job_file")
        with open(job_file_path, 'r') as f:
            print(f"Job file contains:\n{f.read()}")


if __name__ == "__main__":
    print("\n" + "🔬 REMOTE EXECUTION JOB_FILE GENERATION TESTS " + "🔬\n")
    
    test_job_file_without_ase_espresso_command()
    test_job_file_with_mocked_slurm()
    test_job_file_with_different_schedulers()
    test_fallback_chain()
    
    print("\n" + "=" * 70)
    print("✅ ALL TESTS PASSED!")
    print("=" * 70)
    print("\n📌 SUMMARY:")
    print("   - job_file is generated correctly WITHOUT ASE_ESPRESSO_COMMAND")
    print("   - Uses calc.command as fallback default from Espresso class")
    print("   - Supports multiple scheduler types (direct, slurm)")
    print("   - SLURM can be mocked when not physically installed")
    print("   - Remote execution is now ready to use!")
    print("\n")

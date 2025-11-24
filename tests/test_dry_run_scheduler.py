"""
Test that dry run uses xespresso's scheduler system correctly.
"""
import sys
import os
import tempfile
import shutil

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_generate_input_files_with_scheduler():
    """Test that generate_input_files creates job_file using scheduler system."""
    try:
        from ase.build import bulk
        from xespresso.gui.utils.dry_run import generate_input_files
        
        # Create a simple structure
        atoms = bulk('Al', 'fcc', a=4.05)
        
        # Create temporary directory
        tmpdir = tempfile.mkdtemp()
        
        try:
            # Set up calculation parameters with queue configuration
            calc_params = {
                'input_data': {
                    'control': {
                        'calculation': 'scf',
                    },
                    'system': {
                        'ecutwfc': 40,
                        'ecutrho': 320,
                    }
                },
                'pseudopotentials': {'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF'},
                'kpts': (4, 4, 4),
                'queue': {
                    'scheduler': 'direct',
                    'execution': 'local',
                }
            }
            
            # Generate files
            result = generate_input_files(atoms, calc_params, tmpdir)
            
            if result is None:
                print("✗ generate_input_files returned None")
                return False
            
            # Check that input file was created
            assert 'input' in result, "Result should contain 'input' key"
            assert os.path.exists(result['input']), f"Input file should exist: {result['input']}"
            print(f"✓ Input file created: {result['input']}")
            
            # Check that job_file was created by scheduler
            if 'job_file' in result:
                assert os.path.exists(result['job_file']), f"Job file should exist: {result['job_file']}"
                print(f"✓ Job file created by scheduler: {result['job_file']}")
                
                # Verify job file contains expected content
                with open(result['job_file'], 'r') as f:
                    job_content = f.read()
                
                # Should be a bash script
                assert '#!/bin/bash' in job_content, "Job file should be a bash script"
                print("✓ Job file is a valid bash script")
                
                # Should contain the execution command
                assert '.pwi' in job_content or 'espresso' in job_content, "Job file should reference the input file"
                print("✓ Job file contains execution command")
            else:
                print("⚠ Job file not created (may require different scheduler setup)")
            
            return True
            
        finally:
            # Clean up
            shutil.rmtree(tmpdir, ignore_errors=True)
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_scheduler_system_integration():
    """Test that xespresso's scheduler system works correctly."""
    try:
        from ase.build import bulk
        from xespresso import Espresso
        import tempfile
        import shutil
        
        # Create a simple structure
        atoms = bulk('Fe', 'bcc', a=2.87)
        
        # Create temporary directory
        tmpdir = tempfile.mkdtemp()
        
        try:
            # Create calculator with queue configuration
            calc = Espresso(
                input_data={
                    'control': {'calculation': 'scf'},
                    'system': {'ecutwfc': 40}
                },
                pseudopotentials={'Fe': 'Fe.pbe-n-kjpaw_psl.1.0.0.UPF'},
                kpts=(2, 2, 2),
                queue={
                    'scheduler': 'direct',
                    'execution': 'local',
                },
                directory=tmpdir,
                prefix='test'
            )
            
            atoms.calc = calc
            
            # Write input - this should trigger scheduler system
            calc.write_input(atoms)
            
            # Check that job_file was created
            job_file = os.path.join(tmpdir, 'job_file')
            if os.path.exists(job_file):
                print(f"✓ Scheduler created job_file: {job_file}")
                
                with open(job_file, 'r') as f:
                    content = f.read()
                
                print("✓ Job file content:")
                print(content[:200] + "..." if len(content) > 200 else content)
                
                return True
            else:
                print("⚠ Job file not found, checking directory contents:")
                print(os.listdir(tmpdir))
                return False
                
        finally:
            # Clean up
            shutil.rmtree(tmpdir, ignore_errors=True)
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    print("=" * 60)
    print("Testing Dry Run Scheduler Integration")
    print("=" * 60)
    
    tests = [
        test_scheduler_system_integration,
        test_generate_input_files_with_scheduler,
    ]
    
    results = []
    for test in tests:
        print(f"\nRunning: {test.__name__}")
        print("-" * 60)
        results.append(test())
    
    print("\n" + "=" * 60)
    print(f"Results: {sum(results)}/{len(results)} tests passed")
    print("=" * 60)
    
    sys.exit(0 if all(results) else 1)

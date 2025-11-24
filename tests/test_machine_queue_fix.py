"""
Test that Machine objects are properly converted to queue dicts.

This test validates the fix for the issue where Machine objects were passed
to set_queue() causing "argument of type 'Machine' is not iterable" error.
"""
import sys
import os
import tempfile
import shutil

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_machine_object_in_scheduler():
    """Test that set_queue handles Machine objects correctly."""
    try:
        from xespresso.machines.machine import Machine
        from xespresso.scheduler import set_queue
        from unittest.mock import Mock
        
        print("Testing set_queue with Machine object...")
        
        # Create a Machine object
        machine = Machine(
            name='test_machine',
            execution='local',
            scheduler='direct',
            workdir='/tmp/test'
        )
        
        # Create a mock calculator with Machine as queue
        calc = Mock()
        calc.queue = machine
        calc.package = 'pw'
        calc.parallel = ''
        calc.prefix = 'test'
        calc.label = 'test'
        calc.directory = '.'
        
        # This should not raise "argument of type 'Machine' is not iterable"
        set_queue(calc)
        
        # Verify queue was converted to dict
        assert isinstance(calc.queue, dict), "queue should be converted to dict"
        assert calc.queue.get('scheduler') == 'direct', "scheduler should be preserved"
        assert calc.queue.get('execution') == 'local', "execution should be preserved"
        
        print("✓ Machine object successfully converted to queue dict")
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_espresso_calculator_with_machine():
    """Test that Espresso calculator works with Machine object as queue."""
    try:
        from ase.build import bulk
        from xespresso import Espresso
        from xespresso.machines.machine import Machine
        
        print("\nTesting Espresso calculator with Machine object...")
        
        # Create a Machine object
        machine = Machine(
            name='local_test',
            execution='local',
            scheduler='direct',
            workdir='/tmp/test'
        )
        
        # Create a test structure
        atoms = bulk('Si', 'diamond', a=5.43)
        
        # Create temporary directory
        tmpdir = tempfile.mkdtemp()
        
        try:
            # Create calculator with Machine object
            calc = Espresso(
                input_data={
                    'control': {'calculation': 'scf'},
                    'system': {'ecutwfc': 30}
                },
                pseudopotentials={'Si': 'Si.pbe-n-kjpaw_psl.1.0.0.UPF'},
                kpts=(2, 2, 2),
                queue=machine,  # Machine object, not dict
                directory=tmpdir,
                prefix='test_si'
            )
            
            atoms.calc = calc
            
            # Write input - this should trigger set_queue with Machine object
            calc.write_input(atoms)
            
            # Verify queue was converted
            assert isinstance(calc.queue, dict), "calc.queue should be dict after write_input"
            assert calc.queue.get('scheduler') == 'direct'
            
            # Verify input file was created
            pwi_file = os.path.join(tmpdir, 'test_si.pwi')
            assert os.path.exists(pwi_file), f"Input file should exist: {pwi_file}"
            
            print("✓ Espresso calculator works with Machine object")
            return True
            
        finally:
            # Clean up
            shutil.rmtree(tmpdir, ignore_errors=True)
            
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_dict_queue_still_works():
    """Test that dict queue parameter still works (backwards compatibility)."""
    try:
        from xespresso.scheduler import set_queue
        from unittest.mock import Mock
        
        print("\nTesting backwards compatibility with dict queue...")
        
        # Create a mock calculator with dict queue
        calc = Mock()
        calc.queue = {
            'scheduler': 'slurm',
            'execution': 'remote',
            'remote_host': 'cluster.edu'
        }
        calc.package = 'pw'
        calc.parallel = 'mpirun -np 4'
        calc.prefix = 'test'
        calc.label = 'test'
        calc.directory = '.'
        
        # This should still work
        set_queue(calc)
        
        # Verify queue is still a dict
        assert isinstance(calc.queue, dict), "dict queue should remain dict"
        assert calc.queue.get('scheduler') == 'slurm'
        
        print("✓ Dict queue parameter still works")
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    print("=" * 60)
    print("Testing Machine Object to Queue Dict Conversion")
    print("=" * 60)
    
    tests = [
        test_machine_object_in_scheduler,
        test_espresso_calculator_with_machine,
        test_dict_queue_still_works,
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print("\n" + "=" * 60)
    print(f"Results: {sum(results)}/{len(results)} tests passed")
    print("=" * 60)
    
    sys.exit(0 if all(results) else 1)

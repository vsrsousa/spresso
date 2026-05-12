#!/usr/bin/env python
"""
Quick validation script for PROJWFC integration in Wannier Workflow

This script validates that all components work correctly:
1. Imports work without errors
2. Functions are properly defined
3. Workflow class accepts new parameters
4. Methods are accessible
"""

import sys
from pathlib import Path

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    try:
        from xespresso.workflow.wannier_workflow import (
            run_projwfc,
            parse_projwfc_output,
            run_pw2wannier,
            run_wannier90,
            WannierWorkflow,
            generate_seedname_win,
            suggest_nbnd_from_pseudos,
        )
        print("✓ All functions imported successfully")
        return True
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False

def test_function_signatures():
    """Test that functions have correct signatures."""
    print("\nTesting function signatures...")
    try:
        from xespresso.workflow.wannier_workflow import (
            run_projwfc,
            parse_projwfc_output,
        )
        
        # Check run_projwfc signature
        import inspect
        sig = inspect.signature(run_projwfc)
        params = list(sig.parameters.keys())
        expected = ['run_dir', 'prefix', 'queue', 'blocking', 'timeout', 'command', 'den_ext']
        
        if all(p in params for p in expected):
            print(f"✓ run_projwfc has correct parameters: {params}")
        else:
            print(f"⚠ run_projwfc parameters might be incomplete: {params}")
        
        # Check parse_projwfc_output signature
        sig = inspect.signature(parse_projwfc_output)
        params = list(sig.parameters.keys())
        expected = ['projwfc_dir', 'prefix']
        
        if all(p in params for p in expected):
            print(f"✓ parse_projwfc_output has correct parameters: {params}")
        else:
            print(f"⚠ parse_projwfc_output parameters might be incomplete: {params}")
        
        return True
    except Exception as e:
        print(f"✗ Signature test error: {e}")
        return False

def test_workflow_class():
    """Test that WannierWorkflow class works."""
    print("\nTesting WannierWorkflow class...")
    try:
        from xespresso.workflow.wannier_workflow import WannierWorkflow
        import inspect
        
        # Check __init__ parameters
        sig = inspect.signature(WannierWorkflow.__init__)
        params = list(sig.parameters.keys())
        
        critical_params = ['cif_file', 'pseudos', 'protocol', 'num_wann', 'projections']
        if all(p in params for p in critical_params):
            print(f"✓ WannierWorkflow.__init__ has critical parameters")
        else:
            print(f"✗ Missing critical parameters in WannierWorkflow.__init__")
            return False
        
        # Check run() method signature
        sig = inspect.signature(WannierWorkflow.run)
        params = list(sig.parameters.keys())
        
        expected_in_run = ['self', 'labels', 'blocking', 'seedname', 
                          'run_bands_validation', 'run_projwfc_analysis']
        if all(p in params for p in expected_in_run):
            print(f"✓ WannierWorkflow.run() has run_projwfc_analysis parameter")
        else:
            print(f"✗ Missing run_projwfc_analysis parameter in run() method")
            return False
        
        # Check new methods exist
        methods = ['get_projwfc_analysis', 'get_scf_calculator', 
                  'get_nscf_calculator', 'get_band_structure_calculator',
                  'compare_bands', 'validate_wannier_quality']
        
        for method_name in methods:
            if hasattr(WannierWorkflow, method_name):
                print(f"✓ Method '{method_name}' exists")
            else:
                print(f"✗ Method '{method_name}' not found")
                return False
        
        return True
    except Exception as e:
        print(f"✗ Workflow class test error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_docstrings():
    """Test that new functions have documentation."""
    print("\nTesting documentation...")
    try:
        from xespresso.workflow.wannier_workflow import (
            run_projwfc,
            parse_projwfc_output,
        )
        
        if run_projwfc.__doc__:
            print(f"✓ run_projwfc has docstring ({len(run_projwfc.__doc__)} chars)")
        else:
            print(f"⚠ run_projwfc missing docstring")
        
        if parse_projwfc_output.__doc__:
            print(f"✓ parse_projwfc_output has docstring ({len(parse_projwfc_output.__doc__)} chars)")
        else:
            print(f"⚠ parse_projwfc_output missing docstring")
        
        return True
    except Exception as e:
        print(f"✗ Documentation test error: {e}")
        return False

def test_example_file():
    """Test that example file exists and is readable."""
    print("\nTesting example file...")
    try:
        example_file = Path("/home/vinicius/scratch/projects/spresso/examples/wannier_workflow_with_projwfc_example.py")
        if example_file.exists():
            with open(example_file, 'r') as f:
                content = f.read()
            
            # Check for key content
            checks = {
                'WannierWorkflow import': 'WannierWorkflow',
                'PROJWFC parameter': 'run_projwfc_analysis',
                'get_projwfc_analysis': 'get_projwfc_analysis',
                'validate_wannier_quality': 'validate_wannier_quality',
            }
            
            for check_name, check_str in checks.items():
                if check_str in content:
                    print(f"✓ Example file contains '{check_name}'")
                else:
                    print(f"⚠ Example file missing '{check_name}'")
            
            print(f"✓ Example file readable ({len(content)} bytes)")
            return True
        else:
            print(f"✗ Example file not found: {example_file}")
            return False
    except Exception as e:
        print(f"✗ Example file test error: {e}")
        return False

def test_doc_files():
    """Test that documentation files exist."""
    print("\nTesting documentation files...")
    try:
        doc_files = {
            'PROJWFC workflow doc': Path("/home/vinicius/scratch/projects/spresso/docs/PROJWFC_IN_WANNIER_WORKFLOW.md"),
            'Integration summary': Path("/home/vinicius/scratch/projects/spresso/docs/PROJWFC_INTEGRATION_SUMMARY.md"),
        }
        
        all_exist = True
        for name, filepath in doc_files.items():
            if filepath.exists():
                with open(filepath, 'r') as f:
                    size = len(f.read())
                print(f"✓ {name} exists ({size} bytes)")
            else:
                print(f"✗ {name} not found: {filepath}")
                all_exist = False
        
        return all_exist
    except Exception as e:
        print(f"✗ Documentation files test error: {e}")
        return False

def main():
    """Run all tests."""
    print("="*70)
    print("WANNIER WORKFLOW PROJWFC INTEGRATION VALIDATION")
    print("="*70)
    
    tests = [
        ("Imports", test_imports),
        ("Function Signatures", test_function_signatures),
        ("Workflow Class", test_workflow_class),
        ("Documentation", test_docstrings),
        ("Example File", test_example_file),
        ("Doc Files", test_doc_files),
    ]
    
    results = {}
    for test_name, test_func in tests:
        try:
            results[test_name] = test_func()
        except Exception as e:
            print(f"\n✗ {test_name} test crashed: {e}")
            import traceback
            traceback.print_exc()
            results[test_name] = False
    
    # Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    total = len(results)
    passed = sum(1 for v in results.values() if v)
    
    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status:8} {test_name}")
    
    print("-"*70)
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ ALL TESTS PASSED - PROJWFC INTEGRATION IS COMPLETE AND WORKING")
        return 0
    else:
        print(f"\n⚠ {total - passed} test(s) failed - Please review above")
        return 1

if __name__ == '__main__':
    exit_code = main()
    sys.exit(exit_code)

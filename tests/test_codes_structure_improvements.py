"""
Tests for codes configuration structure improvements.

These tests validate the fixes for the issues reported in the problem statement:
1. No redundant top-level modules when versions structure exists
2. Auto-detection of common qe_prefix for codes in same directory
3. Proper storage of version-specific modules and qe_prefix

All tests ensure backward compatibility with existing single-version configs.
"""

import pytest
import os
import tempfile
import json
from pathlib import Path

from xespresso.codes.config import Code, CodesConfig
from xespresso.codes.manager import CodesManager


def test_no_redundant_top_level_modules_with_versions():
    """
    Test that when versions structure exists, modules are NOT duplicated at top level.
    This addresses the main issue in the problem statement.
    """
    config = CodesManager.create_config(
        machine_name="snake5",
        detected_codes={
            'pw': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/pw.x',
            'ph': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/ph.x',
            'pp': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/pp.x',
        },
        qe_version='7.5',
        modules=['qe/7.5-intel_mpi_mkl_scalapack']
    )
    
    config_dict = config.to_dict()
    
    # Should NOT have top-level modules (redundancy removed)
    assert 'modules' not in config_dict or config_dict['modules'] is None, \
        "Top-level modules should not be present when versions structure exists"
    
    # Should have qe_version at top level (indicates default version)
    assert config_dict.get('qe_version') == '7.5', \
        "Top-level qe_version should indicate default version"
    
    # Modules should be in version-specific structure
    assert '7.5' in config_dict['versions'], \
        "Version 7.5 should be in versions structure"
    assert config_dict['versions']['7.5']['modules'] == ['qe/7.5-intel_mpi_mkl_scalapack'], \
        "Modules should be in version-specific structure"


def test_auto_detect_common_qe_prefix():
    """
    Test that common directory prefix is automatically detected when all codes
    are in the same directory.
    """
    config = CodesManager.create_config(
        machine_name="snake5",
        detected_codes={
            'pw': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/pw.x',
            'ph': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/ph.x',
            'pp': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/pp.x',
            'projwfc': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/projwfc.x',
        },
        qe_version='7.5',
        modules=['qe/7.5-intel_mpi_mkl_scalapack']
    )
    
    config_dict = config.to_dict()
    
    # Should have qe_prefix at top level
    assert config_dict.get('qe_prefix') == '/opt/qe/intel/oneapi-2021.4.0/7.5/bin', \
        "Common directory prefix should be auto-detected"
    
    # Should also have qe_prefix in version structure
    assert config_dict['versions']['7.5']['qe_prefix'] == '/opt/qe/intel/oneapi-2021.4.0/7.5/bin', \
        "qe_prefix should also be in version structure"


def test_qe_prefix_stored_in_version_blocks():
    """
    Test that qe_prefix is stored within version blocks for each version.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create config for version 7.5
        config1 = CodesManager.create_config(
            machine_name="cluster",
            detected_codes={
                'pw': '/opt/qe-7.5/bin/pw.x',
                'ph': '/opt/qe-7.5/bin/ph.x',
            },
            qe_version='7.5',
            modules=['qe/7.5']
        )
        filepath = CodesManager.save_config(config1, output_dir=tmpdir, interactive=False)
        
        # Add version 7.4
        config2 = CodesManager.create_config(
            machine_name="cluster",
            detected_codes={
                'pw': '/opt/qe-7.4/bin/pw.x',
                'ph': '/opt/qe-7.4/bin/ph.x',
            },
            qe_version='7.4',
            modules=['qe/7.4']
        )
        CodesManager.save_config(config2, output_dir=tmpdir, merge=True, interactive=False)
        
        # Load and verify
        with open(filepath) as f:
            saved_data = json.load(f)
        
        # Each version should have its own qe_prefix
        assert saved_data['versions']['7.5']['qe_prefix'] == '/opt/qe-7.5/bin', \
            "Version 7.5 should have its own qe_prefix"
        assert saved_data['versions']['7.4']['qe_prefix'] == '/opt/qe-7.4/bin', \
            "Version 7.4 should have its own qe_prefix"
        
        # Each version should have its own modules
        assert saved_data['versions']['7.5']['modules'] == ['qe/7.5'], \
            "Version 7.5 should have its own modules"
        assert saved_data['versions']['7.4']['modules'] == ['qe/7.4'], \
            "Version 7.4 should have its own modules"


def test_detect_common_prefix_helper():
    """Test the detect_common_prefix helper function."""
    # All codes in same directory
    codes1 = {
        'pw': '/opt/qe/7.5/bin/pw.x',
        'ph': '/opt/qe/7.5/bin/ph.x',
        'pp': '/opt/qe/7.5/bin/pp.x',
    }
    prefix1 = CodesManager.detect_common_prefix(codes1)
    assert prefix1 == '/opt/qe/7.5/bin', \
        "Should detect common directory for codes in same path"
    
    # Codes in different directories
    codes2 = {
        'pw': '/opt/qe/7.5/bin/pw.x',
        'ph': '/usr/local/bin/ph.x',
    }
    prefix2 = CodesManager.detect_common_prefix(codes2)
    assert prefix2 is None, \
        "Should return None when codes are in different directories"
    
    # Single code
    codes3 = {
        'pw': '/opt/qe/7.5/bin/pw.x',
    }
    prefix3 = CodesManager.detect_common_prefix(codes3)
    assert prefix3 == '/opt/qe/7.5/bin', \
        "Should return directory for single code"


def test_backward_compatibility_without_versions():
    """
    Test that configs without versions structure still work correctly.
    This ensures backward compatibility.
    """
    config = CodesManager.create_config(
        machine_name="local",
        detected_codes={
            'pw': '/usr/bin/pw.x',
            'ph': '/usr/bin/ph.x',
        },
        qe_version=None  # No version specified
    )
    
    config_dict = config.to_dict()
    
    # Should have codes in main codes dict
    assert 'pw' in config_dict['codes'], \
        "Codes should be in main dict when no version specified"
    assert 'ph' in config_dict['codes'], \
        "All codes should be in main dict"
    
    # Should have qe_prefix at top level
    assert config_dict.get('qe_prefix') == '/usr/bin', \
        "qe_prefix should be at top level for non-versioned config"
    
    # Should not have versions structure
    assert config_dict.get('versions') is None or len(config_dict['versions']) == 0, \
        "Should not have versions structure for non-versioned config"


def test_real_world_scenario_snake5():
    """
    Test the exact scenario from the problem statement with snake5 machine.
    This validates the complete fix.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create config exactly as in problem statement
        config = CodesManager.create_config(
            machine_name="snake5",
            detected_codes={
                'pw': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/pw.x',
                'ph': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/ph.x',
                'pp': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/pp.x',
                'projwfc': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/projwfc.x',
                'dos': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/dos.x',
                'bands': '/opt/qe/intel/oneapi-2021.4.0/7.5/bin/bands.x',
            },
            qe_version='7.5',
            modules=['qe/7.5-intel_mpi_mkl_scalapack']
        )
        
        filepath = CodesManager.save_config(config, output_dir=tmpdir, interactive=False)
        
        # Load and verify JSON structure
        with open(filepath) as f:
            saved_data = json.load(f)
        
        # Verify the structure matches expected improvements
        assert saved_data['machine_name'] == 'snake5'
        assert saved_data['qe_version'] == '7.5'  # OK: indicates default
        assert 'modules' not in saved_data or saved_data.get('modules') is None  # NOT redundant
        assert saved_data['codes'] == {}  # Empty: codes are in versions
        
        # Version-specific data
        assert '7.5' in saved_data['versions']
        version_75 = saved_data['versions']['7.5']
        assert 'codes' in version_75
        assert 'modules' in version_75
        assert version_75['modules'] == ['qe/7.5-intel_mpi_mkl_scalapack']
        assert 'qe_prefix' in version_75
        assert version_75['qe_prefix'] == '/opt/qe/intel/oneapi-2021.4.0/7.5/bin'
        
        print("\n" + "="*60)
        print("✅ Structure improvements validated:")
        print("  - No redundant top-level modules")
        print("  - qe_prefix auto-detected and stored in version")
        print("  - Version-specific modules properly organized")
        print("="*60)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

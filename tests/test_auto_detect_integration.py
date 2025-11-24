"""
Integration test demonstrating the Auto Detect Codes fix end-to-end.

This test simulates the exact scenario that was broken:
1. User auto-detects codes with a QE version specified
2. GUI should properly display the detected codes
3. User should be able to save the configuration
"""

import tempfile
import pytest
from xespresso.codes.manager import detect_qe_codes, CodesManager, DEFAULT_CODES_DIR


def test_auto_detect_with_version_integration():
    """
    Integration test: Auto-detect codes with version specified.
    
    This simulates the exact user workflow that was broken:
    1. Detect codes with qe_version specified
    2. Check that codes are detected (has_any_codes)
    3. Display codes (get_all_codes)
    4. Save codes successfully
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Step 1: User detects codes with version specified
        # (This would normally find real codes, but we'll use mock data)
        from unittest.mock import patch
        
        mock_detected = {
            'pw': '/usr/bin/pw.x',
            'ph': '/usr/bin/ph.x',
            'pp': '/usr/bin/pp.x'
        }
        
        with patch('xespresso.codes.manager.CodesManager.detect_codes', return_value=mock_detected):
            codes_config = detect_qe_codes(
                machine_name='test_cluster',
                qe_version='7.2',
                label='production'
            )
        
        # Step 2: GUI checks if codes were detected
        # Before fix: codes_config.codes would be empty, this would fail
        # After fix: has_any_codes() returns True
        assert codes_config.has_any_codes(), "Should detect codes in versions structure"
        
        # Verify main codes is empty (codes are in versions)
        assert len(codes_config.codes) == 0, "Main codes should be empty when version specified"
        
        # Step 3: GUI gets all codes to display
        # Before fix: codes_config.codes.items() would iterate over nothing
        # After fix: get_all_codes() returns the codes from versions structure
        all_codes = codes_config.get_all_codes()
        assert len(all_codes) == 3, "Should get all 3 detected codes"
        assert 'pw' in all_codes
        assert 'ph' in all_codes
        assert 'pp' in all_codes
        
        # Verify the codes have correct attributes
        assert all_codes['pw'].path == '/usr/bin/pw.x'
        assert all_codes['pw'].version == '7.2'
        
        # Step 4: User saves the configuration
        filepath = CodesManager.save_config(
            codes_config,
            output_dir=tmpdir,
            interactive=False
        )
        
        # Verify save was successful
        assert filepath.endswith('test_cluster.json')
        
        # Step 5: Verify we can load it back
        loaded_config = CodesManager.load_config('test_cluster', tmpdir)
        assert loaded_config is not None
        assert loaded_config.has_any_codes()
        
        # Verify loaded codes match
        loaded_codes = loaded_config.get_all_codes()
        assert len(loaded_codes) == 3
        assert loaded_codes['pw'].path == '/usr/bin/pw.x'


def test_auto_detect_without_version_still_works():
    """
    Integration test: Verify backward compatibility.
    
    Auto-detect without version should still work as before,
    storing codes in main codes dict.
    """
    from unittest.mock import patch
    
    mock_detected = {
        'pw': '/usr/bin/pw.x',
        'ph': '/usr/bin/ph.x',
    }
    
    with patch('xespresso.codes.manager.CodesManager.detect_codes', return_value=mock_detected):
        codes_config = detect_qe_codes(
            machine_name='test_cluster'
            # No version or label specified
        )
    
    # Should detect codes
    assert codes_config.has_any_codes()
    
    # Codes should be in main dict (backward compatible behavior)
    assert len(codes_config.codes) == 2
    
    # get_all_codes should still work
    all_codes = codes_config.get_all_codes()
    assert len(all_codes) == 2
    assert 'pw' in all_codes
    assert 'ph' in all_codes


def test_multiple_versions_in_same_config():
    """
    Integration test: Multiple versions in same machine config.
    
    Users should be able to detect and save multiple QE versions
    for the same machine, with the GUI properly displaying each.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        from unittest.mock import patch
        
        # Detect version 7.2
        mock_detected_72 = {'pw': '/opt/qe-7.2/bin/pw.x', 'ph': '/opt/qe-7.2/bin/ph.x'}
        with patch('xespresso.codes.manager.CodesManager.detect_codes', return_value=mock_detected_72):
            config_72 = detect_qe_codes(
                machine_name='cluster',
                qe_version='7.2',
                label='stable'
            )
        
        # Save first version
        CodesManager.save_config(config_72, output_dir=tmpdir, interactive=False)
        
        # Detect version 7.3
        mock_detected_73 = {'pw': '/opt/qe-7.3/bin/pw.x', 'ph': '/opt/qe-7.3/bin/ph.x', 'hp': '/opt/qe-7.3/bin/hp.x'}
        with patch('xespresso.codes.manager.CodesManager.detect_codes', return_value=mock_detected_73):
            config_73 = detect_qe_codes(
                machine_name='cluster',
                qe_version='7.3',
                label='latest'
            )
        
        # Save second version (merge mode)
        CodesManager.save_config(config_73, output_dir=tmpdir, merge=True, interactive=False)
        
        # Load the merged config
        loaded = CodesManager.load_config('cluster', tmpdir)
        
        # Should have both versions
        assert '7.2' in loaded.list_versions()
        assert '7.3' in loaded.list_versions()
        
        # Should be able to get codes for each version separately
        codes_72 = loaded.get_all_codes(version='7.2')
        assert len(codes_72) == 2
        assert codes_72['pw'].path == '/opt/qe-7.2/bin/pw.x'
        
        codes_73 = loaded.get_all_codes(version='7.3')
        assert len(codes_73) == 3
        assert codes_73['pw'].path == '/opt/qe-7.3/bin/pw.x'
        assert 'hp' in codes_73


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

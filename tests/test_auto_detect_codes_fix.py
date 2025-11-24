"""
Tests for Auto Detect Codes fix.

These tests verify that the auto-detect codes functionality works correctly
when QE version or label is specified, and codes are stored in the versions structure.
"""

import pytest
from xespresso.codes.config import Code, CodesConfig
from xespresso.codes.manager import CodesManager


class TestAutoDetectCodesFix:
    """Tests for the auto-detect codes functionality fix."""
    
    def test_has_any_codes_with_main_codes(self):
        """Test has_any_codes returns True when codes are in main dict."""
        config = CodesConfig(
            machine_name='test',
            codes={'pw': Code(name='pw', path='/usr/bin/pw.x')}
        )
        assert config.has_any_codes() is True
    
    def test_has_any_codes_with_version_codes(self):
        """Test has_any_codes returns True when codes are in versions structure."""
        config = CodesConfig(machine_name='test')
        config.add_code(Code(name='pw', path='/usr/bin/pw.x', version='7.2'), version='7.2')
        assert config.has_any_codes() is True
    
    def test_has_any_codes_empty(self):
        """Test has_any_codes returns False when no codes configured."""
        config = CodesConfig(machine_name='test')
        assert config.has_any_codes() is False
    
    def test_get_all_codes_from_main(self):
        """Test get_all_codes returns codes from main dict."""
        config = CodesConfig(
            machine_name='test',
            codes={
                'pw': Code(name='pw', path='/usr/bin/pw.x'),
                'ph': Code(name='ph', path='/usr/bin/ph.x')
            }
        )
        all_codes = config.get_all_codes()
        assert len(all_codes) == 2
        assert 'pw' in all_codes
        assert 'ph' in all_codes
        assert all_codes['pw'].path == '/usr/bin/pw.x'
    
    def test_get_all_codes_from_version(self):
        """Test get_all_codes returns codes from version structure."""
        config = CodesConfig(machine_name='test', qe_version='7.2')
        config.add_code(Code(name='pw', path='/usr/bin/pw.x', version='7.2'), version='7.2')
        config.add_code(Code(name='ph', path='/usr/bin/ph.x', version='7.2'), version='7.2')
        
        all_codes = config.get_all_codes()
        assert len(all_codes) == 2
        assert 'pw' in all_codes
        assert 'ph' in all_codes
        assert all_codes['pw'].path == '/usr/bin/pw.x'
    
    def test_get_all_codes_with_specific_version(self):
        """Test get_all_codes can get codes from a specific version."""
        config = CodesConfig(machine_name='test')
        config.add_code(Code(name='pw', path='/usr/bin/pw-7.2.x', version='7.2'), version='7.2')
        config.add_code(Code(name='pw', path='/usr/bin/pw-7.3.x', version='7.3'), version='7.3')
        
        codes_7_2 = config.get_all_codes(version='7.2')
        assert len(codes_7_2) == 1
        assert codes_7_2['pw'].path == '/usr/bin/pw-7.2.x'
        
        codes_7_3 = config.get_all_codes(version='7.3')
        assert len(codes_7_3) == 1
        assert codes_7_3['pw'].path == '/usr/bin/pw-7.3.x'
    
    def test_create_config_with_version_stores_in_versions(self):
        """Test that create_config stores codes in versions structure when version is specified."""
        detected_codes = {
            'pw': '/usr/bin/pw.x',
            'ph': '/usr/bin/ph.x',
        }
        
        config = CodesManager.create_config(
            machine_name='test',
            detected_codes=detected_codes,
            qe_version='7.2',
            label='production'
        )
        
        # Codes should be in versions structure, not main codes
        assert len(config.codes) == 0
        assert config.has_any_codes() is True
        
        # Should be able to get codes via get_all_codes
        all_codes = config.get_all_codes()
        assert len(all_codes) == 2
        assert 'pw' in all_codes
        assert 'ph' in all_codes
    
    def test_create_config_with_label_stores_in_versions(self):
        """Test that create_config stores codes in versions structure when label is specified."""
        detected_codes = {
            'pw': '/usr/bin/pw.x',
        }
        
        config = CodesManager.create_config(
            machine_name='test',
            detected_codes=detected_codes,
            qe_version='7.2',
            label='dev'
        )
        
        # Codes should be in versions structure
        assert len(config.codes) == 0
        assert config.has_any_codes() is True
        
        all_codes = config.get_all_codes()
        assert len(all_codes) == 1
        assert 'pw' in all_codes
    
    def test_create_config_without_version_stores_in_main(self):
        """Test that create_config stores codes in main dict when no version specified."""
        detected_codes = {
            'pw': '/usr/bin/pw.x',
            'ph': '/usr/bin/ph.x',
        }
        
        config = CodesManager.create_config(
            machine_name='test',
            detected_codes=detected_codes
        )
        
        # Codes should be in main codes dict
        assert len(config.codes) == 2
        assert config.has_any_codes() is True
        
        all_codes = config.get_all_codes()
        assert len(all_codes) == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

"""
Tests for codes configuration version merging behavior.

These tests verify that when saving codes with different versions/labels,
they are properly stored in the versions structure within a single JSON file
per machine, rather than overwriting each other.
"""

import pytest
import os
import tempfile
from pathlib import Path


def test_save_config_with_version_creates_versions_structure():
    """Test that save_config creates versions structure when version is specified."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    codes = {
        'pw': Code(name='pw', path='/usr/bin/pw.x', version='7.2'),
    }
    
    config = CodesConfig(
        machine_name='test_machine',
        codes=codes,
        qe_version='7.2'
    )
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create config using create_config which should put codes in versions
        from xespresso.codes.manager import CodesManager
        config = CodesManager.create_config(
            machine_name='test_machine',
            detected_codes={'pw': '/usr/bin/pw.x'},
            qe_version='7.2'
        )
        
        filepath = CodesManager.save_config(
            config,
            output_dir=tmpdir,
            interactive=False
        )
        
        # Should save as test_machine.json
        assert os.path.basename(filepath) == 'test_machine.json'
        
        # Load and verify structure
        loaded = CodesConfig.from_json(filepath)
        assert loaded.versions is not None
        assert '7.2' in loaded.versions
        assert 'codes' in loaded.versions['7.2']
        assert 'pw' in loaded.versions['7.2']['codes']


def test_save_multiple_versions_to_same_file():
    """Test that multiple versions are saved to the same machine file."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Save version 7.2
        config1 = CodesManager.create_config(
            machine_name='cluster',
            detected_codes={'pw': '/opt/qe-7.2/pw.x'},
            qe_version='7.2',
            label='production'
        )
        filepath1 = CodesManager.save_config(config1, output_dir=tmpdir, interactive=False)
        
        # Save version 7.3 with merge=True (should add to same file)
        config2 = CodesManager.create_config(
            machine_name='cluster',
            detected_codes={'pw': '/opt/qe-7.3/pw.x'},
            qe_version='7.3',
            label='dev'
        )
        filepath2 = CodesManager.save_config(config2, output_dir=tmpdir, merge=True, interactive=False)
        
        # Both should point to same file
        assert filepath1 == filepath2
        assert os.path.basename(filepath1) == 'cluster.json'
        
        # Load and verify both versions exist
        loaded = CodesConfig.from_json(filepath1)
        assert loaded.versions is not None
        assert '7.2' in loaded.versions
        assert '7.3' in loaded.versions
        assert loaded.versions['7.2']['label'] == 'production'
        assert loaded.versions['7.3']['label'] == 'dev'


def test_merge_preserves_existing_versions():
    """Test that merging a new version doesn't overwrite existing versions."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create initial config with version 7.2
        config1 = CodesManager.create_config(
            machine_name='test',
            detected_codes={'pw': '/opt/qe-7.2/pw.x', 'ph': '/opt/qe-7.2/ph.x'},
            qe_version='7.2',
            label='stable'
        )
        filepath = CodesManager.save_config(config1, output_dir=tmpdir, interactive=False)
        
        # Add version 7.3
        config2 = CodesManager.create_config(
            machine_name='test',
            detected_codes={'pw': '/opt/qe-7.3/pw.x'},
            qe_version='7.3',
            label='latest'
        )
        CodesManager.save_config(config2, output_dir=tmpdir, merge=True, interactive=False)
        
        # Load and verify both versions are preserved
        loaded = CodesConfig.from_json(filepath)
        assert '7.2' in loaded.versions
        assert '7.3' in loaded.versions
        
        # Version 7.2 should still have both codes
        assert 'pw' in loaded.versions['7.2']['codes']
        assert 'ph' in loaded.versions['7.2']['codes']
        assert loaded.versions['7.2']['label'] == 'stable'
        
        # Version 7.3 should have its code
        assert 'pw' in loaded.versions['7.3']['codes']
        assert loaded.versions['7.3']['label'] == 'latest'


def test_save_without_version_uses_main_codes():
    """Test that configs without version go to main codes dict (backward compatible)."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        config = CodesManager.create_config(
            machine_name='test',
            detected_codes={'pw': '/usr/bin/pw.x'},
            qe_version=None  # No version specified
        )
        
        filepath = CodesManager.save_config(config, output_dir=tmpdir, interactive=False)
        
        # Load and verify it's in main codes, not versions
        loaded = CodesConfig.from_json(filepath)
        assert 'pw' in loaded.codes
        assert loaded.versions is None or len(loaded.versions) == 0


def test_merge_same_version_updates_codes():
    """Test that merging same version updates its codes."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Initial save with just pw
        config1 = CodesManager.create_config(
            machine_name='test',
            detected_codes={'pw': '/opt/qe-7.2/pw.x'},
            qe_version='7.2'
        )
        filepath = CodesManager.save_config(config1, output_dir=tmpdir, interactive=False)
        
        # Merge with ph in same version
        config2 = CodesManager.create_config(
            machine_name='test',
            detected_codes={'ph': '/opt/qe-7.2/ph.x'},
            qe_version='7.2'
        )
        CodesManager.save_config(config2, output_dir=tmpdir, merge=True, interactive=False)
        
        # Load and verify both codes in version 7.2
        loaded = CodesConfig.from_json(filepath)
        assert '7.2' in loaded.versions
        assert 'pw' in loaded.versions['7.2']['codes']
        assert 'ph' in loaded.versions['7.2']['codes']


def test_label_stored_in_version_structure():
    """Test that labels are stored within version structures."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        config = CodesManager.create_config(
            machine_name='test',
            detected_codes={'pw': '/opt/qe-7.2/pw.x'},
            qe_version='7.2',
            label='production',
            modules=['qe/7.2']
        )
        
        filepath = CodesManager.save_config(config, output_dir=tmpdir, interactive=False)
        
        # Load and verify label is in version structure
        loaded = CodesConfig.from_json(filepath)
        assert '7.2' in loaded.versions
        assert 'label' in loaded.versions['7.2']
        assert loaded.versions['7.2']['label'] == 'production'
        assert 'modules' in loaded.versions['7.2']
        assert loaded.versions['7.2']['modules'] == ['qe/7.2']


def test_save_without_merge_raises_error_if_file_exists():
    """Test that saving without merge raises error when file exists."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # First save
        config1 = CodesManager.create_config(
            machine_name='test',
            detected_codes={'pw': '/usr/bin/pw.x'},
            qe_version='7.2'
        )
        CodesManager.save_config(config1, output_dir=tmpdir, interactive=False)
        
        # Second save without merge should raise error
        config2 = CodesManager.create_config(
            machine_name='test',
            detected_codes={'ph': '/usr/bin/ph.x'},
            qe_version='7.3'
        )
        
        with pytest.raises(FileExistsError):
            CodesManager.save_config(config2, output_dir=tmpdir, merge=False, interactive=False)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

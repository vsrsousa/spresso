"""
Tests for GUI codes configuration save functionality.

These tests verify that the GUI can properly save codes configuration
using xespresso's CodesManager with interactive=False to prevent hanging.
"""

import pytest
import os
import tempfile
from pathlib import Path


def test_save_config_interactive_false():
    """Test that save_config works with interactive=False (for GUI usage)."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    # Create a test codes configuration
    codes = {
        'pw': Code(name='pw', path='/usr/bin/pw.x', version='7.2'),
        'ph': Code(name='ph', path='/usr/bin/ph.x', version='7.2'),
    }
    
    config = CodesConfig(
        machine_name='test_machine',
        codes=codes,
        qe_version='7.2',
        label='test_label'
    )
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test save with interactive=False (should not hang or prompt)
        filepath = CodesManager.save_config(
            config,
            output_dir=tmpdir,
            overwrite=False,
            merge=False,
            interactive=False  # Critical for GUI usage
        )
        
        assert os.path.exists(filepath)
        
        # Verify the file contains expected data
        loaded_config = CodesConfig.from_json(filepath)
        assert loaded_config.machine_name == 'test_machine'
        assert loaded_config.qe_version == '7.2'
        assert loaded_config.label == 'test_label'
        assert len(loaded_config.codes) == 2


def test_save_config_merge_with_interactive_false():
    """Test that save_config merge mode works with interactive=False."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # First save with version 7.2
        config1 = CodesManager.create_config(
            machine_name='test_machine',
            detected_codes={'pw': '/usr/bin/pw.x'},
            qe_version='7.2',
            label='production'
        )
        
        filepath = CodesManager.save_config(
            config1,
            output_dir=tmpdir,
            interactive=False
        )
        
        # Second save with version 7.3 and merge=True (should add new version)
        config2 = CodesManager.create_config(
            machine_name='test_machine',
            detected_codes={'ph': '/usr/bin/ph.x'},
            qe_version='7.3',
            label='dev'
        )
        
        filepath = CodesManager.save_config(
            config2,
            output_dir=tmpdir,
            merge=True,
            interactive=False
        )
        
        # Verify merge worked - should have both versions now
        loaded_config = CodesConfig.from_json(filepath)
        assert loaded_config.versions is not None
        assert '7.2' in loaded_config.versions
        assert '7.3' in loaded_config.versions
        assert loaded_config.versions['7.2']['label'] == 'production'
        assert loaded_config.versions['7.3']['label'] == 'dev'


def test_save_config_same_machine_different_versions():
    """Test that different versions for same machine are stored in same file."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    import os
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # First save with version 7.2
        config1 = CodesManager.create_config(
            machine_name='test_machine',
            detected_codes={'pw': '/usr/bin/pw-7.2.x'},
            qe_version='7.2',
            label='stable'
        )
        
        filepath1 = CodesManager.save_config(
            config1,
            output_dir=tmpdir,
            interactive=False
        )
        
        # Second save with version 7.3 - should go to same file
        config2 = CodesManager.create_config(
            machine_name='test_machine',
            detected_codes={'pw': '/usr/bin/pw-7.3.x'},
            qe_version='7.3',
            label='latest'
        )
        
        filepath2 = CodesManager.save_config(
            config2,
            output_dir=tmpdir,
            merge=True,
            interactive=False
        )
        
        # Both should point to the same file
        assert filepath1 == filepath2
        assert os.path.basename(filepath1) == 'test_machine.json'
        
        # Verify both versions exist in the file
        loaded_config = CodesConfig.from_json(filepath1)
        assert loaded_config.versions is not None
        assert '7.2' in loaded_config.versions
        assert '7.3' in loaded_config.versions

def test_save_config_file_exists_error_with_interactive_false():
    """Test that save_config raises error when file exists without merge/overwrite."""
    from xespresso.codes.config import Code, CodesConfig
    from xespresso.codes.manager import CodesManager
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # First save
        config1 = CodesConfig(
            machine_name='test_machine',
            codes={'pw': Code(name='pw', path='/usr/bin/pw.x', version='7.2')},
        )
        
        filepath = CodesManager.save_config(
            config1,
            output_dir=tmpdir,
            interactive=False
        )
        
        # Second save without merge or overwrite should raise FileExistsError
        config2 = CodesConfig(
            machine_name='test_machine',
            codes={'ph': Code(name='ph', path='/usr/bin/ph.x', version='7.2')},
        )
        
        with pytest.raises(FileExistsError):
            CodesManager.save_config(
                config2,
                output_dir=tmpdir,
                overwrite=False,
                merge=False,
                interactive=False
            )


def test_codes_config_page_uses_interactive_false():
    """Test that the GUI codes_config page uses interactive=False correctly."""
    # This test verifies the actual GUI code
    import re
    from pathlib import Path
    
    # Read the codes_config.py file
    gui_file = Path(__file__).parent.parent / 'xespresso' / 'gui' / 'pages' / 'codes_config.py'
    
    with open(gui_file, 'r') as f:
        content = f.read()
    
    # Check that CodesManager.save_config is called
    assert 'CodesManager.save_config' in content, "CodesManager.save_config not found in GUI"
    
    # Check that interactive=False is set
    # Look for the pattern: save_config(..., interactive=False)
    pattern = r'CodesManager\.save_config\([^)]*interactive\s*=\s*False'
    matches = re.search(pattern, content, re.DOTALL)
    
    assert matches is not None, "interactive=False not found in CodesManager.save_config call"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

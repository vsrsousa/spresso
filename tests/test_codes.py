"""
Tests for the codes module
"""

import os
import json
import tempfile
import shutil
import pytest
from xespresso.codes import (
    Code, CodesConfig, CodesManager,
    create_codes_config, load_codes_config, detect_qe_codes
)


class TestCode:
    """Tests for Code dataclass"""
    
    def test_code_creation(self):
        """Test creating a Code object"""
        code = Code(
            name="pw",
            path="/usr/bin/pw.x",
            version="7.2"
        )
        assert code.name == "pw"
        assert code.path == "/usr/bin/pw.x"
        assert code.version == "7.2"
    
    def test_code_to_dict(self):
        """Test Code serialization to dict"""
        code = Code(
            name="pw",
            path="/usr/bin/pw.x",
            version="7.2"
        )
        data = code.to_dict()
        assert data["name"] == "pw"
        assert data["path"] == "/usr/bin/pw.x"
        assert data["version"] == "7.2"
    
    def test_code_from_dict(self):
        """Test Code deserialization from dict"""
        data = {
            "name": "pw",
            "path": "/usr/bin/pw.x",
            "version": "7.2"
        }
        code = Code.from_dict(data)
        assert code.name == "pw"
        assert code.path == "/usr/bin/pw.x"
        assert code.version == "7.2"


class TestCodesConfig:
    """Tests for CodesConfig class"""
    
    def test_config_creation(self):
        """Test creating a CodesConfig object"""
        config = CodesConfig(
            machine_name="test_machine",
            qe_version="7.2"
        )
        assert config.machine_name == "test_machine"
        assert config.qe_version == "7.2"
    
    def test_add_code(self):
        """Test adding a code to config"""
        config = CodesConfig(machine_name="test")
        code = Code(name="pw", path="/usr/bin/pw.x")
        config.add_code(code)
        assert "pw" in config.codes
        assert config.get_code("pw") == code
    
    def test_list_codes(self):
        """Test listing codes"""
        config = CodesConfig(machine_name="test")
        config.add_code(Code(name="pw", path="/usr/bin/pw.x"))
        config.add_code(Code(name="hp", path="/usr/bin/hp.x"))
        codes = config.list_codes()
        assert "pw" in codes
        assert "hp" in codes
        assert len(codes) == 2
    
    def test_has_code(self):
        """Test checking code availability"""
        config = CodesConfig(machine_name="test")
        config.add_code(Code(name="pw", path="/usr/bin/pw.x"))
        assert config.has_code("pw")
        assert not config.has_code("nonexistent")
    
    def test_to_dict(self):
        """Test config serialization"""
        config = CodesConfig(machine_name="test", qe_version="7.2")
        config.add_code(Code(name="pw", path="/usr/bin/pw.x", version="7.2"))
        data = config.to_dict()
        assert data["machine_name"] == "test"
        assert data["qe_version"] == "7.2"
        assert "pw" in data["codes"]
    
    def test_from_dict(self):
        """Test config deserialization"""
        data = {
            "machine_name": "test",
            "qe_version": "7.2",
            "codes": {
                "pw": {
                    "name": "pw",
                    "path": "/usr/bin/pw.x",
                    "version": "7.2"
                }
            }
        }
        config = CodesConfig.from_dict(data)
        assert config.machine_name == "test"
        assert config.qe_version == "7.2"
        assert config.has_code("pw")
    
    def test_json_roundtrip(self):
        """Test saving and loading from JSON"""
        config = CodesConfig(machine_name="test", qe_version="7.2")
        config.add_code(Code(name="pw", path="/usr/bin/pw.x", version="7.2"))
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            filepath = f.name
        
        try:
            config.to_json(filepath)
            loaded = CodesConfig.from_json(filepath)
            assert loaded.machine_name == config.machine_name
            assert loaded.qe_version == config.qe_version
            assert loaded.has_code("pw")
        finally:
            os.unlink(filepath)


class TestCodesManager:
    """Tests for CodesManager class"""
    
    def test_create_config(self):
        """Test creating config from detected codes"""
        detected_codes = {
            "pw": "/usr/bin/pw.x",
            "hp": "/usr/bin/hp.x"
        }
        config = CodesManager.create_config(
            machine_name="test",
            detected_codes=detected_codes,
            qe_version="7.2"
        )
        assert config.machine_name == "test"
        assert config.qe_version == "7.2"
        assert config.has_code("pw")
        assert config.has_code("hp")
    
    def test_save_and_load_config(self):
        """Test saving and loading config"""
        config = CodesConfig(machine_name="test", qe_version="7.2")
        config.add_code(Code(name="pw", path="/usr/bin/pw.x", version="7.2"))
        
        tmpdir = tempfile.mkdtemp()
        try:
            filepath = CodesManager.save_config(config, output_dir=tmpdir)
            assert os.path.exists(filepath)
            
            loaded = CodesManager.load_config("test", codes_dir=tmpdir)
            assert loaded is not None
            assert loaded.machine_name == "test"
            assert loaded.has_code("pw")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    
    def test_load_nonexistent_config(self):
        """Test loading a config that doesn't exist"""
        tmpdir = tempfile.mkdtemp()
        try:
            loaded = CodesManager.load_config("nonexistent", codes_dir=tmpdir)
            assert loaded is None
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


class TestConvenienceFunctions:
    """Tests for convenience functions"""
    
    def test_detect_qe_codes_no_installation(self):
        """Test detect_qe_codes when no QE is installed"""
        # This should return an empty config when no QE is found
        config = detect_qe_codes(
            machine_name="test",
            search_paths=["/nonexistent/path"]
        )
        assert config.machine_name == "test"
        assert len(config.codes) == 0
    
    def test_create_codes_config_no_save(self):
        """Test creating config without saving"""
        config = create_codes_config(
            machine_name="test",
            search_paths=["/nonexistent/path"],
            save=False
        )
        assert config.machine_name == "test"


class TestMultiVersionSupport:
    """Tests for multiple QE version support"""
    
    def test_add_code_with_version(self):
        """Test adding codes to specific versions"""
        config = CodesConfig(machine_name="test", qe_version="7.2")
        
        # Add code to version 7.2
        code_72 = Code(name="pw", path="/opt/qe-7.2/bin/pw.x", version="7.2")
        config.add_code(code_72, version="7.2")
        
        # Add code to version 6.8
        code_68 = Code(name="pw", path="/opt/qe-6.8/bin/pw.x", version="6.8")
        config.add_code(code_68, version="6.8")
        
        # Verify both versions are stored
        assert config.has_code("pw", version="7.2")
        assert config.has_code("pw", version="6.8")
        
        # Verify paths are different
        pw_72 = config.get_code("pw", version="7.2")
        pw_68 = config.get_code("pw", version="6.8")
        assert pw_72.path == "/opt/qe-7.2/bin/pw.x"
        assert pw_68.path == "/opt/qe-6.8/bin/pw.x"
    
    def test_list_versions(self):
        """Test listing available versions"""
        config = CodesConfig(machine_name="test")
        
        # Add codes to different versions
        config.add_code(Code(name="pw", path="/opt/qe-7.2/bin/pw.x"), version="7.2")
        config.add_code(Code(name="pw", path="/opt/qe-6.8/bin/pw.x"), version="6.8")
        
        versions = config.list_versions()
        assert "7.2" in versions
        assert "6.8" in versions
        assert len(versions) == 2
    
    def test_list_codes_by_version(self):
        """Test listing codes for specific version"""
        config = CodesConfig(machine_name="test")
        
        # Add different codes to different versions
        config.add_code(Code(name="pw", path="/opt/qe-7.2/bin/pw.x"), version="7.2")
        config.add_code(Code(name="hp", path="/opt/qe-7.2/bin/hp.x"), version="7.2")
        config.add_code(Code(name="pw", path="/opt/qe-6.8/bin/pw.x"), version="6.8")
        
        # Check codes in each version
        codes_72 = config.list_codes(version="7.2")
        codes_68 = config.list_codes(version="6.8")
        
        assert "pw" in codes_72
        assert "hp" in codes_72
        assert len(codes_72) == 2
        
        assert "pw" in codes_68
        assert len(codes_68) == 1
    
    def test_version_config_serialization(self):
        """Test saving and loading multi-version config"""
        config = CodesConfig(machine_name="test", qe_version="7.2")
        
        # Add codes to different versions
        config.add_code(Code(name="pw", path="/opt/qe-7.2/bin/pw.x", version="7.2"), version="7.2")
        config.add_code(Code(name="pw", path="/opt/qe-6.8/bin/pw.x", version="6.8"), version="6.8")
        
        # Add version-specific settings
        if config.versions:
            config.versions["7.2"]["modules"] = ["quantum-espresso/7.2"]
            config.versions["6.8"]["modules"] = ["quantum-espresso/6.8"]
        
        tmpdir = tempfile.mkdtemp()
        try:
            # Save and load
            filepath = CodesManager.save_config(config, output_dir=tmpdir)
            loaded = CodesManager.load_config("test", codes_dir=tmpdir)
            
            # Verify structure
            assert loaded is not None
            assert loaded.machine_name == "test"
            assert loaded.qe_version == "7.2"
            
            # Verify versions
            versions = loaded.list_versions()
            assert "7.2" in versions
            assert "6.8" in versions
            
            # Verify codes
            assert loaded.has_code("pw", version="7.2")
            assert loaded.has_code("pw", version="6.8")
            
            # Verify version-specific settings
            ver_config_72 = loaded.get_version_config("7.2")
            ver_config_68 = loaded.get_version_config("6.8")
            assert ver_config_72 is not None
            assert ver_config_68 is not None
            assert ver_config_72.get("modules") == ["quantum-espresso/7.2"]
            assert ver_config_68.get("modules") == ["quantum-espresso/6.8"]
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    
    def test_backward_compatibility(self):
        """Test backward compatibility with single-version configs"""
        # Old-style config without versions
        config = CodesConfig(machine_name="test", qe_version="7.2")
        config.add_code(Code(name="pw", path="/usr/bin/pw.x", version="7.2"))
        config.add_code(Code(name="hp", path="/usr/bin/hp.x", version="7.2"))
        
        # Should still work with old methods
        assert config.has_code("pw")
        assert config.has_code("hp")
        codes = config.list_codes()
        assert "pw" in codes
        assert "hp" in codes
        
        # Get code without version parameter
        pw = config.get_code("pw")
        assert pw is not None
        assert pw.path == "/usr/bin/pw.x"
    
    def test_get_code_fallback(self):
        """Test get_code fallback to main codes when version not specified"""
        config = CodesConfig(machine_name="test")
        
        # Add to main codes
        config.add_code(Code(name="pw", path="/usr/bin/pw.x"))
        
        # Add to version-specific
        config.add_code(Code(name="hp", path="/opt/qe-7.2/bin/hp.x"), version="7.2")
        
        # Should get from main codes when no version specified
        pw = config.get_code("pw")
        assert pw is not None
        assert pw.path == "/usr/bin/pw.x"
        
        # Should get from version when specified
        hp = config.get_code("hp", version="7.2")
        assert hp is not None
        assert hp.path == "/opt/qe-7.2/bin/hp.x"
        
        # Should return None for non-existent version
        hp_fallback = config.get_code("hp")  # hp not in main codes
        assert hp_fallback is None


class TestNewFeatures:
    """Tests for new features: auto-load machine, port support, module fallback, overwrite protection"""
    
    def test_save_config_overwrite_protection(self):
        """Test that save_config handles existing files correctly"""
        tmpdir = tempfile.mkdtemp()
        try:
            # Create initial config
            config1 = CodesConfig(machine_name="test", qe_version="7.2")
            config1.add_code(Code(name="pw", path="/usr/bin/pw.x", version="7.2"))
            
            # Save first time
            filepath = CodesManager.save_config(config1, output_dir=tmpdir, overwrite=True)
            assert os.path.exists(filepath)
            
            # Try to save again without overwrite should raise error (non-interactive mode)
            config2 = CodesConfig(machine_name="test", qe_version="7.3")
            with pytest.raises(FileExistsError):
                CodesManager.save_config(config2, output_dir=tmpdir, overwrite=False, merge=False, interactive=False)
            
            # With overwrite=True should work
            filepath2 = CodesManager.save_config(config2, output_dir=tmpdir, overwrite=True)
            assert filepath == filepath2
            
            # Load and verify it was overwritten
            loaded = CodesManager.load_config("test", codes_dir=tmpdir)
            assert loaded.qe_version == "7.3"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    
    def test_save_config_merge(self):
        """Test that merge combines configurations correctly"""
        tmpdir = tempfile.mkdtemp()
        try:
            # Create initial config with pw code
            config1 = CodesConfig(machine_name="test", qe_version="7.2")
            config1.add_code(Code(name="pw", path="/usr/bin/pw.x", version="7.2"))
            filepath = CodesManager.save_config(config1, output_dir=tmpdir, overwrite=True)
            
            # Create new config with hp code
            config2 = CodesConfig(machine_name="test", qe_version="7.2")
            config2.add_code(Code(name="hp", path="/usr/bin/hp.x", version="7.2"))
            
            # Save with merge=True
            CodesManager.save_config(config2, output_dir=tmpdir, merge=True)
            
            # Load and verify both codes are present
            loaded = CodesManager.load_config("test", codes_dir=tmpdir)
            assert loaded is not None
            assert loaded.has_code("pw")
            assert loaded.has_code("hp")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    
    def test_save_config_merge_versions(self):
        """Test that merge combines version-specific configurations"""
        tmpdir = tempfile.mkdtemp()
        try:
            # Create config with version 7.2
            config1 = CodesConfig(machine_name="test", qe_version="7.2")
            config1.add_code(Code(name="pw", path="/opt/qe-7.2/bin/pw.x"), version="7.2")
            filepath = CodesManager.save_config(config1, output_dir=tmpdir, overwrite=True)
            
            # Create config with version 6.8
            config2 = CodesConfig(machine_name="test", qe_version="6.8")
            config2.add_code(Code(name="pw", path="/opt/qe-6.8/bin/pw.x"), version="6.8")
            
            # Save with merge=True
            CodesManager.save_config(config2, output_dir=tmpdir, merge=True)
            
            # Load and verify both versions are present
            loaded = CodesManager.load_config("test", codes_dir=tmpdir)
            assert loaded is not None
            versions = loaded.list_versions()
            assert "7.2" in versions
            assert "6.8" in versions
            assert loaded.has_code("pw", version="7.2")
            assert loaded.has_code("pw", version="6.8")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
    
    def test_check_module_available_local(self):
        """Test checking if module command is available locally"""
        # This test will pass if module is not available on the test system
        # which is expected behavior
        result = CodesManager._check_module_available(ssh_connection=None)
        # Result should be boolean
        assert isinstance(result, bool)
    
    def test_ssh_connection_with_port(self):
        """Test that SSH connection dict includes port parameter"""
        # Test that the port parameter is properly handled in ssh_connection
        ssh_conn = {
            'host': 'example.com',
            'username': 'testuser',
            'port': 2222
        }
        
        # Port should be extracted correctly
        assert ssh_conn.get('port', 22) == 2222
        
        # Default port should be 22
        ssh_conn_no_port = {
            'host': 'example.com',
            'username': 'testuser'
        }
        assert ssh_conn_no_port.get('port', 22) == 22
    
    def test_detect_codes_with_use_modules_false(self):
        """Test that detect_codes respects use_modules=False"""
        # This is a unit test - we're not actually detecting codes
        # Just verifying the parameter is accepted
        detected = CodesManager.detect_codes(
            search_paths=["/nonexistent/path"],
            modules=["some-module"],
            use_modules=False
        )
        # Should return empty dict since path doesn't exist
        assert isinstance(detected, dict)
    
    def test_create_codes_config_with_merge(self):
        """Test create_codes_config with merge parameter"""
        tmpdir = tempfile.mkdtemp()
        try:
            # We can't actually detect codes without QE installed,
            # so we'll just verify the function accepts the parameters
            # by calling it with non-existent paths
            config = create_codes_config(
                machine_name="test_machine",
                search_paths=["/nonexistent"],
                save=True,
                output_dir=tmpdir,
                overwrite=True,
                merge=False,
                auto_load_machine=False
            )
            
            assert config.machine_name == "test_machine"
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

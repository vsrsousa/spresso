"""
Tests for the pseudopotential configuration management.
"""

import pytest
import tempfile
import shutil
import json
from pathlib import Path
from xespresso.utils.pseudo_config import (
    get_config_dir,
    ensure_config_dir,
    save_pseudo_config,
    load_pseudo_config,
    list_pseudo_configs,
    delete_pseudo_config,
    get_pseudo_info,
)


@pytest.fixture
def temp_config_dir(monkeypatch, tmp_path):
    """Create a temporary config directory for testing."""
    # Monkeypatch the home directory to use tmp_path
    monkeypatch.setenv("HOME", str(tmp_path))
    
    # Ensure the config directory exists
    config_dir = ensure_config_dir()
    
    yield config_dir
    
    # Cleanup is automatic with tmp_path


def test_get_config_dir():
    """Test that config dir path is correct."""
    config_dir = get_config_dir()
    assert config_dir.name == ".xespresso"
    assert config_dir.parent == Path.home()


def test_ensure_config_dir(temp_config_dir):
    """Test that ensure_config_dir creates the directory."""
    assert temp_config_dir.exists()
    assert temp_config_dir.is_dir()


def test_save_pseudo_config(temp_config_dir):
    """Test saving a pseudopotential configuration."""
    config_data = {
        "name": "test_config",
        "pseudopotentials": {
            "Si": "Si.pbe.UPF",
            "O": "O.pbe.UPF"
        }
    }
    
    save_pseudo_config("test_config", config_data)
    
    config_file = temp_config_dir / "test_config.json"
    assert config_file.exists()
    
    # Verify the content
    with open(config_file, 'r') as f:
        loaded = json.load(f)
    
    assert loaded == config_data


def test_save_pseudo_config_overwrite(temp_config_dir):
    """Test overwrite protection and behavior."""
    config_data = {
        "name": "test_config",
        "pseudopotentials": {"Si": "Si.pbe.UPF"}
    }
    
    # First save should succeed
    save_pseudo_config("test_config", config_data)
    
    # Second save without overwrite should fail
    with pytest.raises(FileExistsError):
        save_pseudo_config("test_config", config_data, overwrite=False)
    
    # Save with overwrite should succeed
    updated_data = {
        "name": "test_config",
        "pseudopotentials": {"Si": "Si.pbe.v2.UPF"}
    }
    save_pseudo_config("test_config", updated_data, overwrite=True)
    
    # Verify it was updated
    loaded = load_pseudo_config("test_config")
    assert loaded["pseudopotentials"]["Si"] == "Si.pbe.v2.UPF"


def test_load_pseudo_config(temp_config_dir):
    """Test loading a pseudopotential configuration."""
    config_data = {
        "name": "test_config",
        "description": "Test configuration",
        "pseudopotentials": {
            "Si": "Si.pbe.UPF",
            "O": "O.pbe.UPF"
        }
    }
    
    save_pseudo_config("test_config", config_data)
    loaded = load_pseudo_config("test_config")
    
    assert loaded == config_data
    assert loaded["name"] == "test_config"
    assert loaded["description"] == "Test configuration"
    assert "Si" in loaded["pseudopotentials"]


def test_load_pseudo_config_not_found(temp_config_dir):
    """Test loading a non-existent configuration."""
    with pytest.raises(FileNotFoundError, match="not found"):
        load_pseudo_config("nonexistent_config")


def test_list_pseudo_configs(temp_config_dir):
    """Test listing all configurations."""
    # Initially empty
    configs = list_pseudo_configs()
    assert configs == []
    
    # Save some configurations
    save_pseudo_config("config1", {"name": "config1", "pseudopotentials": {}})
    save_pseudo_config("config2", {"name": "config2", "pseudopotentials": {}})
    save_pseudo_config("config3", {"name": "config3", "pseudopotentials": {}})
    
    configs = list_pseudo_configs()
    assert len(configs) == 3
    assert "config1" in configs
    assert "config2" in configs
    assert "config3" in configs


def test_list_pseudo_configs_no_directory(monkeypatch, tmp_path):
    """Test listing configs when directory doesn't exist."""
    # Point to a non-existent directory
    monkeypatch.setenv("HOME", str(tmp_path))
    
    configs = list_pseudo_configs()
    assert configs == []


def test_delete_pseudo_config(temp_config_dir):
    """Test deleting a configuration."""
    config_data = {"name": "test_config", "pseudopotentials": {"Si": "Si.pbe.UPF"}}
    
    save_pseudo_config("test_config", config_data)
    
    # Verify it exists
    assert "test_config" in list_pseudo_configs()
    
    # Delete it
    delete_pseudo_config("test_config")
    
    # Verify it's gone
    assert "test_config" not in list_pseudo_configs()


def test_delete_pseudo_config_not_found(temp_config_dir):
    """Test deleting a non-existent configuration."""
    with pytest.raises(FileNotFoundError, match="not found"):
        delete_pseudo_config("nonexistent_config")


def test_get_pseudo_info(temp_config_dir):
    """Test getting pseudopotential info for a specific element."""
    config_data = {
        "name": "test_config",
        "pseudopotentials": {
            "Si": "Si.pbe.UPF",
            "O": "O.pbe.UPF",
            "Fe": "Fe.pbe.UPF"
        }
    }
    
    save_pseudo_config("test_config", config_data)
    
    # Test getting existing element
    pseudo = get_pseudo_info("test_config", "Si")
    assert pseudo == "Si.pbe.UPF"
    
    # Test getting another element
    pseudo = get_pseudo_info("test_config", "Fe")
    assert pseudo == "Fe.pbe.UPF"
    
    # Test getting non-existent element
    pseudo = get_pseudo_info("test_config", "Au")
    assert pseudo is None


def test_get_pseudo_info_no_config(temp_config_dir):
    """Test getting pseudo info from non-existent config."""
    pseudo = get_pseudo_info("nonexistent_config", "Si")
    assert pseudo is None


def test_complex_config_structure(temp_config_dir):
    """Test saving and loading a complex configuration."""
    config_data = {
        "name": "complex_config",
        "description": "A complex configuration with metadata",
        "functional": "PBE",
        "spin_polarized": True,
        "recommended_ecutwfc": 80.0,
        "recommended_ecutrho": 640.0,
        "pseudopotentials": {
            "Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF",
            "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
        },
        "notes": "This is a test configuration with various metadata"
    }
    
    save_pseudo_config("complex_config", config_data)
    loaded = load_pseudo_config("complex_config")
    
    assert loaded == config_data
    assert loaded["functional"] == "PBE"
    assert loaded["spin_polarized"] is True
    assert loaded["recommended_ecutwfc"] == 80.0
    assert "notes" in loaded


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])

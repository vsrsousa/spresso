"""
Tests for the Hubbard parameters module
"""

import pytest
from xespresso.hubbard import HubbardConfig, build_hubbard_str


class TestHubbardConfig:
    """Tests for HubbardConfig class"""
    
    def test_create_empty_config(self):
        """Test creating an empty config"""
        config = HubbardConfig()
        assert config.projector == 'atomic'
        assert len(config.u_params) == 0
        assert len(config.v_params) == 0
    
    def test_add_u_parameter(self):
        """Test adding U parameter"""
        config = HubbardConfig(use_new_format=True)
        config.add_u("Fe", 4.3, orbital="3d")
        assert "Fe-3d" in config.u_params
        assert config.u_params["Fe-3d"] == 4.3
    
    def test_add_u_parameter_old_format(self):
        """Test adding U parameter in old format"""
        config = HubbardConfig(use_new_format=False)
        config.add_u("Fe", 4.3)
        assert "Fe" in config.u_params
        assert config.u_params["Fe"] == 4.3
    
    def test_add_v_parameter(self):
        """Test adding V parameter"""
        config = HubbardConfig(use_new_format=True)
        config.add_v("Fe", "O", 0.5, orbital1="3d", orbital2="2p")
        assert len(config.v_params) == 1
        assert config.v_params[0][0] == "Fe-3d"
        assert config.v_params[0][1] == "O-2p"
        assert config.v_params[0][4] == 0.5
    
    def test_should_use_new_format_explicit(self):
        """Test explicit format selection"""
        config_new = HubbardConfig(use_new_format=True)
        assert config_new.should_use_new_format() == True
        
        config_old = HubbardConfig(use_new_format=False)
        assert config_old.should_use_new_format() == False
    
    def test_should_use_new_format_auto_detect(self):
        """Test automatic format detection"""
        config = HubbardConfig()
        config.add_u("Fe", 4.3, orbital="3d")
        assert config.should_use_new_format() == True
        
        config2 = HubbardConfig()
        config2.add_u("Fe", 4.3)  # No orbital
        assert config2.should_use_new_format() == False
    
    def test_to_new_format_card(self):
        """Test generating new format HUBBARD card"""
        config = HubbardConfig(use_new_format=True)
        config.add_u("Fe", 4.3, orbital="3d")
        config.add_u("O", 3.0, orbital="2p")
        
        lines = config.to_new_format_card()
        assert len(lines) > 0
        assert "HUBBARD" in lines[0]
        assert any("Fe-3d" in line for line in lines)
        assert any("O-2p" in line for line in lines)
    
    def test_from_input_data_new_format(self):
        """Test creating config from input_data with new format"""
        input_data = {
            "qe_version": "7.2",
            "hubbard": {
                "projector": "atomic",
                "u": {
                    "Fe-3d": 4.3,
                    "O-2p": 3.0
                }
            }
        }
        config = HubbardConfig.from_input_data(input_data)
        assert config.should_use_new_format() == True
        assert "Fe-3d" in config.u_params
        assert "O-2p" in config.u_params
    
    def test_from_input_data_old_format(self):
        """Test creating config from input_data with old format"""
        input_data = {
            "input_ntyp": {
                "Hubbard_U": {
                    "Fe": 4.3,
                    "O": 3.0
                }
            }
        }
        config = HubbardConfig.from_input_data(input_data)
        assert "Fe" in config.u_params
        assert "O" in config.u_params
    
    def test_from_input_data_with_hubbard_v(self):
        """Test creating config with old format V parameters"""
        input_data = {
            "input_ntyp": {
                "Hubbard_U": {"Fe": 4.3}
            },
            "hubbard_v": {
                "(1,1,1)": 0.5,
                "(3,3,1)": 1.0
            }
        }
        config = HubbardConfig.from_input_data(input_data)
        assert len(config.v_params) == 2


class TestBuildHubbardStr:
    """Tests for build_hubbard_str function"""
    
    def test_no_hubbard_parameters(self):
        """Test with no Hubbard parameters"""
        input_data = {
            "ecutwfc": 30.0
        }
        species_info = {}
        result = build_hubbard_str(input_data, species_info)
        assert len(result) == 0
    
    def test_with_new_format(self):
        """Test generating HUBBARD card"""
        input_data = {
            "qe_version": "7.2",
            "hubbard": {
                "projector": "atomic",
                "u": {
                    "Fe-3d": 4.3
                }
            }
        }
        species_info = {"Fe": {"index": 1}}
        result = build_hubbard_str(input_data, species_info, qe_version="7.2")
        assert len(result) > 0
        assert any("HUBBARD" in line for line in result)
    
    def test_with_old_format(self):
        """Test with old format (should return empty as it goes to SYSTEM)"""
        input_data = {
            "input_ntyp": {
                "Hubbard_U": {"Fe": 4.3}
            }
        }
        species_info = {"Fe": {"index": 1}}
        result = build_hubbard_str(input_data, species_info)
        # Old format goes to SYSTEM namelist, so this should be empty
        assert len(result) == 0


class TestBackwardCompatibility:
    """Tests for backward compatibility with old format"""
    
    def test_old_format_still_works(self):
        """Test that old format input_data still works"""
        input_data = {
            "lda_plus_u": True,
            "input_ntyp": {
                "Hubbard_U": {
                    "Fe": 4.3,
                    "O": 3.0
                },
                "starting_magnetization": {
                    "Fe": 0.5
                }
            }
        }
        config = HubbardConfig.from_input_data(input_data)
        assert "Fe" in config.u_params
        assert "O" in config.u_params
        assert config.u_params["Fe"] == 4.3
        assert config.u_params["O"] == 3.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

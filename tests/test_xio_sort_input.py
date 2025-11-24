"""
Tests for sort_qe_input function to ensure special parameters are preserved
"""

import pytest
from xespresso.xio import sort_qe_input, check_qe_input


class TestSortQEInput:
    """Tests for sort_qe_input function"""
    
    def test_preserves_hubbard_parameters(self):
        """Test that Hubbard parameters are preserved"""
        input_data = {
            "ecutwfc": 30.0,
            "hubbard": {
                "projector": "atomic",
                "u": {"Fe-3d": 4.3}
            }
        }
        kwargs = {"input_data": input_data, "kpts": (4, 4, 4)}
        
        sorted_params, unused = sort_qe_input(kwargs)
        
        assert "hubbard" in sorted_params["input_data"]
        assert sorted_params["input_data"]["hubbard"]["projector"] == "atomic"
        assert sorted_params["input_data"]["hubbard"]["u"]["Fe-3d"] == 4.3
        assert "hubbard" not in unused
    
    def test_preserves_qe_version(self):
        """Test that qe_version is preserved"""
        input_data = {
            "ecutwfc": 30.0,
            "qe_version": "7.2"
        }
        kwargs = {"input_data": input_data}
        
        sorted_params, unused = sort_qe_input(kwargs)
        
        assert "qe_version" in sorted_params["input_data"]
        assert sorted_params["input_data"]["qe_version"] == "7.2"
        assert "qe_version" not in unused
    
    def test_preserves_hubbard_v(self):
        """Test that hubbard_v is preserved"""
        input_data = {
            "ecutwfc": 30.0,
            "hubbard_v": {"(1,1,1)": 0.5}
        }
        kwargs = {"input_data": input_data}
        
        sorted_params, unused = sort_qe_input(kwargs)
        
        assert "hubbard_v" in sorted_params["input_data"]
        assert sorted_params["input_data"]["hubbard_v"]["(1,1,1)"] == 0.5
        assert "hubbard_v" not in unused
    
    def test_normal_parameters_still_sorted(self):
        """Test that normal parameters are still sorted correctly"""
        input_data = {
            "ecutwfc": 30.0,
            "ecutrho": 240.0,
            "occupations": "smearing",
            "nspin": 2
        }
        kwargs = {"input_data": input_data}
        
        sorted_params, unused = sort_qe_input(kwargs)
        
        # These should be in SYSTEM section
        assert "ecutwfc" in sorted_params["input_data"]["SYSTEM"]
        assert "ecutrho" in sorted_params["input_data"]["SYSTEM"]
        assert "occupations" in sorted_params["input_data"]["SYSTEM"]
        assert "nspin" in sorted_params["input_data"]["SYSTEM"]
    
    def test_combined_parameters(self):
        """Test that both special and normal parameters work together"""
        input_data = {
            "ecutwfc": 30.0,
            "nspin": 2,
            "qe_version": "7.2",
            "hubbard": {"projector": "atomic", "u": {"Fe-3d": 4.3}}
        }
        kwargs = {"input_data": input_data, "kpts": (4, 4, 4)}
        
        sorted_params, unused = sort_qe_input(kwargs)
        
        # Normal parameters in SYSTEM
        assert "ecutwfc" in sorted_params["input_data"]["SYSTEM"]
        assert "nspin" in sorted_params["input_data"]["SYSTEM"]
        
        # Special parameters preserved
        assert "qe_version" in sorted_params["input_data"]
        assert "hubbard" in sorted_params["input_data"]


class TestCheckQEInput:
    """Tests for check_qe_input function"""
    
    def test_handles_special_parameters(self):
        """Test that check_qe_input doesn't fail on special parameters"""
        input_parameters = {
            "CONTROL": {},
            "SYSTEM": {"ecutwfc": 30.0, "nspin": 2},
            "ELECTRONS": {},
            "IONS": {},
            "CELL": {},
            "INPUT_NTYP": {},
            "qe_version": "7.2",
            "hubbard": {"projector": "atomic", "u": {"Fe-3d": 4.3}},
            "hubbard_v": {"(1,1,1)": 0.5}
        }
        
        # This should not raise an exception
        check_qe_input(input_parameters)
    
    def test_checks_normal_parameters(self):
        """Test that check_qe_input still checks normal parameters"""
        input_parameters = {
            "CONTROL": {},
            "SYSTEM": {"ecutwfc": 30.0, "nspin": 2},
            "ELECTRONS": {},
            "IONS": {},
            "CELL": {},
            "INPUT_NTYP": {}
        }
        
        # This should not raise an exception
        check_qe_input(input_parameters)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

"""
Test auto-detection of ecutrho_ratio based on pseudopotential type
"""

from unittest.mock import patch
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from xespresso.pseudopotentials.detector import get_ecutrho_ratio_from_pseudos


class TestEcutrhoRatioAutodetection:
    """Verify ecutrho_ratio is auto-detected from pseudopotential types"""
    
    def test_us_paw_pseudo_returns_8(self):
        """Ultrasoft/PAW pseudos should have ecutrho_ratio = 8.0"""
        # Test the detector function directly
        # Mock parse_upf_header to return PAW type
        with patch('xespresso.pseudopotentials.detector.parse_upf_header') as mock_parse:
            mock_parse.return_value = {
                'type': 'PAW',
                'element': 'Gd',
                'suggested_ecutwfc': 69.0
            }
            
            with patch('xespresso.pseudopotentials.detector.os.path.exists', return_value=True):
                ratio = get_ecutrho_ratio_from_pseudos(
                    {'Gd': 'Gd.paw.UPF'},
                    pseudopotentials_base_path='/fake'
                )
            
            assert ratio == 8.0, f"PAW pseudo should have ratio=8.0, got {ratio}"
            print(f"✓ PAW pseudopotential: ecutrho_ratio = {ratio}")
    
    
    def test_nc_pseudo_returns_4(self):
        """Norm-Conserving pseudos should have ecutrho_ratio = 4.0"""
        with patch('xespresso.pseudopotentials.detector.parse_upf_header') as mock_parse:
            mock_parse.return_value = {
                'type': 'Norm-Conserving',
                'element': 'Si',
                'suggested_ecutwfc': 25.0
            }
            
            with patch('xespresso.pseudopotentials.detector.os.path.exists', return_value=True):
                ratio = get_ecutrho_ratio_from_pseudos(
                    {'Si': 'Si.nc.UPF'},
                    pseudopotentials_base_path='/fake'
                )
            
            assert ratio == 4.0, f"NC pseudo should have ratio=4.0, got {ratio}"
            print(f"✓ Norm-Conserving pseudopotential: ecutrho_ratio = {ratio}")
    
    
    def test_mixed_pseudos_returns_8(self):
        """Mixed types (NC + PAW) should use higher ratio = 8.0"""
        with patch('xespresso.pseudopotentials.detector.parse_upf_header') as mock_parse:
            def side_effect(path):
                if 'Si' in path:
                    return {'type': 'Norm-Conserving', 'element': 'Si', 'suggested_ecutwfc': 25.0}
                else:  # Gd
                    return {'type': 'PAW', 'element': 'Gd', 'suggested_ecutwfc': 69.0}
            
            mock_parse.side_effect = side_effect
            
            with patch('xespresso.pseudopotentials.detector.os.path.exists', return_value=True):
                ratio = get_ecutrho_ratio_from_pseudos(
                    {'Si': 'Si.nc.UPF', 'Gd': 'Gd.paw.UPF'},
                    pseudopotentials_base_path='/fake'
                )
            
            assert ratio == 8.0, f"Mixed types should have ratio=8.0 (higher), got {ratio}"
            print(f"✓ Mixed pseudopotentials (NC + PAW): ecutrho_ratio = {ratio} (uses higher ratio)")
    
    
    def test_unknown_pseudo_returns_8(self):
        """Unknown/unreadable pseudos should default to safe value 8.0"""
        with patch('xespresso.pseudopotentials.detector.parse_upf_header') as mock_parse:
            mock_parse.side_effect = Exception("Can't parse")  # Simulates parse error
            
            with patch('xespresso.pseudopotentials.detector.os.path.exists', return_value=True):
                ratio = get_ecutrho_ratio_from_pseudos(
                    {'X': 'X.unknown.UPF'},
                    pseudopotentials_base_path='/fake'
                )
            
            assert ratio == 8.0, f"Unknown pseudo should default to 8.0, got {ratio}"
            print(f"✓ Unknown pseudopotential type: ecutrho_ratio = {ratio} (safe default)")
    
    
    @patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory')
    def test_workflow_auto_detects_ecutrho(self, mock_discover):
        """Workflow should auto-detect ecutrho_ratio from dict pseudopotentials"""
        atoms = bulk('Gd', 'hcp')
        
        # Mock the discovery to return fake Gd pseudo
        mock_discover.return_value = (
            {'Gd': '/fake/Gd.paw.UPF'},
            '/fake'
        )
        
        # Mock parse_upf_header to return PAW type
        with patch('xespresso.pseudopotentials.detector.parse_upf_header') as mock_parse:
            mock_parse.return_value = {
                'type': 'PAW',
                'element': 'Gd',
                'suggested_ecutwfc': 69.0
            }
            
            with patch('xespresso.pseudopotentials.detector.os.path.exists', return_value=True):
                with patch('xespresso.workflow.convergence_workflow.os.environ', {}):
                    workflow = ConvergenceWorkflow(
                        atoms=atoms,
                        pseudopotentials={'Gd': '/fake/Gd.paw.UPF'},
                        protocol='moderate',
                        precision='low'
                    )
            
            assert workflow.ecutrho_ratio == 8.0, f"Detected ratio should be 8.0, got {workflow.ecutrho_ratio}"
            print(f"✓ ConvergenceWorkflow auto-detected: ecutrho_ratio = {workflow.ecutrho_ratio}")


if __name__ == '__main__':
    test = TestEcutrhoRatioAutodetection()
    
    print("\n" + "="*80)
    print("AUTO-DETECTION OF ECUTRHO_RATIO")
    print("="*80)
    
    print("\nTEST 1: Ultrasoft/PAW pseu returns 8.0")
    test.test_us_paw_pseudo_returns_8()
    
    print("\nTEST 2: Norm-Conserving pseudo returns 4.0")
    test.test_nc_pseudo_returns_4()
    
    print("\nTEST 3: Mixed types (NC + PAW) returns 8.0 (uses higher)")
    test.test_mixed_pseudos_returns_8()
    
    print("\nTEST 4: Unknown pseudo defaults to 8.0")
    test.test_unknown_pseudo_returns_8()
    
    print("\nTEST 5: ConvergenceWorkflow auto-detects from dict")
    test.test_workflow_auto_detects_ecutrho()
    
    print("\n" + "="*80)
    print("✓ ALL ECUTRHO_RATIO AUTO-DETECTION TESTS PASSED")
    print("="*80)
    print("\nSummary:")
    print("  - Norm-Conserving (NC): ecutrho_ratio = 4.0")
    print("  - Ultrasoft (US): ecutrho_ratio = 8.0 (default)")
    print("  - PAW: ecutrho_ratio = 8.0 (default)")
    print("  - Mixed types: uses higher ratio (8.0) to be safe")
    print("  - Unknown types: defaults to 8.0 (safe)")

"""
Test visualization utilities with multiple viewer backends.
"""
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def test_imports():
    """Test that visualization module can be imported."""
    try:
        from xespresso.gui.utils import visualization
        assert hasattr(visualization, 'create_3d_structure_plot')
        assert hasattr(visualization, 'create_x3d_viewer')
        assert hasattr(visualization, 'create_jmol_viewer')
        assert hasattr(visualization, 'create_py3dmol_viewer')
        assert hasattr(visualization, 'launch_ase_viewer')
        assert hasattr(visualization, 'render_structure_viewer')
        print("✓ All visualization functions are defined")
        return True
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False


def test_jmol_viewer_creation():
    """Test JMol viewer HTML generation."""
    try:
        from ase.build import bulk
        from xespresso.gui.utils.visualization import create_jmol_viewer
        
        # Create a simple structure
        atoms = bulk('Fe', 'bcc', a=2.87)
        
        # Generate JMol viewer HTML
        html_content = create_jmol_viewer(atoms)
        
        if html_content is None:
            print("✓ JMol viewer returned None (expected if dependencies missing)")
            return True
        
        # Check that HTML contains expected JMol elements
        assert 'JSmol' in html_content or 'jmol' in html_content.lower(), "JMol HTML should contain JSmol references"
        assert 'script' in html_content.lower(), "JMol HTML should contain script tags"
        print("✓ JMol viewer HTML generated successfully")
        return True
    except Exception as e:
        print(f"✗ JMol viewer test failed: {e}")
        return False


def test_py3dmol_viewer_availability():
    """Test py3Dmol viewer availability check."""
    try:
        from xespresso.gui.utils.visualization import PY3DMOL_AVAILABLE
        print(f"✓ py3Dmol availability: {PY3DMOL_AVAILABLE}")
        return True
    except Exception as e:
        print(f"✗ py3Dmol availability check failed: {e}")
        return False


def test_ase_viewer_availability():
    """Test ASE viewer availability check."""
    try:
        from xespresso.gui.utils.visualization import ASE_VIEWER_AVAILABLE
        print(f"✓ ASE viewer availability: {ASE_VIEWER_AVAILABLE}")
        return True
    except Exception as e:
        print(f"✗ ASE viewer availability check failed: {e}")
        return False


def test_plotly_viewer_availability():
    """Test Plotly viewer availability check."""
    try:
        from xespresso.gui.utils.visualization import PLOTLY_AVAILABLE
        print(f"✓ Plotly availability: {PLOTLY_AVAILABLE}")
        return True
    except Exception as e:
        print(f"✗ Plotly availability check failed: {e}")
        return False


def test_x3d_viewer_creation():
    """Test X3D viewer HTML generation."""
    try:
        from ase.build import bulk
        from xespresso.gui.utils.visualization import create_x3d_viewer
        
        # Create a simple structure
        atoms = bulk('Al', 'fcc', a=4.05)
        
        # Generate X3D viewer HTML
        html_content = create_x3d_viewer(atoms)
        
        if html_content is None:
            print("✓ X3D viewer returned None (expected if dependencies missing)")
            return True
        
        # Check that HTML contains X3D elements
        assert 'x3d' in html_content.lower() or 'X3D' in html_content, "X3D HTML should contain X3D elements"
        print("✓ X3D viewer HTML generated successfully")
        return True
    except Exception as e:
        print(f"✗ X3D viewer test failed: {e}")
        return False


if __name__ == '__main__':
    print("=" * 60)
    print("Testing Visualization Utilities")
    print("=" * 60)
    
    tests = [
        test_imports,
        test_ase_viewer_availability,
        test_plotly_viewer_availability,
        test_py3dmol_viewer_availability,
        test_jmol_viewer_creation,
        test_x3d_viewer_creation,
    ]
    
    results = []
    for test in tests:
        print(f"\nRunning: {test.__name__}")
        print("-" * 60)
        results.append(test())
    
    print("\n" + "=" * 60)
    print(f"Results: {sum(results)}/{len(results)} tests passed")
    print("=" * 60)
    
    sys.exit(0 if all(results) else 1)

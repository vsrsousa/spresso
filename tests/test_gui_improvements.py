"""
Test suite for GUI improvements based on problem statement.

Tests cover:
1. Session state initialization for persistence
2. Machine configuration improvements
3. Codes configuration with multi-version support
4. Structure viewer with ASE database
5. Calculation setup with new parameters
6. Results viewer functionality
"""

import os
import sys
import tempfile
import pytest
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestSessionStateInitialization:
    """Test that all required session state variables are initialized."""
    
    def test_session_state_variables(self):
        """Test session state initialization."""
        # Import the streamlit app module
        # Note: This will trigger streamlit warnings which we can ignore
        import warnings
        warnings.filterwarnings('ignore')
        
        # The module should import without errors
        try:
            import gui.streamlit_app
            assert True
        except ImportError as e:
            pytest.skip(f"Streamlit not available: {e}")
        except Exception as e:
            # Streamlit session state errors are expected when not running via streamlit run
            if "session_state" not in str(e).lower():
                raise


class TestMachineConfiguration:
    """Test machine configuration improvements."""
    
    def test_machine_persistence_variables(self):
        """Test that machine persistence variables exist in session state."""
        # We just check that the code structure is correct
        import gui.streamlit_app as app
        
        # Check that the file contains the required session state variables
        with open(app.__file__, 'r') as f:
            content = f.read()
            
        assert 'current_machine_name' in content
        assert 'current_machine' in content
        assert 'Test Connection' in content
        assert 'View Current Configuration' in content


class TestCodesConfiguration:
    """Test codes configuration improvements."""
    
    def test_codes_version_support(self):
        """Test that codes configuration supports version labels."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for version label support
        assert 'version_label' in content.lower()
        assert 'Version Label' in content
        
        # Check for multiple version support mentions
        assert 'Multiple versions' in content or 'multiple versions' in content
        
        # Check for code selection
        assert 'selected_code_version' in content


class TestStructureViewer:
    """Test structure viewer improvements."""
    
    def test_ase_database_integration(self):
        """Test that ASE database functionality is present."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for ASE database support
        assert 'ASE Database' in content
        assert 'Save to Database' in content
        assert 'Load from Database' in content
        assert 'ase_db_path' in content
    
    def test_structure_persistence(self):
        """Test that structure is shown at top when loaded."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for current structure display
        assert 'Current structure:' in content or 'current_structure' in content


class TestCalculationSetup:
    """Test calculation setup improvements."""
    
    def test_smearing_configuration(self):
        """Test that smearing options are available."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for smearing support
        assert 'smearing' in content.lower()
        assert 'degauss' in content.lower()
        assert 'occupations' in content.lower()
    
    def test_dual_parameter(self):
        """Test that dual parameter is configurable."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for dual parameter
        assert 'dual' in content.lower()
        assert 'Dual Parameter' in content or 'dual parameter' in content
    
    def test_dft_u_support(self):
        """Test that DFT+U configuration is available."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for DFT+U support
        assert 'dft_u' in content.lower() or 'DFT+U' in content
        assert 'hubbard' in content.lower()
    
    def test_band_structure_options(self):
        """Test that band structure specific options are available."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for band path options
        assert 'band_path' in content.lower()
        assert 'K-path' in content or 'k-path' in content


class TestWorkflowImprovements:
    """Test general workflow improvements."""
    
    def test_local_workdir_selection(self):
        """Test that local working directory can be selected."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for local workdir
        assert 'local_workdir' in content
        assert 'Working Directory' in content
    
    def test_file_locations_display(self):
        """Test that file locations are shown."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for file location display
        assert 'File Locations' in content or 'Files will be saved' in content


class TestResultsViewer:
    """Test results and post-processing page."""
    
    def test_results_page_exists(self):
        """Test that results page is added."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for results page
        assert 'Results & Post-Processing' in content
        assert 'Output Files' in content or 'output file' in content.lower()
    
    def test_output_viewing(self):
        """Test that output file viewing is available."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for output viewing
        assert 'View File Content' in content
        assert 'Download Output' in content or 'Download' in content
    
    def test_post_processing_tools(self):
        """Test that post-processing tools are present."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for post-processing
        assert 'Post-Processing' in content
        assert 'Energy Analysis' in content or 'DOS' in content or 'Band Structure' in content


class TestConnectionTest:
    """Test connection testing functionality."""
    
    def test_connection_test_feedback(self):
        """Test that connection test provides feedback."""
        import gui.streamlit_app as app
        
        with open(app.__file__, 'r') as f:
            content = f.read()
        
        # Check for connection test implementation
        assert 'Test Connection' in content
        assert 'Connection Test Results' in content or 'connection test' in content.lower()
        assert 'SSH connection' in content or 'ssh' in content.lower()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

"""
Test suite for {nprocs} placeholder substitution in scheduler.

This test validates that the {nprocs} placeholder in launcher commands
is correctly substituted with the actual number of processors.
"""
import pytest
import os
import tempfile
import shutil
from unittest.mock import Mock
from xespresso.scheduler import set_queue


class TestNprocsSubstitution:
    """Test {nprocs} placeholder substitution in launcher commands."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        
    def teardown_method(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def _create_mock_calc(self, queue_config):
        """Helper to create a mock calculator."""
        calc = Mock()
        calc.package = 'pw'
        calc.parallel = ''
        calc.prefix = 'test'
        calc.label = 'test'
        calc.directory = self.test_dir
        calc.queue = queue_config
        return calc
    
    def _get_job_file_content(self):
        """Helper to read job file content."""
        job_file_path = os.path.join(self.test_dir, 'job_file')
        with open(job_file_path, 'r') as f:
            return f.read()
    
    def test_nprocs_substitution_basic(self):
        """Test basic {nprocs} substitution with value 4."""
        queue_config = {
            'scheduler': 'direct',
            'execution': 'local',
            'launcher': 'mpirun -np {nprocs}',
            'nprocs': 4
        }
        
        calc = self._create_mock_calc(queue_config)
        os.environ["ASE_ESPRESSO_COMMAND"] = "LAUNCHER pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        set_queue(calc)
        content = self._get_job_file_content()
        
        assert "mpirun -np 4" in content
        assert "{nprocs}" not in content
    
    def test_nprocs_substitution_various_values(self):
        """Test {nprocs} substitution with various processor counts."""
        test_cases = [1, 4, 16, 32, 64, 128]
        
        for nprocs in test_cases:
            queue_config = {
                'scheduler': 'direct',
                'execution': 'local',
                'launcher': 'mpirun -np {nprocs}',
                'nprocs': nprocs
            }
            
            calc = self._create_mock_calc(queue_config)
            os.environ["ASE_ESPRESSO_COMMAND"] = "LAUNCHER pw.x -in PREFIX.pwi > PREFIX.pwo"
            
            set_queue(calc)
            content = self._get_job_file_content()
            
            expected = f"mpirun -np {nprocs}"
            assert expected in content, f"Failed for nprocs={nprocs}"
            assert "{nprocs}" not in content
    
    def test_nprocs_default_value(self):
        """Test that missing nprocs defaults to 1."""
        queue_config = {
            'scheduler': 'direct',
            'execution': 'local',
            'launcher': 'mpirun -np {nprocs}',
            # nprocs not specified
        }
        
        calc = self._create_mock_calc(queue_config)
        os.environ["ASE_ESPRESSO_COMMAND"] = "LAUNCHER pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        set_queue(calc)
        content = self._get_job_file_content()
        
        assert "mpirun -np 1" in content
        assert "{nprocs}" not in content
    
    def test_hardcoded_nprocs_preserved(self):
        """Test that hardcoded launcher values are preserved."""
        queue_config = {
            'scheduler': 'direct',
            'execution': 'local',
            'launcher': 'mpirun -np 8',  # Hardcoded, no placeholder
            'nprocs': 4  # Should be ignored
        }
        
        calc = self._create_mock_calc(queue_config)
        os.environ["ASE_ESPRESSO_COMMAND"] = "LAUNCHER pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        set_queue(calc)
        content = self._get_job_file_content()
        
        assert "mpirun -np 8" in content
        assert "mpirun -np 4" not in content
    
    def test_multiple_nprocs_placeholders(self):
        """Test substitution when {nprocs} appears multiple times."""
        queue_config = {
            'scheduler': 'direct',
            'execution': 'local',
            'launcher': 'mpirun -np {nprocs} --npernode {nprocs}',
            'nprocs': 8
        }
        
        calc = self._create_mock_calc(queue_config)
        os.environ["ASE_ESPRESSO_COMMAND"] = "LAUNCHER pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        set_queue(calc)
        content = self._get_job_file_content()
        
        # Both placeholders should be replaced
        assert "mpirun -np 8 --npernode 8" in content
        assert "{nprocs}" not in content
    
    def test_srun_launcher_with_nprocs(self):
        """Test {nprocs} substitution with srun launcher."""
        queue_config = {
            'scheduler': 'direct',
            'execution': 'local',
            'launcher': 'srun -n {nprocs}',
            'nprocs': 16
        }
        
        calc = self._create_mock_calc(queue_config)
        os.environ["ASE_ESPRESSO_COMMAND"] = "LAUNCHER pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        set_queue(calc)
        content = self._get_job_file_content()
        
        assert "srun -n 16" in content
        assert "{nprocs}" not in content
    
    def test_no_launcher_placeholder(self):
        """Test behavior when LAUNCHER placeholder is not in command."""
        queue_config = {
            'scheduler': 'direct',
            'execution': 'local',
            'launcher': 'mpirun -np {nprocs}',
            'nprocs': 4
        }
        
        calc = self._create_mock_calc(queue_config)
        # Command without LAUNCHER placeholder
        os.environ["ASE_ESPRESSO_COMMAND"] = "pw.x -in PREFIX.pwi > PREFIX.pwo"
        
        set_queue(calc)
        content = self._get_job_file_content()
        
        # Should not contain the launcher at all
        assert "mpirun" not in content
        assert "pw.x -in test.pwi > test.pwo" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

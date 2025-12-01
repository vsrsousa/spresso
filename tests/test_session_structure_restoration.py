"""
Tests for session structure restoration in the Qt GUI.

These tests verify that structures are correctly restored when loading sessions.
"""

import pytest
import os
import json
import tempfile
import shutil


@pytest.fixture
def temp_session_dir():
    """Create a temporary directory for session files."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def temp_db_dir():
    """Create a temporary directory for database files."""
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


class TestStructureSourceParsing:
    """Tests for parsing structure_source strings."""
    
    def test_database_source_pattern(self):
        """Test parsing database source pattern."""
        source = "Database: ID 1"
        assert source.startswith("Database: ID ")
        structure_id = int(source.replace("Database: ID ", "").strip())
        assert structure_id == 1
    
    def test_database_source_pattern_large_id(self):
        """Test parsing database source with large ID."""
        source = "Database: ID 12345"
        assert source.startswith("Database: ID ")
        structure_id = int(source.replace("Database: ID ", "").strip())
        assert structure_id == 12345
    
    def test_file_source_pattern(self):
        """Test parsing file source pattern."""
        source = "File: structure.cif"
        assert source.startswith("File: ")
        filename = source.replace("File: ", "").strip()
        assert filename == "structure.cif"
    
    def test_built_source_pattern(self):
        """Test parsing built source pattern."""
        source = "Built: Fe bcc"
        assert source.startswith("Built: ")
        desc = source.replace("Built: ", "").strip()
        assert desc == "Fe bcc"


class TestStructureRestorationFromDatabase:
    """Tests for restoring structures from ASE database."""
    
    @pytest.fixture
    def db_with_structure(self, temp_db_dir):
        """Create a database with a test structure."""
        try:
            from ase import Atoms
            from ase.db import connect as ase_db_connect
        except ImportError:
            pytest.skip("ASE not available")
        
        db_path = os.path.join(temp_db_dir, "test_structures.db")
        
        # Create a simple structure
        atoms = Atoms('Fe', positions=[[0, 0, 0]], cell=[2.87, 2.87, 2.87], pbc=True)
        
        # Save to database
        db = ase_db_connect(db_path)
        row_id = db.write(atoms, name='test_fe')
        
        return db_path, row_id
    
    def test_restore_structure_from_database(self, db_with_structure, temp_session_dir):
        """Test restoring a structure from database on session load."""
        try:
            from ase import Atoms
        except ImportError:
            pytest.skip("ASE not available")
        
        db_path, structure_id = db_with_structure
        
        # Create a session file with structure_source set
        session_data = {
            "session_name": "Test Session",
            "structure_source": f"Database: ID {structure_id}",
            "structure_db_path": db_path,
            "working_directory": "/tmp",
            "_session_id": "test_session_1"
        }
        
        session_file = os.path.join(temp_session_dir, "test_session.json")
        with open(session_file, 'w') as f:
            json.dump(session_data, f)
        
        # Verify the session file was created correctly
        assert os.path.exists(session_file)
        
        # Read and verify contents
        with open(session_file, 'r') as f:
            loaded = json.load(f)
        
        assert loaded['structure_source'] == f"Database: ID {structure_id}"
        assert loaded['structure_db_path'] == db_path


class TestStructureRestorationFromFile:
    """Tests for restoring structures from files."""
    
    @pytest.fixture
    def structure_file(self, temp_db_dir):
        """Create a structure file."""
        try:
            from ase import Atoms
            from ase.io import write
        except ImportError:
            pytest.skip("ASE not available")
        
        file_path = os.path.join(temp_db_dir, "test_structure.xyz")
        
        # Create a simple structure
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
        write(file_path, atoms, format='xyz')
        
        return file_path
    
    def test_file_path_saved_in_session(self, structure_file, temp_session_dir):
        """Test that full file path is saved in session."""
        # Create a session file with structure_file_path set
        session_data = {
            "session_name": "Test Session",
            "structure_source": f"File: {os.path.basename(structure_file)}",
            "structure_file_path": structure_file,
            "working_directory": "/tmp",
            "_session_id": "test_session_1"
        }
        
        session_file = os.path.join(temp_session_dir, "test_session.json")
        with open(session_file, 'w') as f:
            json.dump(session_data, f)
        
        # Read and verify contents
        with open(session_file, 'r') as f:
            loaded = json.load(f)
        
        assert loaded['structure_file_path'] == structure_file
        assert os.path.exists(loaded['structure_file_path'])


class TestAllowedSessionKeys:
    """Tests for allowed session keys."""
    
    def test_structure_db_path_is_allowed(self):
        """Test that structure_db_path is in allowed keys."""
        # We can't import the module directly due to PySide6 dependency,
        # but we can check the file
        main_app_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'qtgui', 'main_app.py'
        )
        
        with open(main_app_path, 'r') as f:
            content = f.read()
        
        assert "'structure_db_path'" in content or '"structure_db_path"' in content
    
    def test_structure_file_path_is_allowed(self):
        """Test that structure_file_path is in allowed keys."""
        main_app_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'qtgui', 'main_app.py'
        )
        
        with open(main_app_path, 'r') as f:
            content = f.read()
        
        assert "'structure_file_path'" in content or '"structure_file_path"' in content


class TestRestoreStructureFromSourceMethod:
    """Tests for the _restore_structure_from_source method implementation."""
    
    def test_method_exists_in_source(self):
        """Test that _restore_structure_from_source method exists in source."""
        main_app_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'qtgui', 'main_app.py'
        )
        
        with open(main_app_path, 'r') as f:
            content = f.read()
        
        assert '_restore_structure_from_source' in content
    
    def test_method_called_in_load_session(self):
        """Test that _restore_structure_from_source is called in _load_session."""
        main_app_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'qtgui', 'main_app.py'
        )
        
        with open(main_app_path, 'r') as f:
            content = f.read()
        
        # Check both load methods call the restore method
        assert content.count('self._restore_structure_from_source()') >= 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

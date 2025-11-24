"""
Tests for magnetic and Hubbard UI persistence in sessions.
"""

import pytest
import tempfile
import os
import json
from xespresso.gui.utils.session_manager import save_session, load_session, restore_session
from unittest.mock import Mock, patch
import streamlit as st


class TestMagneticHubbardPersistence:
    """Test that magnetic and Hubbard configurations persist in sessions."""
    
    @patch('streamlit.session_state', {})
    def test_magnetic_config_saved_in_session(self):
        """
        Test that magnetic configuration is saved when session is saved.
        """
        # Set up workflow_config with magnetic settings
        st.session_state['workflow_config'] = {
            'enable_magnetism': True,
            'magnetic_config': {
                'Fe': [1.0, -1.0],
                'O': [0]
            },
            'expand_cell': True,
            'label': 'scf/Fe2O3',
            'calc_type': 'scf'
        }
        st.session_state['selected_machine'] = 'local'
        st.session_state['working_directory'] = '/tmp/test'
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save session
            session_file = save_session(
                filename='test_magnetic.json',
                session_dir=tmpdir,
                session_name='Test Magnetic Session'
            )
            
            # Verify file was created
            assert os.path.exists(session_file)
            
            # Load and verify contents
            with open(session_file, 'r') as f:
                session_data = json.load(f)
            
            # Check that workflow_config is in the saved state
            assert 'state' in session_data
            assert 'workflow_config' in session_data['state']
            
            # Verify magnetic configuration was saved
            saved_config = session_data['state']['workflow_config']
            assert saved_config['enable_magnetism'] is True
            assert 'magnetic_config' in saved_config
            assert saved_config['magnetic_config']['Fe'] == [1.0, -1.0]
    
    @patch('streamlit.session_state', {})
    def test_hubbard_config_saved_in_session(self):
        """
        Test that Hubbard configuration is saved when session is saved.
        """
        # Set up workflow_config with Hubbard settings
        st.session_state['workflow_config'] = {
            'enable_hubbard': True,
            'hubbard_u': {
                'Fe': 4.3,
                'Ni': 6.0
            },
            'hubbard_format': 'new',
            'hubbard_projector': 'atomic',
            'label': 'scf/FeNi',
            'calc_type': 'scf'
        }
        st.session_state['selected_machine'] = 'local'
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save session
            session_file = save_session(
                filename='test_hubbard.json',
                session_dir=tmpdir,
                session_name='Test Hubbard Session'
            )
            
            # Load and verify contents
            with open(session_file, 'r') as f:
                session_data = json.load(f)
            
            # Verify Hubbard configuration was saved
            saved_config = session_data['state']['workflow_config']
            assert saved_config['enable_hubbard'] is True
            assert 'hubbard_u' in saved_config
            assert saved_config['hubbard_u']['Fe'] == 4.3
            assert saved_config['hubbard_u']['Ni'] == 6.0
            assert saved_config['hubbard_format'] == 'new'
    
    @patch('streamlit.session_state', {})
    def test_combined_magnetic_hubbard_saved(self):
        """
        Test that both magnetic and Hubbard configurations are saved together.
        """
        # Set up workflow_config with both magnetic and Hubbard settings
        st.session_state['workflow_config'] = {
            'enable_magnetism': True,
            'magnetic_config': {
                'Fe': [2.0, -2.0]
            },
            'enable_hubbard': True,
            'hubbard_u': {
                'Fe': 4.3
            },
            'hubbard_format': 'new',
            'label': 'scf/Fe',
            'calc_type': 'scf'
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save session
            session_file = save_session(
                filename='test_combined.json',
                session_dir=tmpdir
            )
            
            # Load and verify both configurations are saved
            with open(session_file, 'r') as f:
                session_data = json.load(f)
            
            saved_config = session_data['state']['workflow_config']
            
            # Verify magnetic config
            assert saved_config['enable_magnetism'] is True
            assert saved_config['magnetic_config']['Fe'] == [2.0, -2.0]
            
            # Verify Hubbard config
            assert saved_config['enable_hubbard'] is True
            assert saved_config['hubbard_u']['Fe'] == 4.3
    
    @patch('streamlit.session_state', {})
    def test_session_restore_preserves_magnetic_hubbard(self):
        """
        Test that loading a session restores magnetic and Hubbard configurations.
        """
        # Create initial state
        initial_config = {
            'enable_magnetism': True,
            'magnetic_config': {'Fe': [1.5]},
            'enable_hubbard': True,
            'hubbard_u': {'Fe': 5.0},
            'label': 'scf/Fe'
        }
        st.session_state['workflow_config'] = initial_config
        st.session_state['test_value'] = 'preserved'
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save session
            session_file = save_session(
                filename='test_restore.json',
                session_dir=tmpdir
            )
            
            # Clear session state (simulate new session)
            st.session_state.clear()
            
            # Load session
            state, session_name = load_session(session_file)
            
            # Restore to session_state
            for key, value in state.items():
                st.session_state[key] = value
            
            # Verify workflow_config was restored
            assert 'workflow_config' in st.session_state
            restored_config = st.session_state['workflow_config']
            
            # Check magnetic config
            assert restored_config['enable_magnetism'] is True
            assert restored_config['magnetic_config']['Fe'] == [1.5]
            
            # Check Hubbard config
            assert restored_config['enable_hubbard'] is True
            assert restored_config['hubbard_u']['Fe'] == 5.0


class TestUIAutoExpand:
    """Test that UI elements auto-expand when configurations are present."""
    
    def test_magnetic_expander_logic(self):
        """
        Test the logic for determining if magnetic expanders should be expanded.
        """
        config = {
            'magnetic_config': {
                'Fe': [1.0, -1.0],  # Configured - should expand
                'O': [0],           # Default - should not expand
                'Al': []            # Empty - should not expand
            }
        }
        
        # Test Fe (configured with non-zero values)
        current_mag_fe = config['magnetic_config'].get('Fe', [0])
        is_configured_fe = bool(
            current_mag_fe 
            and current_mag_fe != [0] 
            and (len(current_mag_fe) > 1 or (len(current_mag_fe) == 1 and current_mag_fe[0] != 0))
        )
        assert is_configured_fe is True
        
        # Test O (default [0])
        current_mag_o = config['magnetic_config'].get('O', [0])
        is_configured_o = bool(
            current_mag_o 
            and current_mag_o != [0] 
            and (len(current_mag_o) > 1 or (len(current_mag_o) == 1 and current_mag_o[0] != 0))
        )
        assert is_configured_o is False
        
        # Test Al (empty list)
        current_mag_al = config['magnetic_config'].get('Al', [0])
        is_configured_al = bool(
            current_mag_al 
            and current_mag_al != [0] 
            and (len(current_mag_al) > 1 or (len(current_mag_al) == 1 and current_mag_al[0] != 0))
        )
        assert is_configured_al is False
        
        # Test non-existent element (defaults to [0])
        current_mag_cu = config['magnetic_config'].get('Cu', [0])
        is_configured_cu = bool(
            current_mag_cu 
            and current_mag_cu != [0] 
            and (len(current_mag_cu) > 1 or (len(current_mag_cu) == 1 and current_mag_cu[0] != 0))
        )
        assert is_configured_cu is False
    
    def test_hubbard_expander_logic(self):
        """
        Test the logic for determining if Hubbard expanders should be expanded.
        """
        config = {
            'hubbard_u': {
                'Fe': 4.3,  # Configured - should expand
                'Ni': 0.0,  # Zero - should not expand
            }
        }
        
        # Test Fe (configured)
        is_configured_fe = 'Fe' in config['hubbard_u'] and config['hubbard_u']['Fe'] != 0.0
        assert is_configured_fe is True
        
        # Test Ni (zero)
        is_configured_ni = 'Ni' in config['hubbard_u'] and config['hubbard_u']['Ni'] != 0.0
        assert is_configured_ni is False
        
        # Test non-existent element
        is_configured_cu = 'Cu' in config['hubbard_u'] and config['hubbard_u']['Cu'] != 0.0
        assert is_configured_cu is False


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

"""
Database interface for querying provenance and structure databases.
"""

import streamlit as st
from pathlib import Path
from .config import STRUCTURES_DB, PROVENANCE_DB
import pandas as pd
import logging

logger = logging.getLogger(__name__)


@st.cache_resource
def get_databases():
    """Load database connections (cached)."""
    try:
        import ase.db
        from xespresso.db import DatabaseWorkflow
        
        wf = DatabaseWorkflow(
            db_path=str(STRUCTURES_DB),
            provenance_path=str(PROVENANCE_DB)
        )
        return wf
    except Exception as e:
        logger.error(f"Failed to load databases: {e}")
        return None


def get_all_structures():
    """Get all structures from ASE database."""
    wf = get_databases()
    if not wf:
        return pd.DataFrame()
    
    try:
        structures = []
        for row in wf.db.select():
            structures.append({
                'id': row.id,
                'formula': row.get('formula', 'N/A'),
                'natoms': row.natoms,
                'energy': row.get('energy'),
                'calculator': row.get('calculator'),
                'timestamp': row.get('timestamp'),
            })
        return pd.DataFrame(structures)
    except Exception as e:
        logger.error(f"Error fetching structures: {e}")
        return pd.DataFrame()


def get_structure_details(structure_id: int):
    """Get detailed information about a structure."""
    wf = get_databases()
    if not wf:
        return None
    
    try:
        return wf.db.get_atoms(structure_id)
    except Exception as e:
        logger.error(f"Error fetching structure {structure_id}: {e}")
        return None


def get_all_calculations():
    """Get all calculations from provenance database."""
    wf = get_databases()
    if not wf:
        return pd.DataFrame()
    
    try:
        cursor = wf.provenance.conn.execute('''
            SELECT 
                id, calculation_hash, calculation_method, energy, 
                convergence_steps, success, timestamp
            FROM calculations
            ORDER BY timestamp DESC
        ''')
        
        calculations = []
        for row in cursor:
            calculations.append({
                'id': row[0],
                'hash': row[1][:8] + '...',
                'method': row[2],
                'energy': row[3],
                'steps': row[4],
                'success': row[5],
                'timestamp': row[6],
            })
        return pd.DataFrame(calculations)
    except Exception as e:
        logger.error(f"Error fetching calculations: {e}")
        return pd.DataFrame()


def get_derivations(source_structure_id: int):
    """Get all derivations from a source structure."""
    wf = get_databases()
    if not wf:
        return pd.DataFrame()
    
    try:
        cursor = wf.provenance.conn.execute('''
            SELECT source_structure_id, derived_structure_id, 
                   derivation_method, energy_change, derivation_timestamp
            FROM derivations
            WHERE source_structure_id = ?
            ORDER BY derivation_timestamp DESC
        ''', (source_structure_id,))
        
        derivations = []
        for row in cursor:
            derivations.append({
                'source': row[0],
                'derived': row[1],
                'method': row[2],
                'energy_change': row[3],
                'timestamp': row[4],
            })
        return pd.DataFrame(derivations)
    except Exception as e:
        logger.error(f"Error fetching derivations: {e}")
        return pd.DataFrame()


def get_structure_creator(structure_id: int):
    """Get the calculation that created a structure (if it's a derived structure)."""
    wf = get_databases()
    if not wf:
        return None
    
    try:
        cursor = wf.provenance.conn.execute('''
            SELECT c.id, c.calculation_hash, c.calculation_method, c.energy, 
                   c.convergence_steps, c.success, c.timestamp, c.input_structure_id,
                   d.source_structure_id, d.energy_change
            FROM calculations c
            JOIN derivations d ON c.output_structure_id = d.derived_structure_id
            WHERE d.derived_structure_id = ?
            LIMIT 1
        ''', (structure_id,))
        
        row = cursor.fetchone()
        if row:
            return {
                'calc_id': row[0],
                'hash': row[1][:8] + '...',
                'method': row[2],
                'energy': row[3],
                'steps': row[4],
                'success': row[5],
                'timestamp': row[6],
                'input_structure_id': row[7],
                'source_structure_id': row[8],
                'energy_change': row[9],
            }
        return None
    except Exception as e:
        logger.error(f"Error fetching structure creator: {e}")
        return None


def get_execution_history(calculation_id: int):
    """Get execution history for a calculation."""
    wf = get_databases()
    if not wf:
        return pd.DataFrame()
    
    try:
        cursor = wf.provenance.conn.execute('''
            SELECT machine, wall_time, xespresso_version, 
                   qe_version, success, execution_timestamp
            FROM execution_history
            WHERE calculation_id = ?
            ORDER BY execution_timestamp DESC
        ''', (calculation_id,))
        
        executions = []
        for row in cursor:
            executions.append({
                'machine': row[0],
                'wall_time': row[1],
                'xespresso_version': row[2],
                'qe_version': row[3],
                'success': row[4],
                'timestamp': row[5],
            })
        return pd.DataFrame(executions)
    except Exception as e:
        logger.error(f"Error fetching execution history: {e}")
        return pd.DataFrame()


def search_structures(formula: str = None, energy_min: float = None, 
                     energy_max: float = None):
    """Search structures by criteria."""
    all_structs = get_all_structures()
    
    if all_structs.empty:
        return pd.DataFrame()
    
    result = all_structs.copy()
    
    if formula:
        result = result[result['formula'].str.contains(formula, case=False, na=False)]
    
    if energy_min is not None:
        result = result[result['energy'] >= energy_min]
    
    if energy_max is not None:
        result = result[result['energy'] <= energy_max]
    
    return result


def get_statistics():
    """Get overall statistics."""
    wf = get_databases()
    if not wf:
        return {}
    
    try:
        structs = wf.db.select()
        calcs = wf.provenance.conn.execute('SELECT COUNT(*) as cnt FROM calculations')
        derivs = wf.provenance.conn.execute('SELECT COUNT(*) as cnt FROM derivations')
        exec_hist = wf.provenance.conn.execute('SELECT COUNT(*) as cnt FROM execution_history')
        
        return {
            'total_structures': len(list(structs)),
            'total_calculations': calcs.fetchone()[0],
            'total_derivations': derivs.fetchone()[0],
            'total_executions': exec_hist.fetchone()[0],
        }
    except Exception as e:
        logger.error(f"Error fetching statistics: {e}")
        return {}

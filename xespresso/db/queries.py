"""
Query and analysis functions for ASE Database and Provenance.

Provides utilities for data retrieval, filtering, and large-scale analysis.
"""

import logging
import json
from typing import List, Dict, Optional

import pandas as pd

logger = logging.getLogger(__name__)


def get_structure_history(provenance_db, structure_id: int) -> List[Dict]:
    """
    Get complete lineage of a structure (original → current).
    
    Args:
        provenance_db: ProvenanceDB instance
        structure_id: ASE database structure ID
    
    Returns:
        List of dicts showing structure evolution with energies and methods
    """
    return provenance_db.get_calculation_history(structure_id)


def get_all_derivatives(provenance_db, structure_id: int) -> List[Dict]:
    """
    Get all structures derived from a given structure.
    
    Args:
        provenance_db: ProvenanceDB instance
        structure_id: ASE database structure ID
    
    Returns:
        List of derived structure information
    """
    derivations = provenance_db.query_derivations(structure_id)
    
    results = []
    for deriv in derivations:
        results.append({
            'id': deriv['id'],
            'derived_structure_id': deriv['derived_structure_id'],
            'method': deriv['derivation_method'],
            'energy_change': deriv['energy_change'],
            'timestamp': deriv['derivation_timestamp']
        })
    
    return results


def validate_consistency(provenance_db, structure_id: int) -> bool:
    """
    Verify all calculations on a structure used consistent parameters.
    
    Args:
        provenance_db: ProvenanceDB instance
        structure_id: ASE database structure ID
    
    Returns:
        bool: True if consistent, False otherwise
    """
    calcs = provenance_db.query_by_structure(structure_id)
    
    if not calcs:
        logger.info("No calculations found for this structure")
        return True
    
    first_params = json.loads(calcs[0]['convergence_params']) if calcs[0]['convergence_params'] else None
    
    inconsistencies = []
    for calc in calcs[1:]:
        calc_params = json.loads(calc['convergence_params']) if calc['convergence_params'] else None
        if calc_params != first_params:
            inconsistencies.append({
                'calculation_id': calc['id'],
                'expected': first_params,
                'actual': calc_params
            })
    
    if inconsistencies:
        logger.warning(f"⚠️  Found {len(inconsistencies)} inconsistent calculations")
        for inc in inconsistencies:
            logger.warning(f"  Calculation {inc['calculation_id']}: different parameters")
        return False
    
    logger.info(f"✅ All calculations on structure {structure_id} are consistent")
    return True


def query_by_protocol(provenance_db, protocol: str) -> List[Dict]:
    """
    Get all calculations using a specific protocol.
    
    Args:
        provenance_db: ProvenanceDB instance
        protocol: Protocol name (fast, moderate, accurate)
    
    Returns:
        List of calculation records
    """
    calcs = provenance_db.query_all_by_protocol(protocol)
    
    results = []
    for calc in calcs:
        results.append({
            'id': calc['id'],
            'hash': calc['calculation_hash'][:8],
            'energy': calc['energy'],
            'method': calc['calculation_method'],
            'convergence_steps': calc['convergence_steps'],
            'timestamp': calc['timestamp']
        })
    
    logger.info(f"Found {len(results)} calculations with protocol '{protocol}'")
    return results


def query_by_element(ase_db, provenance_db, element: str) -> List[Dict]:
    """
    Get all calculations containing a specific element.
    
    Args:
        ase_db: ASE Database connection
        provenance_db: ProvenanceDB instance
        element: Element symbol (e.g., 'Fe', 'O')
    
    Returns:
        List of calculation records with element
    """
    # Query ASE database for structures containing element
    rows = ase_db.select(f'molecule==False')  # Get all structures
    
    matching = []
    for row in rows:
        atoms = row.toatoms()
        if element in atoms.get_chemical_symbols():
            calc = provenance_db.query_by_structure(row.id)
            if calc:
                for c in calc:
                    matching.append({
                        'structure_id': row.id,
                        'element': element,
                        'formula': atoms.get_chemical_formula(),
                        'energy': c['energy'],
                        'method': c['calculation_method'],
                        'timestamp': c['timestamp']
                    })
    
    logger.info(f"Found {len(matching)} calculations containing '{element}'")
    return matching


def export_to_dataframe(ase_db, provenance_db, filter_protocol: Optional[str] = None) -> pd.DataFrame:
    """
    Export calculations to pandas DataFrame for analysis.
    
    Args:
        ase_db: ASE Database connection
        provenance_db: ProvenanceDB instance
        filter_protocol: Optional protocol filter
    
    Returns:
        pandas DataFrame with calculation data
    """
    data = []
    
    # Get all calculations
    if filter_protocol:
        calcs = provenance_db.query_all_by_protocol(filter_protocol)
    else:
        calcs = provenance_db.conn.execute('SELECT * FROM calculations').fetchall()
    
    for calc in calcs:
        conv_params = json.loads(calc['convergence_params']) if calc['convergence_params'] else {}
        
        # Get structure formula if available
        formula = None
        if calc['input_structure_id']:
            try:
                row = ase_db.get_atoms(calc['input_structure_id'])
                formula = row.toatoms().get_chemical_formula()
            except:
                pass
        
        # Get execution history
        executions = provenance_db.query_executions(calc['id'])
        machines = [e['machine'] for e in executions if e['machine']]
        wall_times = [e['wall_time'] for e in executions if e['wall_time']]
        
        data.append({
            'calculation_id': calc['id'],
            'hash': calc['calculation_hash'][:8],
            'formula': formula,
            'energy': calc['energy'],
            'method': calc['calculation_method'],
            'convergence_steps': calc['convergence_steps'],
            'protocol': conv_params.get('protocol'),
            'ecutwfc': conv_params.get('ecutwfc'),
            'kspacing': conv_params.get('kspacing'),
            'conv_thr': conv_params.get('conv_thr'),
            'machines': ','.join(machines) if machines else None,
            'wall_time_max': max(wall_times) if wall_times else None,
            'timestamp': calc['timestamp']
        })
    
    df = pd.DataFrame(data)
    logger.info(f"Exported {len(df)} calculations to DataFrame")
    return df


def get_calculation_stats(provenance_db) -> Dict:
    """
    Get overall statistics about calculations.
    
    Args:
        provenance_db: ProvenanceDB instance
    
    Returns:
        Dict with statistics
    """
    calcs = provenance_db.conn.execute(
        'SELECT COUNT(*) as total, calculation_method FROM calculations GROUP BY calculation_method'
    ).fetchall()
    
    methods = {row['calculation_method']: row['total'] for row in calcs}
    
    all_calcs = provenance_db.conn.execute(
        'SELECT * FROM calculations'
    ).fetchall()
    
    energies = [c['energy'] for c in all_calcs if c['energy'] is not None]
    convergence_steps = [c['convergence_steps'] for c in all_calcs if c['convergence_steps'] is not None]
    
    return {
        'total_calculations': len(all_calcs),
        'by_method': methods,
        'energy_min': min(energies) if energies else None,
        'energy_max': max(energies) if energies else None,
        'energy_mean': sum(energies) / len(energies) if energies else None,
        'avg_convergence_steps': sum(convergence_steps) / len(convergence_steps) if convergence_steps else None,
    }


def compare_structures(ase_db, struct_id_1: int, struct_id_2: int) -> Dict:
    """
    Compare two structures and their calculations.
    
    Args:
        ase_db: ASE Database connection
        struct_id_1: First structure ID
        struct_id_2: Second structure ID
    
    Returns:
        Dict with comparison information
    """
    row1 = ase_db.get_atoms(struct_id_1)
    row2 = ase_db.get_atoms(struct_id_2)
    
    atoms1 = row1.toatoms()
    atoms2 = row2.toatoms()
    
    return {
        'structure_1': {
            'id': struct_id_1,
            'formula': atoms1.get_chemical_formula(),
            'energy': row1.get('energy'),
            'calculation_method': row1.get('calculation_method')
        },
        'structure_2': {
            'id': struct_id_2,
            'formula': atoms2.get_chemical_formula(),
            'energy': row2.get('energy'),
            'calculation_method': row2.get('calculation_method')
        },
        'energy_difference': (row2.get('energy') or 0) - (row1.get('energy') or 0),
    }

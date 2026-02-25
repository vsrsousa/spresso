"""
Provenance tracking using SQLite database.

Tracks calculation history, dependencies, and execution metadata.
"""

import sqlite3
import logging
from datetime import datetime
from pathlib import Path
import json

logger = logging.getLogger(__name__)


class ProvenanceDB:
    """
    SQLite-based provenance database for tracking calculations.
    
    Stores:
    - Calculation metadata (hash, energy, convergence)
    - Execution history (machine, wall_time, versions)
    - Derivation relationships (parent structure → derived structure)
    - Dependencies between calculations
    """
    
    def __init__(self, db_path='~/.xespresso/provenance.db'):
        """
        Initialize provenance database.
        
        Args:
            db_path: Path to SQLite database file (default: ~/.xespresso/provenance.db)
        """
        self.db_path = Path(db_path).expanduser()
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row
        self._create_tables()
        logger.info(f"ProvenanceDB initialized at {self.db_path}")
    
    def _create_tables(self):
        """Create database schema."""
        # Main calculations table
        self.conn.execute('''
            CREATE TABLE IF NOT EXISTS calculations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                calculation_hash TEXT UNIQUE NOT NULL,
                input_structure_id INTEGER,
                output_structure_id INTEGER,
                energy REAL,
                convergence_steps INTEGER,
                convergence_params TEXT,
                calculation_method TEXT,
                timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
                
                INDEX idx_hash (calculation_hash),
                INDEX idx_input_structure (input_structure_id),
                INDEX idx_output_structure (output_structure_id),
                INDEX idx_timestamp (timestamp)
            )
        ''')
        
        # Execution history table
        self.conn.execute('''
            CREATE TABLE IF NOT EXISTS execution_history (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                calculation_id INTEGER NOT NULL,
                machine TEXT,
                wall_time REAL,
                xespresso_version TEXT,
                qe_version TEXT,
                execution_timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
                
                FOREIGN KEY(calculation_id) REFERENCES calculations(id),
                INDEX idx_calculation (calculation_id),
                INDEX idx_machine (machine)
            )
        ''')
        
        # Structure derivation table
        self.conn.execute('''
            CREATE TABLE IF NOT EXISTS derivations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                source_structure_id INTEGER NOT NULL,
                derived_structure_id INTEGER NOT NULL,
                derivation_method TEXT,
                energy_change REAL,
                derivation_timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
                
                FOREIGN KEY(source_structure_id) REFERENCES calculations(input_structure_id),
                FOREIGN KEY(derived_structure_id) REFERENCES calculations(output_structure_id),
                INDEX idx_source (source_structure_id),
                INDEX idx_derived (derived_structure_id)
            )
        ''')
        
        # Dependencies table
        self.conn.execute('''
            CREATE TABLE IF NOT EXISTS dependencies (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                calculation_id INTEGER NOT NULL,
                depends_on INTEGER NOT NULL,
                dependency_type TEXT,
                reason TEXT,
                
                FOREIGN KEY(calculation_id) REFERENCES calculations(id),
                FOREIGN KEY(depends_on) REFERENCES calculations(id),
                INDEX idx_calculation (calculation_id),
                INDEX idx_depends_on (depends_on)
            )
        ''')
        
        self.conn.commit()
    
    def log_calculation(self, calc_hash, input_structure_id=None, output_structure_id=None,
                       energy=None, convergence_steps=None, convergence_params=None,
                       calculation_method=None):
        """
        Log a calculation to the database.
        
        Args:
            calc_hash: Unique hash of calculation inputs
            input_structure_id: ASE database ID of input structure
            output_structure_id: ASE database ID of output structure
            energy: Total energy in eV
            convergence_steps: Number of SCF/structural iterations
            convergence_params: Dict of convergence parameters used
            calculation_method: Type of calculation (scf, relax, phonon, etc)
        
        Returns:
            int: Calculation ID in provenance database
        """
        params_json = json.dumps(convergence_params) if convergence_params else None
        
        cursor = self.conn.execute(
            '''INSERT INTO calculations 
               (calculation_hash, input_structure_id, output_structure_id, 
                energy, convergence_steps, convergence_params, calculation_method)
               VALUES (?, ?, ?, ?, ?, ?, ?)''',
            (calc_hash, input_structure_id, output_structure_id,
             energy, convergence_steps, params_json, calculation_method)
        )
        self.conn.commit()
        
        calc_id = cursor.lastrowid
        logger.info(f"Logged calculation {calc_id} with hash {calc_hash[:8]}...")
        return calc_id
    
    def log_execution(self, calculation_id, machine=None, wall_time=None,
                     xespresso_version=None, qe_version=None):
        """
        Log execution details for a calculation (can be called multiple times per calculation).
        
        Args:
            calculation_id: ID from log_calculation()
            machine: Machine where calculation was executed
            wall_time: Wall clock time in seconds
            xespresso_version: Version of xespresso used
            qe_version: Version of Quantum ESPRESSO used
        
        Returns:
            int: Execution history ID
        """
        cursor = self.conn.execute(
            '''INSERT INTO execution_history
               (calculation_id, machine, wall_time, xespresso_version, qe_version)
               VALUES (?, ?, ?, ?, ?)''',
            (calculation_id, machine, wall_time, xespresso_version, qe_version)
        )
        self.conn.commit()
        
        logger.info(f"Logged execution of calculation {calculation_id} on {machine}")
        return cursor.lastrowid
    
    def log_derivation(self, source_structure_id, derived_structure_id,
                      derivation_method=None, energy_change=None):
        """
        Log structural derivation (e.g., relaxation, adsorbate addition).
        
        Args:
            source_structure_id: ASE database ID of parent structure
            derived_structure_id: ASE database ID of derived structure
            derivation_method: Type of derivation (vc-relax, relax, adsorbate, etc)
            energy_change: Energy difference (derived - source) in eV
        
        Returns:
            int: Derivation ID
        """
        cursor = self.conn.execute(
            '''INSERT INTO derivations
               (source_structure_id, derived_structure_id, derivation_method, energy_change)
               VALUES (?, ?, ?, ?)''',
            (source_structure_id, derived_structure_id, derivation_method, energy_change)
        )
        self.conn.commit()
        
        logger.info(f"Logged derivation: {source_structure_id} → {derived_structure_id}")
        return cursor.lastrowid
    
    def log_dependency(self, calculation_id, depends_on, dependency_type=None, reason=None):
        """
        Log dependency between calculations.
        
        Args:
            calculation_id: ID of dependent calculation
            depends_on: ID of required calculation
            dependency_type: Type of dependency (convergence-params, relaxed-structure, etc)
            reason: Human-readable reason for dependency
        """
        self.conn.execute(
            '''INSERT INTO dependencies
               (calculation_id, depends_on, dependency_type, reason)
               VALUES (?, ?, ?, ?)''',
            (calculation_id, depends_on, dependency_type, reason)
        )
        self.conn.commit()
        logger.info(f"Logged dependency: {calculation_id} depends on {depends_on}")
    
    def query_by_hash(self, calc_hash):
        """
        Find calculation by hash.
        
        Args:
            calc_hash: Calculation hash string
        
        Returns:
            Row object or None if not found
        """
        return self.conn.execute(
            'SELECT * FROM calculations WHERE calculation_hash = ?',
            (calc_hash,)
        ).fetchone()
    
    def query_by_structure(self, structure_id):
        """
        Find all calculations related to a structure (as input or output).
        
        Args:
            structure_id: ASE database structure ID
        
        Returns:
            List of Row objects
        """
        return self.conn.execute(
            '''SELECT * FROM calculations 
               WHERE input_structure_id = ? OR output_structure_id = ?
               ORDER BY timestamp DESC''',
            (structure_id, structure_id)
        ).fetchall()
    
    def query_executions(self, calculation_id):
        """
        Get all execution records for a calculation.
        
        Args:
            calculation_id: Calculation ID
        
        Returns:
            List of Row objects
        """
        return self.conn.execute(
            'SELECT * FROM execution_history WHERE calculation_id = ? ORDER BY execution_timestamp',
            (calculation_id,)
        ).fetchall()
    
    def query_derivations(self, structure_id):
        """
        Get all structures derived from a given structure.
        
        Args:
            structure_id: ASE database structure ID
        
        Returns:
            List of Row objects
        """
        return self.conn.execute(
            'SELECT * FROM derivations WHERE source_structure_id = ? ORDER BY derivation_timestamp',
            (structure_id,)
        ).fetchall()
    
    def query_all_by_protocol(self, protocol):
        """
        Get all calculations using a specific protocol.
        
        Args:
            protocol: Protocol name (fast, moderate, accurate)
        
        Returns:
            List of Row objects
        """
        return self.conn.execute(
            '''SELECT * FROM calculations 
               WHERE convergence_params LIKE ?
               ORDER BY timestamp DESC''',
            (f'%"protocol": "{protocol}"%',)
        ).fetchall()
    
    def query_all_by_method(self, method):
        """
        Get all calculations of a specific type.
        
        Args:
            method: Calculation method (scf, relax, phonon, band, etc)
        
        Returns:
            List of Row objects
        """
        return self.conn.execute(
            'SELECT * FROM calculations WHERE calculation_method = ? ORDER BY timestamp DESC',
            (method,)
        ).fetchall()
    
    def get_calculation_history(self, structure_id):
        """
        Get complete lineage of a structure (parent → current).
        
        Args:
            structure_id: ASE database structure ID
        
        Returns:
            List of dicts containing structure evolution
        """
        # This is a simplified version; in practice would need recursive query
        derivations = self.conn.execute(
            '''SELECT source_structure_id, derived_structure_id, derivation_method, 
                      energy_change, derivation_timestamp
               FROM derivations 
               WHERE derived_structure_id = ?
               ORDER BY derivation_timestamp DESC''',
            (structure_id,)
        ).fetchall()
        
        history = []
        current = structure_id
        
        while current:
            calc = self.conn.execute(
                'SELECT * FROM calculations WHERE output_structure_id = ?',
                (current,)
            ).fetchone()
            
            if calc:
                history.append({
                    'id': current,
                    'calculation_id': calc['id'],
                    'energy': calc['energy'],
                    'method': calc['calculation_method'],
                    'convergence_params': json.loads(calc['convergence_params']) if calc['convergence_params'] else None,
                    'timestamp': calc['timestamp']
                })
                
                # Find parent
                parent = self.conn.execute(
                    'SELECT source_structure_id FROM derivations WHERE derived_structure_id = ?',
                    (current,)
                ).fetchone()
                
                current = parent['source_structure_id'] if parent else None
            else:
                break
        
        return history[::-1]  # Reverse to show parent → child order
    
    def close(self):
        """Close database connection."""
        self.conn.close()
        logger.info("ProvenanceDB connection closed")

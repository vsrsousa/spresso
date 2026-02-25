"""
Structure and Directory Overview for Database + Provenance Integration
========================================================================

New files added to xespresso:

xespresso/db/
├── __init__.py                 - Package initialization, exports main classes
├── provenance.py               - ProvenanceDB class (SQLite backend)
├── database_workflow.py         - DatabaseWorkflow class (integration layer)
└── queries.py                  - Query and analysis functions

examples/
└── database_workflow_example.py - Full usage example


FILE DESCRIPTIONS
=================

1. xespresso/db/__init__.py
   - Exports: ProvenanceDB, DatabaseWorkflow, query functions
   - Single import point for users

2. xespresso/db/provenance.py
   - ProvenanceDB class
   - Creates SQLite database with 4 tables:
     * calculations - Main results and energies
     * execution_history - Machine/version info
     * derivations - Structure parent-child relationships
     * dependencies - Calculation dependencies
   - Methods:
     * log_calculation() - Store new calculation
     * log_execution() - Record where it ran
     * log_derivation() - Track structure evolution
     * query_by_hash() - Find by calculation hash
     * query_by_structure() - Find related calculations
     * get_calculation_history() - Get full lineage

3. xespresso/db/database_workflow.py
   - DatabaseWorkflow class (wraps CalculationWorkflow)
   - Features:
     * get_or_calculate() - Cache + execute logic
     * _compute_calculation_hash() - Based on QE params + structure only
     * get_structure_with_params() - Retrieve + convergence params
     * run_properties_on_structure() - Inherit params for derived calcs
     * log_structure_derivation() - Track relaxations/optimizations
     * validate_consistency() - Check parameter agreement
   - Key insight: hash based on QE params, not machine

4. xespresso/db/queries.py
   - Standalone query functions (don't require DatabaseWorkflow instance)
   - Functions:
     * get_structure_history() - Full lineage
     * get_all_derivatives() - Child structures
     * validate_consistency() - Parameter checking
     * query_by_protocol() - Filter by fast/moderate/accurate
     * query_by_element() - Find structures containing element
     * export_to_dataframe() - Convert to pandas for analysis
     * get_calculation_stats() - Overall statistics
     * compare_structures() - Compare two structures

5. examples/database_workflow_example.py
   - Real workflow example:
     1. Initialize DatabaseWorkflow
     2. Run SCF on original structure (calculates)
     3. Run relax (calculates)
     4. Run phonon on relaxed (inherits params)
     5. Run SCF again (cached! no recalculation)
     6. Analyze evolution
     7. Export to DataFrame


USAGE PATTERN
=============

Basic workflow:

    from xespresso.db import DatabaseWorkflow
    
    # Initialize
    db = DatabaseWorkflow()
    
    # First time: calculates
    result1, from_cache = db.get_or_calculate(
        atoms,
        {'protocol': 'moderate', 'machine': 'medusa'},
        calculation_method='scf'
    )
    
    # Second time: returns from cache
    result2, from_cache = db.get_or_calculate(
        atoms,  # Same structure
        {'protocol': 'moderate', 'machine': 'medusa'},  # Same params
        calculation_method='scf'
    )
    # from_cache = True, result2 = result1


Structure derivation:

    # Calculate relaxation
    relax_result, _ = db.get_or_calculate(
        atoms,
        params,
        calculation_method='relax',
        input_structure_id=original_id
    )
    relaxed_id = len(db.db) - 1
    
    # Log derivation
    db.log_structure_derivation(
        source_structure_id=original_id,
        derived_structure_id=relaxed_id,
        derivation_method='vc-relax',
        energy_change=energy_diff
    )
    
    # Properties on relaxed structure
    phonon, _ = db.run_properties_on_structure(
        structure_id=relaxed_id,
        calc_type='phonon'
        # Automatically inherits convergence params from relaxation
    )


Database querying:

    from xespresso.db.queries import (
        get_structure_history,
        export_to_dataframe,
        get_calculation_stats
    )
    
    # Get evolution
    history = get_structure_history(db.provenance, struct_id)
    for step in history:
        print(f"{step['method']}: {step['energy']} eV")
    
    # Export for analysis
    df = export_to_dataframe(db.db, db.provenance, filter_protocol='moderate')
    print(df.describe())
    
    # Statistics
    stats = get_calculation_stats(db.provenance)
    print(f"Total calcs: {stats['total_calculations']}")


DATABASE FILES
==============

~/.xespresso/
├── database.db          ← ASE Database (binary)
│                          Stores: Atoms, energies, calculation metadata
│
└── provenance.db        ← SQLite (binary)
                           Stores: Calculation hashes, history, dependencies


HASH COMPUTATION
================

Calculation hash includes ONLY:
  - Structure (composition + atomic positions)
  - Quantum mechanics parameters:
    * protocol (fast/moderate/accurate)
    * ecutwfc, ecutrho
    * kspacing
    * conv_thr
    * mixing_beta
    * electron_maxstep
    * pseudopotentials
    * magnetic_config

Does NOT include:
  - Machine (medusa/local/cluster)
  - Timestamp
  - Execution version

This means: same calculation computed on different machines uses same cache.


CONVERGENCE PARAMETERS INHERITANCE
===================================

Workflow:

1. User runs relaxation with protocol='moderate'
   - Stores convergence_params in database
   - enc={ecutwfc=50, kspacing=0.3, conv_thr=1e-8, ...}

2. User runs phonon on relaxed structure
   - DatabaseWorkflow automatically retrieves saved params
   - run_properties_on_structure() uses inherited params
   - Phonon calculation guaranteed to use same parameters

3. User runs band structure on relaxed structure
   - Same inherited parameters apply
   - Consistency validated automatically


SCIENTIFIC BENEFITS
====================

✅ No Redundancy
   - Same calculation never computed twice
   - Saves computational resources

✅ Consistency
   - Related calculations use consistent parameters
   - No accidental parameter changes

✅ Traceability
   - Complete history of structure evolution
   - Know which machine was used for each calculation
   - Reproducibility guaranteed

✅ Large-Scale Analysis
   - Export to pandas DataFrame
   - Filter by element, protocol, method
   - Statistical analysis on thousands of calculations

✅ Scientific Integrity
   - Automatic validation of parameter consistency
   - Provenance trail for publications
   - Audit trail of all modifications


EXAMPLES IN CODE
================

See: examples/database_workflow_example.py

Demonstrates:
1. Initialize DatabaseWorkflow
2. Store and cache SCF calculations
3. Track structure evolution through relaxation
4. Inherit convergence parameters
5. Query calculation history
6. Export to DataFrame for analysis
7. Validate consistency
"""

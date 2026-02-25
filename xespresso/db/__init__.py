"""
Database and Provenance tracking for xespresso.

This module provides:
- ASE Database integration for storing structures and calculation results
- ProvenanceDB for tracking calculation history and dependencies
- DatabaseWorkflow for automated caching and result retrieval
- Query functions for analysis and provenance tracking
"""

from xespresso.db.provenance import ProvenanceDB
from xespresso.db.database_workflow import DatabaseWorkflow
from xespresso.db.queries import (
    get_structure_history,
    get_all_derivatives,
    validate_consistency,
    query_by_protocol,
    query_by_element,
    export_to_dataframe,
)

__all__ = [
    'ProvenanceDB',
    'DatabaseWorkflow',
    'get_structure_history',
    'get_all_derivatives',
    'validate_consistency',
    'query_by_protocol',
    'query_by_element',
    'export_to_dataframe',
]

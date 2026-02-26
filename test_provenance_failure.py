#!/usr/bin/env python3
"""
Test provenance failure logging
"""

from xespresso.db import DatabaseWorkflow
import tempfile
import os

# Create temporary databases for testing
with tempfile.TemporaryDirectory() as tmpdir:
    db_path = os.path.join(tmpdir, 'test.db')
    ase_db_path = os.path.join(tmpdir, 'ase.db')

    wf = DatabaseWorkflow(db_path=db_path, ase_db_path=ase_db_path)

    # Test logging a failed calculation
    calc_hash = "test_failed_hash_123"
    calc_id = wf.provenance.log_calculation(
        calc_hash=calc_hash,
        calculation_method='scf',
        success=False,
        error_message="Test failure: job failed on remote cluster"
    )

    print(f"Logged failed calculation with ID: {calc_id}")

    # Test that query_by_hash doesn't return failed calculations
    result = wf.provenance.query_by_hash(calc_hash)
    if result is None:
        print("✅ Correctly: Failed calculation not returned by query_by_hash")
    else:
        print("❌ Error: Failed calculation was returned by query_by_hash")

    # Test logging successful calculation
    calc_hash_success = "test_success_hash_456"
    calc_id_success = wf.provenance.log_calculation(
        calc_hash=calc_hash_success,
        calculation_method='scf',
        success=True
    )

    print(f"Logged successful calculation with ID: {calc_id_success}")

    # Test that successful calculation is returned
    result_success = wf.provenance.query_by_hash(calc_hash_success)
    if result_success:
        print("✅ Correctly: Successful calculation returned by query_by_hash")
        print(f"   Success status: {result_success['success']}")
    else:
        print("❌ Error: Successful calculation not returned by query_by_hash")

print("\n🎉 Provenance failure logging test completed!")
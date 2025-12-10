"""
Tests for job monitor session state access.

These tests verify that the job monitor correctly accesses
working_directory from session state using dictionary-style access.
"""

import os


def test_job_monitor_methods_use_correct_access():
    """Test that job monitor helper methods use correct session state access."""
    # Read the source file directly to avoid import issues
    job_submission_file = os.path.join(
        os.path.dirname(__file__), 
        "..", 
        "qtgui", 
        "pages", 
        "job_submission.py"
    )
    
    with open(job_submission_file, 'r') as f:
        source = f.read()
    
    # Find the _open_job_monitor method
    open_monitor_start = source.find("def _open_job_monitor(self):")
    assert open_monitor_start != -1, "_open_job_monitor method not found"
    
    # Get a chunk of code that includes the method (next 500 chars should be enough)
    open_monitor_chunk = source[open_monitor_start:open_monitor_start + 500]
    
    # Find the _add_job_to_monitor method
    add_job_start = source.find("def _add_job_to_monitor(self, job_info):")
    assert add_job_start != -1, "_add_job_to_monitor method not found"
    
    # Get a chunk of code that includes the method
    add_job_chunk = source[add_job_start:add_job_start + 500]
    
    # Verify they use .get() method, not attribute access
    assert 'session_state.get("working_directory"' in open_monitor_chunk, \
        "_open_job_monitor should use session_state.get() method"
    assert 'session_state.get("working_directory"' in add_job_chunk, \
        "_add_job_to_monitor should use session_state.get() method"
    
    # Verify they don't use incorrect attribute access
    assert 'session_state.working_directory' not in open_monitor_chunk, \
        "_open_job_monitor should not use attribute access"
    assert 'session_state.working_directory' not in add_job_chunk, \
        "_add_job_to_monitor should not use attribute access"

"""Pytest setup helpers for running tests from alternative working dirs.

Ensure the project root (this file's directory) is on sys.path so that
`import xespresso` works even when pytest is invoked from a different
path such as `/scratch/projects/...`.
"""
import os
import sys

# Insert project root at front of sys.path
_ROOT = os.path.abspath(os.path.dirname(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

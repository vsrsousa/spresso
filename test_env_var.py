#!/usr/bin/env python3
import sys
import os

print(f"BEFORE import xespresso: ASE_ESPRESSO_COMMAND = {os.environ.get('ASE_ESPRESSO_COMMAND', 'NOT SET')}")

from ase.build import bulk
from xespresso import CalculationWorkflow

print(f"AFTER import xespresso: ASE_ESPRESSO_COMMAND = {os.environ.get('ASE_ESPRESSO_COMMAND', 'NOT SET')}")

atoms = bulk("Si", cubic=True)

workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials_config="default",
    protocol="moderate",
    machine="snake5"
)

print(f"AFTER creating workflow: ASE_ESPRESSO_COMMAND = {os.environ.get('ASE_ESPRESSO_COMMAND', 'NOT SET')}")

workflow.write_input(label="test_debug/si-test")

print(f"AFTER write_input: ASE_ESPRESSO_COMMAND = {os.environ.get('ASE_ESPRESSO_COMMAND', 'NOT SET')}")

# Check job file
job_file = 'test_debug/si-test/job_file'
if os.path.exists(job_file):
    with open(job_file) as f:
        print("\n===== JOB FILE =====")
        print(f.read())

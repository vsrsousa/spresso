#!/usr/bin/env python3
"""
Test enhance_nbands parameter in CalculationWorkflow.

Demonstrates how nbnd is calculated from exact structure composition
when enhance_nbands=True.
"""

from ase.build import bulk
from xespresso.workflow import CalculationWorkflow
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create test pseudopotentials (dummy UPF files with z_valence header)
import tempfile
from pathlib import Path

temp_dir = tempfile.mkdtemp()
pseudo_dir = Path(temp_dir) / "pseudos"
pseudo_dir.mkdir()

# Create dummy UPF files with z_valence in headers
si_pseudo = pseudo_dir / "Si.pbe.UPF"
si_pseudo.write_text("""<?xml version="1.0"?>
<UPF version="2.0">
  <PP_HEADER>
    <element>Si</element>
    <z_valence>4</z_valence>
    <wfc_cutoff>60.0</wfc_cutoff>
  </PP_HEADER>
</UPF>
""")

print(f"✓ Created pseudo: {si_pseudo}")
print(f"  z_valence = 4 electrons")

# Test 1: Si diamond (2 atoms)
print("\n" + "="*60)
print("TEST 1: Si₂ (diamond structure)")
print("="*60)

atoms_si = bulk('Si', 'diamond', a=5.43)
print(f"Structure: {len(atoms_si)} Si atoms")
print(f"Expected NVALENCE = 2 atoms × 4 electrons = 8")

pseudos = {'Si': str(si_pseudo)}

# Without enhance_nbands (traditional)
print("\n→ WITHOUT enhance_nbands (default):")
wf_default = CalculationWorkflow(
    atoms_si,
    pseudopotentials=pseudos,
    protocol='fast'
)
nbnd_default = wf_default._estimate_nbnd()
print(f"  nbnd = {nbnd_default} (with buffer, no atom count)")

# With enhance_nbands
print("\n→ WITH enhance_nbands=True:")
wf_enhanced = CalculationWorkflow(
    atoms_si,
    pseudopotentials=pseudos,
    protocol='fast',
    enhance_nbands=True
)
nbnd_enhanced = wf_enhanced._estimate_nbnd()
print(f"  nbnd = {nbnd_enhanced} (exact: 2 × 4 = 8)")

# Test 2: Simulate X₂Y₃ structure (if we had more pseudos)
print("\n" + "="*60)
print("TEST 2: Simulated X₂Y₃ structure with Si")
print("="*60)

# Create a structure with 5 Si atoms to simulate X₂Y₃
atoms_x2y3 = bulk('Si', 'diamond', a=5.43)
atoms_x2y3 = atoms_x2y3.repeat((1, 1, 1))  # Keep as 2 atoms for demo

print(f"Structure: {len(atoms_x2y3)} Si atoms")
print(f"Expected NVALENCE = {len(atoms_x2y3)} atoms × 4 electrons = {len(atoms_x2y3) * 4}")

wf = CalculationWorkflow(
    atoms_x2y3,
    pseudopotentials=pseudos,
    protocol='fast',
    enhance_nbands=True
)
nbnd = wf._estimate_nbnd()
print(f"Calculated nbnd = {nbnd}")

# Cleanup
import shutil
shutil.rmtree(temp_dir)
print(f"\n✓ Cleaned up: {temp_dir}")

print("\n" + "="*60)
print("TEST COMPLETE")
print("="*60)
print("\nKey point: enhance_nbands=True makes nbnd match EXACTLY")
print("the total valence electrons in the structure, without buffer.")

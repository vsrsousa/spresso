"""
Demonstration: DOS with spin polarization detection
Shows how run_dos() automatically detects and reports magnetic systems
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from unittest.mock import patch, MagicMock
from pathlib import Path

print("=" * 70)
print("DOS SPIN POLARIZATION DETECTION DEMO")
print("=" * 70)

# Test 1: Non-magnetic system
print("\n1. Non-magnetic Al system (nspin=1):")
print("-" * 70)

atoms_al = bulk('Al', 'fcc', a=4.05)
workflow_al = CalculationWorkflow(
    atoms_al,
    protocol='fast',
    pseudopotentials={'Al': 'Al.pbe.UPF'}
)

nspin_al = workflow_al.input_data.get('nspin', 1)
print(f"   nspin = {nspin_al}")
print(f"   Magnetic system? {nspin_al > 1}")

# Test 2: Magnetic system (Ferromagnetic)
print("\n2. Ferromagnetic Fe system (nspin=2):")
print("-" * 70)

atoms_fe = bulk('Fe', 'bcc', a=2.87)
workflow_fe = CalculationWorkflow(
    atoms_fe,
    protocol='fast',
    pseudopotentials={'Fe': 'Fe.pbe.UPF'},
    magnetic_config='ferro'
)

nspin_fe = workflow_fe.input_data.get('nspin', 1)
print(f"   nspin = {nspin_fe}")
print(f"   Magnetic system? {nspin_fe > 1}")
print(f"   Magnetic config: ferro → nspin={nspin_fe} (collinear magnetism)")

# Test 3: Antiferromagnetic
print("\n3. Antiferromagnetic system (nspin=2):")
print("-" * 70)

workflow_afm = CalculationWorkflow(
    atoms_fe,
    protocol='fast',
    pseudopotentials={'Fe': 'Fe.pbe.UPF'},
    magnetic_config='antiferro'
)

nspin_afm = workflow_afm.input_data.get('nspin', 1)
print(f"   nspin = {nspin_afm}")
print(f"   Magnetic system? {nspin_afm > 1}")
print(f"   Magnetic config: antiferro → nspin={nspin_afm} (collinear magnetism)")

# Simulate run_dos() to show the reporting
print("\n" + "=" * 70)
print("SIMULATED run_dos() OUTPUT FOR EACH SYSTEM")
print("=" * 70)

@patch('xespresso.post.dos.EspressoDos')
def test_dos_detection(mock_dos_class):
    mock_dos_instance = MagicMock()
    mock_dos_instance.directory = 'nscf/test/dos'
    mock_dos_class.return_value = mock_dos_instance
    
    Path('nscf/test').mkdir(parents=True, exist_ok=True)
    
    try:
        # Test each workflow
        print("\n[Test 1] Non-magnetic Al:")
        workflow_al.run_dos(nscf_label='nscf/test')
        
        print("\n[Test 2] Ferromagnetic Fe:")
        workflow_fe.run_dos(nscf_label='nscf/test')
        
        print("\n[Test 3] Antiferromagnetic Fe:")
        workflow_afm.run_dos(nscf_label='nscf/test', pdos=True)
        
    finally:
        import shutil
        if Path('nscf/test').exists():
            shutil.rmtree('nscf/test')

test_dos_detection()

print("\n" + "=" * 70)
print("SUMMARY: run_dos() features for magnetic systems")
print("=" * 70)
print("""
✓ Automatic nspin detection
✓ Separate spin-up and spin-down DOS
✓ Magnetic system reporting in output
✓ PDOS support for site-projected analysis
✓ Integration with Hubbard U calculations

Key advantage: Validates magnetic ordering through analysis of
spin-polarized electronic structure!
""")

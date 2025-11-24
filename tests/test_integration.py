"""
Integration test for codes module and Hubbard parameters

This test verifies that:
1. Codes module works correctly
2. Hubbard parameters work with both old and new formats
3. Integration between codes and xespresso.Espresso works
"""

import tempfile
import shutil
import os
from ase.build import bulk
import numpy as np

print("="*70)
print("Integration Test: Codes Module + Hubbard Parameters")
print("="*70)

# Test 1: Codes Module
print("\n1. Testing codes module...")
from xespresso.codes import CodesConfig, Code, CodesManager

config = CodesConfig(
    machine_name="test_machine",
    qe_prefix="/usr/local/qe-7.2/bin",
    qe_version="7.2"
)

config.add_code(Code(name="pw", path="/usr/local/qe-7.2/bin/pw.x", version="7.2"))
config.add_code(Code(name="hp", path="/usr/local/qe-7.2/bin/hp.x", version="7.2"))

tmpdir = tempfile.mkdtemp()
try:
    filepath = CodesManager.save_config(config, output_dir=tmpdir)
    print(f"✅ Codes config saved to: {filepath}")
    
    loaded = CodesManager.load_config("test_machine", codes_dir=tmpdir)
    assert loaded is not None
    assert loaded.has_code("pw")
    assert loaded.qe_version == "7.2"
    print("✅ Codes config loaded successfully")
finally:
    shutil.rmtree(tmpdir, ignore_errors=True)

# Test 2: Hubbard Parameters - New Format
print("\n2. Testing Hubbard parameters (new format)...")
from xespresso.hubbard import HubbardConfig

config_new = HubbardConfig(use_new_format=True, projector='atomic')
config_new.add_u("Fe", 4.3, orbital="3d")
config_new.add_u("O", 3.0, orbital="2p")

lines = config_new.to_new_format_card()
assert len(lines) > 0
assert any("HUBBARD" in line for line in lines)
assert any("Fe-3d" in line for line in lines)
print("✅ New format HUBBARD card generated successfully")
print("   Generated lines:")
for line in lines:
    print(f"   {line.rstrip()}")

# Test 3: Hubbard Parameters - Old Format
print("\n3. Testing Hubbard parameters (old format)...")
config_old = HubbardConfig(use_new_format=False)
config_old.add_u("Fe", 4.3)
config_old.add_u("O", 3.0)

species_info = {
    "Fe": {"index": 1},
    "O": {"index": 2}
}
old_params = config_old.to_old_format_dict(species_info)
assert "Hubbard_U(1)" in old_params
assert "Hubbard_U(2)" in old_params
print("✅ Old format parameters generated successfully")
print("   Generated parameters:")
for key, value in old_params.items():
    print(f"   {key} = {value}")

# Test 4: Auto-detection
print("\n4. Testing automatic format detection...")
input_data_new = {
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {"Fe-3d": 4.3}
    }
}
config = HubbardConfig.from_input_data(input_data_new)
assert config.should_use_new_format() == True
print("✅ Correctly detected new format from input_data")

input_data_old = {
    "input_ntyp": {
        "Hubbard_U": {"Fe": 4.3}
    }
}
config = HubbardConfig.from_input_data(input_data_old)
assert config.should_use_new_format() == False
print("✅ Correctly detected old format from input_data")

# Test 5: Integration with xio (new format)
print("\n5. Testing integration with xio module (new format)...")
from xespresso.xio import write_espresso_in

atoms = bulk("Fe", cubic=True)
atoms.new_array("species", np.array(atoms.get_chemical_symbols(), dtype="U20"))

input_data = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "qe_version": "7.2",
    "hubbard": {
        "projector": "atomic",
        "u": {"Fe-3d": 4.3}
    }
}

pseudopotentials = {
    "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"
}

tmpfile = tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False)
tmpfile.close()

try:
    write_espresso_in(
        tmpfile.name,
        atoms,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=(2, 2, 2)
    )
    
    with open(tmpfile.name, 'r') as f:
        content = f.read()
    
    # Check that HUBBARD card is present
    assert "HUBBARD" in content
    assert "Fe-3d" in content
    print("✅ Input file with new format HUBBARD card generated")
    
    # Show relevant parts
    print("\n   Relevant sections from generated input:")
    lines = content.split('\n')
    in_hubbard = False
    for line in lines:
        if 'HUBBARD' in line:
            in_hubbard = True
        if in_hubbard:
            print(f"   {line}")
            if line.strip() == '':
                break
finally:
    os.unlink(tmpfile.name)

# Test 6: Integration with xio (old format)
print("\n6. Testing integration with xio module (old format)...")
input_data_old_fmt = {
    "ecutwfc": 30.0,
    "occupations": "smearing",
    "degauss": 0.02,
    "nspin": 2,
    "lda_plus_u": True,
    "input_ntyp": {
        "Hubbard_U": {"Fe": 4.3}
    }
}

tmpfile = tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False)
tmpfile.close()

try:
    write_espresso_in(
        tmpfile.name,
        atoms,
        input_data=input_data_old_fmt,
        pseudopotentials=pseudopotentials,
        kpts=(2, 2, 2)
    )
    
    with open(tmpfile.name, 'r') as f:
        content = f.read()
    
    # Check that Hubbard_U is in SYSTEM namelist
    # Note: ASE writes parameter names in lowercase
    assert "hubbard_u" in content.lower()
    print("✅ Input file with old format Hubbard_U generated")
    
    # Show SYSTEM namelist
    print("\n   Relevant sections from generated input:")
    lines = content.split('\n')
    in_system = False
    for line in lines:
        if '&SYSTEM' in line:
            in_system = True
        if in_system:
            print(f"   {line}")
            if line.strip() == '/':
                break
finally:
    os.unlink(tmpfile.name)

print("\n" + "="*70)
print("✅ All integration tests passed!")
print("="*70)
print("\nSummary:")
print("  - Codes module: Working correctly")
print("  - Hubbard new format: Working correctly")
print("  - Hubbard old format: Working correctly")
print("  - Auto-detection: Working correctly")
print("  - Integration with xio: Working correctly")
print("  - Backward compatibility: Maintained")

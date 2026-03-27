"""
Test: Verify _get_default_orbital function works correctly
"""

from xespresso.hubbard import _get_default_orbital

print("="*80)
print("ORBITAL DEFAULT DETECTION FUNCTION TEST")
print("="*80)

test_cases = [
    # (element, expected_orbital)
    ("Fe", "3d"),      # Transition metal
    ("Mn", "3d"),      # Transition metal
    ("Co", "3d"),      # Transition metal
    ("Zr", "4d"),      # 4d transition metal
    ("Ru", "4d"),      # 4d transition metal
    ("Hf", "5d"),      # 5d transition metal
    ("Gd", "4f"),      # Lanthanide
    ("Nd", "4f"),      # Lanthanide
    ("Ce", "4f"),      # Lanthanide
    ("U", "5f"),       # Actinide
    ("Pu", "5f"),      # Actinide
    ("O", "2p"),       # Main group (fallback)
]

print("\nTesting orbital detection for various elements:")
print("-" * 80)

all_pass = True
for element, expected in test_cases:
    result = _get_default_orbital(element)
    status = "✅" if result == expected else "❌"
    if result != expected:
        all_pass = False
    print(f"{status} {element:3s} → {result:3s} (expected: {expected})")

print("\n" + "="*80)
if all_pass:
    print("✅ ALL TESTS PASSED!")
else:
    print("❌ Some tests failed")
print("="*80)

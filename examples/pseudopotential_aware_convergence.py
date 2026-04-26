#!/usr/bin/env python3
"""
Example demonstrating automatic pseudopotential-aware range adjustment.

This example shows how the ConvergenceWorkflow automatically adjusts
ecutwfc ranges based on pseudopotential requirements.
"""

import tempfile
import os
from ase.build import bulk
from xespresso.workflow import ConvergenceWorkflow

def create_mock_upf(element, suggested_ecutwfc, filename):
    """Create a mock UPF file with specified suggested ecutwfc."""
    content = f'''<PP_INFO>
Element: {element}
suggested_ecutwfc="{suggested_ecutwfc}"
</PP_INFO>'''
    with open(filename, 'w') as f:
        f.write(content)

def main():
    print("=== Pseudopotential-Aware Range Adjustment Example ===\n")

    # Create a simple silicon structure
    atoms = bulk('Si', 'diamond', a=5.43)

    # Create temporary mock pseudopotential files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Mock Si pseudopotential with high ecutwfc requirement
        si_upf = os.path.join(temp_dir, 'Si_high.UPF')
        create_mock_upf('Si', 100.0, si_upf)  # Requires 100 Ry

        pseudopotentials = {'Si': si_upf}

        print("Mock Si pseudopotential requires ecutwfc = 100 Ry")
        print("Testing different precision levels:\n")

        for precision in ['low', 'medium', 'high', 'ultra']:
            print(f"Precision '{precision}':")

            workflow = ConvergenceWorkflow(
                atoms=atoms,
                pseudopotentials=pseudopotentials,
                precision=precision
            )

            base_ranges = {
                'low': [30, 40, 50],
                'medium': [40, 50, 60, 70],
                'high': [50, 60, 70, 80, 90],
                'ultra': [60, 80, 100, 120, 140]
            }

            print(f"  Base range: {base_ranges[precision]}")
            print(f"  Adjusted range: {workflow.ecut_range}")
            print(f"  Max value: {max(workflow.ecut_range)} Ry")
            print(f"  Covers requirement: {max(workflow.ecut_range) >= 100.0}")
            print()

    print("=== Key Benefits ===")
    print("✓ Automatic detection of pseudopotential requirements")
    print("✓ Ranges extended when needed (120% of suggested value)")
    print("✓ No manual parameter tuning required")
    print("✓ Works with any pseudopotential library")
    print("✓ Maintains precision level characteristics")

    print("\n=== For Pseudopotentials Without Suggested Values ===")
    print("If a pseudopotential doesn't specify suggested_ecutwfc,")
    print("the system falls back to default precision ranges.")
    print("Users can still override with custom ranges if needed.")

if __name__ == '__main__':
    main()
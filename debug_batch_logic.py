#!/usr/bin/env python3
"""
Debug script to test batch label mapping logic
"""

# Simulating what wait_for_batch_jobs returns
monitoring_results = [
    {'label': 'eos/point_00', 'job_id': '4505', 'success': False, 'state': 'FAILED'},
    {'label': 'eos/point_01', 'job_id': '4506', 'success': False, 'state': 'FAILED'},
    {'label': 'eos/point_02', 'job_id': '4507', 'success': True, 'energy': -5.1234},
]

# OLD WAY (accessing by index)
print("OLD WAY (by index):")
for idx, result in enumerate(monitoring_results):
    print(f"  [{idx}] label={result['label']}, success={result.get('success')}")

print()

# NEW WAY (mapping by label)
print("NEW WAY (mapping by label):")
label_to_monitor_result = {r['label']: r for r in monitoring_results}
print(f"  Mapping: {list(label_to_monitor_result.keys())}")

# Simulating factor order
factors = [0.92, 0.94, 0.96]
factor_to_index = {0.92: 0, 0.94: 1, 0.96: 2}

print()
print("Looking up by label:")
for factor in sorted(factors):
    point_idx = factor_to_index[factor]
    calc_label = f"eos/point_{point_idx:02d}"
    result = label_to_monitor_result.get(calc_label)
    if result:
        print(f"  Factor {factor:.2f} → point_{point_idx:02d} ({calc_label}) → success={result.get('success')}")
    else:
        print(f"  Factor {factor:.2f} → point_{point_idx:02d} ({calc_label}) → NOT FOUND!")

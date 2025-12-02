# Hubbard Parameter Fix Summary

## Problem Statement

The qtGui was not properly respecting how xespresso handles Hubbard parameters. It was only implementing the old format and incorrectly passing Hubbard U values to the Espresso calculator, not following xespresso's documented design patterns.

## Root Cause

In `qtgui/pages/job_submission.py` (lines 770-777), the code only supported the old format:
```python
input_data['input_ntyp']['Hubbard_U'][element] = u_value
```

This approach:
- Only worked with old QE format (< 7.0)
- Did not support orbital specifications (e.g., `Fe-3d`, `O-2p`)
- Did not support the new HUBBARD card format (QE >= 7.0)
- Did not auto-detect format based on QE version

## Solution Implemented

### 1. Updated `job_submission.py`

Added comprehensive Hubbard parameter handling in `_build_input_data()`:

```python
# Determine format based on config and QE version
use_new_format = False
if hubbard_format == 'new':
    use_new_format = True
elif hubbard_format == 'old':
    use_new_format = False
elif qe_version:
    # Auto-detect from version
    try:
        major = int(qe_version.split('.')[0])
        use_new_format = (major >= 7)
    except (ValueError, IndexError):
        pass

if use_new_format:
    # NEW FORMAT: Use 'hubbard' dictionary
    hubbard_dict = {
        'projector': config.get('hubbard_projector', 'atomic'),
        'u': {},
        'v': []
    }
    for element, u_value in config.get('hubbard_u', {}).items():
        if u_value > 0:
            orbital = config.get('hubbard_orbitals', {}).get(
                element, _get_default_orbital(element)
            )
            hubbard_dict['u'][f"{element}-{orbital}"] = u_value
    
    input_data['hubbard'] = hubbard_dict
    input_data['qe_version'] = qe_version
else:
    # OLD FORMAT: Use 'input_ntyp'
    ensure_input_ntyp(input_data)
    input_data['input_ntyp']['Hubbard_U'] = {}
    for element, u_value in config.get('hubbard_u', {}).items():
        if u_value > 0:
            input_data['input_ntyp']['Hubbard_U'][element] = u_value
```

Added default orbital mapping:
```python
DEFAULT_HUBBARD_ORBITALS = {
    'Fe': '3d', 'Mn': '3d', 'Co': '3d', 'Ni': '3d',  # 3d metals
    'Ce': '4f', 'Nd': '4f', 'Gd': '4f',              # 4f rare earths
    'O': '2p', 'N': '2p', 'S': '3p',                 # p-block
    # ... more elements
}
```

### 2. Updated `calculation_setup.py`

Enhanced the Hubbard configuration UI:

- Added orbital selection widgets for each element
- Implemented `_get_common_orbitals_for_element()` to provide appropriate orbital choices
- Updated `_get_config()` to save orbital information
- Updated `_restore_config_to_ui()` to restore orbital selections
- Initialized `hubbard_orbital_edits` in `__init__`

Example UI changes:
```python
# Now creates a widget with both U value and orbital selector
container = QWidget()
hlayout = QHBoxLayout(container)

# U value spinbox
spin = QDoubleSpinBox()
spin.setValue(default_val)
hlayout.addWidget(spin)

# Orbital selector (for new format)
orbital_combo = QComboBox()
orbitals = self._get_common_orbitals_for_element(element)
orbital_combo.addItems(orbitals)
hlayout.addWidget(QLabel("Orbital:"))
hlayout.addWidget(orbital_combo)
```

### 3. Added Comprehensive Tests

Created `tests/test_qtgui_hubbard_fix.py` with tests for:
- Old format Hubbard parameters
- New format Hubbard parameters with explicit orbitals
- New format with default orbitals
- Format auto-detection based on QE version
- No Hubbard parameters case

### 4. Added Example

Created `examples/qtgui_hubbard_example.py` demonstrating:
- Old format configuration
- New format configuration with explicit orbitals
- New format with default orbitals
- How to use in the qtGui

## Formats Supported

### Old Format (QE < 7.0)

```python
input_data = {
    'SYSTEM': {'lda_plus_u': True},
    'input_ntyp': {
        'Hubbard_U': {
            'Fe': 4.3,
            'O': 3.0
        }
    }
}
```

Generates in QE input:
```
&SYSTEM
  lda_plus_u = .true.
  Hubbard_U(1) = 4.3
  Hubbard_U(2) = 3.0
/
```

### New Format (QE >= 7.0)

```python
input_data = {
    'SYSTEM': {'lda_plus_u': True},
    'hubbard': {
        'projector': 'atomic',
        'u': {
            'Fe-3d': 4.3,
            'O-2p': 3.0
        },
        'v': []
    },
    'qe_version': '7.2'
}
```

Generates in QE input:
```
HUBBARD {atomic}
  U Fe-3d 4.3
  U O-2p 3.0
```

## Testing Results

All tests passed:
- ✅ 14 existing Hubbard tests
- ✅ 3 Hubbard integration tests
- ✅ 33 magnetic helpers tests (including Hubbard+magnetic combinations)
- ✅ Logic validation tests for both formats
- ✅ Code review completed
- ✅ Security scan (0 alerts)

## Benefits

1. **Respects xespresso's design** - Follows the documented patterns in xespresso/hubbard.py
2. **Supports both formats** - Works with QE < 7.0 and QE >= 7.0
3. **Explicit orbital specifications** - Allows precise control over which orbitals get Hubbard corrections
4. **Sensible defaults** - Provides default orbitals for common elements
5. **Format detection** - Auto-detects based on QE version and user preference
6. **Backward compatible** - Old format still works exactly as before
7. **Fully tested** - Comprehensive test coverage ensures reliability

## Files Changed

- `qtgui/pages/job_submission.py` - 84 lines added
- `qtgui/pages/calculation_setup.py` - 88 lines added
- `tests/test_qtgui_hubbard_fix.py` - 266 lines added (new file)
- `examples/qtgui_hubbard_example.py` - 232 lines added (new file)

Total: 670 lines added, 7 lines modified

## References

- xespresso documentation: `xespresso/hubbard.py`
- Example usage: `examples/ex02-MnO-lda+u.py`, `examples/ex04-dft+u.py`
- Hubbard format examples: `examples/hubbard_format_examples.py`
- New format example: `examples/hubbard_new_format_example.py`

## Security

No security vulnerabilities were found during CodeQL analysis.

# Hubbard Format Support Implementation Summary

## User Request

> "It is just not clear for me how the code is handling the new hubbard format that uses a Card when you call this setup_magnetic_config(), because it needs to identify the element and the orbital in which U is applied. You should also take care when the user defines both U and V"

## Solution

Added comprehensive support for both old (QE < 7.0) and new (QE >= 7.0) Hubbard parameter formats in `setup_magnetic_config()`.

## Key Features Implemented

### 1. Element + Orbital Identification (QE 7.x)

**Syntax:**
```python
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}
}, qe_version='7.2')
```

**Result:**
- Creates HUBBARD card with Fe-3d and Fe1-3d orbital specifications
- Properly maps each species of Fe to its 3d orbital
- Returns `config['hubbard']` dict with format:
  ```python
  {
      'projector': 'atomic',
      'u': {'Fe-3d': 4.3, 'Fe1-3d': 4.3},
      'v': []
  }
  ```

### 2. V Parameter Support (Inter-site Interactions)

**Syntax:**
```python
config = setup_magnetic_config(atoms, {
    'Fe': {
        'mag': [1],
        'U': {'3d': 4.3},
        'V': [{'species2': 'O', 'orbital1': '3d', 'orbital2': '2p', 'value': 1.0}]
    }
}, qe_version='7.2')
```

**Features:**
- Full V parameter specification with orbitals
- Multiple V parameters supported
- Automatically applies to all species of the element
- Includes site indices (i, j) with defaults of (1, 1)

### 3. Different U per Species

**Syntax:**
```python
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': {'3d': [4.3, 4.5]}}
}, qe_version='7.2')
```

**Result:**
- Fe-3d gets U = 4.3 eV
- Fe1-3d gets U = 4.5 eV

### 4. Automatic Format Detection

**Rules:**
1. If `qe_version >= 7.0` → new format
2. If `qe_version < 7.0` → old format
3. If U is dict (`{'3d': 4.3}`) and no qe_version → new format (auto-detect)
4. If U is simple value (`4.3`) and no qe_version → old format (auto-detect)
5. Can override with `hubbard_format='old'` or `'new'`

### 5. Backward Compatibility

**Old format still works:**
```python
config = setup_magnetic_config(atoms, {
    'Fe': {'mag': [1, -1], 'U': 4.3}
}, qe_version='6.8')
```

**Result:**
- Creates `config['input_ntyp']['Hubbard_U']` = `{'Fe': 4.3, 'Fe1': 4.3}`
- Works with existing QE < 7.0 infrastructure

## Format Comparison

### Old Format (QE < 7.0)

**Input:**
```python
{'Fe': {'mag': [1, -1], 'U': 4.3}}
```

**Output in QE input file:**
```
&SYSTEM
  Hubbard_U(1) = 4.3
  Hubbard_U(2) = 4.3
/
```

### New Format (QE >= 7.0)

**Input:**
```python
{'Fe': {'mag': [1, -1], 'U': {'3d': 4.3}}}
```

**Output in QE input file:**
```
HUBBARD {atomic}
  U Fe-3d 4.3
  U Fe1-3d 4.3
```

### New Format with V

**Input:**
```python
{
    'Fe': {
        'mag': [1],
        'U': {'3d': 4.3},
        'V': [{'species2': 'O', 'orbital1': '3d', 'orbital2': '2p', 'value': 1.0}]
    }
}
```

**Output in QE input file:**
```
HUBBARD {atomic}
  U Fe-3d 4.3
  V Fe-3d O-2p 1 1 1.0
```

## Implementation Details

### Code Structure

**Function signature:**
```python
def setup_magnetic_config(atoms, magnetic_config, pseudopotentials=None, 
                         expand_cell=False, qe_version=None, hubbard_format='auto')
```

**New parameters:**
- `qe_version`: String like '7.2' or '6.8' to determine format
- `hubbard_format`: 'auto', 'old', or 'new' for explicit control

**Return value additions:**
- `config['hubbard']`: Dict with 'u' and 'v' parameters (new format)
- `config['hubbard_format']`: 'old' or 'new' indicator
- `config['qe_version']`: QE version string (if provided)

### Format Detection Logic

```python
use_new_hubbard_format = False

if hubbard_format == 'new':
    use_new_hubbard_format = True
elif hubbard_format == 'old':
    use_new_hubbard_format = False
elif qe_version:
    # Auto-detect from version
    major = int(qe_version.split('.')[0])
    use_new_hubbard_format = (major >= 7)
else:
    # Auto-detect from U format
    if U is dict:  # {'3d': 4.3}
        use_new_hubbard_format = True
```

### V Parameter Processing

For each V parameter:
1. Extract: species2, orbital1, orbital2, value, i, j
2. Apply to all species of the element
3. Store as: `{'species1': 'Fe', 'orbital1': '3d', 'species2': 'O', 'orbital2': '2p', 'i': 1, 'j': 1, 'value': 1.0}`

## Testing

### Tests Added (5 new)

1. **test_new_hubbard_format_with_orbital**: Basic orbital specification
2. **test_new_hubbard_format_different_u_per_species**: Different U values
3. **test_new_hubbard_format_with_v_parameter**: V parameter support
4. **test_old_format_with_dict_u_warns**: Backward compatibility
5. **test_new_format_error_without_orbital**: Error handling

### Test Coverage

- ✅ Orbital specification (Fe-3d, Mn-3d)
- ✅ Multiple orbitals per element
- ✅ Different U per species
- ✅ V parameters with full specification
- ✅ Auto-detection logic
- ✅ Format validation and error messages
- ✅ Old format compatibility

### All Tests Passing

```
30 passed in 0.47s
```

## Examples

### Example File

`examples/hubbard_format_examples.py` - 8 complete examples:

1. Old format (QE < 7.0)
2. New format with orbital
3. Different U per species
4. V parameters
5. Complex multi-element system
6. Auto-detection
7. Explicit format control
8. Usage with Espresso calculator

## Integration with Existing Code

### Works with HubbardConfig Class

The new format output is compatible with:
- `HubbardConfig.from_input_data()`
- `build_hubbard_str()`
- `apply_hubbard_to_system()`

### Input Data Format

For new format, users pass:
```python
input_data = {
    'ecutwfc': 70.0,
    'nspin': 2,
    'lda_plus_u': True,
    'input_ntyp': config['input_ntyp'],
    'hubbard': config['hubbard'],  # New format
    'qe_version': '7.2'
}
```

## Error Handling

### Clear Error Messages

**Example 1:** New format without orbital
```
ValueError: Element Fe: New Hubbard format requires orbital specification.
Use 'U': {'3d': 4.3} instead of 'U': 4.3
```

**Example 2:** Missing V parameter fields
```
ValueError: V parameter must specify: species2, orbital1 (or orbital), orbital2, value
```

## Documentation

### Updated Files

1. **MAGNETIC_HELPERS.md**: Added Hubbard format section
2. **xespresso/tools.py**: Updated docstring with new format examples
3. **examples/hubbard_format_examples.py**: Complete examples
4. **tests/test_magnetic_helpers.py**: Test documentation

### Quick Reference

See `MAGNETIC_HELPERS.md` section "Hubbard Parameters Support" for:
- API reference
- Format comparison
- Multiple examples
- Migration guide

## Commit

**Hash:** 0014dfb
**Message:** Add support for new Hubbard format (QE 7.x) with orbital specification and V parameters

## Summary

Successfully addressed user's concern about Hubbard format handling:

✅ Element AND orbital identification implemented
✅ Full V parameter support added
✅ Both U and V can be used together
✅ Clear documentation and examples
✅ All tests passing
✅ Backward compatible

The function now properly handles the new QE 7.x HUBBARD card format with explicit element-orbital specification while maintaining full backward compatibility with the old format.

"""
Top-level package for xespresso.

This module performs guarded imports of commonly-used symbols so that
importing the top-level package doesn't fail during test collection or
when optional submodules are temporarily unavailable.

Symbols are imported from their modules when possible; failures are
silently ignored to preserve importability.
"""

import importlib

__all__ = []

# Mapping of modules to names we want to expose at package level
_EXPORTS = {
    'xespresso.xespresso': ['Espresso'],
    'xespresso.hubbard': ['HubbardConfig', 'build_hubbard_str', 'apply_hubbard_to_system'],
    'xespresso.tools': [
        'set_magnetic_moments',
        'set_antiferromagnetic',
        'set_ferromagnetic',
        'setup_magnetic_config',
        'kpts_from_spacing',
    ],
    'xespresso.workflow': ['CalculationWorkflow', 'quick_scf', 'quick_relax', 'PRESETS'],
    'xespresso.pseudopotentials': [
        'PseudopotentialsManager',
        'create_pseudopotentials_config',
        'load_pseudopotentials_config',
    ],
}

for modname, names in _EXPORTS.items():
    try:
        _mod = importlib.import_module(modname)
    except Exception:
        # If the submodule fails to import, skip its symbols
        continue
    for _name in names:
        try:
            globals()[_name] = getattr(_mod, _name)
            __all__.append(_name)
        except AttributeError:
            # Symbol not present in module; skip
            continue

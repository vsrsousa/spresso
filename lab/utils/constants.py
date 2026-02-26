"""
Constants for XESPRESSO Lab UI.
"""

# Page titles and descriptions
PAGE_TITLES = {
    'structures': 'Structures Explorer',
    'calculations': 'Calculations Timeline',
    'lineage': 'Derivation Lineage',
    'search': 'Advanced Search',
}

PAGE_DESCRIPTIONS = {
    'structures': 'Browse all structures in the database',
    'calculations': 'View calculation history and execution details',
    'lineage': 'Explore derivation relationships between structures',
    'search': 'Search and filter by multiple criteria',
}

# Calculation methods
CALCULATION_METHODS = {
    'scf': 'Self-Consistent Field',
    'relax': 'Ionic Relaxation',
    'vc-relax': 'Cell + Ionic Relaxation',
    'phonon': 'Phonon Calculation',
    'band': 'Band Structure',
    'dos': 'Density of States',
}

# Status symbols and text
STATUS_TEXT = {
    True: 'Success',
    False: 'Failed',
}

# Default values
DEFAULT_ENERGY_MIN = None
DEFAULT_ENERGY_MAX = None
DEFAULT_ITEMS_PER_PAGE = 20

"""Helper to load default element color palette from JSON resource.

Provides DEFAULT_ELEMENT_COLORS dict and a loader function.
"""
import json
import os

_THIS_DIR = os.path.dirname(__file__)
_JSON_PATH = os.path.join(_THIS_DIR, 'element_colors.json')

try:
    with open(_JSON_PATH, 'r', encoding='utf-8') as f:
        DEFAULT_ELEMENT_COLORS = json.load(f)
except Exception:
    # Fallback minimal mapping
    DEFAULT_ELEMENT_COLORS = {
        'H': '#FFFFFF', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D',
        'S': '#FFFF30', 'P': '#FF8000', 'Fe': '#E06633'
    }

def get_default_element_colors():
    return dict(DEFAULT_ELEMENT_COLORS)

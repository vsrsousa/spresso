"""
Configuration for XESPRESSO Lab.
"""

from pathlib import Path
import os

# Database paths
DATABASES_DIR = Path.home() / '.xespresso'
STRUCTURES_DB = DATABASES_DIR / 'structures.db'
PROVENANCE_DB = DATABASES_DIR / 'provenance.db'

# Ensure directory exists
DATABASES_DIR.mkdir(parents=True, exist_ok=True)

# Streamlit configuration
STREAMLIT_CONFIG = {
    'theme': {
        'primaryColor': '#4A90E2',      # Neutral blue
        'backgroundColor': '#FFFFFF',
        'secondaryBackgroundColor': '#F5F5F5',
        'textColor': '#1F1F1F',
        'font': 'sans serif',
    },
    'layout': 'wide',
}

# Cache settings
CACHE_TTL = 3600  # 1 hour

# Pagination
ITEMS_PER_PAGE = 20

# Plot colors (neutral palette)
COLORS = {
    'primary': '#4A90E2',
    'success': '#5FB878',
    'warning': '#F5A623',
    'error': '#D0021B',
    'info': '#50E3C2',
    'background_light': '#F5F5F5',
    'background_dark': '#EEEEEE',
}

# Status mapping
STATUS_COLORS = {
    True: COLORS['success'],
    False: COLORS['error'],
}

STATUS_SYMBOLS = {
    True: '✅',
    False: '❌',
}

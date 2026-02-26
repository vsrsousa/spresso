# XESPRESSO Lab

Interactive Provenance Explorer for XESPRESSO calculations.

A Streamlit-based interface for exploring calculation history, derivations, structures, and execution details stored in the XESPRESSO provenance database.

## Features

- **Structures Explorer**: Browse atomic structures, view properties, compare geometries
- **Calculations Timeline**: View calculation history with execution details
- **Derivation Lineage**: Explore relationships between structures and calculations
- **Advanced Search**: Search across structures and calculations with multiple criteria
- **Execution History**: Track where each calculation was executed

## Installation

Install dependencies:

```bash
pip install -r lab/requirements.txt
```

## Usage

### Run the interface

```bash
cd /home/vinicius/projects/spresso
python -m streamlit run lab/app.py
```

The interface will open at `http://localhost:8501`

### Database Locations

The interface expects databases in:

- `~/.xespresso/structures.db` - ASE database with atomic structures
- `~/.xespresso/provenance.db` - Provenance database with calculation metadata

These are created automatically by `DatabaseWorkflow`.

## Architecture

```
lab/
├── app.py                    # Main Streamlit application
├── config.py                 # Configuration (db paths, colors, cache)
├── database_interface.py     # Database queries and caching
├── pages/                    # Streamlit pages
│   ├── 01_structures.py     # Structures explorer
│   ├── 02_calculations.py   # Calculations timeline
│   ├── 03_lineage.py        # Derivation lineage
│   └── 04_search.py         # Advanced search
├── components/               # Reusable UI components
│   ├── structure_viewer.py  # Structure visualization
│   ├── derivation_graph.py  # Derivation graphs
│   ├── metrics.py           # Metric cards
│   └── filters.py           # Filter widgets
└── utils/                    # Utilities
    ├── formatters.py        # Data formatting
    ├── constants.py         # UI constants
    └── cache.py             # Caching utilities
```

## Quick Reference

### View all structures
1. Go to "Structures" page
2. Click "Search" to display all structures
3. Click on a structure ID for detailed view

### Explore derivations
1. Go to "Lineage" page
2. Select a source structure
3. View the derivation tree and energy evolution
4. Compare structures side-by-side

### Track calculations
1. Go to "Calculations" page
2. View timeline of all calculations
3. Click on a calculation to see execution history
4. Check which machine executed each calculation

### Custom search
1. Go to "Search" page
2. Configure criteria (formula, energy range, calculation type)
3. View results with detailed information

## Customization

### Change color palette

Edit `lab/config.py`:

```python
COLORS = {
    'primary': '#4A90E2',      # Change these values
    'success': '#5FB878',
    'warning': '#F5A623',
    'error': '#D0021B',
    ...
}
```

### Modify cache TTL

Edit `lab/config.py`:

```python
CACHE_TTL = 3600  # seconds
```

### Add custom pages

Create a new file in `lab/pages/` and import it in `lab/app.py`.

## Notes

- Neutral color palette (blue, green, orange, red) for professional appearance
- Responsive design with Streamlit columns
- Automatic caching of database queries
- Query results cached for 1 hour by default
- Zero Node.js dependencies - pure Python

## Future Enhancements

- Real-time calculation monitoring
- Export results (CSV, JSON)
- Custom visualization plugins
- Performance profiling
- Batch operations

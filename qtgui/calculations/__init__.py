"""
Calculations module for xespresso GUI.

This module provides calculation preparation functionality for the GUI,
encapsulating the creation of atoms and Espresso calculator objects
from GUI configuration.

The key principle: calculation modules prepare objects, job submission executes them.
"""

from qtgui.calculations.base import BaseCalculationPreparation
from qtgui.calculations.preparation import (
    CalculationPreparation,
    prepare_calculation_from_gui,
    dry_run_calculation
)

__all__ = [
    'BaseCalculationPreparation',
    'CalculationPreparation',
    'prepare_calculation_from_gui',
    'dry_run_calculation',
]

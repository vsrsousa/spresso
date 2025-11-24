"""
machines package initializer.

Exposes top-level config and resource utilities.
"""

from .config import load_machine, create_machine, edit_machine, migrate_machines, rollback_migration
from .machine import Machine
#from .resources import set_partition, set_serial_partition, set_gpu, set_debug_mode

# xespresso/utils/slurm.py

import shutil
import subprocess
import os
from xespresso.config import VERBOSE_ERRORS

def check_slurm_available():
    """
    Validates that the SLURM job scheduler is installed and operational on the system.

    This function performs two checks:
    1. Verifies that the 'sbatch' command is available in the system's PATH.
    2. Confirms that the SLURM controller daemon ('slurmctld') is responsive via 'scontrol ping'.

    Behavior:
    - If the environment variable XESPRESSO_FORCE_SCHEDULER is set to '1', all checks are skipped.
    - If VERBOSE_ERRORS is True (via XESPRESSO_VERBOSE_ERRORS), full Python tracebacks will be shown.
      Otherwise, errors are raised cleanly without traceback clutter.

    Raises:
        RuntimeError: If SLURM is not installed or not responding, unless overridden by environment.

    Example:
        >>> check_slurm_available()
        # Raises RuntimeError with a clean message if SLURM is unavailable

    Environment Variables:
        XESPRESSO_FORCE_SCHEDULER: Set to '1' to bypass SLURM checks entirely.
        XESPRESSO_VERBOSE_ERRORS: Set to '1', 'true', or 'yes' to enable full tracebacks on failure.
    """
    if os.getenv("XESPRESSO_FORCE_SCHEDULER") == "1":
        return

    if shutil.which("sbatch") is None:
        msg = (
            "SLURM scheduler requested but 'sbatch' command not found.\n"
            "Please install SLURM or use 'direct' as the scheduler for local execution."
        )
        raise RuntimeError(msg) if VERBOSE_ERRORS else RuntimeError(msg) from None

    try:
        subprocess.run(
            ["scontrol", "ping"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    except subprocess.CalledProcessError:
        msg = (
            "SLURM is installed but the controller is not responding.\n"
            "Make sure slurmctld is running and accessible."
        )
        raise RuntimeError(msg) if VERBOSE_ERRORS else RuntimeError(msg) from None

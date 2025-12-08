import os
import logging
from xespresso.schedulers.factory import get_scheduler
from xespresso.config import VERBOSE_ERRORS
from xespresso.utils.slurm import check_slurm_available

def set_queue(calc, package=None, parallel=None, queue=None, command=None):
    """
    Configures the calculator for job submission in xespresso.

    This function prepares the job execution environment by:
        - Replacing placeholders in the command string using calc attributes.
        - Defaulting to local execution with the 'direct' scheduler if none is specified.
        - Initializing the appropriate scheduler via get_scheduler().
        - Writing the job script using scheduler.write_script().
        - Executing the job using scheduler.run(), which supports local or remote execution.
        - Storing the final command string in calc.command and optionally in calc.profile.

    SLURM availability is validated only for local jobs using the SLURM scheduler.

    Args:
        calc: ASE calculator object.
        package (str, optional): Executable name (e.g., 'pw').
        parallel (str, optional): Parallel execution command (e.g., 'mpirun -np 4').
        queue (dict or Machine, optional): Scheduler configuration dictionary or Machine object.
        command (str, optional): Command template string with placeholders.

    Raises:
        RuntimeError: If SLURM is selected for local execution but not available.
        ValueError: If scheduler initialization fails.
    """
    logger = logging.getLogger(__name__)
    queue = queue or calc.queue
    
    # Convert Machine object to queue dict if necessary
    from xespresso.machines.machine import Machine
    if isinstance(queue, Machine):
        queue = queue.to_queue()
    
    calc.queue = queue
    package = package or calc.package
    parallel = parallel or calc.parallel
    command = command or os.environ.get("ASE_ESPRESSO_COMMAND", "")

    # Replace placeholders
    if "LAUNCHER" in command:
        command = command.replace("LAUNCHER", queue.get("launcher", ""))
    if "PACKAGE" in command:
        command = command.replace("PACKAGE", package)
    if "PREFIX" in command:
        command = command.replace("PREFIX", calc.prefix)
    if "PARALLEL" in command:
        command = command.replace("PARALLEL", parallel)

    logger.debug(f"Espresso command: {command}")

    # Default to local direct scheduler if none is defined
    if not queue or "scheduler" not in queue:
        queue = {
            "execution": "local",
            "scheduler": "direct"
        }
        calc.queue = queue
        logger.debug("No scheduler defined. Defaulting to local direct execution.")

    # Validate SLURM only for local execution
    if queue.get("scheduler") == "slurm" and queue.get("execution", "local") == "local":
        try:
            check_slurm_available()
        except RuntimeError as e:
            msg = f"SLURM scheduler is not available locally: {e}"
            raise RuntimeError(msg) if VERBOSE_ERRORS else RuntimeError(msg) from None

    # Initialize scheduler
    try:
        scheduler = get_scheduler(calc, queue, command)
    except Exception as e:
        msg = f"Failed to initialize scheduler '{queue.get('scheduler')}': {e}"
        raise ValueError(msg) if VERBOSE_ERRORS else ValueError(msg) from None

    # Write job script only — defer execution to scheduler
    scheduler.write_script()

    # Store scheduler for later use
    calc.scheduler = scheduler
    
    # Store the command - it will be used by scheduler.run() for both local and remote
    # For remote execution: used in the remote job script
    # For local execution: used by subprocess in base Scheduler.run()
    calc.command = scheduler.submit_command()
    
    if hasattr(calc, "profile"):
        calc.profile.command = calc.command
    
    logger.debug(f"Scheduler configured for {queue.get('execution', 'local')} execution with command: {calc.command}")

from .slurm import SlurmScheduler
from .direct import DirectScheduler

def get_scheduler(calc, queue, command):
    """
    Factory function that returns the appropriate Scheduler instance.

    This function inspects the 'scheduler' key in the queue dictionary and
    instantiates the corresponding scheduler class. It abstracts away the
    selection logic, allowing the main interface to remain clean and modular.

    Supported schedulers:
        - "slurm": Uses SlurmScheduler (submits via sbatch)
        - "direct": Uses DirectScheduler (runs via bash script)

    Args:
        calc (Calculator): The calculator instance to be scheduled.
        queue (dict): Dictionary containing scheduler configuration.
        command (str): The command to execute the job.

    Returns:
        Scheduler: An instance of the appropriate subclass of Scheduler.

    Raises:
        ValueError: If the scheduler type is unsupported or missing.
    """
    scheduler_type = queue.get("scheduler")
    
    if not scheduler_type:
        raise ValueError("Scheduler type must be specified in queue configuration")
    
    scheduler_type = scheduler_type.lower()

    if scheduler_type == "slurm":
        return SlurmScheduler(calc, queue, command)
    elif scheduler_type == "direct":
        return DirectScheduler(calc, queue, command)
    else:
        raise ValueError(f"Unsupported scheduler: {scheduler_type}")

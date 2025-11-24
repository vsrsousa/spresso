from .base import Scheduler
from .remote_mixin import RemoteExecutionMixin
import os

class SlurmScheduler(RemoteExecutionMixin, Scheduler):
    """
    SLURM-compatible job scheduler with remote execution support.

    Inherits from:
        - Scheduler: Abstract base class for job script generation and submission.
        - RemoteExecutionMixin: Adds remote SSH execution, file transfer, and pseudopotential syncing.

    Features:
        - Generates SLURM job scripts using queue parameters.
        - Supports remote execution via SSH.
        - Loads environment modules and shell setup from `config_script`.
        - Appends cleanup or post-processing commands via `post_script`.
        - Submits jobs using `sbatch`.

    Attributes:
        calc (Calculator): The calculator instance (e.g., Espresso).
        queue (dict): Scheduler and remote execution configuration.
        command (str): The command to run (e.g., 'pw.x -in scf.pwi').
        job_file (str): Name of the job script file.
        script_dir (str): Directory where the job script is saved.
        config_script (str): Pre-execution shell setup block.
        post_script (str): Post-execution shell block.
    """

    def write_script(self):
        """
        Generates a SLURM job script and writes it to disk.

        The script includes:
            - SBATCH directives from queue["resources"]
            - Environment setup via config_script
            - Execution command
            - Optional post-processing via post_script
        """
        lines = ["#!/bin/bash", ""]

        # SLURM directives
        sbatch_keys = {
            "job-name": self.queue.get("job_name", self.calc.prefix),
            "output": self.queue.get("output", f"{self.calc.prefix}.out"),
            "error": self.queue.get("error", f"{self.calc.prefix}.err")
        }

        for key, value in sbatch_keys.items():
            lines.append(f"#SBATCH --{key}={value}")

        # Additional SBATCH options
        for key, value in self.queue.get("resources", {}).items():
            lines.append(f"#SBATCH --{key}={value}")

        lines.append("")  # Blank line after SBATCH block

        # Environment setup
        if self.config_script:
            lines.append(self.config_script)
            lines.append("")

        # Main execution command
        lines.append(self.command)
        lines.append("")

        # Post-execution block
        if self.post_script:
            lines.append(self.post_script)

        # Write to job file
        with open(os.path.join(self.script_dir, self.job_file), "w") as f:
            f.write("\n".join(lines))

    def submit_command(self):
        """
        Returns the SLURM submission command.

        Returns:
            str: Command to submit the job script via sbatch.
        """
        return f"sbatch {self.job_file}"

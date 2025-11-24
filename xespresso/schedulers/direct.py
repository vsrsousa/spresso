import os
from .base import Scheduler
from .remote_mixin import RemoteExecutionMixin

class DirectScheduler(RemoteExecutionMixin, Scheduler):
    """
    Scheduler for direct bash execution (no job manager), supporting both
    local and remote execution via RemoteExecutionMixin.

    Attributes:
        calc (Calculator): The calculator instance.
        queue (dict): Scheduler configuration parameters.
        command (str): Execution command with placeholders resolved.
        job_file (str): Name of the job script file.
        script_dir (str): Directory where the job script is written.
        config_script (str): Optional environment setup script content.
        post_script (str): Optional post-execution commands.

    Methods:
        write_script(): Writes a bash script (.job_file) with structured layout.
        submit_command(): Returns the bash execution command.
        run(): Executes the job locally or remotely (via mixin).
    """

    def write_script(self):
        lines = ["#!/bin/bash", ""]

        # Prepend block (includes modules and config)
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
        return f"bash {self.job_file}"

import os
import subprocess

class Scheduler:
    """
    Abstract base class for job schedulers used in xespresso.

    This class defines the interface and shared logic for all scheduler types
    (e.g., SLURM, Direct). It handles job script generation, submission, and
    environment setup via the `queue` dictionary.

    Supported keys in `queue`:
        - "prepend": Shell commands to run before the main job (e.g. source env.sh)
        - "modules": List of modules to load (e.g. ["quantum-espresso"])
        - "use_modules": Boolean flag to enable module loading
        - "postpend": Shell commands to run after the main job
        - "xespressorc": Optional shell config file in $HOME
        - "execution": Either "local" or "remote" (default: "local")

    Attributes:
        calc (Calculator): The calculator instance associated with the job.
        queue (dict): Dictionary containing scheduler parameters and metadata.
        command (str): The command to execute the job.
        job_file (str): Name of the job script file to be written.
        script_dir (str): Directory where the job script will be saved.
        config_script (str): Pre-execution environment setup commands.
        post_script (str): Post-execution cleanup or final commands.

    Methods:
        _load_config_script(): Builds the environment setup block.
        _load_post_script(): Loads post-execution commands.
        write_script(): Abstract method to generate the job script.
        submit_command(): Abstract method to return the job submission command.
        run(): Executes the job locally (default). Can be overridden for remote execution.
    """

    def __init__(self, calc, queue, command):
        self.calc = calc
        self.queue = queue
        self.command = command
        self.job_file = "job_file"
        self.script_dir = calc.directory
        self.config_script = self._load_config_script()
        self.post_script = self._load_post_script()

    def _load_config_script(self):
        """
        Builds the environment setup block to be placed before the main execution line.

        Priority:
            1. queue["prepend"]: Custom shell commands
            2. queue["use_modules"] + queue["modules"]: Load modules
            3. queue["xespressorc"]: Optional shell config file in $HOME

        Returns:
            str: Combined shell commands for environment setup.
        """
        lines = []

        if "prepend" in self.queue:
            lines.append(self.queue["prepend"])
        lines.append("")

        if self.queue.get("use_modules", False):
            modules = self.queue.get("modules", [])
            if isinstance(modules, list) and modules:
                lines.append("module purge")
                for mod in modules:
                    lines.append(f"module load {mod}")
        lines.append("")

        config_name = self.queue.get("xespressorc")
        if config_name:
            home_path = os.path.join(os.environ.get("HOME", ""), config_name)
            if os.path.exists(home_path):
                with open(home_path, "r") as f:
                    lines.append(f.read())

        return "\n".join(lines)

    def _load_post_script(self):
        """
        Loads post-execution commands to be placed after the main execution line.

        Returns:
            str: Shell commands for cleanup or final steps.
        """
        return self.queue.get("postpend", "")

    def write_script(self):
        """
        Abstract method to write the job script.

        Subclasses must implement this method to generate the appropriate job file
        based on the scheduler type and queue parameters.
        """
        raise NotImplementedError

    def submit_command(self):
        """
        Abstract method to return the job submission command.

        Subclasses must implement this method to return the correct command
        for submitting the job script (e.g., 'sbatch', 'bash').
        """
        raise NotImplementedError

    def run(self):
        """
        Executes the job locally using os.system.

        This method can be overridden by subclasses to support remote execution
        via SSH or other mechanisms.

        Returns:
            tuple: (stdout, stderr) if applicable, else (None, None)
        """
#        os.system(self.submit_command())
        subprocess.run(self.submit_command(), shell=True, cwd=self.script_dir, check=True)
        return None, None

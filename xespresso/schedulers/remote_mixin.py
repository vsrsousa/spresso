import os
import hashlib
import time
from xespresso.utils.auth import RemoteAuth
from xespresso.utils import warnings as warnings  # Custom warning system

# Apply custom formatting globally (if your module supports it)
warnings.formatwarning = lambda msg, cat, fname, lineno, *_: (
    f"\n⚠️ {cat.__name__} in {fname}:{lineno}\n→ {msg}\n"
)


class RemoteExecutionMixin:
    """
    Mixin class that adds remote execution capabilities to any Scheduler.

    Features:
    - Reuses SSH connection across multiple calculations on the same server
    - Automatically opens a new connection if the server or user changes
    - Dynamically computes unique remote working directory based on calc.directory
    - Uses hash of local path to ensure uniqueness across different local directories
    - Avoids redundant remote_path setup if calc.directory hasn't changed
    - Transfers required pseudopotentials to remote ./pseudo directory
    - Logs file transfers and warnings
    - Validates file integrity using SHA256 checksums

    Assumes:
    - self.calc: an Espresso calculator with .prefix, .package, and .directory
    - self.queue: dictionary with keys:
        - 'execution': "remote"
        - 'remote_host': hostname or IP
        - 'remote_user': SSH username
        - 'remote_auth': dict with 'method', 'ssh_key' or 'password'
        - 'remote_dir': base path on remote machine
        - 'env_setup': (optional) shell commands to set up environment
                      (e.g., "source /etc/profile"). Defaults to "source /etc/profile"
                      if not provided. Useful for making module command available
                      in non-interactive SSH sessions.
    - self.job_file: name of the job script
    - self.submit_command(): method that returns the job submission command
    - self.logger: optional logger object with .info() and .warning()
    """

    _remote_sessions = {}
    _last_remote_path = None

    def _setup_remote(self):
        key = (self.queue["remote_host"], self.queue["remote_user"])
        if key not in self._remote_sessions:
            remote = RemoteAuth(
                username=self.queue["remote_user"],
                host=self.queue["remote_host"],
                auth_config=self.queue["remote_auth"],
            )
            remote.connect()
            self._remote_sessions[key] = remote
        self.remote = self._remote_sessions[key]

        # Generate unique remote directory name to avoid collisions
        # when multiple jobs with the same label run from different local directories
        # Format: {basename}_{hash} where hash is derived from full local path
        dir_basename = os.path.basename(self.calc.directory.rstrip(os.sep))
        
        # Create a hash of the full local directory path for uniqueness
        # Use first 8 characters of MD5 hash for brevity while maintaining uniqueness
        local_dir_hash = hashlib.md5(self.calc.directory.encode()).hexdigest()[:8]
        unique_dirname = f"{dir_basename}_{local_dir_hash}"
        
        current_path = os.path.join(self.queue["remote_dir"], unique_dirname)
        if current_path != self._last_remote_path:
            self.remote_path = current_path
            self._last_remote_path = current_path

    def _sha256(self, filepath):
        """Returns SHA256 checksum of a file."""
        h = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()

    def _transfer_pseudopotentials(self, max_retries=1):
        pseudopotentials = self.calc.parameters.get("pseudopotentials", {})
        remote_pseudo_dir = os.path.join(self.remote_path, "pseudo")
        self.remote.run_command(f"mkdir -p {remote_pseudo_dir}")

        search_dirs = []
        control = self.calc.parameters.get("input_data", {}).get("CONTROL", {})
        if "pseudo_dir" in control:
            search_dirs.append(control["pseudo_dir"])
        if "ESPRESSO_PSEUDO" in os.environ:
            search_dirs.append(os.path.join(os.environ["ESPRESSO_PSEUDO"]))
        search_dirs.append(os.path.expanduser("~/espresso/pseudo/"))

        missing_pseudos = []
        for symbol, pseudo_file in pseudopotentials.items():
            found = False
            for attempt in range(max_retries + 1):
                for pseudo_dir in search_dirs:
                    local_path = os.path.join(pseudo_dir, pseudo_file)
                    if os.path.exists(local_path):
                        remote_path = os.path.join(remote_pseudo_dir, pseudo_file)
                        self.remote.send_file(local_path, remote_path)

                        # Verify checksum
                        local_hash = self._sha256(local_path)
                        remote_hash = self.remote.sha256(remote_path)
                        if local_hash != remote_hash:
                            warnings.warn(
                                f"Checksum mismatch for {pseudo_file} after transfer."
                            )
                            if hasattr(self, "logger"):
                                self.logger.warning(f"Checksum mismatch: {pseudo_file}")
                        else:
                            if hasattr(self, "logger"):
                                self.logger.info(
                                    f"Transferred {pseudo_file} for {symbol} with verified checksum."
                                )
                        found = True
                        break
                if found:
                    break
            if not found:
                missing_pseudos.append((symbol, pseudo_file))
                warnings.warn(
                    f"Pseudopotential '{pseudo_file}' not found in any known directory."
                )
                if hasattr(self, "logger"):
                    self.logger.warning(
                        f"Missing pseudopotential: {pseudo_file} for {symbol}"
                    )

        # Raise exception if any pseudopotentials are missing
        if missing_pseudos:
            missing_list = ", ".join(
                [f"{symbol}: {pseudo_file}" for symbol, pseudo_file in missing_pseudos]
            )
            error_msg = f"Cannot proceed with calculation. Missing pseudopotentials: {missing_list}"
            if hasattr(self, "logger"):
                self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        self.calc.parameters["input_data"]["CONTROL"]["pseudo_dir"] = "./pseudo"
        self.calc.write_input(self.calc.atoms)

    def run(self):
        """
        Executes the calculation remotely if 'execution' is set to 'remote' in the queue.

        Behavior:
            - For remote execution:
                - Sets up SSH connection and remote working directory
                - Transfers input and job script files
                - Submits job using the scheduler's submit_command()
                - If scheduler is SLURM:
                    - Extracts job ID from submission output
                    - If wait_for_completion=False (default): returns immediately with job ID
                    - If wait_for_completion=True: polls job status until completion
                - If scheduler is direct:
                    - If wait_for_completion=False (default): runs job in background, captures PID, returns immediately
                    - If wait_for_completion=True: runs job in foreground, blocks until completion
                - Retrieves output file only in blocking mode after job finishes
            - For direct/local execution:
                - Falls back to the base Scheduler.run() method

        Queue Parameters:
            - wait_for_completion (bool): If False (default), submits job and returns immediately.
                                         If True, blocks until job completes.
                                         - For SLURM: job ID stored in calc.last_job_id
                                         - For direct: PID stored in calc.last_job_id as "PID:{pid}"
                                         Default is False to prevent GUI freezing.
            - job_timeout (int): Maximum time in seconds to wait for job completion (default: 3600).
                                Only used when wait_for_completion=True with SLURM scheduler.

        Returns:
            tuple: (stdout, stderr) from the job submission command
        """
        if self.queue.get("execution") != "remote":
            return super().run()

        self._setup_remote()

        input_file = f"{self.calc.prefix}.{self.calc.package}i"
        output_file = f"{self.calc.prefix}.{self.calc.package}o"
        job_file = self.job_file

        local_input = os.path.join(self.calc.directory, input_file)
        local_output = os.path.join(self.calc.directory, output_file)
        local_job = os.path.join(self.calc.directory, job_file)

        self._transfer_pseudopotentials()

        self.remote.send_file(local_input, f"{self.remote_path}/{input_file}")
        self.remote.send_file(local_job, f"{self.remote_path}/{job_file}")

        if hasattr(self, "logger"):
            self.logger.info(f"Submitting job via: {self.submit_command()}")

        # Use env_setup from queue if available, otherwise default to sourcing /etc/profile
        env_setup = self.queue.get("env_setup", "source /etc/profile")
        
        # Check if non-blocking mode is requested (default is False for non-blocking)
        wait_for_completion = self.queue.get("wait_for_completion", False)
        
        # Handle non-blocking execution for both SLURM and direct schedulers
        if not wait_for_completion:
            # For non-blocking mode, run job in background
            if self.queue.get("scheduler") == "slurm":
                # SLURM: submit job and get job ID
                command = f"cd {self.remote_path} && {env_setup} && {self.submit_command()}"
                stdout, stderr = self.remote.run_command(command)
                
                import re
                match = re.search(r"Submitted batch job (\d+)", stdout)
                job_id = match.group(1) if match else None
                
                if job_id:
                    if hasattr(self, "logger"):
                        self.logger.info(f"SLURM job {job_id} submitted (non-blocking mode).")
                    self.calc.last_job_id = job_id
                    # Store the actual remote path for monitoring
                    self.calc.last_remote_path = self.remote_path
                    return stdout, stderr
                else:
                    if hasattr(self, "logger"):
                        self.logger.warning("Could not extract SLURM job ID from sbatch output")
            else:
                # Direct scheduler: run in background and capture PID
                # Redirect output to log files and run in background
                bg_command = f"cd {self.remote_path} && {env_setup} && {self.submit_command()} > {output_file} 2>&1 & echo $!"
                stdout, stderr = self.remote.run_command(bg_command)
                
                # Extract PID from stdout
                pid = stdout.strip()
                if pid and pid.isdigit():
                    if hasattr(self, "logger"):
                        self.logger.info(f"Direct job started with PID: {pid} (non-blocking mode).")
                    self.calc.last_job_id = f"PID:{pid}"
                    # Store the actual remote path for monitoring
                    self.calc.last_remote_path = self.remote_path
                    return stdout, stderr
                else:
                    if hasattr(self, "logger"):
                        self.logger.warning(f"Could not extract PID from background job. Output: {stdout}")
        else:
            # Blocking mode: wait for completion
            command = f"cd {self.remote_path} && {env_setup} && {self.submit_command()}"
            stdout, stderr = self.remote.run_command(command)
            
            # If SLURM, extract job ID and wait for completion
            if self.queue.get("scheduler") == "slurm":
                import re
                
                match = re.search(r"Submitted batch job (\d+)", stdout)
                job_id = match.group(1) if match else None
                
                if job_id:
                    if hasattr(self, "logger"):
                        self.logger.info(f"Job submitted with ID: {job_id}")
                    self._wait_for_slurm_completion(job_id)
                else:
                    if hasattr(self, "logger"):
                        self.logger.warning("Could not extract job ID from sbatch output")
            
            # For both SLURM (after completion) and direct (already completed),
            # verify and retrieve output file
            self._verify_and_retrieve_output_file(output_file, local_output)

        return stdout, stderr

    def _wait_for_slurm_completion(self, job_id):
        """
        Wait for SLURM job completion with timeout and proper status checking.

        Args:
            job_id (str): The SLURM job ID to monitor
        """
        if hasattr(self, "logger"):
            self.logger.info(f"Waiting for SLURM job {job_id} to complete...")

        # Wait for job completion with timeout and better status checking
        timeout = self.queue.get("job_timeout", 3600)  # Default 1 hour timeout
        start_time = time.time()

        while True:
            if time.time() - start_time > timeout:
                error_msg = f"Job {job_id} timed out after {timeout} seconds"
                if hasattr(self, "logger"):
                    self.logger.error(error_msg)
                raise RuntimeError(error_msg)

            status_stdout, status_stderr = self.remote.run_command(
                f"squeue -j {job_id} -h -o '%T'"
            )

            # If squeue returns no output, job is no longer in queue (completed or failed)
            if not status_stdout.strip():
                # Check if job completed successfully using sacct
                sacct_stdout, sacct_stderr = self.remote.run_command(
                    f"sacct -j {job_id} -n -o State --parsable2"
                )
                if sacct_stdout.strip():
                    job_state = sacct_stdout.strip().split("\n")[0]
                    if hasattr(self, "logger"):
                        self.logger.info(
                            f"Job {job_id} finished with state: {job_state}"
                        )

                    if job_state not in ["COMPLETED", "COMPLETING"]:
                        error_msg = f"Job {job_id} failed with state: {job_state}"
                        if hasattr(self, "logger"):
                            self.logger.error(error_msg)
                        raise RuntimeError(error_msg)
                break

            time.sleep(10)

    def _verify_and_retrieve_output_file(
        self, output_file, local_output, max_retries=3
    ):
        """
        Verifies that the output file exists on the remote system and retrieves it with retry logic.

        Args:
            output_file (str): Name of the output file
            local_output (str): Local path where the output file should be saved
            max_retries (int): Maximum number of retry attempts for file retrieval
        """
        remote_output_path = f"{self.remote_path}/{output_file}"

        # Check if output file exists on remote system
        check_cmd = f"[ -f {remote_output_path} ] && echo 'exists' || echo 'missing'"
        stdout, stderr = self.remote.run_command(check_cmd)

        if "missing" in stdout:
            error_msg = (
                f"Output file {remote_output_path} does not exist on remote system"
            )
            if hasattr(self, "logger"):
                self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)

        # Retrieve file with retry logic
        for attempt in range(max_retries):
            try:
                self.remote.retrieve_file(remote_output_path, local_output)
                if hasattr(self, "logger"):
                    self.logger.info(
                        f"Successfully retrieved {output_file} to {local_output}"
                    )
                return
            except Exception as e:
                if attempt < max_retries - 1:
                    if hasattr(self, "logger"):
                        self.logger.warning(
                            f"Attempt {attempt + 1} failed to retrieve {output_file}: {e}. Retrying..."
                        )
                    time.sleep(2**attempt)  # Exponential backoff
                else:
                    error_msg = f"Failed to retrieve {output_file} after {max_retries} attempts: {e}"
                    if hasattr(self, "logger"):
                        self.logger.error(error_msg)
                    raise RuntimeError(error_msg)

    @classmethod
    def close_all_connections(cls):
        """Closes all remote connections and clears the session cache."""
        for remote in cls._remote_sessions.values():
            remote.close()
        cls._remote_sessions.clear()
        cls._last_remote_path = None

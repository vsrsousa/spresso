"""
Monitor for remote job execution (SLURM, direct scheduler, etc).

Provides a unified interface to check job status, retrieve output, and monitor progress.
"""

import time
import logging
from typing import Optional, Dict, Tuple

logger = logging.getLogger(__name__)


class RemoteJobMonitor:
    """
    Monitor for tracking remote job execution.
    
    Supports SLURM and direct (bash) job execution with unified interface.
    
    Examples:
        >>> # After running non-blocking remote job
        >>> calc = workflow.run_scf(label='scf/si-test')
        >>> monitor = RemoteJobMonitor(calc)
        >>> print(monitor.status())  # 'running', 'completed', 'failed'
        >>> if monitor.wait(timeout=3600):
        ...     output = monitor.retrieve_output()
    """
    
    def __init__(self, calc, remote_connection=None):
        """
        Initialize job monitor.
        
        Args:
            calc: Espresso calculator object with last_job_id and last_remote_path
            remote_connection: Optional remote connection object (auto-detected if None)
        """
        self.calc = calc
        self.job_id = getattr(calc, 'last_job_id', None)
        self.remote_path = getattr(calc, 'last_remote_path', None)
        
        # Auto-detect remote connection from calc if available
        if remote_connection is None:
            # Try calc.remote first
            remote_connection = getattr(calc, 'remote', None)
            # If not found, try calc.scheduler.remote (for Espresso calculator)
            if remote_connection is None and hasattr(calc, 'scheduler'):
                remote_connection = getattr(calc.scheduler, 'remote', None)
        
        self.remote = remote_connection
        
        if not self.job_id or not self.remote_path:
            raise ValueError("Calculator must have last_job_id and last_remote_path set (from non-blocking remote execution)")
        
        if not self.remote:
            raise ValueError("No remote connection available. Please provide remote_connection parameter or ensure calc.remote is set.")
        
        # Detect job type from job_id format
        if isinstance(self.job_id, str) and self.job_id.startswith('PID:'):
            self.job_type = 'direct'
            self.pid = self.job_id.replace('PID:', '')
        else:
            self.job_type = 'slurm'
            self.slurm_job_id = self.job_id
    
    def status(self) -> str:
        """
        Get current job status.
        
        Returns:
            str: 'running', 'completed', 'failed', 'unknown'
        """
        if not self.remote:
            logger.warning("No remote connection available for status check")
            return 'unknown'
        
        try:
            if self.job_type == 'slurm':
                return self._check_slurm_status()
            else:  # direct scheduler
                return self._check_direct_status()
        except Exception as e:
            logger.error(f"Error checking job status: {e}")
            return 'unknown'
    
    def _check_slurm_status(self) -> str:
        """Check SLURM job status via squeue/sacct."""
        # First try squeue (running jobs)
        stdout, stderr = self.remote.run_command(
            f"squeue -j {self.slurm_job_id} -h -o '%T'"
        )
        
        if stdout.strip():
            state = stdout.strip()
            if state in ['RUNNING', 'PENDING', 'CONFIGURING']:
                return 'running'
            elif state in ['COMPLETED', 'COMPLETING']:
                return 'completed'
            else:
                return 'failed'
        
        # Job not in squeue, check sacct (completed/failed jobs)
        stdout, stderr = self.remote.run_command(
            f"sacct -j {self.slurm_job_id} -n -o State --parsable2"
        )
        
        if stdout.strip():
            state = stdout.strip().split('\n')[0]
            if state in ['COMPLETED', 'COMPLETING']:
                return 'completed'
            else:
                return 'failed'
        
        return 'unknown'
    
    def _check_direct_status(self) -> str:
        """Check direct job status via ps."""
        stdout, stderr = self.remote.run_command(f"ps -p {self.pid} -o pid=")
        
        if stdout.strip():
            return 'running'
        
        # Process not found, check if output was generated
        output_file = f"{self.remote_path}/{self.calc.prefix}.{self.calc.package}o"
        stdout, stderr = self.remote.run_command(f"[ -f {output_file} ] && echo 'exists' || echo 'missing'")
        
        if 'exists' in stdout:
            return 'completed'
        else:
            return 'failed'
    
    def wait(self, timeout: int = 3600, poll_interval: int = 10) -> bool:
        """
        Wait for job completion.
        
        Args:
            timeout: Maximum time to wait in seconds (default 1 hour)
            poll_interval: How often to check status in seconds (default 10)
            
        Returns:
            bool: True if completed, False if timeout
        """
        start_time = time.time()
        
        while time.time() - start_time < timeout:
            status = self.status()
            
            if status == 'completed':
                logger.info(f"Job {self.job_id} completed")
                return True
            elif status == 'failed':
                logger.warning(f"Job {self.job_id} failed")
                return False
            
            logger.info(f"Job {self.job_id} still running... ({int(time.time() - start_time)}s)")
            time.sleep(poll_interval)
        
        logger.error(f"Job {self.job_id} timed out after {timeout}s")
        return False
    
    def retrieve_output(self) -> Tuple[str, str]:
        """
        Retrieve job output files from remote.
        
        Returns:
            tuple: (local_output_path, output_file_content)
        """
        output_file = f"{self.calc.prefix}.{self.calc.package}o"
        remote_output_path = f"{self.remote_path}/{output_file}"
        local_output = f"{self.calc.directory}/{output_file}"
        
        try:
            self.remote.retrieve_file(remote_output_path, local_output)
            
            with open(local_output, 'r') as f:
                content = f.read()
            
            logger.info(f"Retrieved output: {local_output}")
            return local_output, content
        except Exception as e:
            logger.error(f"Failed to retrieve output: {e}")
            raise
    
    def cancel(self) -> bool:
        """
        Cancel the running job.
        
        Returns:
            bool: True if cancellation succeeded
        """
        try:
            if self.job_type == 'slurm':
                self.remote.run_command(f"scancel {self.slurm_job_id}")
            else:  # direct scheduler
                self.remote.run_command(f"kill {self.pid}")
            
            logger.info(f"Job {self.job_id} cancelled")
            return True
        except Exception as e:
            logger.error(f"Failed to cancel job: {e}")
            return False
    
    def info(self) -> Dict:
        """
        Get detailed job information.
        
        Returns:
            dict: Status, resource info, etc.
        """
        info = {
            'job_id': self.job_id,
            'job_type': self.job_type,
            'status': self.status(),
            'remote_path': self.remote_path,
        }
        
        if self.job_type == 'slurm':
            try:
                stdout, _ = self.remote.run_command(
                    f"squeue -j {self.slurm_job_id} -o '%j %T %M %l' -h"
                )
                if stdout.strip():
                    parts = stdout.strip().split()
                    info['time_used'] = parts[2] if len(parts) > 2 else 'N/A'
                    info['time_limit'] = parts[3] if len(parts) > 3 else 'N/A'
            except Exception:
                pass
        
        return info

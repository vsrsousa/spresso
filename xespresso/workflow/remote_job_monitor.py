"""
Remote Job Monitoring for xespresso Calculations

Monitors the status of remote SLURM jobs and retrieves output when complete.
"""

import time
import logging
from typing import Optional

logger = logging.getLogger(__name__)


class RemoteJobMonitor:
    """
    Monitor remote SLURM job status and retrieve output.
    
    Used to track parallel job submissions and retrieve results
    from the remote scheduler.
    """
    
    def __init__(self, calc):
        """
        Initialize monitor for a calculator's remote job.
        
        Args:
            calc: Espresso calculator object with remote execution info
        """
        self.calc = calc
        self.scheduler = getattr(calc, 'scheduler', None)
        self.job_id = getattr(calc, 'last_job_id', None)
        self.directory = getattr(calc, 'directory', None)
        
        # Try to get remote connection from calc, with fallback to scheduler
        self.remote = getattr(calc, 'remote', None)
        if self.remote is None and self.scheduler is not None:
            self.remote = getattr(self.scheduler, 'remote', None)
        
        if not self.job_id:
            logger.warning(f"No job_id found in calculator (needed for monitoring)")
        if not self.remote:
            logger.warning("No remote connection found in calculator or scheduler")
    
    def wait(self, timeout: int = 3600, poll_interval: int = 30) -> bool:
        """
        Wait for remote job to complete.
        
        Args:
            timeout: Maximum time to wait (seconds)
            poll_interval: Time between status checks (seconds)
            
        Returns:
            True if job completed successfully, False if timeout
        """
        if not self.remote or not self.job_id:
            logger.warning(f"Cannot monitor job: remote={self.remote is not None}, job_id={self.job_id}")
            return False
        
        start_time = time.time()
        
        logger.info(f"\nMonitoring remote job {self.job_id}...")
        logger.info(f"Timeout: {timeout}s | Poll interval: {poll_interval}s")
        logger.info("-" * 70)
        
        last_status = None
        unchanged_count = 0
        
        while True:
            elapsed = time.time() - start_time
            
            # Check timeout
            if elapsed > timeout:
                logger.warning(f"[{elapsed:.0f}s] Timeout reached ({timeout}s). Stopping monitoring.")
                logger.warning(f"Job {self.job_id} may still be running on remote.")
                logger.warning(f"Check manually with: squeue -j {self.job_id} or sacct -j {self.job_id}")
                return False
            
            # Get job status
            try:
                status = self._get_job_status()
                
                # Log status changes
                if status != last_status:
                    logger.info(f"[{elapsed:.0f}s] State: {status}")
                    last_status = status
                    unchanged_count = 0
                else:
                    unchanged_count += 1
                    # Only log periodically if status hasn't changed
                    if unchanged_count % 5 == 0:  # Log every 5 polls
                        logger.debug(f"[{elapsed:.0f}s] State: {status} (unchanged)")
                
                if status == 'COMPLETED':
                    logger.info(f"[{elapsed:.0f}s] ✓ Job COMPLETED successfully")
                    return True
                elif status == 'FAILED':
                    logger.error(f"[{elapsed:.0f}s] ✗ Job FAILED")
                    return False
                # else: RUNNING, PENDING, or UNKNOWN - continue polling
                    
            except Exception as e:
                logger.error(f"Error checking job status: {e}")
                logger.debug(f"Exception details: {type(e).__name__}")
            
            # Wait before next poll
            time.sleep(poll_interval)
    
    def _get_job_status(self) -> str:
        """
        Get current job status using squeue (running) and sacct (completed/failed).
        
        Returns:
            'PENDING', 'RUNNING', 'COMPLETED', 'FAILED', or 'UNKNOWN'
        """
        if not self.remote or not self.job_id:
            return 'UNKNOWN'
        
        try:
            # First, try squeue (for running/pending jobs)
            stdout, stderr = self.remote.run_command(
                f"squeue -j {self.job_id} -h -o '%T'"
            )
            
            if stdout.strip():
                # Job is in queue
                state = stdout.strip().split('\n')[0]  # Get first line
                if state in ['RUNNING', 'PENDING', 'CONFIGURING']:
                    return state
                elif state in ['COMPLETED', 'COMPLETING']:
                    return 'COMPLETED'
                else:
                    return 'FAILED'
            
            # Job not in squeue - check sacct for completed/failed status
            logger.debug(f"Job {self.job_id} not in squeue, checking sacct...")
            
            stdout_sacct, stderr_sacct = self.remote.run_command(
                f"sacct -j {self.job_id} -n -o State --parsable2"
            )
            
            if stdout_sacct.strip():
                # Get last (most recent) state line
                lines = [l.strip() for l in stdout_sacct.strip().split('\n') if l.strip()]
                state = lines[-1] if lines else 'UNKNOWN'
                
                logger.debug(f"Job {self.job_id} sacct state: {state}")
                
                if state in ['COMPLETED', 'COMPLETING']:
                    return 'COMPLETED'
                elif state in ['FAILED', 'TIMEOUT', 'CANCELLED', 'OUT_OF_MEMORY']:
                    return 'FAILED'
                else:
                    # Unknown state from sacct
                    return 'UNKNOWN'
            
            # No info from either squeue or sacct
            logger.warning(f"Job {self.job_id} not found in squeue or sacct")
            return 'UNKNOWN'
            
        except Exception as e:
            logger.debug(f"Error checking job status: {e}")
            return 'UNKNOWN'
    
    def retrieve_output(self) -> bool:
        """
        Retrieve output files from remote job.
        
        Returns:
            True if successful, False otherwise
        """
        if not self.directory:
            logger.warning("No directory information available")
            return False
        
        try:
            # Try to read results from remote directory
            if hasattr(self.calc, 'read'):
                # ASE calculator's read method retrieves results
                self.calc.read(self.directory)
                logger.info(f"Retrieved output from {self.directory}")
                return True
            else:
                logger.warning("Calculator has no read() method")
                return False
                
        except Exception as e:
            logger.error(f"Failed to retrieve output: {e}")
            return False

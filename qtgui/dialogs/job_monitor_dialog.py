"""
Job Monitor Dialog for tracking remote job submissions.

This module provides a non-blocking dialog window that displays all submitted
remote jobs, their status, and allows checking status and retrieving results.
Job information persists across sessions using a JSON file.
"""

import json
import os
from datetime import datetime
from pathlib import Path

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidget,
    QTableWidgetItem, QLabel, QMessageBox, QHeaderView, QWidget,
    QTextEdit
)
from PySide6.QtCore import Qt, Signal, QTimer, QThread
from PySide6.QtGui import QColor


class JobStatusWorker(QThread):
    """
    Worker thread for checking job status without blocking the UI.
    
    Signals:
        status_updated: Emitted when status check completes (row, new_status, error_msg)
    """
    status_updated = Signal(int, str, str)  # row, status, error_message
    
    def __init__(self, row, job, parent=None):
        super().__init__(parent)
        self.row = row
        self.job = job
    
    def run(self):
        """Check job status in background thread."""
        job_id = self.job.get('job_id')
        scheduler = self.job.get('scheduler')
        queue = self.job.get('queue', {})
        
        # Skip if job_id is unknown or None
        if not job_id or job_id == 'Unknown':
            self.status_updated.emit(self.row, 'unknown', 'Job ID is unknown')
            return
        
        try:
            # Import here to avoid circular dependencies
            from xespresso.schedulers.remote_connection import RemoteConnection
            
            # Create remote connection
            remote = RemoteConnection(queue)
            
            # Check status based on scheduler type
            if scheduler == 'slurm':
                # Check SLURM job status
                cmd = f"squeue -j {job_id} -h -o '%T' 2>/dev/null || echo 'NOT_FOUND'"
                result = remote.run_command(cmd)
                status_output = result.strip()
                
                if status_output == 'NOT_FOUND' or not status_output:
                    # Job not in queue, check if completed
                    cmd_history = f"sacct -j {job_id} -n -o State -X 2>/dev/null | head -1"
                    result_history = remote.run_command(cmd_history)
                    history_status = result_history.strip()
                    
                    if 'COMPLETED' in history_status:
                        new_status = 'completed'
                    elif 'FAILED' in history_status or 'CANCELLED' in history_status:
                        new_status = 'failed'
                    else:
                        new_status = 'unknown'
                elif 'RUNNING' in status_output:
                    new_status = 'running'
                elif 'PENDING' in status_output:
                    new_status = 'pending'
                elif 'COMPLETED' in status_output:
                    new_status = 'completed'
                else:
                    new_status = 'unknown'
                    
            elif scheduler == 'direct':
                # Check process status using PID
                pid = job_id.replace('PID:', '') if 'PID:' in str(job_id) else job_id
                cmd = f"ps -p {pid} -o state= 2>/dev/null || echo 'NOT_FOUND'"
                result = remote.run_command(cmd)
                
                if 'NOT_FOUND' in result or not result.strip():
                    # Process not found, assume completed
                    new_status = 'completed'
                else:
                    new_status = 'running'
            else:
                new_status = 'unknown'
            
            remote.disconnect()
            self.status_updated.emit(self.row, new_status, '')
            
        except Exception as e:
            self.status_updated.emit(self.row, '', str(e))


class JobMonitorDialog(QDialog):
    """
    Non-blocking job monitor dialog for tracking remote job submissions.
    
    Displays all submitted jobs with their status and provides buttons to:
    - Refresh job status from remote server
    - Retrieve results when jobs complete
    - Cancel running jobs
    - View job details
    
    Job information is persisted to .spresso_jobs.json in the working directory.
    
    Signals:
        job_completed: Emitted when a job completes successfully (job_id, label)
    """
    
    job_completed = Signal(str, str)
    
    def __init__(self, working_dir=None, parent=None):
        """
        Initialize the job monitor dialog.
        
        Args:
            working_dir: Directory where jobs file is stored (defaults to current directory)
            parent: Parent widget (main window)
        """
        super().__init__(parent)
        self.working_dir = working_dir or os.getcwd()
        self.jobs_file = os.path.join(self.working_dir, ".spresso_jobs.json")
        self.jobs = self._load_jobs()
        
        self.setWindowTitle("🔍 Job Monitor")
        self.setMinimumSize(1000, 600)
        
        # Make the dialog non-modal (doesn't block main window)
        self.setModal(False)
        
        # Set window flags to keep dialog on top but allow interaction with main window
        self.setWindowFlags(
            Qt.Window | Qt.WindowCloseButtonHint | Qt.WindowMinimizeButtonHint
        )
        
        self._init_ui()
        self._refresh_table()
        
        # Auto-refresh timer (optional, disabled by default)
        self.auto_refresh_timer = QTimer()
        self.auto_refresh_timer.timeout.connect(self._refresh_all_status)
        
    def _init_ui(self):
        """Initialize the user interface."""
        layout = QVBoxLayout()
        
        # Header
        header_label = QLabel("<h2>🔍 Remote Job Monitor</h2>")
        header_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(header_label)
        
        info_label = QLabel(
            "Track submitted remote jobs, check their status, and retrieve results."
        )
        info_label.setAlignment(Qt.AlignCenter)
        info_label.setStyleSheet("color: gray; margin-bottom: 10px;")
        layout.addWidget(info_label)
        
        # Jobs table
        self.table = QTableWidget()
        self.table.setColumnCount(7)
        self.table.setHorizontalHeaderLabels([
            "Label", "Job ID", "Status", "Scheduler", "Submitted Time",
            "Remote Host", "Actions"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        
        # Set column widths
        self.table.setColumnWidth(0, 150)  # Label
        self.table.setColumnWidth(1, 120)  # Job ID
        self.table.setColumnWidth(2, 100)  # Status
        self.table.setColumnWidth(3, 80)   # Scheduler
        self.table.setColumnWidth(4, 150)  # Time
        self.table.setColumnWidth(5, 200)  # Host
        self.table.setColumnWidth(6, 200)  # Actions
        
        layout.addWidget(self.table)
        
        # Button panel
        button_layout = QHBoxLayout()
        
        self.refresh_all_btn = QPushButton("🔄 Refresh All")
        self.refresh_all_btn.clicked.connect(self._refresh_all_status)
        button_layout.addWidget(self.refresh_all_btn)
        
        self.clear_completed_btn = QPushButton("🗑️ Clear Completed")
        self.clear_completed_btn.clicked.connect(self._clear_completed_jobs)
        button_layout.addWidget(self.clear_completed_btn)
        
        button_layout.addStretch()
        
        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.close)
        button_layout.addWidget(self.close_btn)
        
        layout.addLayout(button_layout)
        
        # Status bar
        self.status_label = QLabel(f"Total jobs: {len(self.jobs)}")
        self.status_label.setStyleSheet("padding: 5px; background-color: #f0f0f0;")
        layout.addWidget(self.status_label)
        
        self.setLayout(layout)
    
    def _load_jobs(self):
        """Load jobs from JSON file."""
        if os.path.exists(self.jobs_file):
            try:
                with open(self.jobs_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading jobs file: {e}")
                return []
        return []
    
    def _save_jobs(self):
        """Save jobs to JSON file."""
        try:
            with open(self.jobs_file, 'w') as f:
                json.dump(self.jobs, f, indent=2)
        except Exception as e:
            print(f"Error saving jobs file: {e}")
    
    def add_job(self, job_info):
        """
        Add a new job to the monitor.
        
        Args:
            job_info: Dictionary containing job information:
                - label: Job label/name
                - job_id: Job ID or PID
                - scheduler: 'slurm' or 'direct'
                - remote_host: Remote hostname
                - remote_user: Remote username
                - remote_dir: Remote working directory
                - local_dir: Local working directory
                - queue: Full queue configuration dict
                - submitted_time: ISO timestamp (optional, auto-generated)
                - status: Initial status (optional, defaults to 'submitted')
        """
        if 'submitted_time' not in job_info:
            job_info['submitted_time'] = datetime.now().isoformat()
        if 'status' not in job_info:
            job_info['status'] = 'submitted'
        
        self.jobs.append(job_info)
        self._save_jobs()
        self._refresh_table()
        
        # Show the dialog if it's hidden
        if not self.isVisible():
            self.show()
    
    def _refresh_table(self):
        """Refresh the jobs table display."""
        self.table.setRowCount(len(self.jobs))
        
        for row, job in enumerate(self.jobs):
            # Label
            self.table.setItem(row, 0, QTableWidgetItem(job.get('label', 'N/A')))
            
            # Job ID
            job_id = job.get('job_id', 'N/A')
            self.table.setItem(row, 1, QTableWidgetItem(str(job_id)))
            
            # Status with color coding
            status = job.get('status', 'unknown')
            status_item = QTableWidgetItem(status.upper())
            if status == 'completed':
                status_item.setBackground(QColor(200, 255, 200))  # Light green
            elif status == 'failed':
                status_item.setBackground(QColor(255, 200, 200))  # Light red
            elif status == 'running':
                status_item.setBackground(QColor(200, 230, 255))  # Light blue
            elif status == 'pending':
                status_item.setBackground(QColor(255, 255, 200))  # Light yellow
            self.table.setItem(row, 2, status_item)
            
            # Scheduler
            scheduler = job.get('scheduler', 'N/A')
            self.table.setItem(row, 3, QTableWidgetItem(scheduler))
            
            # Submitted time
            submitted_time = job.get('submitted_time', 'N/A')
            if submitted_time != 'N/A':
                try:
                    dt = datetime.fromisoformat(submitted_time)
                    time_str = dt.strftime("%Y-%m-%d %H:%M:%S")
                except:
                    time_str = submitted_time
            else:
                time_str = 'N/A'
            self.table.setItem(row, 4, QTableWidgetItem(time_str))
            
            # Remote host
            remote_host = job.get('remote_host', 'N/A')
            self.table.setItem(row, 5, QTableWidgetItem(remote_host))
            
            # Actions buttons
            actions_widget = self._create_action_buttons(row, job)
            self.table.setCellWidget(row, 6, actions_widget)
        
        # Update status label
        running_count = sum(1 for j in self.jobs if j.get('status') in ['submitted', 'running', 'pending'])
        self.status_label.setText(
            f"Total jobs: {len(self.jobs)} | Running/Pending: {running_count}"
        )
    
    def _create_action_buttons(self, row, job):
        """Create action buttons for a job row."""
        widget = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(2, 2, 2, 2)
        
        # Check Status button
        check_btn = QPushButton("Check")
        check_btn.setMaximumWidth(60)
        check_btn.clicked.connect(lambda: self._check_job_status(row))
        layout.addWidget(check_btn)
        
        # Retrieve Results button (only if completed)
        if job.get('status') == 'completed':
            retrieve_btn = QPushButton("Retrieve")
            retrieve_btn.setMaximumWidth(70)
            retrieve_btn.clicked.connect(lambda: self._retrieve_results(row))
            layout.addWidget(retrieve_btn)
        
        # Details button
        details_btn = QPushButton("Details")
        details_btn.setMaximumWidth(70)
        details_btn.clicked.connect(lambda: self._show_job_details(row))
        layout.addWidget(details_btn)
        
        layout.addStretch()
        widget.setLayout(layout)
        return widget
    
    def _check_job_status(self, row):
        """Check status of a specific job."""
        if row >= len(self.jobs):
            return
        
        job = self.jobs[row]
        job_id = job.get('job_id')
        scheduler = job.get('scheduler')
        queue = job.get('queue', {})
        
        try:
            # Import here to avoid circular dependencies
            from xespresso.schedulers.remote_connection import RemoteConnection
            
            # Create remote connection
            remote = RemoteConnection(queue)
            
            # Check status based on scheduler type
            if scheduler == 'slurm':
                # Check SLURM job status
                cmd = f"squeue -j {job_id} -h -o '%T' 2>/dev/null || echo 'NOT_FOUND'"
                result = remote.run_command(cmd)
                status_output = result.strip()
                
                if status_output == 'NOT_FOUND' or not status_output:
                    # Job not in queue, check if completed
                    cmd_history = f"sacct -j {job_id} -n -o State -X 2>/dev/null | head -1"
                    result_history = remote.run_command(cmd_history)
                    history_status = result_history.strip()
                    
                    if 'COMPLETED' in history_status:
                        job['status'] = 'completed'
                    elif 'FAILED' in history_status or 'CANCELLED' in history_status:
                        job['status'] = 'failed'
                    else:
                        job['status'] = 'unknown'
                elif 'RUNNING' in status_output:
                    job['status'] = 'running'
                elif 'PENDING' in status_output:
                    job['status'] = 'pending'
                elif 'COMPLETED' in status_output:
                    job['status'] = 'completed'
                else:
                    job['status'] = 'unknown'
                    
            elif scheduler == 'direct':
                # Check process status using PID
                pid = job_id.replace('PID:', '') if 'PID:' in job_id else job_id
                cmd = f"ps -p {pid} -o state= 2>/dev/null || echo 'NOT_FOUND'"
                result = remote.run_command(cmd)
                
                if 'NOT_FOUND' in result or not result.strip():
                    # Process not found, assume completed
                    job['status'] = 'completed'
                else:
                    job['status'] = 'running'
            
            remote.disconnect()
            
            # Save updated status
            self._save_jobs()
            self._refresh_table()
            
            QMessageBox.information(
                self,
                "Status Updated",
                f"Job '{job.get('label')}' status: {job['status'].upper()}"
            )
            
        except Exception as e:
            QMessageBox.warning(
                self,
                "Status Check Failed",
                f"Failed to check job status:\n{str(e)}"
            )
    
    def _refresh_all_status(self):
        """Refresh status of all non-completed jobs."""
        self.refresh_all_btn.setEnabled(False)
        self.refresh_all_btn.setText("Refreshing...")
        
        updated_count = 0
        for row in range(len(self.jobs)):
            job = self.jobs[row]
            if job.get('status') not in ['completed', 'failed']:
                self._check_job_status(row)
                updated_count += 1
        
        self.refresh_all_btn.setEnabled(True)
        self.refresh_all_btn.setText("🔄 Refresh All")
        
        if updated_count > 0:
            QMessageBox.information(
                self,
                "Refresh Complete",
                f"Checked status of {updated_count} job(s)"
            )
    
    def _retrieve_results(self, row):
        """Retrieve results from a completed job."""
        if row >= len(self.jobs):
            return
        
        job = self.jobs[row]
        queue = job.get('queue', {})
        
        try:
            from xespresso.schedulers.remote_connection import RemoteConnection
            
            remote = RemoteConnection(queue)
            
            # Determine output file name
            label = job.get('label', 'calculation')
            prefix = label.split('/')[-1]  # Get last part if it's a path
            remote_dir = job.get('remote_dir', '')
            local_dir = job.get('local_dir', os.getcwd())
            
            # Create local directory if it doesn't exist
            os.makedirs(local_dir, exist_ok=True)
            
            # Download output files
            output_files = [f"{prefix}.pwo", f"{prefix}.pwi", "job_file"]
            downloaded = []
            
            for filename in output_files:
                remote_file = f"{remote_dir}/{label}/{filename}"
                local_file = os.path.join(local_dir, label, filename)
                
                try:
                    # Create local subdirectory
                    os.makedirs(os.path.dirname(local_file), exist_ok=True)
                    remote.get_file(remote_file, local_file)
                    downloaded.append(filename)
                except Exception as e:
                    print(f"Could not download {filename}: {e}")
            
            remote.disconnect()
            
            if downloaded:
                QMessageBox.information(
                    self,
                    "Results Retrieved",
                    f"Downloaded {len(downloaded)} file(s) to:\n{local_dir}/{label}\n\n"
                    f"Files: {', '.join(downloaded)}"
                )
            else:
                QMessageBox.warning(
                    self,
                    "No Files Retrieved",
                    "Could not download any output files. They may not exist yet."
                )
                
        except Exception as e:
            QMessageBox.warning(
                self,
                "Retrieval Failed",
                f"Failed to retrieve results:\n{str(e)}"
            )
    
    def _show_job_details(self, row):
        """Show detailed information about a job."""
        if row >= len(self.jobs):
            return
        
        job = self.jobs[row]
        
        details = f"""
<h3>Job Details</h3>
<table>
<tr><td><b>Label:</b></td><td>{job.get('label', 'N/A')}</td></tr>
<tr><td><b>Job ID:</b></td><td>{job.get('job_id', 'N/A')}</td></tr>
<tr><td><b>Status:</b></td><td>{job.get('status', 'N/A').upper()}</td></tr>
<tr><td><b>Scheduler:</b></td><td>{job.get('scheduler', 'N/A')}</td></tr>
<tr><td><b>Remote Host:</b></td><td>{job.get('remote_host', 'N/A')}</td></tr>
<tr><td><b>Remote User:</b></td><td>{job.get('remote_user', 'N/A')}</td></tr>
<tr><td><b>Remote Directory:</b></td><td>{job.get('remote_dir', 'N/A')}/{job.get('label', '')}</td></tr>
<tr><td><b>Local Directory:</b></td><td>{job.get('local_dir', 'N/A')}/{job.get('label', '')}</td></tr>
<tr><td><b>Submitted:</b></td><td>{job.get('submitted_time', 'N/A')}</td></tr>
</table>
"""
        
        # Create a dialog to show details
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Job Details - {job.get('label')}")
        dialog.setMinimumSize(500, 400)
        
        layout = QVBoxLayout()
        text_edit = QTextEdit()
        text_edit.setHtml(details)
        text_edit.setReadOnly(True)
        layout.addWidget(text_edit)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)
        
        dialog.setLayout(layout)
        dialog.exec()
    
    def _clear_completed_jobs(self):
        """Remove completed and failed jobs from the list."""
        reply = QMessageBox.question(
            self,
            "Clear Completed Jobs",
            "Remove all completed and failed jobs from the monitor?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            initial_count = len(self.jobs)
            self.jobs = [j for j in self.jobs if j.get('status') not in ['completed', 'failed']]
            removed_count = initial_count - len(self.jobs)
            
            self._save_jobs()
            self._refresh_table()
            
            QMessageBox.information(
                self,
                "Jobs Cleared",
                f"Removed {removed_count} completed/failed job(s)"
            )

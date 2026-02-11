from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QTabWidget, QFormLayout, QComboBox,
    QTextEdit, QPushButton, QHBoxLayout, QApplication, QMessageBox
)
from qtpy.QtCore import Qt
import os

# Lightweight minimal CalculationWindow with Test Connection copied from machine_config.py
try:
    from xespresso.machines.config.loader import list_machines, load_machine
except Exception:
    list_machines = lambda *a, **k: []
    load_machine = lambda *a, **k: None
try:
    from xespresso.codes.manager import load_codes_config
    CODES_AVAILABLE = True
except Exception:
    load_codes_config = lambda *a, **k: None
    CODES_AVAILABLE = False


class CalculationWindow(QWidget):
    def __init__(self, calc_name: str, session_state=None, parent=None):
        super().__init__(parent)
        self.calc_name = calc_name
        self.session_state = session_state or {}
        self.setWindowTitle(f"Calculation: {calc_name}")
        self.resize(800, 520)

        layout = QVBoxLayout(self)
        tabs = QTabWidget()
        layout.addWidget(tabs)

        # Machine tab
        machine = QWidget()
        mform = QFormLayout(machine)

        self.machine_combo = QComboBox()
        try:
            for m in (list_machines() or []):
                self.machine_combo.addItem(str(m))
        except Exception:
            pass
        mform.addRow('Machine:', self.machine_combo)

        self.test_conn_btn = QPushButton('Test Connection')
        self.results_label = QLabel('')
        self.results_label.setWordWrap(True)

        def _test_conn():
            mn = self.machine_combo.currentText() or ''
            if not mn and self.session_state:
                mn = self.session_state.get('current_machine_name') or ''
            machine_obj = None
            try:
                machine_obj = load_machine(mn)
            except Exception:
                try:
                    machine_obj = load_machine(None, mn, None, return_object=True)
                except Exception:
                    machine_obj = None

            execution = None
            try:
                execution = getattr(machine_obj, 'execution', None)
            except Exception:
                execution = None
            if not execution and self.session_state:
                execution = self.session_state.get('machines', {}).get(mn, {}).get('execution')

            if execution == 'local':
                workdir = getattr(machine_obj, 'workdir', None) or (self.session_state.get('machines', {}).get(mn, {}).get('workdir') if self.session_state else None) or './calculations'
                user = os.environ.get('USER', 'unknown')
                self.results_label.setText(f"✅ Local machine - connection OK\nWorking directory: {workdir}\nCurrent user: {user}")
                self.results_label.setStyleSheet('color: green;')
                return

            host = None; username = None; port = 22; ssh_key = None
            try:
                if machine_obj is not None:
                    host = getattr(machine_obj, 'host', None)
                    username = getattr(machine_obj, 'username', None)
                    port = getattr(machine_obj, 'port', 22) or 22
                    auth = getattr(machine_obj, 'auth', None)
                    if isinstance(auth, dict):
                        ssh_key = auth.get('ssh_key') or auth.get('key')
                if host is None and self.session_state:
                    host = self.session_state.get('machines', {}).get(mn, {}).get('host')
                if username is None and self.session_state:
                    username = self.session_state.get('machines', {}).get(mn, {}).get('username')
                if ssh_key is None and self.session_state:
                    ssh_key = self.session_state.get('machines', {}).get(mn, {}).get('ssh_key')
            except Exception:
                pass

            if not host or not username:
                QMessageBox.warning(self, "Warning", "Please enter host and username")
                return

            key_path = os.path.expanduser(ssh_key) if ssh_key else None
            if not key_path or not os.path.isfile(key_path):
                self.results_label.setText(f"❌ SSH key not found: {key_path}\n💡 Check the SSH key path")
                self.results_label.setStyleSheet('color: red;')
                return

            self.results_label.setText("Testing SSH connection...")
            self.results_label.setStyleSheet('color: blue;')
            QApplication.processEvents()
            try:
                from xespresso.utils.auth import test_ssh_connection
            except Exception:
                def test_ssh_connection(username, host, key_path=None, port=22):
                    import subprocess
                    key_path = os.path.expanduser(key_path) if key_path else None
                    cmd = ["ssh", "-p", str(port), "-o", "PasswordAuthentication=no", "-o", "BatchMode=yes", "-o", "ConnectTimeout=5"]
                    if key_path:
                        cmd += ["-i", key_path]
                    cmd += [f"{username}@{host}", "echo 'Connection successful'"]
                    try:
                        subprocess.run(cmd, check=True, capture_output=True, text=True)
                        return True
                    except Exception:
                        return False

            try:
                success = test_ssh_connection(username, host, key_path, port)
                if success:
                    self.results_label.setText(f"✅ SSH connection successful!\nConnected to: {username}@{host}:{port}")
                    self.results_label.setStyleSheet('color: green;')
                else:
                    self.results_label.setText("❌ SSH connection failed.\n💡 Check your credentials and SSH key configuration")
                    self.results_label.setStyleSheet('color: red;')
            except Exception as e:
                self.results_label.setText(f"❌ Test failed: {e}")
                self.results_label.setStyleSheet('color: red;')

        try:
            self.test_conn_btn.clicked.connect(_test_conn)
        except Exception:
            pass

        mform.addRow(self.test_conn_btn)
        mform.addRow(self.results_label)
        tabs.addTab(machine, 'Machine')

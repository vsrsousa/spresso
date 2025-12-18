from qtpy.QtWidgets import QDialog, QVBoxLayout, QLabel, QTextEdit, QPushButton
from qtpy.QtCore import Qt, Signal
import threading
import tempfile

from xespresso.workflow.tasks import ScfTask, RelaxTask, WorkflowRunner
from qtgui.calculations.preparation import prepare_calculation_from_gui


class WorkflowRunDialog(QDialog):
    """Run a workflow (prototype). Runs tasks in a background thread and shows log/status.

    UI updates are marshalled to the Qt main thread via signals to avoid
    unsafe cross-thread widget access which can cause crashes.
    """

    log_signal = Signal(str)
    status_signal = Signal(str)

    def __init__(self, session, preset_name: str = "SCF", run_label: str = None, parent=None):
        super().__init__(parent)
        self.session = session
        self.preset_name = preset_name
        self.run_label = run_label or f"{preset_name.lower()}-run"
        self.setWindowTitle(f"Run: {self.run_label}")
        self.resize(600, 360)
        self._init_ui()
        # connect signals for thread-safe UI updates
        self.log_signal.connect(self._append_from_signal)
        self.status_signal.connect(self._set_status_from_signal)

        self._start_run()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        self.status = QLabel("Starting...")
        layout.addWidget(self.status)

        self.log = QTextEdit()
        self.log.setReadOnly(True)
        layout.addWidget(self.log)

        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.accept)
        layout.addWidget(self.close_btn, alignment=Qt.AlignRight)

    def _append(self, text: str):
        # Internal helper used on main thread; prefer using signals from worker
        self.log.append(text)

    def _append_from_signal(self, text: str):
        self.log.append(text)

    def _set_status_from_signal(self, text: str):
        self.status.setText(text)

    def _build_tasks(self):
        # Obtain atoms and workflow configuration from session state
        atoms = None
        try:
            atoms = self.session.get("current_structure")
        except Exception:
            atoms = None

        config = {}
        try:
            config = self.session.get("workflow_config") or {}
        except Exception:
            config = {}

        # Merge session-level provenance and queue settings into runner context
        queue = config.get("queue") or self.session.get("workflow_machine") or self.session.get("current_machine")

        # Common params applied to tasks
        common = {
            "label": self.run_label,
            "pseudopotentials": config.get("pseudopotentials", {}),
            "kpts": config.get("kpts"),
            "input_data": config.get("input_data", {}),
            "queue": queue,
            "debug": True,
        }

        tasks = []
        # Map human-friendly preset names to task sequences
        if self.preset_name == "SCF":
            params = dict(common)
            params.update({"atoms": atoms})
            tasks.append(ScfTask(name="scf", params=params))
        elif self.preset_name == "Relax" or self.preset_name == "Geometry Optimization":
            params = dict(common)
            params.update({"atoms": atoms})
            tasks.append(RelaxTask(name="relax", params=params))
        elif self.preset_name == "SCF+Relax":
            params1 = dict(common)
            params1.update({"label": self.run_label + "_scf", "atoms": atoms})
            params2 = dict(common)
            params2.update({"label": self.run_label + "_relax", "atoms": atoms})
            tasks.extend([ScfTask(name="scf", params=params1), RelaxTask(name="relax", params=params2)])
        elif self.preset_name == "Convergence Test":
            # Single convergence task which internally will run SCF steps
            from xespresso.workflow.tasks import ConvergenceTask
            params = dict(common)
            params.update({"atoms": atoms})
            tasks.append(ConvergenceTask(name="convergence", params=params))
        elif self.preset_name == "NEB":
            from xespresso.workflow.tasks import NebTask
            params = dict(common)
            params.update({"atoms": atoms})
            tasks.append(NebTask(name="neb", params=params))
        elif self.preset_name == "Post-Processing":
            from xespresso.workflow.tasks import PpTask
            params = dict(common)
            params.update({"atoms": atoms})
            tasks.append(PpTask(name="pp", params=params))
        else:
            # Unknown preset: fallback to a single SCF
            params = dict(common)
            params.update({"atoms": atoms})
            tasks.append(ScfTask(name="scf", params=params))

        return tasks

    def _start_run(self):
        self._append(f"Launching workflow preset: {self.preset_name}")
        tasks = self._build_tasks()
        runner = WorkflowRunner(tasks=tasks, context={})

        prov_dir = tempfile.mkdtemp(prefix="xespresso_prov_")
        runner.context["debug"] = True
        runner.context["provenance_dir"] = prov_dir

        def _target():
            try:
                self.log_signal.emit("Running tasks...")
                results = runner.run_all()
                self.log_signal.emit("Run finished.")
                for r in results:
                    self.log_signal.emit(str(r))
                # store run summary in session
                runs = None
                try:
                    runs = self.session.get("workflow_runs") or []
                except Exception:
                    runs = getattr(self.session, "workflow_runs", []) or []

                runs.append({"label": self.run_label, "preset": self.preset_name, "results": results})
                try:
                    self.session["workflow_runs"] = runs
                except Exception:
                    try:
                        setattr(self.session, "workflow_runs", runs)
                    except Exception:
                        pass
                self.status_signal.emit("Finished")
            except Exception as exc:
                self.log_signal.emit(f"Run failed: {exc}")
                self.status_signal.emit("Failed")

        th = threading.Thread(target=_target, daemon=True)
        th.start()

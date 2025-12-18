from qtpy.QtWidgets import QDialog, QVBoxLayout, QLabel, QComboBox, QLineEdit, QPushButton
from qtpy.QtCore import Qt

from .workflow_run import WorkflowRunDialog


class WorkflowLauncherDialog(QDialog):
    """Dialog to choose a workflow preset and launch it."""

    def __init__(self, session, workflow_name=None, parent=None):
        super().__init__(parent)
        self.session = session
        self.setWindowTitle("Launch Workflow")
        self.resize(420, 160)
        self._init_ui(workflow_name)

    def _init_ui(self, workflow_name):
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Workflow"))
        self.combo = QComboBox()
        presets = [
            "SCF",
            "Relax",
            "SCF+Relax",
            "Convergence Test",
            "Geometry Optimization",
            "NEB",
            "Post-Processing",
        ]
        self.combo.addItems(presets)
        if workflow_name and workflow_name in presets:
            self.combo.setCurrentText(workflow_name)
        layout.addWidget(self.combo)

        layout.addWidget(QLabel("Run label"))
        self.label_edit = QLineEdit("run-1")
        layout.addWidget(self.label_edit)

        self.launch_btn = QPushButton("Launch")
        self.launch_btn.clicked.connect(self._on_launch)
        layout.addWidget(self.launch_btn, alignment=Qt.AlignRight)

    def _on_launch(self):
        preset = self.combo.currentText()
        label = self.label_edit.text().strip() or None
        dlg = WorkflowRunDialog(self.session, preset_name=preset, run_label=label, parent=self)
        dlg.exec_()
        self.accept()

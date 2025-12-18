from qtpy import QtWidgets
from qtpy.QtWidgets import QWidget, QVBoxLayout, QListWidget, QPushButton, QLabel

from ..dialogs.workflow_launcher import WorkflowLauncherDialog


class WorkflowsPage(QWidget):
    """Simple workflows page that lists available workflows and lets user launch them."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)

        title = QLabel("Workflows")
        title.setObjectName("pageTitle")
        layout.addWidget(title)

        self.list = QListWidget()
        # For prototype we show a few built-in workflow names. Real app should
        # discover workflows from configured sources.
        presets = [
            "SCF",
            "Relax",
            "SCF+Relax",
            "Convergence Test",
            "Geometry Optimization",
            "NEB",
            "Post-Processing",
        ]
        self.list.addItems(presets)
        layout.addWidget(self.list)

        btn_layout = QtWidgets.QHBoxLayout()
        self.launch_btn = QPushButton("Launch")
        self.launch_btn.clicked.connect(self._on_launch)
        btn_layout.addWidget(self.launch_btn)

        layout.addLayout(btn_layout)
        layout.addStretch()

    def _on_launch(self):
        sel = self.list.currentItem()
        if not sel:
            return
        name = sel.text()
        dlg = WorkflowLauncherDialog(self.session, workflow_name=name, parent=self)
        dlg.exec_()

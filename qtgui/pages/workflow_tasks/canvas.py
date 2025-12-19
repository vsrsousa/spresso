try:
    from qtpy.QtWidgets import QWidget, QVBoxLayout, QListWidget, QPushButton
except Exception:
    QWidget = object


class WorkflowCanvas(QWidget):
    """Very small canvas/list view for workflow tasks.

    This placeholder stores a list of task descriptors and allows adding
    and removing items. Intended as a minimal, testable placeholder for
    a richer drag-and-drop canvas.
    """
    def __init__(self, parent=None):
        if hasattr(self.__class__, "__mro__"):
            super().__init__(parent)
        try:
            self.layout = QVBoxLayout()
            self.list = QListWidget()
            self.add_btn = QPushButton("Add Task")
            self.remove_btn = QPushButton("Remove Selected")
            self.layout.addWidget(self.list)
            self.layout.addWidget(self.add_btn)
            self.layout.addWidget(self.remove_btn)
            if hasattr(self, 'setLayout'):
                self.setLayout(self.layout)
            self.add_btn.clicked.connect(self._on_add)
            self.remove_btn.clicked.connect(self._on_remove)
        except Exception:
            pass

    def _on_add(self):
        try:
            self.list.addItem("New Task")
        except Exception:
            pass

    def _on_remove(self):
        try:
            idx = self.list.currentRow()
            if idx >= 0:
                self.list.takeItem(idx)
        except Exception:
            pass

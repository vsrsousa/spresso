try:
    from qtpy.QtWidgets import QWidget, QVBoxLayout, QListWidget, QPushButton
except Exception:
    QWidget = object


class TaskPalette(QWidget):
    """Simple task palette listing available task types.

    This is a lightweight placeholder used by the GUI. It can be extended
    to provide drag-and-drop creation of tasks onto a canvas.
    """
    def __init__(self, parent=None):
        if hasattr(self.__class__, "__mro__"):
            super().__init__(parent)
        try:
            self.layout = QVBoxLayout()
            self.list = QListWidget()
            self.list.addItems(["SCF", "Relax", "Convergence", "NEB", "PostProcessing"])
            self.new_btn = QPushButton("New Task")
            self.layout.addWidget(self.list)
            self.layout.addWidget(self.new_btn)
            if hasattr(self, 'setLayout'):
                self.setLayout(self.layout)
        except Exception:
            pass


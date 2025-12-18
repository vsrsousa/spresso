try:
    from PySide6.QtWidgets import QWidget, QVBoxLayout, QListWidget, QPushButton
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
try:
    from PySide6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QListWidget
except Exception:
    QWidget = object

class TaskPalette(QWidget):
    """Simple task palette placeholder for the GUI.

    Shows available task types and returns a factory when an item is selected.
    This is a lightweight placeholder so the GUI can import and instantiate
    the palette without requiring the full app wiring.
    """

    def __init__(self, parent=None):
        if hasattr(self.__class__, "__mro__"):
            super().__init__(parent)
            self.layout = QVBoxLayout()
            self.list = QListWidget()
            self.list.addItems(["SCF", "Relax", "Convergence", "NEB", "PP"])
            self.layout.addWidget(self.list)
            self.setLayout(self.layout)
        else:
            pass

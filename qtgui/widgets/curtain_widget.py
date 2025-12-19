from qtpy.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QToolButton, QLabel, QFrame
from qtpy.QtCore import Qt


class CurtainWidget(QWidget):
    """A simple collapsible curtain widget with a header and content area.

    Usage:
        c = CurtainWidget('Header')
        c.content_layout.addWidget(some_widget)
        parent_layout.addWidget(c)
    """

    def __init__(self, title: str = '', parent=None, expanded: bool = True):
        super().__init__(parent)
        self._expanded = bool(expanded)

        self.main = QVBoxLayout(self)
        self.main.setContentsMargins(0, 0, 0, 0)

        # Header
        header = QHBoxLayout()
        header.setContentsMargins(4, 2, 4, 2)
        self.toggle = QToolButton(text='')
        self.toggle.setCheckable(True)
        self.toggle.setChecked(self._expanded)
        self.toggle.setToolButtonStyle(Qt.ToolButtonIconOnly)
        self.toggle.setArrowType(Qt.DownArrow if self._expanded else Qt.RightArrow)
        self.toggle.clicked.connect(self._on_toggled)

        self.title = QLabel(str(title))
        self.title.setObjectName('curtainTitle')

        header.addWidget(self.toggle, 0, Qt.AlignLeft)
        header.addWidget(self.title, 1, Qt.AlignLeft)
        header.addStretch()

        header_widget = QWidget()
        header_widget.setLayout(header)
        header_widget.setObjectName('curtainHeader')

        self.main.addWidget(header_widget)

        # Content area
        self.content = QFrame()
        self.content.setFrameShape(QFrame.StyledPanel)
        self.content_layout = QVBoxLayout(self.content)
        self.content_layout.setContentsMargins(8, 4, 8, 8)
        self.main.addWidget(self.content)

        self._update_visibility()

    def _on_toggled(self):
        self._expanded = not self._expanded
        self.toggle.setArrowType(Qt.DownArrow if self._expanded else Qt.RightArrow)
        self._update_visibility()

    def _update_visibility(self):
        self.content.setVisible(self._expanded)

    def setTitle(self, title: str):
        self.title.setText(str(title))

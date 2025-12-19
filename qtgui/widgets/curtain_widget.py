from qtpy.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QToolButton, QLabel, QFrame
from qtpy.QtCore import Qt
from qtpy.QtGui import QFont


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
        # simple curtain: toggle then title label

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

        # Title label placed beside the toggle (compact, consistent placement)
        display = str(title)
        self.title_label = QLabel(display)
        self.title_label.setObjectName('curtainHeader')
        try:
            f = self.title_label.font()
            f.setBold(True)
            self.title_label.setFont(f)
        except Exception:
            pass

        # No icon label: only show the title text beside the toggle
        self.icon_label = None

        # Place the toggle at the left, then title, with stretch to the right.
        header.addWidget(self.toggle)
        header.addWidget(self.title_label)
        header.addStretch()
        self.main.addLayout(header)

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
        try:
            self.title_label.setText(str(title))
        except Exception:
            pass

    def setIcon(self, icon: str):
        # Icon support removed; no-op kept for API compatibility.
        return

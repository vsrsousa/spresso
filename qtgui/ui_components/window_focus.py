from qtpy.QtCore import QObject, QEvent
from qtpy.QtCore import QObject, QEvent
from qtpy.QtWidgets import QGraphicsDropShadowEffect
from qtpy.QtGui import QColor


class WindowFocusDecorator(QObject):
    """Event filter that applies a drop shadow when a window/dialog gains focus.

    Install on top-level widgets to give a visible highlight when active.
    This implementation avoids drawing any external overlay to maintain
    native window stacking behavior and compatibility.
    """

    def __init__(self, parent=None, active_color=QColor(37, 99, 235, 200)):
        super().__init__(parent)
        self._active_color = active_color
        self._effects = {}

    def eventFilter(self, obj, event):
        # Use WindowActivate and WindowDeactivate as reliable signals for focus
        if event.type() == QEvent.WindowActivate:
            self._apply_effect(obj)
        elif event.type() == QEvent.WindowDeactivate:
            self._remove_effect(obj)
        elif event.type() == QEvent.FocusIn:
            self._apply_effect(obj)
        elif event.type() == QEvent.FocusOut:
            self._remove_effect(obj)
        return super().eventFilter(obj, event)

    def _apply_effect(self, widget):
        try:
            if widget in self._effects:
                return
            effect = QGraphicsDropShadowEffect(widget)
            effect.setBlurRadius(30)
            effect.setOffset(0, 8)
            effect.setColor(self._active_color)
            widget.setGraphicsEffect(effect)
            self._effects[widget] = effect
        except Exception:
            pass

    def _remove_effect(self, widget):
        try:
            eff = self._effects.pop(widget, None)
            if eff is not None:
                widget.setGraphicsEffect(None)
        except Exception:
            pass


def apply_focus_decorator(widget, active_color=None):
    from qtpy.QtGui import QColor
    color = QColor(37, 99, 235, 200) if active_color is None else active_color
    decorator = WindowFocusDecorator(widget, active_color=color)
    widget.installEventFilter(decorator)
    # Keep reference to decorator on widget so it isn't garbage-collected
    setattr(widget, "_focus_decorator", decorator)
    # For already-active windows, apply immediately
    try:
        decorator._apply_effect(widget)
    except Exception:
        pass
    return decorator

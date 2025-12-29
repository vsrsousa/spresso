"""Window focus helper: non-blocking safe no-op decorator.

Older versions installed an application-wide event filter which could
perform work during `eventFilter` and block the Qt event loop. That
led to UI hangs and made debugging difficult. Replace the behavior with
an extremely lightweight helper that does not install filters or do
long-running work.
"""
from qtpy.QtWidgets import QWidget


def apply_focus_decorator(window: QWidget) -> QWidget:
	"""Return the window unchanged; keep this hook so callers can use it.

	This implementation intentionally avoids installing any event
	filters or running synchronous work on window activation.
	"""
	return window

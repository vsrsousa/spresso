"""Compatibility shim: re-export `SessionState` from the new module.

This file is intentionally small and preserves the original import
path `qtgui.session_state` so existing modules keep working while the
implementation lives in `qtgui.session.state_manager`.
"""

from qtgui.session.state_manager import SessionState, session_state

__all__ = ["SessionState", "session_state"]
"""Compatibility shim: re-export `SessionState` from the new module.

This file is intentionally small and preserves the original import
path `qtgui.session_state` so existing modules keep working while the
implementation lives in `qtgui.session.state_manager`.
"""

from qtgui.session.state_manager import SessionState, session_state

__all__ = ["SessionState", "session_state"]

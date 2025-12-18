"""Session state manager moved into a dedicated module.

This preserves the original `SessionState` API but centralizes the
implementation away from the large `qtgui/session_state.py` file so
we can iterate on splitting I/O and UI concerns safely.
"""
import os
import json
import re
from datetime import datetime


# Default session data directory
DEFAULT_SESSIONS_DIR = os.path.expanduser("~/.xespresso/sessions")

# Default database path for structures
DEFAULT_STRUCTURES_DB_PATH = os.path.expanduser("~/.xespresso/structures.db")

# Invalid characters in filenames (cross-platform)
INVALID_FILENAME_CHARS = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']


class SessionState:
    """
    Enhanced session state manager for application state.

    This is a near-exact port of the previous implementation but uses
    absolute imports when necessary so it can live in a subpackage.
    """
    _instance = None

    # Allowed keys for session state validation (security measure)
    ALLOWED_SESSION_KEYS = {
        'current_structure', 'current_machine', 'current_machine_name',
        'current_codes', 'selected_code_version', 'workflow_config',
        'workflow_tabs_config', 'open_workflow_tabs', 'workflow_runs',
        'working_directory', 'session_name', 'session_created',
        'session_modified', 'calc_machine', 'selected_machine',
        'selected_qe_version', 'structure_source', 'workflow_machine',
        'structure_file_path', 'structure_db_path',
        # provenance configuration keys
        'provenance_db_path', 'provenance_dir'
    }

    def __new__(cls, *, isolated=False):
        if isolated:
            instance = super().__new__(cls)
            instance._state = {}
            instance._sessions = {}
            instance._current_session_id = None
            instance._sessions_dir = DEFAULT_SESSIONS_DIR
            instance._listeners = []
            instance._active_session_names = []
            instance._notifying = False
            instance._structure_locked = False
            instance._last_error_message = None
            instance._isolated = True
            instance._initialize_defaults()
            instance._load_sessions_index()
            return instance

        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._state = {}
            cls._instance._sessions = {}
            cls._instance._current_session_id = None
            cls._instance._sessions_dir = DEFAULT_SESSIONS_DIR
            cls._instance._listeners = []
            cls._instance._active_session_names = []
            cls._instance._notifying = False
            cls._instance._structure_locked = False
            cls._instance._last_error_message = None
            cls._instance._isolated = False
            cls._instance._initialize_defaults()
            cls._instance._load_sessions_index()
        return cls._instance

    def _initialize_defaults(self):
        self._state['current_structure'] = None
        self._state['current_machine'] = None
        self._state['current_machine_name'] = None
        self._state['current_codes'] = None
        self._state['selected_code_version'] = None
        self._state['workflow_config'] = {}
        self._state['working_directory'] = os.path.expanduser("~")
        self._state['session_name'] = "Default Session"
        self._state['session_created'] = datetime.now().isoformat()
        self._state['session_modified'] = datetime.now().isoformat()

    def __getitem__(self, key):
        val = self._state.get(key)
        if val is not None:
            return val
        if self.is_isolated():
            return SessionState()._state.get(key)
        return None

    def __setitem__(self, key, value):
        if key == 'current_structure' and getattr(self, '_structure_locked', False):
            print("Warning: Structure is locked for this session; ignoring change")
            return

        self._state[key] = value
        if self.is_isolated() and key in (
            'current_structure', 'structure_source', 'structure_file_path', 'structure_db_path',
            'workflow_config', 'espresso_calculator', 'prepared_atoms'
        ):
            try:
                shared = SessionState()
                shared._state[key] = value
                shared._state['session_modified'] = self._state.get('session_modified', datetime.now().isoformat())
            except Exception:
                pass
        self._state['session_modified'] = datetime.now().isoformat()
        if self.is_isolated() and key in (
            'current_machine', 'current_machine_name',
            'current_codes', 'selected_code_version',
            'selected_qe_version'
        ):
            shared = SessionState()
            shared._state[key] = value
            shared._state['session_modified'] = self._state['session_modified']
        self._notify_listeners()

    def __contains__(self, key):
        return key in self._state

    def get(self, key, default=None):
        val = self._state.get(key)
        if val is not None:
            return self._state.get(key, default)
        if self.is_isolated():
            return SessionState()._state.get(key, default)
        return default

    def get_config_state(self):
        return SessionState()

    def is_isolated(self):
        return getattr(self, "_isolated", False)

    def clear_state(self):
        self._state = {}
        self._current_session_id = None
        self._initialize_defaults()
        self._notify_listeners()

    def keys(self):
        return self._state.keys()

    def items(self):
        return self._state.items()

    def add_listener(self, callback):
        if callback not in self._listeners:
            self._listeners.append(callback)

    def remove_listener(self, callback):
        if callback in self._listeners:
            self._listeners.remove(callback)

    def _notify_listeners(self):
        if self._notifying:
            return
        self._notifying = True
        try:
            for listener in self._listeners:
                try:
                    listener()
                except Exception:
                    pass
        finally:
            self._notifying = False

    def _load_sessions_index(self):
        from qtgui.session import persistence
        persistence.load_sessions_index(self)

    def _save_sessions_index(self):
        from qtgui.session import persistence
        persistence.save_sessions_index(self)

    def list_sessions(self):
        return dict(self._sessions)

    def get_active_session_names(self):
        return list(self._active_session_names)

    def add_active_session(self, session_name):
        if session_name and session_name not in self._active_session_names:
            self._active_session_names.append(session_name)

    def rename_active_session(self, old_name, new_name):
        if old_name in self._active_session_names:
            idx = self._active_session_names.index(old_name)
            self._active_session_names[idx] = new_name

    def get_sessions_dir(self):
        return self._sessions_dir

    def list_session_files(self):
        from qtgui.session import persistence
        return persistence.list_session_files(self)

    def load_session_from_file(self, file_path):
        from qtgui.session import persistence
        return persistence.load_session_from_file(self, file_path)

    def get_current_session_id(self):
        return self._current_session_id or "default"

    def reset(self):
        """Reset session to defaults (test-friendly API)."""
        self._state = {}
        self._current_session_id = "default"
        self._initialize_defaults()
        self._notify_listeners()

    def get_session_name(self):
        return self._state.get('session_name', 'Unnamed Session')

    def create_session(self, name):
        session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
        self.save_session()
        self._current_session_id = session_id
        self._initialize_defaults()
        self._state['session_name'] = name
        self._sessions[session_id] = {
            'name': name,
            'created': self._state['session_created'],
            'modified': self._state['session_modified']
        }
        self._save_sessions_index()
        return session_id

    def switch_session(self, session_id):
        if session_id not in self._sessions and session_id != 'default':
            return False
        self.save_session()
        self._current_session_id = session_id
        self._load_session(session_id)
        self._notify_listeners()
        return True

    def rename_session(self, new_name):
        old_name = self._state.get('session_name', '')
        self._state['session_name'] = new_name
        if self._current_session_id in self._sessions:
            self._sessions[self._current_session_id]['name'] = new_name
            self._save_sessions_index()
        if old_name and old_name != new_name:
            old_safe = old_name
            new_safe = new_name
            for char in INVALID_FILENAME_CHARS:
                old_safe = old_safe.replace(char, '_')
                new_safe = new_safe.replace(char, '_')
            old_path = os.path.join(self._sessions_dir, f"{old_safe}.json")
            new_path = os.path.join(self._sessions_dir, f"{new_safe}.json")
            if os.path.exists(old_path) and old_path != new_path:
                try:
                    with open(old_path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                    data['session_name'] = new_name
                    with open(new_path, 'w', encoding='utf-8') as f:
                        json.dump(data, f, indent=2)
                    os.remove(old_path)
                except Exception as e:
                    print(f"Warning: Could not rename session file: {e}")

    def delete_session(self, session_id):
        if session_id == self._current_session_id:
            return False
        if session_id in self._sessions:
            del self._sessions[session_id]
            self._save_sessions_index()
            session_path = os.path.join(self._sessions_dir, f"{session_id}.json")
            if os.path.exists(session_path):
                os.remove(session_path)
            return True
        return False

    def save_session(self, session_id=None):
        from qtgui.session import persistence
        return persistence.save_session(self, session_id)

    def is_structure_locked(self):
        return getattr(self, '_structure_locked', False)

    def get_last_error_message(self):
        return getattr(self, '_last_error_message', None)

    def set_structure_from_file(self, file_path):
        from qtgui.session import persistence
        return persistence.set_structure_from_file(self, file_path)

    def set_structure_from_database(self, db_id, db_path=None):
        from qtgui.session import persistence
        return persistence.set_structure_from_database(self, db_id, db_path)

    def add_structure_file_to_db(self, file_path, db_path=None):
        from qtgui.session import persistence
        return persistence.add_structure_file_to_db(self, file_path, db_path)

    def _load_session(self, session_id):
        from qtgui.session import persistence
        return persistence.load_session(self, session_id)

    def _restore_structure_from_source(self):
        from qtgui.session import persistence
        return persistence.restore_structure_from_source(self)


# Backwards compatible singleton
session_state = SessionState()

__all__ = ["SessionState", "session_state"]

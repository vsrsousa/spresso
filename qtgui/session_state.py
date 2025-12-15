"""
Session state management for xespresso GUI.

This module contains the `SessionState` class and the shared
`session_state` instance for backward compatibility.
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


_page_modules = {}


class SessionState:
    """
    Enhanced session state manager for application state.

    This provides a central place to store and access state across pages,
    with support for multiple sessions, saving/loading, and session switching.
    """
    _instance = None

    # Allowed keys for session state validation (security measure)
    ALLOWED_SESSION_KEYS = {
        'current_structure', 'current_machine', 'current_machine_name',
        'current_codes', 'selected_code_version', 'workflow_config',
        'working_directory', 'session_name', 'session_created',
        'session_modified', 'calc_machine', 'selected_machine',
        'selected_qe_version', 'structure_source', 'workflow_machine',
        'structure_file_path', 'structure_db_path'
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
        # If local value is present and not None, return it. If it's None and
        # this is an isolated session, fall back to the shared session value.
        val = self._state.get(key)
        # Treat empty mappings as 'not set' for certain keys so isolated
        # sessions can fall back to shared state values.
        if val is not None:
            if isinstance(val, dict) and len(val) == 0 and key in ('workflow_config',):
                pass
            else:
                return val
        if self.is_isolated():
            return SessionState()._state.get(key)
        return None

    def __setitem__(self, key, value):
        # Prevent changing the session structure if it has been locked
        if key == 'current_structure' and getattr(self, '_structure_locked', False):
            # Ignore attempts to change the structure once locked
            print("Warning: Structure is locked for this session; ignoring change")
            return

        self._state[key] = value
        # If this is an isolated session instance, mirror certain keys
        # (structure and prepared/workflow config) to the shared SessionState
        # so other UI instances can see them.
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
        # Prefer local non-None value, otherwise fall back to shared for isolated
        val = self._state.get(key)
        if val is not None:
            if isinstance(val, dict) and len(val) == 0 and key in ('workflow_config',):
                pass
            else:
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
        os.makedirs(self._sessions_dir, exist_ok=True)
        index_path = os.path.join(self._sessions_dir, "sessions_index.json")
        if os.path.exists(index_path):
            try:
                with open(index_path, 'r', encoding='utf-8') as f:
                    loaded_data = json.load(f)
                if isinstance(loaded_data, dict):
                    self._sessions = {k: v for k, v in loaded_data.items() if isinstance(k, str) and isinstance(v, dict)}
                else:
                    self._sessions = {}
            except Exception:
                self._sessions = {}
        else:
            self._sessions = {}

    def _save_sessions_index(self):
        os.makedirs(self._sessions_dir, exist_ok=True)
        index_path = os.path.join(self._sessions_dir, "sessions_index.json")
        try:
            with open(index_path, 'w', encoding='utf-8') as f:
                json.dump(self._sessions, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save sessions index: {e}")

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
        sessions = []
        if not os.path.isdir(self._sessions_dir):
            return sessions
        for filename in os.listdir(self._sessions_dir):
            if filename.endswith('.json') and filename != 'sessions_index.json':
                file_path = os.path.join(self._sessions_dir, filename)
                session_id = filename[:-5]
                session_name = filename[:-5]
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                        if isinstance(data, dict):
                            if '_session_id' in data:
                                session_id = data['_session_id']
                            if 'session_name' in data:
                                session_name = data['session_name']
                except Exception:
                    pass
                sessions.append((session_id, session_name, file_path))
        sessions.sort(key=lambda x: x[1].lower())
        return sessions

    def load_session_from_file(self, file_path):
        if not os.path.exists(file_path):
            return False
        # reset last error message
        self._last_error_message = None
        filename = os.path.basename(file_path)
        if not filename.endswith('.json'):
            return False
        session_name = filename[:-5]
        self.save_session()
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                loaded_state = json.load(f)
            if not isinstance(loaded_state, dict):
                raise ValueError("Invalid session data format")
            # Require session files to include structure metadata so we can
            # reliably restore session state. If none of the expected
            # structure-related keys are present, refuse to load the file.
            if not any(k in loaded_state for k in ('structure_source', 'structure_file_path', 'structure_db_path')):
                msg = "Session file lacks structure metadata; refusing to load."
                # Keep a last-error message for the caller to display in dialogs
                self._last_error_message = msg
                return False
            session_id = loaded_state.get('_session_id', session_name)
            self._current_session_id = session_id
            self._initialize_defaults()
            for key, value in loaded_state.items():
                if isinstance(key, str) and key in self.ALLOWED_SESSION_KEYS:
                    self._state[key] = value
            self._state['session_name'] = session_name
            self._restore_structure_from_source()
            # If a structure was restored from the file, lock it for this session
            if self._state.get('current_structure') is not None:
                self._structure_locked = True
            if session_id not in self._sessions:
                self._sessions[session_id] = {
                    'name': session_name,
                    'created': self._state.get('session_created', ''),
                    'modified': self._state.get('session_modified', '')
                }
                self._save_sessions_index()
            else:
                self._sessions[session_id]['name'] = session_name
                self._save_sessions_index()
            self._notify_listeners()
            return True
        except Exception as e:
            msg = f"Could not load session from file: {e}"
            self._last_error_message = msg
            self._initialize_defaults()
            return False

    def get_current_session_id(self):
        return self._current_session_id

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
        if session_id is None:
            session_id = self._current_session_id
        os.makedirs(self._sessions_dir, exist_ok=True)
        session_name = self._state.get('session_name', 'Unnamed')
        safe_filename = session_name
        for char in INVALID_FILENAME_CHARS:
            safe_filename = safe_filename.replace(char, '_')
        session_path = os.path.join(self._sessions_dir, f"{safe_filename}.json")
        serializable_state = {}
        for key, value in self._state.items():
            if key in ('current_structure', 'current_machine', 'current_codes'):
                continue
            try:
                # Try normal JSON serialization first
                json.dumps(value)
                serializable_state[key] = value
            except (TypeError, ValueError):
                try:
                    # Fallback: serialize using default=str to coerce non-serializables
                    dumped = json.dumps(value, default=str)
                    # Load back to Python object to store as JSON-serializable structure
                    serializable_state[key] = json.loads(dumped)
                except Exception:
                    # Last resort: store string representation
                    try:
                        serializable_state[key] = str(value)
                    except Exception:
                        pass
        serializable_state['_session_id'] = session_id
        try:
            with open(session_path, 'w', encoding='utf-8') as f:
                json.dump(serializable_state, f, indent=2)
            if session_id not in self._sessions:
                self._sessions[session_id] = {}
            self._sessions[session_id]['name'] = session_name
            self._sessions[session_id]['modified'] = datetime.now().isoformat()
            self._save_sessions_index()
        except Exception as e:
            print(f"Warning: Could not save session: {e}")

    def is_structure_locked(self):
        return getattr(self, '_structure_locked', False)

    def get_last_error_message(self):
        return getattr(self, '_last_error_message', None)

    def set_structure_from_file(self, file_path):
        """Load a structure from a file and set it as the session structure if not locked."""
        try:
            from .utils import read_structure
            atoms = read_structure(file_path)
        except Exception as e:
            print(f"Warning: Could not read structure file: {e}")
            return False

        if getattr(self, '_structure_locked', False):
            print("Warning: Structure is locked; cannot set from file")
            return False

        self._state['current_structure'] = atoms
        self._state['structure_source'] = f"File: {file_path}"
        self._state['structure_file_path'] = file_path
        self._structure_locked = True
        self._notify_listeners()
        return True

    def set_structure_from_database(self, db_id, db_path=None):
        """Set structure by database ID (restores atoms) and lock it."""
        if getattr(self, '_structure_locked', False):
            print("Warning: Structure is locked; cannot set from database")
            return False
        try:
            from ase.db import connect as ase_db_connect
            if not db_path:
                db_path = DEFAULT_STRUCTURES_DB_PATH
            db_path = os.path.abspath(os.path.expanduser(db_path))
            db = ase_db_connect(db_path)
            row = db.get(id=int(db_id))
            atoms = row.toatoms()
            self._state['current_structure'] = atoms
            self._state['structure_source'] = f"Database: ID {db_id}"
            self._state['structure_db_path'] = db_path
            self._structure_locked = True
            self._notify_listeners()
            return True
        except Exception as e:
            print(f"Warning: Could not set structure from database: {e}")
            return False

    def add_structure_file_to_db(self, file_path, db_path=None):
        """Read a structure file and save it into the ASE structures DB.

        Returns the inserted database id on success, or None on failure.
        """
        try:
            from .utils import read_structure
            atoms = read_structure(file_path)
        except Exception as e:
            self._last_error_message = f"Could not read structure file: {e}"
            return None

        try:
            from ase.db import connect as ase_db_connect
            if not db_path:
                db_path = DEFAULT_STRUCTURES_DB_PATH
            db_path = os.path.abspath(os.path.expanduser(db_path))
            os.makedirs(os.path.dirname(db_path), exist_ok=True)
            db = ase_db_connect(db_path)
            # write returns the id of the new row in many ASE versions
            row_id = db.write(atoms)
            return int(row_id) if row_id is not None else None
        except Exception as e:
            self._last_error_message = f"Could not save structure to DB: {e}"
            return None

    def _load_session(self, session_id):
        session_path = os.path.join(self._sessions_dir, f"{session_id}.json")
        if os.path.exists(session_path):
            try:
                with open(session_path, 'r', encoding='utf-8') as f:
                    loaded_state = json.load(f)
                if not isinstance(loaded_state, dict):
                    raise ValueError("Invalid session data format")
                self._initialize_defaults()
                for key, value in loaded_state.items():
                    if isinstance(key, str) and key in self.ALLOWED_SESSION_KEYS:
                        self._state[key] = value
                self._restore_structure_from_source()
            except Exception as e:
                print(f"Warning: Could not load session: {e}")
                self._initialize_defaults()
        else:
            self._initialize_defaults()

    def _restore_structure_from_source(self):
        structure_source = self._state.get('structure_source', '')
        if not structure_source:
            return
        try:
            from ase import io as ase_io
        except ImportError:
            print("Warning: ASE not available, cannot restore structure")
            return
        db_match = re.match(r'^Database:\s*ID\s*(\d+)$', structure_source)
        if db_match:
            try:
                from ase.db import connect as ase_db_connect
                structure_id = int(db_match.group(1))
                db_path = self._state.get('structure_db_path')
                if not db_path:
                    db_path = DEFAULT_STRUCTURES_DB_PATH
                db_path = os.path.abspath(os.path.expanduser(db_path))
                if os.path.exists(db_path):
                    db = ase_db_connect(db_path)
                    row = db.get(id=structure_id)
                    atoms = row.toatoms()
                    self._state['current_structure'] = atoms
                    print(f"Restored structure from database: ID {structure_id}")
                else:
                    print(f"Warning: Database not found at {db_path}")
            except Exception as e:
                print(f"Warning: Could not restore structure from database: {e}")
        elif structure_source.startswith("File: "):
            try:
                file_path = self._state.get('structure_file_path')
                if file_path and os.path.exists(file_path):
                    from .utils import read_structure
                    atoms = read_structure(file_path)
                    self._state['current_structure'] = atoms
                    print(f"Restored structure from file: {file_path}")
                else:
                    filename = structure_source.replace("File: ", "").strip()
                    print(f"Warning: Structure file path not saved or file not found. Original filename was: {filename}")
            except Exception as e:
                print(f"Warning: Could not restore structure from file: {e}")
        elif structure_source.startswith("Built: "):
            print(f"Info: Built structure '{structure_source}' cannot be automatically restored. Please rebuild the structure.")


# Global session state instance (backwards compatible)
session_state = SessionState()

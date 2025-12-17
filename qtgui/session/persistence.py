"""Persistence helpers for session state I/O and ASE DB operations.

These functions operate on a `SessionState` instance passed as the
first argument so they can be used from `state_manager` without causing
circular imports.
"""
import os
import json
import re
from datetime import datetime

DEFAULT_SESSIONS_DIR = os.path.expanduser("~/.xespresso/sessions")
DEFAULT_STRUCTURES_DB_PATH = os.path.expanduser("~/.xespresso/structures.db")
INVALID_FILENAME_CHARS = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']


def load_sessions_index(session):
    os.makedirs(session._sessions_dir, exist_ok=True)
    index_path = os.path.join(session._sessions_dir, "sessions_index.json")
    if os.path.exists(index_path):
        try:
            with open(index_path, 'r', encoding='utf-8') as f:
                loaded_data = json.load(f)
            if isinstance(loaded_data, dict):
                session._sessions = {k: v for k, v in loaded_data.items() if isinstance(k, str) and isinstance(v, dict)}
            else:
                session._sessions = {}
        except Exception:
            session._sessions = {}
    else:
        session._sessions = {}


def save_sessions_index(session):
    os.makedirs(session._sessions_dir, exist_ok=True)
    index_path = os.path.join(session._sessions_dir, "sessions_index.json")
    try:
        with open(index_path, 'w', encoding='utf-8') as f:
            json.dump(session._sessions, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not save sessions index: {e}")


def list_session_files(session):
    sessions = []
    if not os.path.isdir(session._sessions_dir):
        return sessions
    for filename in os.listdir(session._sessions_dir):
        if filename.endswith('.json') and filename != 'sessions_index.json':
            file_path = os.path.join(session._sessions_dir, filename)
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


def save_session(session, session_id=None):
    # Ensure we always have a session id when saving (use 'default' as fallback)
    if session_id is None:
        session_id = session._current_session_id or "default"
    os.makedirs(session._sessions_dir, exist_ok=True)
    session_name = session._state.get('session_name', 'Unnamed')
    safe_filename = session_name
    for char in INVALID_FILENAME_CHARS:
        safe_filename = safe_filename.replace(char, '_')
    session_path = os.path.join(session._sessions_dir, f"{safe_filename}.json")
    serializable_state = {}
    for key, value in session._state.items():
        if key in ('current_structure', 'current_machine', 'current_codes'):
            continue
        try:
            json.dumps(value)
            serializable_state[key] = value
        except (TypeError, ValueError):
            try:
                dumped = json.dumps(value, default=str)
                serializable_state[key] = json.loads(dumped)
            except Exception:
                try:
                    serializable_state[key] = str(value)
                except Exception:
                    pass
    serializable_state['_session_id'] = session_id
    try:
        with open(session_path, 'w', encoding='utf-8') as f:
            json.dump(serializable_state, f, indent=2)
        if session_id not in session._sessions:
            session._sessions[session_id] = {}
        session._sessions[session_id]['name'] = session_name
        session._sessions[session_id]['modified'] = datetime.now().isoformat()
        save_sessions_index(session)
    except Exception as e:
        print(f"Warning: Could not save session: {e}")


def load_session(session, session_id):
    # Treat falsy session_id as the default logical session id
    if not session_id:
        session_id = "default"

    # Sessions are stored with filenames based on the user-visible
    # session name, while the file contents carry a `_session_id` that
    # identifies the logical session. Locate the file by scanning the
    # sessions directory for a matching `_session_id` first, falling
    # back to a literal file named "{session_id}.json".
    session_path = None
    if os.path.isdir(session._sessions_dir):
        for filename in os.listdir(session._sessions_dir):
            if not filename.endswith('.json') or filename == 'sessions_index.json':
                continue
            file_path = os.path.join(session._sessions_dir, filename)
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                if isinstance(data, dict) and data.get('_session_id') == session_id:
                    session_path = file_path
                    break
            except Exception:
                continue

    # Fallback to direct file by id (legacy behaviour)
    if session_path is None:
        candidate = os.path.join(session._sessions_dir, f"{session_id}.json")
        if os.path.exists(candidate):
            session_path = candidate

    if session_path and os.path.exists(session_path):
        try:
            with open(session_path, 'r', encoding='utf-8') as f:
                loaded_state = json.load(f)
            if not isinstance(loaded_state, dict):
                raise ValueError("Invalid session data format")
            session._initialize_defaults()
            for key, value in loaded_state.items():
                if isinstance(key, str) and key in session.ALLOWED_SESSION_KEYS:
                    session._state[key] = value
            restore_structure_from_source(session)
        except Exception as e:
            print(f"Warning: Could not load session: {e}")
            session._initialize_defaults()
    else:
        session._initialize_defaults()


def load_session_from_file(session, file_path):
    if not os.path.exists(file_path):
        return False
    session._last_error_message = None
    filename = os.path.basename(file_path)
    if not filename.endswith('.json'):
        return False
    session_name = filename[:-5]
    save_session(session)
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            loaded_state = json.load(f)
        if not isinstance(loaded_state, dict):
            raise ValueError("Invalid session data format")
        # Older session files may not include structure metadata; accept
        # files that only contain basic session state (working directory,
        # name, timestamps, etc.). We don't require structure metadata to
        # proceed with loading the session into an isolated workspace.
        session_id = loaded_state.get('_session_id', session_name)
        session._current_session_id = session_id
        session._initialize_defaults()
        for key, value in loaded_state.items():
            if isinstance(key, str) and key in session.ALLOWED_SESSION_KEYS:
                session._state[key] = value
        session._state['session_name'] = session_name
        restore_structure_from_source(session)
        if session._state.get('current_structure') is not None:
            session._structure_locked = True
        if session_id not in session._sessions:
            session._sessions[session_id] = {
                'name': session_name,
                'created': session._state.get('session_created', ''),
                'modified': session._state.get('session_modified', '')
            }
            save_sessions_index(session)
        else:
            session._sessions[session_id]['name'] = session_name
            save_sessions_index(session)
        session._notify_listeners()
        return True
    except Exception as e:
        msg = f"Could not load session from file: {e}"
        session._last_error_message = msg
        session._initialize_defaults()
        return False


def restore_structure_from_source(session):
    structure_source = session._state.get('structure_source', '')
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
            db_path = session._state.get('structure_db_path')
            if not db_path:
                db_path = DEFAULT_STRUCTURES_DB_PATH
            db_path = os.path.abspath(os.path.expanduser(db_path))
            if os.path.exists(db_path):
                db = ase_db_connect(db_path)
                row = db.get(id=structure_id)
                atoms = row.toatoms()
                session._state['current_structure'] = atoms
                print(f"Restored structure from database: ID {structure_id}")
            else:
                print(f"Warning: Database not found at {db_path}")
        except Exception as e:
            print(f"Warning: Could not restore structure from database: {e}")
    elif structure_source.startswith("File: "):
        try:
            file_path = session._state.get('structure_file_path')
            if file_path and os.path.exists(file_path):
                from qtgui.utils import read_structure
                atoms = read_structure(file_path)
                session._state['current_structure'] = atoms
                print(f"Restored structure from file: {file_path}")
            else:
                filename = structure_source.replace("File: ", "").strip()
                print(f"Warning: Structure file path not saved or file not found. Original filename was: {filename}")
        except Exception as e:
            print(f"Warning: Could not restore structure from file: {e}")
    elif structure_source.startswith("Built: "):
        print(f"Info: Built structure '{structure_source}' cannot be automatically restored. Please rebuild the structure.")


def set_structure_from_file(session, file_path):
    try:
        from qtgui.utils import read_structure
        atoms = read_structure(file_path)
    except Exception as e:
        print(f"Warning: Could not read structure file: {e}")
        return False

    if getattr(session, '_structure_locked', False):
        print("Warning: Structure is locked; cannot set from file")
        return False

    session._state['current_structure'] = atoms
    session._state['structure_source'] = f"File: {file_path}"
    session._state['structure_file_path'] = file_path
    session._structure_locked = True
    session._notify_listeners()
    return True


def set_structure_from_database(session, db_id, db_path=None):
    if getattr(session, '_structure_locked', False):
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
        session._state['current_structure'] = atoms
        session._state['structure_source'] = f"Database: ID {db_id}"
        session._state['structure_db_path'] = db_path
        session._structure_locked = True
        session._notify_listeners()
        return True
    except Exception as e:
        print(f"Warning: Could not set structure from database: {e}")
        return False


def add_structure_file_to_db(session, file_path, db_path=None):
    try:
        from qtgui.utils import read_structure
        atoms = read_structure(file_path)
    except Exception as e:
        session._last_error_message = f"Could not read file: {e}"
        return None

    try:
        from ase.db import connect as ase_db_connect
        if not db_path:
            db_path = DEFAULT_STRUCTURES_DB_PATH
        db_path = os.path.abspath(os.path.expanduser(db_path))
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        db = ase_db_connect(db_path)
        row_id = db.write(atoms)
        return int(row_id) if row_id is not None else None
    except Exception as e:
        session._last_error_message = f"Could not save structure to DB: {e}"
        return None

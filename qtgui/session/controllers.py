"""Non-GUI session controllers: helpers for structure DB/file operations.

These functions encapsulate logic that does not require PySide6, making
`session_window.py` safer to import in headless environments.
"""
from typing import Optional, Tuple


def add_structure_file_and_load(session, file_path: str, db_path: Optional[str] = None) -> Tuple[Optional[int], bool, str]:
    """Try to add a structure file to the DB and load it into `session`.

    Returns a tuple (db_id or None, loaded_bool, message).
    """
    try:
        db_id = session.add_structure_file_to_db(file_path, db_path=db_path)
    except Exception as e:
        db_id = None

    if not db_id:
        try:
            ok = session.set_structure_from_file(file_path)
            if ok:
                return None, True, "Loaded from file"
            else:
                return None, False, session.get_last_error_message() or "Could not set structure from file"
        except Exception as e:
            return None, False, str(e)

    # if saved to DB, try to load from DB into session
    try:
        ok2 = session.set_structure_from_database(db_id)
        if ok2:
            return int(db_id), True, "Saved to DB and loaded"
        else:
            return int(db_id), False, session.get_last_error_message() or "Could not load saved structure from DB"
    except Exception as e:
        return int(db_id), False, str(e)


def load_structure_from_db(session, db_id: int, db_path: Optional[str] = None) -> Tuple[bool, str]:
    """Load a structure from DB id into session. Returns (ok, message)."""
    try:
        ok = session.set_structure_from_database(db_id, db_path=db_path)
        if ok:
            return True, "Loaded"
        return False, session.get_last_error_message() or "Could not load structure"
    except Exception as e:
        return False, str(e)


def list_db_rows(db_path: Optional[str]):
    """Return rows iterator/list from ASE DB or raise ImportError if ASE not present."""
    try:
        from ase.db import connect as ase_db_connect
    except Exception:
        raise
    if not db_path:
        from os import path
        db_path = path.expanduser("~/.xespresso/structures.db")
    db_path = os.path.abspath(os.path.expanduser(db_path))
    db = ase_db_connect(db_path)
    return list(db.select())

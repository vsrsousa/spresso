"""
Session management utilities for xespresso GUI.

Provides functionality to:
- Save current session state to a file
- Load session state from a file
- Reset session state
"""

import streamlit as st
import json
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List

# Default directory for session files
DEFAULT_SESSION_DIR = os.path.expanduser("~/.xespresso/sessions")


def _is_widget_key(key: str) -> bool:
    """
    Check if a session state key appears to be a widget key.
    
    Widget keys are created by Streamlit widgets and should not be
    saved/restored as they conflict with widget rendering.
    
    Args:
        key: Session state key to check
        
    Returns:
        True if the key appears to be a widget key
    """
    # Common widget key patterns used in the GUI
    widget_patterns = [
        '_quick_',      # Quick access buttons in directory browser
        '_up',          # Up/parent directory button
        '_subdir_selector',  # Subfolder selectbox
        '_enter_subdir',     # Enter subfolder button
        '_custom_path',      # Custom path text input
        '_go_custom',        # Go to custom path button
        '_create_dir',       # Create directory button
        '_tkinter_browse',   # Native file dialog browse button
        '_new',              # New session button
        '_save',             # Save session button
        '_rename_',          # Rename session buttons/inputs
        '_switch_',          # Switch session controls
        '_close_',           # Close session button
        '_load_',            # Load session button
        '_renaming',         # Renaming flag
    ]
    
    # Check if key matches any widget pattern
    for pattern in widget_patterns:
        if pattern in key:
            return True
    
    # Check for common widget suffixes
    # Button widgets typically end with _btn
    # Other interactive widgets may have similar patterns
    widget_suffixes = [
        '_btn',         # Button widgets (e.g., build_crystal_btn)
        '_button',      # Alternative button naming
    ]
    
    for suffix in widget_suffixes:
        if key.endswith(suffix):
            return True
    
    return False


def get_serializable_state(exclude_keys: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Extract serializable items from session state.
    
    Excludes machine and code configurations as these are persistent
    and stored separately in ~/.xespresso/
    
    Special handling for ASE Atoms objects - converts to JSON string format.
    
    Args:
        exclude_keys: List of keys to exclude from serialization
        
    Returns:
        Dictionary of serializable session state items
    """
    if exclude_keys is None:
        # Default keys to exclude (non-serializable or internal Streamlit keys)
        exclude_keys = [
            # Streamlit internal keys
            'FormSubmitter',
            'FileUploader',
            # Widget keys that will be recreated
            '_widget_state',
            # Machine and code configurations (persistent, not session-specific)
            'current_machine',  # Machine object - config stored in ~/.xespresso/machines/
            'current_codes',    # Codes config - stored in ~/.xespresso/codes/
            # Calculator objects (will be recreated from config)
            'espresso_calculator',
            'calc_machine',
            'workflow_machine',
        ]
    
    serializable_state = {}
    
    for key, value in st.session_state.items():
        # Skip excluded keys
        if key in exclude_keys:
            continue
            
        # Skip keys starting with underscore (usually internal)
        if key.startswith('_'):
            continue
        
        # Skip widget keys to avoid conflicts when restoring
        if _is_widget_key(key):
            continue
        
        # Special handling for ASE Atoms objects
        is_atoms = False
        try:
            from ase import Atoms
            from ase.io import write
            import io
            
            if isinstance(value, Atoms):
                is_atoms = True
                # Convert Atoms to JSON string
                sio = io.StringIO()
                write(sio, value, format='json')
                serializable_state[key] = {
                    '__type__': 'ase.Atoms',
                    '__data__': sio.getvalue()
                }
                continue
        except ImportError:
            # ASE not available, can't check if it's an Atoms object
            pass
        except Exception as e:
            # Error serializing Atoms - log and skip this key
            print(f"Warning: Could not serialize {key} as Atoms: {e}")
            continue
        
        # Try to serialize the value as regular JSON
        try:
            # Test if value is JSON serializable
            json.dumps(value)
            serializable_state[key] = value
        except (TypeError, ValueError):
            # Skip non-serializable values
            # Don't convert to string - this was causing the bug
            # where Atoms objects would become strings if serialization failed
            print(f"Warning: Skipping non-serializable key '{key}' of type {type(value)}")
            pass
    
    return serializable_state


def save_session(filename: Optional[str] = None, session_dir: Optional[str] = None, session_name: Optional[str] = None) -> str:
    """
    Save current session state to a JSON file.
    
    Args:
        filename: Name of the session file. If None, generates name from session_name or timestamp
        session_dir: Directory to save session files. If None, uses default
        session_name: Name of the session to store in metadata. If provided and filename is None,
                     this will be used as the base filename.
        
    Returns:
        Path to saved session file
        
    Raises:
        IOError: If unable to save session
    """
    if session_dir is None:
        session_dir = DEFAULT_SESSION_DIR
    
    # Create session directory if it doesn't exist
    os.makedirs(session_dir, exist_ok=True)
    
    # Generate filename if not provided
    if filename is None:
        if session_name:
            # Use session name as the base filename
            filename = f"{session_name.replace(' ', '_')}.json"
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"session_{timestamp}.json"
    
    # Ensure .json extension
    if not filename.endswith('.json'):
        filename += '.json'
    
    filepath = os.path.join(session_dir, filename)
    
    # Get serializable state
    state = get_serializable_state()
    
    # Add metadata including session name
    session_data = {
        "metadata": {
            "saved_at": datetime.now().isoformat(),
            "version": "1.0",
            "session_name": session_name or filename.replace('.json', '').replace('_', ' '),
        },
        "state": state
    }
    
    # Save to file
    with open(filepath, 'w') as f:
        json.dump(session_data, f, indent=2)
    
    return filepath


def load_session(filepath: str) -> tuple[Dict[str, Any], Optional[str]]:
    """
    Load session state from a JSON file.
    
    Args:
        filepath: Path to session file
        
    Returns:
        Tuple of (session state dictionary, session name from metadata)
        
    Raises:
        IOError: If unable to read file
        ValueError: If file format is invalid
    """
    if not os.path.exists(filepath):
        raise IOError(f"Session file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        session_data = json.load(f)
    
    # Validate format
    if "state" not in session_data:
        raise ValueError("Invalid session file format: missing 'state' key")
    
    # Extract session name from metadata
    session_name = None
    if "metadata" in session_data:
        session_name = session_data["metadata"].get("session_name")
    
    return session_data["state"], session_name


def restore_session(state: Dict[str, Any], clear_first: bool = True):
    """
    Restore session state from a dictionary.
    
    Handles special deserialization for ASE Atoms objects.
    
    Args:
        state: Dictionary of session state to restore
        clear_first: Whether to clear current session state first
    """
    if clear_first:
        # Clear current session state but keep some internal keys
        keys_to_keep = [k for k in st.session_state.keys() if k.startswith('_')]
        for key in list(st.session_state.keys()):
            if key not in keys_to_keep:
                del st.session_state[key]
    
    # Restore state, with special handling for certain types
    for key, value in state.items():
        # Skip widget keys - they should not be restored
        if _is_widget_key(key):
            continue
        
        # Check for special types that need deserialization
        if isinstance(value, dict) and '__type__' in value:
            if value['__type__'] == 'ase.Atoms':
                # Deserialize ASE Atoms object
                try:
                    from ase.io import read
                    import io
                    
                    sio = io.StringIO(value['__data__'])
                    atoms = read(sio, format='json')
                    st.session_state[key] = atoms
                    print(f"Successfully restored Atoms object for key '{key}': {atoms.get_chemical_formula()}")
                    continue
                except (ImportError, Exception) as e:
                    # If deserialization fails, log error and skip this key
                    print(f"ERROR: Could not deserialize {key} as Atoms object: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
        
        # Regular value - restore directly
        st.session_state[key] = value
    
    # Restore ESPRESSO_PSEUDO environment variable if it was saved
    if "_espresso_pseudo_path" in state:
        import os
        os.environ["ESPRESSO_PSEUDO"] = state["_espresso_pseudo_path"]
    
    # Ensure directory browser syncs with working_directory if it wasn't saved separately
    # This handles sessions saved before the browser key was added
    if 'working_directory' in st.session_state:
        workdir = st.session_state['working_directory']
        browser_key = 'workdir_browser_current_path'
        # If browser key wasn't in the saved state or is empty, sync it
        if browser_key not in state or not st.session_state.get(browser_key):
            st.session_state[browser_key] = workdir

def reset_session(keep_keys: Optional[List[str]] = None):
    """
    Reset session state, optionally keeping specified keys.
    
    Args:
        keep_keys: List of keys to keep. If None, clears everything
    """
    if keep_keys is None:
        keep_keys = []
    
    # Store values to keep
    kept_values = {key: st.session_state[key] for key in keep_keys if key in st.session_state}
    
    # Clear all session state
    for key in list(st.session_state.keys()):
        if not key.startswith('_'):  # Keep internal Streamlit keys
            del st.session_state[key]
    
    # Restore kept values
    for key, value in kept_values.items():
        st.session_state[key] = value


def list_sessions(session_dir: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    List available session files.
    
    Args:
        session_dir: Directory containing session files. If None, uses default
        
    Returns:
        List of dictionaries with session info (filename, path, saved_at)
    """
    if session_dir is None:
        session_dir = DEFAULT_SESSION_DIR
    
    if not os.path.exists(session_dir):
        return []
    
    sessions = []
    for filename in os.listdir(session_dir):
        if filename.endswith('.json'):
            filepath = os.path.join(session_dir, filename)
            try:
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    saved_at = data.get('metadata', {}).get('saved_at', 'Unknown')
                    sessions.append({
                        'filename': filename,
                        'path': filepath,
                        'saved_at': saved_at
                    })
            except:
                # Skip invalid files
                pass
    
    # Sort by saved_at (most recent first)
    sessions.sort(key=lambda x: x['saved_at'], reverse=True)
    
    return sessions


def get_active_sessions() -> Dict[str, Dict[str, Any]]:
    """
    Get all active sessions from session state.
    
    Returns:
        Dictionary mapping session IDs to session data
    """
    if '_active_sessions' not in st.session_state:
        st.session_state._active_sessions = {}
    return st.session_state._active_sessions


def get_current_session_id() -> str:
    """
    Get the current active session ID.
    
    Returns:
        Current session ID
    """
    if '_current_session_id' not in st.session_state:
        # Create first session automatically
        st.session_state._current_session_id = 'session_1'
        st.session_state._session_counter = 1
        
        # Initialize active sessions and add the first session
        if '_active_sessions' not in st.session_state:
            st.session_state._active_sessions = {}
        
        st.session_state._active_sessions['session_1'] = {
            'name': 'Session 1',
            'created_at': datetime.now().isoformat(),
            'state': {}
        }
    
    return st.session_state._current_session_id


def create_new_session() -> str:
    """
    Create a new calculation session.
    
    Returns:
        New session ID
    """
    # Save current session state before creating new one
    if '_active_sessions' in st.session_state and '_current_session_id' in st.session_state:
        current_id = st.session_state._current_session_id
        if current_id in st.session_state._active_sessions:
            st.session_state._active_sessions[current_id]['state'] = get_serializable_state()
    
    # Increment session counter
    if '_session_counter' not in st.session_state:
        st.session_state._session_counter = 1
    else:
        st.session_state._session_counter += 1
    
    # Create new session ID
    new_session_id = f"session_{st.session_state._session_counter}"
    
    # Initialize sessions dict if needed
    if '_active_sessions' not in st.session_state:
        st.session_state._active_sessions = {}
    
    # Create new session with empty state
    st.session_state._active_sessions[new_session_id] = {
        'name': f"Session {st.session_state._session_counter}",
        'created_at': datetime.now().isoformat(),
        'state': {}
    }
    
    # Switch to new session
    st.session_state._current_session_id = new_session_id
    
    # Clear current calculation state (but keep machine/code configs)
    # Reset to defaults for session-specific values
    keys_to_clear = [
        'current_structure', 'workflow_config',
        'selected_code_version', 'current_machine_name',
        'espresso_calculator', 'prepared_atoms'
    ]
    for key in keys_to_clear:
        if key in st.session_state:
            del st.session_state[key]
    
    # Initialize working_directory for the new session to default
    import os
    st.session_state.working_directory = os.path.expanduser("~")
    
    return new_session_id


def switch_session(session_id: str):
    """
    Switch to a different active session.
    
    Args:
        session_id: ID of session to switch to
    """
    if '_active_sessions' not in st.session_state:
        st.session_state._active_sessions = {}
    
    if session_id not in st.session_state._active_sessions:
        return
    
    # Save current session state before switching
    current_id = get_current_session_id()
    if current_id in st.session_state._active_sessions:
        st.session_state._active_sessions[current_id]['state'] = get_serializable_state()
    
    # Switch to new session
    st.session_state._current_session_id = session_id
    
    # Restore new session state
    session_data = st.session_state._active_sessions[session_id]
    if session_data.get('state'):
        restore_session(session_data['state'], clear_first=True)
    else:
        # If no saved state yet, initialize with defaults
        import os
        if 'working_directory' not in st.session_state:
            st.session_state.working_directory = os.path.expanduser("~")


def close_session(session_id: str):
    """
    Close an active session.
    
    Args:
        session_id: ID of session to close
    """
    if '_active_sessions' not in st.session_state:
        return
    
    if session_id in st.session_state._active_sessions:
        del st.session_state._active_sessions[session_id]
    
    # If closing current session, switch to another or create new
    if st.session_state._current_session_id == session_id:
        remaining = list(st.session_state._active_sessions.keys())
        if remaining:
            switch_session(remaining[0])
        else:
            create_new_session()


def rename_session(session_id: str, new_name: str):
    """
    Rename an active session.
    
    Args:
        session_id: ID of session to rename
        new_name: New name for the session
    """
    if '_active_sessions' not in st.session_state:
        return
    
    if session_id in st.session_state._active_sessions:
        st.session_state._active_sessions[session_id]['name'] = new_name


def render_session_manager(key: str = "session_manager"):
    """
    Render multi-session management UI component.
    
    Uses a compact selectbox design that scales well with many sessions.
    Shows current session info prominently, with actions and dropdown for switching.
    
    Args:
        key: Unique key for the component
    """
    st.sidebar.markdown("---")
    st.sidebar.subheader("📑 Sessions")
    
    # Get active sessions
    active_sessions = get_active_sessions()
    current_session_id = get_current_session_id()
    
    # Display current session info prominently
    if active_sessions and current_session_id in active_sessions:
        current_session = active_sessions[current_session_id]
        current_name = current_session.get('name', current_session_id)
        
        st.sidebar.markdown("**Current Session:**")
        st.sidebar.info(f"→ {current_name} ✓")
    
    # Session actions row
    col1, col2, col3 = st.sidebar.columns(3)
    with col1:
        if st.button("➕ New", key=f"{key}_new", use_container_width=True, 
                    help="Create a new calculation session"):
            new_id = create_new_session()
            st.rerun()
    
    with col2:
        if st.button("💾 Save", key=f"{key}_save", use_container_width=True, 
                    help="Save current session to file"):
            try:
                # Save current session state first
                if current_session_id in active_sessions:
                    active_sessions[current_session_id]['state'] = get_serializable_state()
                
                # Get session name
                session_name = active_sessions[current_session_id].get('name', 'Session')
                
                # Save with session name as filename (no timestamp)
                filepath = save_session(session_name=session_name)
                st.sidebar.success(f"✅ Saved to {os.path.basename(filepath)}")
            except Exception as e:
                st.sidebar.error(f"❌ Error: {e}")
    
    with col3:
        # Show rename button for current session
        if st.button("✏️ Rename", key=f"{key}_rename_current", use_container_width=True,
                    help="Rename current session"):
            st.session_state[f'{key}_renaming'] = True
            st.rerun()
    
    # Rename input for current session (if renaming)
    if st.session_state.get(f'{key}_renaming', False):
        current_name = active_sessions[current_session_id].get('name', current_session_id)
        new_name = st.sidebar.text_input(
            "New name:",
            value=current_name,
            key=f"{key}_new_name_input"
        )
        col1, col2 = st.sidebar.columns(2)
        with col1:
            if st.button("✓ OK", key=f"{key}_rename_ok", use_container_width=True):
                rename_session(current_session_id, new_name)
                st.session_state[f'{key}_renaming'] = False
                st.rerun()
        with col2:
            if st.button("✗ Cancel", key=f"{key}_rename_cancel", use_container_width=True):
                st.session_state[f'{key}_renaming'] = False
                st.rerun()
    
    # Switch to another session (only show if there are multiple sessions)
    if len(active_sessions) > 1:
        st.sidebar.markdown("---")
        st.sidebar.markdown("**Switch Session:**")
        
        # Create list of session options (exclude current session)
        session_options = []
        session_ids = []
        for sess_id, sess_data in active_sessions.items():
            if sess_id != current_session_id:
                name = sess_data.get('name', sess_id)
                session_options.append(name)
                session_ids.append(sess_id)
        
        if session_options:
            # Use selectbox for switching - more scalable than buttons
            selected_name = st.sidebar.selectbox(
                "Select session to switch to:",
                options=session_options,
                key=f"{key}_switch_select",
                help=f"Switch from '{current_name}' to another session"
            )
            
            col1, col2 = st.sidebar.columns(2)
            with col1:
                if st.button("🔄 Switch", key=f"{key}_switch_btn", use_container_width=True):
                    # Find the session ID for the selected name
                    idx = session_options.index(selected_name)
                    target_session_id = session_ids[idx]
                    switch_session(target_session_id)
                    st.rerun()
            
            with col2:
                # Close button for selected session
                if st.button("✖ Close", key=f"{key}_close_selected", use_container_width=True,
                           help=f"Close '{selected_name}'"):
                    idx = session_options.index(selected_name)
                    target_session_id = session_ids[idx]
                    close_session(target_session_id)
                    st.rerun()
        
        # Show session count
        st.sidebar.caption(f"📊 {len(active_sessions)} active session(s)")
    else:
        # Only one session
        st.sidebar.caption("💡 Create a new session to work on multiple calculations")
    
    # Load saved sessions
    with st.sidebar.expander("📂 Load Saved Session", expanded=False):
        sessions = list_sessions()
        
        if sessions:
            st.markdown("**Saved sessions:**")
            
            for session in sessions[:5]:  # Show last 5
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.caption(f"📄 {session['filename']}")
                    st.caption(f"🕐 {session['saved_at'][:19]}")
                with col2:
                    if st.button("Load", key=f"{key}_load_{session['filename']}", 
                               use_container_width=True):
                        try:
                            state, session_name = load_session(session['path'])
                            
                            # Create new session for loaded state
                            new_id = create_new_session()
                            
                            # Set session name from metadata (or filename as fallback)
                            if session_name:
                                active_sessions[new_id]['name'] = session_name
                            else:
                                # Fallback: use filename without .json
                                name = session['filename'].replace('.json', '').replace('_', ' ')
                                active_sessions[new_id]['name'] = name
                            
                            # Restore state
                            restore_session(state, clear_first=True)
                            active_sessions[new_id]['state'] = state
                            
                            st.success("✅ Session loaded!")
                            st.rerun()
                        except Exception as e:
                            st.error(f"❌ Error: {e}")
                st.markdown("---")
        else:
            st.info("No saved sessions found")

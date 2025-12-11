"""
PySide6 Main Application for xespresso GUI.

This module provides the main window and navigation for the xespresso
configuration interface using PySide6.

Performance optimizations:
- Lazy imports for page modules (loaded when first accessed)
- Efficient Qt6 event handling
"""

import sys
import os
import json
import re
from pathlib import Path
from datetime import datetime

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QStackedWidget, QListWidget, QListWidgetItem, QLabel, QGroupBox,
    QFileDialog, QMessageBox, QSplitter, QFrame, QPushButton,
    QStatusBar, QMenuBar, QMenu, QToolBar, QComboBox,
    QInputDialog, QLineEdit, QSizePolicy
)
from PySide6.QtCore import Qt, QSize, Signal
from PySide6.QtGui import QIcon, QFont, QAction, QScreen


# Default session data directory
DEFAULT_SESSIONS_DIR = os.path.expanduser("~/.xespresso/sessions")

# Default database path for structures
DEFAULT_STRUCTURES_DB_PATH = os.path.expanduser("~/.xespresso/structures.db")

# Invalid characters in filenames (cross-platform)
INVALID_FILENAME_CHARS = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']


# Lazy import helper for page modules (improves startup time)
_page_modules = {}

def _get_page_class(name):
    """Lazy import page class to improve startup time."""
    if name not in _page_modules:
        try:
            if name == 'StructureViewerPage':
                from qtgui.pages.structure_viewer import StructureViewerPage
                _page_modules[name] = StructureViewerPage
            elif name == 'CalculationSetupPage':
                from qtgui.pages.calculation_setup import CalculationSetupPage
                _page_modules[name] = CalculationSetupPage
            elif name == 'WorkflowBuilderPage':
                from qtgui.pages.workflow_builder import WorkflowBuilderPage
                _page_modules[name] = WorkflowBuilderPage
            elif name == 'JobSubmissionPage':
                from qtgui.pages.job_submission import JobSubmissionPage
                _page_modules[name] = JobSubmissionPage
            elif name == 'ResultsPostprocessingPage':
                from qtgui.pages.results_postprocessing import ResultsPostprocessingPage
                _page_modules[name] = ResultsPostprocessingPage
            else:
                raise ValueError(f"Unknown page class: {name}")
        except ImportError as e:
            raise ImportError(f"Failed to import page module '{name}': {e}") from e
    return _page_modules[name]


class SessionState:
    """
    Enhanced session state manager for application state.
    
    This provides a central place to store and access state across pages,
    with support for multiple sessions, saving/loading, and session switching.
    
    Features:
        - Multiple named sessions
        - Save/load sessions to disk
        - Session switching without data loss
        - Automatic session persistence
    """
    _instance = None
    
    # Allowed keys for session state validation (security measure)
    ALLOWED_SESSION_KEYS = {
        'current_structure', 'current_machine', 'current_machine_name',
        'current_codes', 'selected_code_version', 'workflow_config',
        'working_directory', 'session_name', 'session_created',
        'session_modified', 'calc_machine', 'selected_machine',
        'selected_qe_version', 'structure_source', 'workflow_machine',
        'structure_file_path', 'structure_db_path'  # For restoring structures on session load
    }
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._state = {}
            cls._instance._sessions = {}
            cls._instance._current_session_id = None  # No session active initially
            cls._instance._sessions_dir = DEFAULT_SESSIONS_DIR
            cls._instance._listeners = []
            # Track sessions created/loaded in this app session (for combo box)
            cls._instance._active_session_names = []  # List of session names shown in combo
            # Guard against recursive notifications - this prevents infinite loops
            # when listeners update state during notification, which would otherwise
            # trigger another notification cycle
            cls._instance._notifying = False
            cls._instance._initialize_defaults()
            cls._instance._load_sessions_index()
        return cls._instance
    
    def _initialize_defaults(self):
        """Initialize default state values for current session."""
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
        return self._state.get(key)
    
    def __setitem__(self, key, value):
        self._state[key] = value
        self._state['session_modified'] = datetime.now().isoformat()
        self._notify_listeners()
    
    def __contains__(self, key):
        return key in self._state
    
    def get(self, key, default=None):
        return self._state.get(key, default)
    
    def clear_state(self):
        """Clear all state and reinitialize defaults."""
        self._state = {}
        self._current_session_id = None
        self._initialize_defaults()
        self._notify_listeners()
    
    def keys(self):
        """Return all keys in the state."""
        return self._state.keys()
    
    def items(self):
        """Return all items in the state."""
        return self._state.items()
    
    def add_listener(self, callback):
        """Add a callback to be notified when state changes."""
        if callback not in self._listeners:
            self._listeners.append(callback)
    
    def remove_listener(self, callback):
        """Remove a state change listener."""
        if callback in self._listeners:
            self._listeners.remove(callback)
    
    def _notify_listeners(self):
        """Notify all listeners of state change.
        
        Uses a guard to prevent recursive notifications which could occur
        if a listener updates the state during notification.
        """
        # Prevent recursive notifications
        if self._notifying:
            return
        
        self._notifying = True
        try:
            for listener in self._listeners:
                try:
                    listener()
                except Exception:
                    pass  # Don't let listener errors crash the app
        finally:
            self._notifying = False
    
    def _load_sessions_index(self):
        """Load the index of available sessions."""
        os.makedirs(self._sessions_dir, exist_ok=True)
        index_path = os.path.join(self._sessions_dir, "sessions_index.json")
        
        if os.path.exists(index_path):
            try:
                with open(index_path, 'r', encoding='utf-8') as f:
                    loaded_data = json.load(f)
                # Validate the loaded data is a dictionary with expected structure
                if isinstance(loaded_data, dict):
                    self._sessions = {
                        k: v for k, v in loaded_data.items()
                        if isinstance(k, str) and isinstance(v, dict)
                    }
                else:
                    self._sessions = {}
            except Exception:
                self._sessions = {}
        else:
            self._sessions = {}
    
    def _save_sessions_index(self):
        """Save the sessions index to disk."""
        os.makedirs(self._sessions_dir, exist_ok=True)
        index_path = os.path.join(self._sessions_dir, "sessions_index.json")
        
        try:
            with open(index_path, 'w', encoding='utf-8') as f:
                json.dump(self._sessions, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save sessions index: {e}")
    
    def list_sessions(self):
        """
        List all available sessions.
        
        Returns:
            dict: Mapping of session_id to session metadata
        """
        return dict(self._sessions)
    
    def get_active_session_names(self):
        """
        Get list of session names that are currently active in the combo box.
        
        These are sessions that the user has created or loaded in this app session.
        
        Returns:
            list: List of session names (without .json extension)
        """
        return list(self._active_session_names)
    
    def add_active_session(self, session_name):
        """
        Add a session name to the active sessions list (for combo box).
        
        Args:
            session_name: Name of the session (without .json extension)
        """
        if session_name and session_name not in self._active_session_names:
            self._active_session_names.append(session_name)
    
    def rename_active_session(self, old_name, new_name):
        """
        Rename a session in the active sessions list.
        
        Args:
            old_name: Old session name
            new_name: New session name
        """
        if old_name in self._active_session_names:
            idx = self._active_session_names.index(old_name)
            self._active_session_names[idx] = new_name
    
    def get_sessions_dir(self):
        """
        Get the directory where sessions are stored.
        
        Returns:
            str: Path to the sessions directory
        """
        return self._sessions_dir
    
    def list_session_files(self):
        """
        List all session JSON files in the sessions directory.
        
        This scans the actual filesystem rather than relying on the index,
        ensuring all saved sessions are discoverable.
        
        The session ID is read from inside the JSON file (stored as '_session_id'),
        and the filename is the session name.
        
        Returns:
            list: List of tuples (session_id, session_name, file_path)
        """
        sessions = []
        if not os.path.isdir(self._sessions_dir):
            return sessions
        
        for filename in os.listdir(self._sessions_dir):
            if filename.endswith('.json') and filename != 'sessions_index.json':
                file_path = os.path.join(self._sessions_dir, filename)
                # Default session_id to filename without .json extension
                session_id = filename[:-5]
                # Default session_name to filename without .json extension
                session_name = filename[:-5]
                
                # Try to read the session ID and name from the file
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                        if isinstance(data, dict):
                            # Get session ID from inside the file
                            if '_session_id' in data:
                                session_id = data['_session_id']
                            # Get session name from inside the file
                            if 'session_name' in data:
                                session_name = data['session_name']
                except Exception:
                    pass
                
                sessions.append((session_id, session_name, file_path))
        
        # Sort by session name
        sessions.sort(key=lambda x: x[1].lower())
        return sessions
    
    def load_session_from_file(self, file_path):
        """
        Load a session directly from a file path.
        
        The session name becomes the filename (without .json extension).
        This is what appears in the "Active:" dropdown.
        
        Args:
            file_path (str): Path to the session JSON file
            
        Returns:
            bool: True if loading was successful
        """
        if not os.path.exists(file_path):
            return False
        
        # Extract filename
        filename = os.path.basename(file_path)
        if not filename.endswith('.json'):
            return False
        
        # Session name is the filename without .json extension
        session_name = filename[:-5]
        
        # Save current session first
        self.save_session()
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                loaded_state = json.load(f)
            
            # Validate loaded data is a dictionary
            if not isinstance(loaded_state, dict):
                raise ValueError("Invalid session data format")
            
            # Get session ID from inside the file, or fall back to filename
            session_id = loaded_state.get('_session_id', session_name)
            
            # Set new session ID
            self._current_session_id = session_id
            
            # Initialize defaults first, then override with validated values
            self._initialize_defaults()
            for key, value in loaded_state.items():
                # Only load keys that are strings and in allowed set
                if isinstance(key, str) and key in self.ALLOWED_SESSION_KEYS:
                    self._state[key] = value
            
            # IMPORTANT: Set the session name to the filename (without .json)
            # This ensures the "Active:" dropdown shows the correct name
            self._state['session_name'] = session_name
            
            # Try to restore the structure from the saved source
            self._restore_structure_from_source()
            
            # Register in index if not already
            if session_id not in self._sessions:
                self._sessions[session_id] = {
                    'name': session_name,
                    'created': self._state.get('session_created', ''),
                    'modified': self._state.get('session_modified', '')
                }
                self._save_sessions_index()
            else:
                # Update the name in case it differs
                self._sessions[session_id]['name'] = session_name
                self._save_sessions_index()
            
            self._notify_listeners()
            return True
            
        except Exception as e:
            print(f"Warning: Could not load session from file: {e}")
            self._initialize_defaults()
            return False
    
    def get_current_session_id(self):
        """Get the ID of the current session."""
        return self._current_session_id
    
    def get_session_name(self):
        """Get the name of the current session."""
        return self._state.get('session_name', 'Unnamed Session')
    
    def create_session(self, name):
        """
        Create a new session.
        
        Args:
            name: Name for the new session
            
        Returns:
            str: ID of the new session
        """
        # Generate unique session ID using microseconds for uniqueness
        session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}"
        
        # Save current session first
        self.save_session()
        
        # Create new session
        self._current_session_id = session_id
        self._initialize_defaults()
        self._state['session_name'] = name
        
        # Register in index
        self._sessions[session_id] = {
            'name': name,
            'created': self._state['session_created'],
            'modified': self._state['session_modified']
        }
        self._save_sessions_index()
        
        return session_id
    
    def switch_session(self, session_id):
        """
        Switch to a different session.
        
        Args:
            session_id: ID of the session to switch to
            
        Returns:
            bool: True if switch was successful
        """
        if session_id not in self._sessions and session_id != 'default':
            return False
        
        # Save current session
        self.save_session()
        
        # Load new session
        self._current_session_id = session_id
        self._load_session(session_id)
        self._notify_listeners()
        
        return True
    
    def rename_session(self, new_name):
        """
        Rename the current session and its JSON file.
        
        Args:
            new_name: New name for the session
        """
        old_name = self._state.get('session_name', '')
        
        # Update in-memory state
        self._state['session_name'] = new_name
        if self._current_session_id in self._sessions:
            self._sessions[self._current_session_id]['name'] = new_name
            self._save_sessions_index()
        
        # Rename the JSON file if it exists
        if old_name and old_name != new_name:
            # Sanitize filenames
            old_safe = old_name
            new_safe = new_name
            for char in INVALID_FILENAME_CHARS:
                old_safe = old_safe.replace(char, '_')
                new_safe = new_safe.replace(char, '_')
            
            old_path = os.path.join(self._sessions_dir, f"{old_safe}.json")
            new_path = os.path.join(self._sessions_dir, f"{new_safe}.json")
            
            if os.path.exists(old_path) and old_path != new_path:
                try:
                    # Read the old file
                    with open(old_path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                    
                    # Update session name in the data
                    data['session_name'] = new_name
                    
                    # Write to new file
                    with open(new_path, 'w', encoding='utf-8') as f:
                        json.dump(data, f, indent=2)
                    
                    # Remove old file
                    os.remove(old_path)
                except Exception as e:
                    print(f"Warning: Could not rename session file: {e}")
    
    def delete_session(self, session_id):
        """
        Delete a session.
        
        Args:
            session_id: ID of the session to delete
            
        Returns:
            bool: True if deletion was successful
        """
        if session_id == self._current_session_id:
            return False  # Can't delete current session
        
        if session_id in self._sessions:
            del self._sessions[session_id]
            self._save_sessions_index()
            
            # Delete session file
            session_path = os.path.join(self._sessions_dir, f"{session_id}.json")
            if os.path.exists(session_path):
                os.remove(session_path)
            
            return True
        return False
    
    def save_session(self, session_id=None):
        """
        Save the current session to disk.
        
        The file is saved using the session name (e.g., "My Session.json").
        The session ID is stored inside the JSON file.
        
        Args:
            session_id: Optional specific session ID to save to
        """
        if session_id is None:
            session_id = self._current_session_id
        
        os.makedirs(self._sessions_dir, exist_ok=True)
        
        # Use session name for the filename, not the session ID
        session_name = self._state.get('session_name', 'Unnamed')
        # Sanitize filename: replace characters that are invalid in filenames
        safe_filename = session_name
        for char in INVALID_FILENAME_CHARS:
            safe_filename = safe_filename.replace(char, '_')
        session_path = os.path.join(self._sessions_dir, f"{safe_filename}.json")
        
        # Prepare serializable state (exclude non-serializable objects)
        serializable_state = {}
        for key, value in self._state.items():
            if key in ('current_structure', 'current_machine', 'current_codes'):
                continue  # Skip non-serializable objects
            try:
                json.dumps(value)  # Test if serializable
                serializable_state[key] = value
            except (TypeError, ValueError):
                pass
        
        # Store the session ID inside the JSON file
        serializable_state['_session_id'] = session_id
        
        try:
            with open(session_path, 'w', encoding='utf-8') as f:
                json.dump(serializable_state, f, indent=2)
            
            # Update sessions index
            if session_id not in self._sessions:
                self._sessions[session_id] = {}
            self._sessions[session_id]['name'] = session_name
            self._sessions[session_id]['modified'] = datetime.now().isoformat()
            self._save_sessions_index()
            
        except Exception as e:
            print(f"Warning: Could not save session: {e}")
    
    def _load_session(self, session_id):
        """
        Load a session from disk.
        
        Args:
            session_id: ID of the session to load
        """
        session_path = os.path.join(self._sessions_dir, f"{session_id}.json")
        
        if os.path.exists(session_path):
            try:
                with open(session_path, 'r', encoding='utf-8') as f:
                    loaded_state = json.load(f)
                
                # Validate loaded data is a dictionary
                if not isinstance(loaded_state, dict):
                    raise ValueError("Invalid session data format")
                
                # Initialize defaults first, then override with validated values
                self._initialize_defaults()
                for key, value in loaded_state.items():
                    # Only load keys that are strings and in allowed set
                    if isinstance(key, str) and key in self.ALLOWED_SESSION_KEYS:
                        self._state[key] = value
                
                # Try to restore the structure from the saved source
                self._restore_structure_from_source()
                    
            except Exception as e:
                print(f"Warning: Could not load session: {e}")
                self._initialize_defaults()
        else:
            self._initialize_defaults()
    
    def _restore_structure_from_source(self):
        """
        Restore the structure from the saved structure_source.
        
        This method attempts to reload the structure based on the source type:
        - "Database: ID X" - loads from ASE database
        - "File: <path>" - loads from file if structure_file_path is saved
        - "Built: ..." - cannot be automatically restored
        """
        structure_source = self._state.get('structure_source', '')
        if not structure_source:
            return
        
        try:
            from ase import io as ase_io
        except ImportError:
            print("Warning: ASE not available, cannot restore structure")
            return
        
        # Try to restore from database using regex for robust parsing
        db_match = re.match(r'^Database:\s*ID\s*(\d+)$', structure_source)
        if db_match:
            try:
                from ase.db import connect as ase_db_connect
                
                # Extract the ID using regex match
                structure_id = int(db_match.group(1))
                
                # Use saved db_path or default
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
        
        # Try to restore from file
        elif structure_source.startswith("File: "):
            try:
                # Try to use the saved full path
                file_path = self._state.get('structure_file_path')
                if file_path and os.path.exists(file_path):
                    from .utils import read_structure
                    atoms = read_structure(file_path)
                    self._state['current_structure'] = atoms
                    print(f"Restored structure from file: {file_path}")
                else:
                    # File path not saved or file doesn't exist
                    filename = structure_source.replace("File: ", "").strip()
                    print(f"Warning: Structure file path not saved or file not found. "
                          f"Original filename was: {filename}")
                    
            except Exception as e:
                print(f"Warning: Could not restore structure from file: {e}")
        
        # Built structures cannot be automatically restored
        elif structure_source.startswith("Built: "):
            print(f"Info: Built structure '{structure_source}' cannot be automatically restored. "
                  f"Please rebuild the structure.")
    
    def reset(self):
        """Reset the current session to defaults."""
        self._initialize_defaults()
        self._notify_listeners()


# Global session state instance
session_state = SessionState()


class MainWindow(QMainWindow):
    """Main application window for xespresso PySide6 GUI."""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("⚛️ xespresso - Quantum ESPRESSO Configuration GUI")
        
        # Set window size as a proportion of the screen (80% width, 80% height)
        # with a reasonable minimum for usability
        screen = QApplication.primaryScreen().availableGeometry()
        width = max(int(screen.width() * 0.8), 800)
        height = max(int(screen.height() * 0.8), 600)
        self.resize(width, height)
        
        # Set minimum size to ensure usability on small screens
        self.setMinimumSize(800, 600)
        
        # Center the window on the screen
        self.move(
            (screen.width() - width) // 2,
            (screen.height() - height) // 2
        )
        
        # Initialize session state - use singleton directly for robustness
        self.session_state = SessionState()
        
        # Configuration dialog (created on demand)
        self._config_dialog = None
        
        # Job monitor dialog (created on demand, accessible without session)
        self._job_monitor = None
        
        # Guard to prevent recursive updates during session changes
        self._updating = False
        
        # Setup UI
        self._setup_ui()
        self._setup_menu()
        self._setup_toolbar()
        self._setup_statusbar()
        
        # Connect session state listener for UI updates
        self.session_state.add_listener(self._on_session_changed)
        
        # Select default page (Structure)
        # Find the Structure item in the navigation list
        for i in range(self.nav_list.count()):
            if "Structure" in self.nav_list.item(i).text():
                self.nav_list.setCurrentRow(i)
                break
    
    def _setup_ui(self):
        """Setup the main user interface."""
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout with splitter
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(5, 5, 5, 5)
        
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left sidebar
        sidebar = self._create_sidebar()
        splitter.addWidget(sidebar)
        
        # Right content area
        self.content_stack = QStackedWidget()
        self._create_pages()
        splitter.addWidget(self.content_stack)
        
        # Set splitter sizes (sidebar:content = 1:4)
        splitter.setSizes([250, 950])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
    
    def _create_sidebar(self):
        """Create the navigation sidebar."""
        sidebar = QWidget()
        sidebar.setMaximumWidth(300)
        sidebar.setMinimumWidth(200)
        layout = QVBoxLayout(sidebar)
        layout.setContentsMargins(5, 5, 5, 5)
        
        # Title
        title_label = QLabel("⚛️ xespresso")
        title_font = QFont()
        title_font.setPointSize(16)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(title_label)
        
        # Configuration button - opens non-blocking dialog
        config_group = QGroupBox("⚙️ Configuration")
        config_layout = QVBoxLayout(config_group)
        
        config_info = QLabel("Configure machines, codes, and pseudopotentials")
        config_info.setWordWrap(True)
        config_info.setStyleSheet("color: gray;")
        config_layout.addWidget(config_info)
        
        config_btn = QPushButton("⚙️ Open Configuration...")
        config_btn.setToolTip("Open configuration dialog for machines, codes, and pseudopotentials")
        config_btn.clicked.connect(self._open_config_dialog)
        config_layout.addWidget(config_btn)
        
        layout.addWidget(config_group)
        
        # Session section - improved with session management
        session_group = QGroupBox("📋 Session")
        session_layout = QVBoxLayout(session_group)
        
        # Session selector - starts empty, populated when user creates/loads sessions
        session_selector_layout = QHBoxLayout()
        session_selector_layout.addWidget(QLabel("Active:"))
        self.session_combo = QComboBox()
        self.session_combo.setToolTip("Active session (create new or load saved session)")
        self.session_combo.currentTextChanged.connect(self._on_session_selected)
        session_selector_layout.addWidget(self.session_combo, 1)
        session_layout.addLayout(session_selector_layout)
        
        # Session name display
        self.session_name_label = QLabel("")
        self.session_name_label.setWordWrap(True)
        self._update_session_name_label()
        session_layout.addWidget(self.session_name_label)
        
        # Session buttons - row 1: New, Load
        session_btn_layout1 = QHBoxLayout()
        
        new_session_btn = QPushButton("New")
        new_session_btn.setToolTip("Create a new session")
        new_session_btn.clicked.connect(self._new_session)
        session_btn_layout1.addWidget(new_session_btn)
        
        load_session_btn = QPushButton("Load")
        load_session_btn.setToolTip("Load a saved session")
        load_session_btn.clicked.connect(self._load_session_dialog)
        session_btn_layout1.addWidget(load_session_btn)
        
        session_layout.addLayout(session_btn_layout1)
        
        # Session buttons - row 2: Rename, Save
        session_btn_layout2 = QHBoxLayout()
        
        rename_session_btn = QPushButton("Rename")
        rename_session_btn.setToolTip("Rename current session")
        rename_session_btn.clicked.connect(self._rename_session)
        session_btn_layout2.addWidget(rename_session_btn)
        
        save_session_btn = QPushButton("Save")
        save_session_btn.setToolTip("Save current session")
        save_session_btn.clicked.connect(self._save_session)
        session_btn_layout2.addWidget(save_session_btn)
        
        session_layout.addLayout(session_btn_layout2)
        
        # Session buttons - row 3: Close
        session_btn_layout3 = QHBoxLayout()
        
        close_session_btn = QPushButton("Close")
        close_session_btn.setToolTip("Close current session and start fresh")
        close_session_btn.clicked.connect(self._close_session)
        session_btn_layout3.addWidget(close_session_btn)
        
        session_layout.addLayout(session_btn_layout3)
        
        layout.addWidget(session_group)
        
        # Session combo starts empty - don't populate until user creates/loads session
        # self._refresh_session_list()  # Removed - combo starts empty
        
        # Working directory section
        workdir_group = QGroupBox("📁 Working Directory")
        workdir_layout = QVBoxLayout(workdir_group)
        
        self.workdir_label = QLabel(self.session_state.get('working_directory', '~'))
        self.workdir_label.setWordWrap(True)
        workdir_layout.addWidget(self.workdir_label)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_workdir)
        workdir_layout.addWidget(browse_btn)
        layout.addWidget(workdir_group)
        
        # Workflow section - only workflow pages (no config pages)
        workflow_group = QGroupBox("🔬 Workflow")
        workflow_layout = QVBoxLayout(workflow_group)
        
        self.nav_list = QListWidget()
        nav_items = [
            "🔬 Structure",
            "📊 Calculation Setup",
            "🔄 Workflow Builder",
            "🚀 Job Submission",
            "📈 Results & Post-Processing"
        ]
        for item in nav_items:
            self.nav_list.addItem(QListWidgetItem(item))
        self.nav_list.currentRowChanged.connect(self._on_nav_changed)
        workflow_layout.addWidget(self.nav_list)
        layout.addWidget(workflow_group)
        
        # Spacer
        layout.addStretch()
        
        # About section
        about_label = QLabel("""
<b>About</b><br>
<b>xespresso GUI</b> - PySide6 interface for Quantum ESPRESSO calculations<br>
<br>
Version: 1.2.0<br>
<a href="https://github.com/vsrsousa/spresso">Documentation</a> | 
<a href="https://github.com/vsrsousa/spresso/issues">Report Issue</a>
""")
        about_label.setTextFormat(Qt.RichText)
        about_label.setOpenExternalLinks(True)
        about_label.setWordWrap(True)
        layout.addWidget(about_label)
        
        return sidebar
    
    def _create_pages(self):
        """Create all page widgets (workflow pages only) using lazy imports."""
        # Create pages in order matching navigation list
        # Configuration pages are now in the dialog
        # Use lazy imports for faster startup
        self.pages = [
            _get_page_class('StructureViewerPage')(self.session_state),
            _get_page_class('CalculationSetupPage')(self.session_state),
            _get_page_class('WorkflowBuilderPage')(self.session_state),
            _get_page_class('JobSubmissionPage')(self.session_state),
            _get_page_class('ResultsPostprocessingPage')(self.session_state)
        ]
        
        # Store reference to JobSubmissionPage for setting job monitor later
        self._job_submission_page = self.pages[3]  # 4th page in the list
        
        for page in self.pages:
            self.content_stack.addWidget(page)
    
    def _setup_menu(self):
        """Setup the menu bar."""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("&File")
        
        new_session_action = QAction("New Session", self)
        new_session_action.setShortcut("Ctrl+N")
        new_session_action.triggered.connect(self._new_session)
        file_menu.addAction(new_session_action)
        
        save_session_action = QAction("Save Session", self)
        save_session_action.setShortcut("Ctrl+S")
        save_session_action.triggered.connect(self._save_session)
        file_menu.addAction(save_session_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Edit menu
        edit_menu = menubar.addMenu("&Edit")
        
        config_action = QAction("Configuration...", self)
        config_action.setShortcut("Ctrl+,")
        config_action.triggered.connect(self._open_config_dialog)
        edit_menu.addAction(config_action)
        
        # View menu
        view_menu = menubar.addMenu("&View")
        
        # Navigate menu items
        nav_actions = [
            ("Structure", "Ctrl+1"),
            ("Calculation Setup", "Ctrl+2"),
            ("Workflow Builder", "Ctrl+3"),
            ("Job Submission", "Ctrl+4"),
            ("Results", "Ctrl+5"),
        ]
        
        for i, (name, shortcut) in enumerate(nav_actions):
            action = QAction(name, self)
            action.setShortcut(shortcut)
            action.triggered.connect(lambda checked, idx=i: self._navigate_to(idx))
            view_menu.addAction(action)
        
        # Help menu
        help_menu = menubar.addMenu("&Help")
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)
    
    def _setup_toolbar(self):
        """Setup the toolbar."""
        toolbar = QToolBar("Main Toolbar")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        
        # Configuration button
        config_action = QAction("⚙️ Configuration", self)
        config_action.setToolTip("Open configuration dialog")
        config_action.triggered.connect(self._open_config_dialog)
        toolbar.addAction(config_action)
        
        toolbar.addSeparator()
        
        # Session actions
        new_action = QAction("📋 New Session", self)
        new_action.setToolTip("Create a new session")
        new_action.triggered.connect(self._new_session)
        toolbar.addAction(new_action)
        
        save_action = QAction("💾 Save Session", self)
        save_action.setToolTip("Save current session")
        save_action.triggered.connect(self._save_session)
        toolbar.addAction(save_action)
        
        # Add a load session button
        load_action = QAction("📂 Load Session", self)
        load_action.setToolTip("Load a saved session")
        load_action.triggered.connect(self._load_session_dialog)
        toolbar.addAction(load_action)
        
        toolbar.addSeparator()
        
        # Job Monitor button - accessible without session
        job_monitor_action = QAction("🔍 Job Monitor", self)
        job_monitor_action.setToolTip("Track remote job submissions")
        job_monitor_action.triggered.connect(self._open_job_monitor)
        toolbar.addAction(job_monitor_action)
        
        # Spacer to push quit button to the right
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        toolbar.addWidget(spacer)
        
        # Quit button at the end
        quit_action = QAction("🚪 Quit", self)
        quit_action.setToolTip("Quit the application")
        quit_action.triggered.connect(self.close)
        toolbar.addAction(quit_action)
    
    def _setup_statusbar(self):
        """Setup the status bar."""
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        self._update_statusbar()
    
    def _update_statusbar(self):
        """Update status bar with current session info."""
        session_name = self.session_state.get_session_name()
        self.statusbar.showMessage(f"Session: {session_name} | Ready")
    
    def _update_session_name_label(self):
        """Update the session name label in the sidebar."""
        session_name = self.session_state.get_session_name()
        session_id = self.session_state.get_current_session_id()
        self.session_name_label.setText(f"<b>{session_name}</b><br><small>ID: {session_id}</small>")
        self.session_name_label.setTextFormat(Qt.RichText)
    
    def _refresh_session_list(self):
        """Refresh the session list in the combo box.
        
        Only shows sessions that have been created or loaded by the user in this app session.
        The combo box starts empty until the user creates a new session or loads a saved one.
        """
        self.session_combo.blockSignals(True)
        self.session_combo.clear()
        
        # Get list of active session names (created or loaded in this app session)
        active_names = self.session_state.get_active_session_names()
        current_name = self.session_state.get_session_name()
        
        # Add each active session to the combo box
        for name in active_names:
            self.session_combo.addItem(name, name)
        
        # Select current session if it's in the list
        if current_name in active_names:
            index = active_names.index(current_name)
            self.session_combo.setCurrentIndex(index)
        
        self.session_combo.blockSignals(False)
    
    def _on_session_selected(self, text):
        """Handle session selection from combo box."""
        # Guard against recursive calls during updates
        if self._updating:
            return
        
        session_id = self.session_combo.currentData()
        if session_id and session_id != self.session_state.get_current_session_id():
            if self.session_state.switch_session(session_id):
                self._on_session_changed()
                self.statusbar.showMessage(f"Switched to session: {text}")
    
    def _on_session_changed(self):
        """Handle session state changes.
        
        This method is called when the session state changes, either from
        the session listener or when switching sessions. Uses a guard to
        prevent recursive updates.
        """
        # Guard against recursive calls
        if self._updating:
            return
        
        self._updating = True
        try:
            self._update_session_name_label()
            self._update_statusbar()
            self.workdir_label.setText(self.session_state.get('working_directory', '~'))
            
            # Refresh all pages
            for page in self.pages:
                if hasattr(page, 'refresh'):
                    try:
                        page.refresh()
                    except Exception:
                        pass  # Don't let page refresh errors crash the app
        finally:
            self._updating = False
    
    def _open_config_dialog(self):
        """Open the configuration dialog."""
        # Import here to avoid circular imports
        from qtgui.dialogs import ConfigurationDialog
        
        if self._config_dialog is None:
            self._config_dialog = ConfigurationDialog(self.session_state, self)
            self._config_dialog.configuration_changed.connect(self._on_config_changed)
        
        self._config_dialog.show()
        self._config_dialog.raise_()
        self._config_dialog.activateWindow()
    
    def _on_config_changed(self):
        """Handle configuration changes from the dialog."""
        # Refresh pages that depend on configuration
        for page in self.pages:
            if hasattr(page, 'refresh'):
                page.refresh()
        self.statusbar.showMessage("Configuration updated")
    
    def _on_nav_changed(self, index):
        """Handle navigation change."""
        self.content_stack.setCurrentIndex(index)
        
        # Update status bar
        page_names = [
            "Structure",
            "Calculation Setup",
            "Workflow Builder",
            "Job Submission",
            "Results & Post-Processing"
        ]
        if 0 <= index < len(page_names):
            session_name = self.session_state.get_session_name()
            self.statusbar.showMessage(f"Session: {session_name} | Page: {page_names[index]}")
    
    def _navigate_to(self, index):
        """Navigate to a specific page."""
        self.nav_list.setCurrentRow(index)
    
    def _browse_workdir(self):
        """Open directory browser for working directory."""
        current_dir = self.session_state.get('working_directory', os.path.expanduser("~"))
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Working Directory",
            current_dir,
            QFileDialog.ShowDirsOnly
        )
        if directory:
            self.session_state['working_directory'] = directory
            self.workdir_label.setText(directory)
            self.statusbar.showMessage(f"Working directory set to: {directory}")
    
    def _new_session(self):
        """Create a new session."""
        name, ok = QInputDialog.getText(
            self,
            "New Session",
            "Enter a name for the new session:",
            QLineEdit.Normal,
            f"Session {len(self.session_state.list_sessions()) + 1}"
        )
        
        if ok and name:
            session_id = self.session_state.create_session(name)
            # Add to active sessions list (for combo box)
            self.session_state.add_active_session(name)
            self._refresh_session_list()
            self._on_session_changed()
            self.statusbar.showMessage(f"Created new session: {name}")
    
    def _rename_session(self):
        """Rename the current session."""
        current_name = self.session_state.get_session_name()
        name, ok = QInputDialog.getText(
            self,
            "Rename Session",
            "Enter a new name for the session:",
            QLineEdit.Normal,
            current_name
        )
        
        if ok and name:
            # Update the active sessions list
            self.session_state.rename_active_session(current_name, name)
            # Rename the session (in-memory and JSON file)
            self.session_state.rename_session(name)
            self._refresh_session_list()
            self._update_session_name_label()
            self.statusbar.showMessage(f"Session renamed to: {name}")
    
    def _save_session(self):
        """Save the current session.
        
        First collects current state from all pages, then saves to disk.
        """
        # Set updating flag to prevent refresh during save
        # This prevents the session state listener from triggering page refreshes
        # while we're collecting state from pages, which would restore old UI values
        self._updating = True
        try:
            # Collect current state from all pages before saving
            for page in self.pages:
                if hasattr(page, 'save_state'):
                    try:
                        page.save_state()
                    except Exception as e:
                        print(f"Warning: Could not save state from page: {e}")
            
            self.session_state.save_session()
            self.statusbar.showMessage("Session saved")
        finally:
            self._updating = False
    
    def _close_session(self):
        """Close the current session and clear all data."""
        reply = QMessageBox.question(
            self,
            "Close Session",
            "Close the current session and clear all data?\n\n"
            "Make sure to save your session first if you want to keep it.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            # Clear session state
            self.session_state.clear_state()
            
            # Clear the session combo box (no active session)
            self.session_combo.blockSignals(True)
            self.session_combo.clear()
            self.session_combo.blockSignals(False)
            
            # Process events to keep UI responsive
            from PySide6.QtWidgets import QApplication
            QApplication.processEvents()
            
            # Refresh all pages to show empty state
            for page in self.pages:
                if hasattr(page, 'refresh'):
                    try:
                        page.refresh()
                        # Process events after each page refresh for responsiveness
                        QApplication.processEvents()
                    except Exception as e:
                        print(f"Warning: Could not refresh page: {e}")
            
            # Update UI
            self._update_session_name_label()
            self.workdir_label.setText(self.session_state.get('working_directory', '~'))
            
            self.statusbar.showMessage("Session closed")
    
    def _load_session_dialog(self):
        """Open a dialog to load a saved session from the sessions directory.
        
        Lists the JSON filenames (without .json extension) from ~/.xespresso/sessions.
        When user selects one, it becomes the active session.
        """
        sessions_dir = self.session_state.get_sessions_dir()
        
        # Ensure sessions directory exists
        os.makedirs(sessions_dir, exist_ok=True)
        
        # Get list of session files from the actual filesystem
        session_files = self.session_state.list_session_files()
        
        if not session_files:
            # If no sessions found, offer to open file dialog anyway
            reply = QMessageBox.question(
                self,
                "No Saved Sessions",
                f"No saved sessions found in:\n{sessions_dir}\n\n"
                "Would you like to browse for a session file?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes
            )
            if reply == QMessageBox.No:
                return
            # Fall through to file dialog
        else:
            # Build list of session names for selection
            # Display name is the filename without .json extension
            session_items = []
            session_paths = {}
            for session_id, session_name, file_path in session_files:
                # Use filename without .json as the display name
                filename = os.path.basename(file_path)
                display_name = filename[:-5] if filename.endswith('.json') else filename
                session_items.append(display_name)
                session_paths[display_name] = file_path
            
            # Show selection dialog with option to browse
            session_items.append("--- Browse for other file... ---")
            
            from PySide6.QtWidgets import QInputDialog
            selected, ok = QInputDialog.getItem(
                self,
                "Load Session",
                f"Sessions directory: {sessions_dir}\n\nSelect a session to load:",
                session_items,
                0,  # Default index
                False  # Not editable
            )
            
            if not ok:
                return
            
            if selected != "--- Browse for other file... ---":
                # Load the selected session file
                file_path = session_paths[selected]
                if self.session_state.load_session_from_file(file_path):
                    # Add to active sessions list (for combo box)
                    self.session_state.add_active_session(selected)
                    self.workdir_label.setText(self.session_state.get('working_directory', '~'))
                    self._refresh_session_list()
                    self._on_session_changed()
                    self.statusbar.showMessage(f"Loaded session: {selected}")
                else:
                    QMessageBox.warning(
                        self,
                        "Error",
                        f"Could not load session: {selected}"
                    )
                return
        
        # Open file dialog starting in sessions directory
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Session File",
            sessions_dir,
            "Session Files (*.json);;All Files (*)"
        )
        
        if file_path:
            if self.session_state.load_session_from_file(file_path):
                # Add to active sessions list (for combo box) using filename without .json
                filename = os.path.basename(file_path)
                session_name = filename[:-5] if filename.endswith('.json') else filename
                self.session_state.add_active_session(session_name)
                self.workdir_label.setText(self.session_state.get('working_directory', '~'))
                self._refresh_session_list()
                self._on_session_changed()
                session_name = self.session_state.get_session_name()
                self.statusbar.showMessage(f"Loaded session: {session_name}")
            else:
                QMessageBox.warning(
                    self,
                    "Error",
                    f"Could not load session from file:\n{file_path}"
                )
    
    def _show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About xespresso GUI",
            """<h2>xespresso GUI</h2>
<p><b>Version:</b> 1.2.0</p>
<p>PySide6 interface for Quantum ESPRESSO calculations.</p>
<p>This application provides a user-friendly interface for:
<ul>
<li>Configuring machines (local/remote execution environments)</li>
<li>Setting up Quantum ESPRESSO codes</li>
<li>Viewing and selecting molecular structures</li>
<li>Configuring calculations and workflows</li>
<li>Submitting computational jobs</li>
<li>Tracking remote jobs with the Job Monitor</li>
</ul></p>
<p><b>New in 1.2.0:</b>
<ul>
<li>Migrated from PyQt5 to PySide6 for faster startup</li>
<li>Improved session management with save/load</li>
<li>Multiple session support</li>
</ul></p>
<p><a href="https://github.com/vsrsousa/spresso">GitHub Repository</a></p>
"""
        )
    
    def _get_job_monitor(self):
        """
        Get or create the Job Monitor dialog instance.
        
        Returns:
            JobMonitorDialog: The singleton Job Monitor instance
        """
        if self._job_monitor is None:
            from qtgui.dialogs.job_monitor_dialog import JobMonitorDialog
            # Use ~/.xespresso as the base directory for jobs file
            xespresso_dir = os.path.expanduser("~/.xespresso")
            self._job_monitor = JobMonitorDialog(config_dir=xespresso_dir, parent=self)
            
            # Set the job monitor reference in the Job Submission page
            # Use the stored reference (safer than searching by class name)
            if hasattr(self, '_job_submission_page'):
                self._job_submission_page.set_job_monitor(self._job_monitor)
        
        return self._job_monitor
    
    def _open_job_monitor(self):
        """Open or show the Job Monitor dialog."""
        job_monitor = self._get_job_monitor()
        
        # Show and raise the dialog
        job_monitor.show()
        job_monitor.raise_()
        job_monitor.activateWindow()
    
    def closeEvent(self, event):
        """Handle window close event."""
        # Save session before closing
        self.session_state.save_session()
        
        # Close config dialog if open
        if self._config_dialog is not None:
            self._config_dialog.close()
        
        event.accept()


def main():
    """Main entry point for the application."""
    app = QApplication(sys.argv)
    app.setApplicationName("xespresso GUI")
    app.setOrganizationName("xespresso")
    
    # Set application style
    app.setStyle("Fusion")
    
    # Set global stylesheet for consistent dropdown menu styling
    # This ensures all comboboxes use black text with good contrast
    app.setStyleSheet("""
        QComboBox QAbstractItemView {
            selection-background-color: #B3D9FF;
            selection-color: black;
            background-color: #F5F5F5;
            color: black;
        }
        QComboBox QAbstractItemView::item:hover {
            background-color: #D6EBFF;
            color: black;
        }
    """)
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()

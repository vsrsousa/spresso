"""
PyQt5 Main Application for xespresso GUI.

This module provides the main window and navigation for the xespresso
configuration interface using PyQt5.
"""

import sys
import os
import json
from pathlib import Path
from datetime import datetime

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QStackedWidget, QListWidget, QListWidgetItem, QLabel, QGroupBox,
    QFileDialog, QMessageBox, QSplitter, QFrame, QPushButton,
    QStatusBar, QMenuBar, QMenu, QAction, QToolBar, QComboBox,
    QInputDialog, QLineEdit
)
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtGui import QIcon, QFont

# Import page modules
from qtgui.pages import (
    StructureViewerPage,
    CalculationSetupPage,
    WorkflowBuilderPage,
    JobSubmissionPage,
    ResultsPostprocessingPage
)


# Default session data directory
DEFAULT_SESSIONS_DIR = os.path.expanduser("~/.xespresso/sessions")


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
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._state = {}
            cls._instance._sessions = {}
            cls._instance._current_session_id = "default"
            cls._instance._sessions_dir = DEFAULT_SESSIONS_DIR
            cls._instance._listeners = []
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
        """Notify all listeners of state change."""
        for listener in self._listeners:
            try:
                listener()
            except Exception:
                pass  # Don't let listener errors crash the app
    
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
        # Generate unique session ID
        session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
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
        Rename the current session.
        
        Args:
            new_name: New name for the session
        """
        self._state['session_name'] = new_name
        if self._current_session_id in self._sessions:
            self._sessions[self._current_session_id]['name'] = new_name
            self._save_sessions_index()
    
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
        
        Args:
            session_id: Optional specific session ID to save to
        """
        if session_id is None:
            session_id = self._current_session_id
        
        os.makedirs(self._sessions_dir, exist_ok=True)
        session_path = os.path.join(self._sessions_dir, f"{session_id}.json")
        
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
        
        try:
            with open(session_path, 'w', encoding='utf-8') as f:
                json.dump(serializable_state, f, indent=2)
            
            # Update sessions index
            if session_id not in self._sessions:
                self._sessions[session_id] = {}
            self._sessions[session_id]['name'] = self._state.get('session_name', 'Unnamed')
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
        
        # Define allowed keys for security validation
        allowed_keys = {
            'current_structure', 'current_machine', 'current_machine_name',
            'current_codes', 'selected_code_version', 'workflow_config',
            'working_directory', 'session_name', 'session_created',
            'session_modified', 'calc_machine', 'selected_machine',
            'selected_qe_version', 'structure_source'
        }
        
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
                    if isinstance(key, str) and key in allowed_keys:
                        self._state[key] = value
                    
            except Exception as e:
                print(f"Warning: Could not load session: {e}")
                self._initialize_defaults()
        else:
            self._initialize_defaults()
    
    def reset(self):
        """Reset the current session to defaults."""
        self._initialize_defaults()
        self._notify_listeners()


# Global session state instance
session_state = SessionState()


class MainWindow(QMainWindow):
    """Main application window for xespresso PyQt GUI."""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle("⚛️ xespresso - Quantum ESPRESSO Configuration GUI")
        self.setMinimumSize(1200, 800)
        
        # Initialize session state
        self.session_state = session_state
        
        # Configuration dialog (created on demand)
        self._config_dialog = None
        
        # Setup UI
        self._setup_ui()
        self._setup_menu()
        self._setup_toolbar()
        self._setup_statusbar()
        
        # Connect session state listener for UI updates
        self.session_state.add_listener(self._on_session_changed)
        
        # Select default page (Structure Viewer)
        self.nav_list.setCurrentRow(0)  # Now index 0 is Structure Viewer
    
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
        
        # Session selector
        session_selector_layout = QHBoxLayout()
        session_selector_layout.addWidget(QLabel("Active:"))
        self.session_combo = QComboBox()
        self.session_combo.setToolTip("Select active session")
        self.session_combo.currentTextChanged.connect(self._on_session_selected)
        session_selector_layout.addWidget(self.session_combo, 1)
        session_layout.addLayout(session_selector_layout)
        
        # Session name display
        self.session_name_label = QLabel("")
        self.session_name_label.setWordWrap(True)
        self._update_session_name_label()
        session_layout.addWidget(self.session_name_label)
        
        # Session buttons
        session_btn_layout = QHBoxLayout()
        
        new_session_btn = QPushButton("New")
        new_session_btn.setToolTip("Create a new session")
        new_session_btn.clicked.connect(self._new_session)
        session_btn_layout.addWidget(new_session_btn)
        
        rename_session_btn = QPushButton("Rename")
        rename_session_btn.setToolTip("Rename current session")
        rename_session_btn.clicked.connect(self._rename_session)
        session_btn_layout.addWidget(rename_session_btn)
        
        save_session_btn = QPushButton("Save")
        save_session_btn.setToolTip("Save current session")
        save_session_btn.clicked.connect(self._save_session)
        session_btn_layout.addWidget(save_session_btn)
        
        session_layout.addLayout(session_btn_layout)
        
        layout.addWidget(session_group)
        
        # Populate session combo
        self._refresh_session_list()
        
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
            "🔬 Structure Viewer",
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
<b>xespresso GUI</b> - PyQt interface for Quantum ESPRESSO calculations<br>
<br>
Version: 1.1.0<br>
<a href="https://github.com/vsrsousa/spresso">Documentation</a> | 
<a href="https://github.com/vsrsousa/spresso/issues">Report Issue</a>
""")
        about_label.setTextFormat(Qt.RichText)
        about_label.setOpenExternalLinks(True)
        about_label.setWordWrap(True)
        layout.addWidget(about_label)
        
        return sidebar
    
    def _create_pages(self):
        """Create all page widgets (workflow pages only)."""
        # Create pages in order matching navigation list
        # Configuration pages are now in the dialog
        self.pages = [
            StructureViewerPage(self.session_state),
            CalculationSetupPage(self.session_state),
            WorkflowBuilderPage(self.session_state),
            JobSubmissionPage(self.session_state),
            ResultsPostprocessingPage(self.session_state)
        ]
        
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
            ("Structure Viewer", "Ctrl+1"),
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
        """Refresh the session list in the combo box."""
        self.session_combo.blockSignals(True)
        self.session_combo.clear()
        
        sessions = self.session_state.list_sessions()
        current_id = self.session_state.get_current_session_id()
        
        # Add default session if not in list
        if 'default' not in sessions:
            self.session_combo.addItem("Default Session", "default")
        
        for session_id, info in sessions.items():
            name = info.get('name', session_id)
            self.session_combo.addItem(name, session_id)
        
        # Select current session
        for i in range(self.session_combo.count()):
            if self.session_combo.itemData(i) == current_id:
                self.session_combo.setCurrentIndex(i)
                break
        
        self.session_combo.blockSignals(False)
    
    def _on_session_selected(self, text):
        """Handle session selection from combo box."""
        session_id = self.session_combo.currentData()
        if session_id and session_id != self.session_state.get_current_session_id():
            if self.session_state.switch_session(session_id):
                self._on_session_changed()
                self.statusbar.showMessage(f"Switched to session: {text}")
    
    def _on_session_changed(self):
        """Handle session state changes."""
        self._update_session_name_label()
        self._update_statusbar()
        self.workdir_label.setText(self.session_state.get('working_directory', '~'))
        
        # Refresh all pages
        for page in self.pages:
            if hasattr(page, 'refresh'):
                page.refresh()
    
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
            "Structure Viewer",
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
            self.session_state.rename_session(name)
            self._refresh_session_list()
            self._update_session_name_label()
            self.statusbar.showMessage(f"Session renamed to: {name}")
    
    def _save_session(self):
        """Save the current session."""
        self.session_state.save_session()
        self.statusbar.showMessage("Session saved")
    
    def _show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About xespresso GUI",
            """<h2>xespresso GUI</h2>
<p><b>Version:</b> 1.1.0</p>
<p>PyQt5 interface for Quantum ESPRESSO calculations.</p>
<p>This application provides a user-friendly interface for:
<ul>
<li>Configuring machines (local/remote execution environments)</li>
<li>Setting up Quantum ESPRESSO codes</li>
<li>Viewing and selecting molecular structures</li>
<li>Configuring calculations and workflows</li>
<li>Submitting computational jobs</li>
</ul></p>
<p><b>New in 1.1.0:</b>
<ul>
<li>Non-blocking configuration dialog</li>
<li>Improved session management with save/load</li>
<li>Multiple session support</li>
</ul></p>
<p><a href="https://github.com/vsrsousa/spresso">GitHub Repository</a></p>
"""
        )
    
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
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

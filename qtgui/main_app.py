"""
PyQt5 Main Application for xespresso GUI.

This module provides the main window and navigation for the xespresso
configuration interface using PyQt5.
"""

import sys
import os
from pathlib import Path

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QStackedWidget, QListWidget, QListWidgetItem, QLabel, QGroupBox,
    QFileDialog, QMessageBox, QSplitter, QFrame, QPushButton,
    QStatusBar, QMenuBar, QMenu, QAction, QToolBar
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon, QFont

# Import page modules
from qtgui.pages import (
    MachineConfigPage,
    CodesConfigPage,
    PseudopotentialsConfigPage,
    StructureViewerPage,
    CalculationSetupPage,
    WorkflowBuilderPage,
    JobSubmissionPage,
    ResultsPostprocessingPage
)


class SessionState:
    """
    Singleton class to manage application state (similar to st.session_state in Streamlit).
    
    This provides a central place to store and access state across pages.
    """
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._state = {}
            cls._instance._initialize_defaults()
        return cls._instance
    
    def _initialize_defaults(self):
        """Initialize default state values."""
        self._state['current_structure'] = None
        self._state['current_machine'] = None
        self._state['current_machine_name'] = None
        self._state['current_codes'] = None
        self._state['selected_code_version'] = None
        self._state['workflow_config'] = {}
        self._state['working_directory'] = os.path.expanduser("~")
    
    def __getitem__(self, key):
        return self._state.get(key)
    
    def __setitem__(self, key, value):
        self._state[key] = value
    
    def __contains__(self, key):
        return key in self._state
    
    def get(self, key, default=None):
        return self._state.get(key, default)


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
        
        # Setup UI
        self._setup_ui()
        self._setup_menu()
        self._setup_statusbar()
        
        # Select default page (Structure Viewer)
        self.nav_list.setCurrentRow(3)  # Index of Structure Viewer
    
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
        
        # Configuration section (collapsible group)
        config_group = QGroupBox("⚙️ Configuration")
        config_layout = QVBoxLayout(config_group)
        
        self.config_list = QListWidget()
        self.config_list.setMaximumHeight(120)
        config_items = [
            "🖥️ Machine Configuration",
            "⚙️ Codes Configuration",
            "🧪 Pseudopotentials"
        ]
        for item in config_items:
            self.config_list.addItem(QListWidgetItem(item))
        self.config_list.itemClicked.connect(self._on_config_item_clicked)
        config_layout.addWidget(self.config_list)
        layout.addWidget(config_group)
        
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
        
        # Workflow section
        workflow_group = QGroupBox("🔬 Session & Workflow")
        workflow_layout = QVBoxLayout(workflow_group)
        
        self.nav_list = QListWidget()
        nav_items = [
            "🖥️ Machine Configuration",
            "⚙️ Codes Configuration",
            "🧪 Pseudopotentials",
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
Version: 1.0.0<br>
<a href="https://github.com/vsrsousa/spresso">Documentation</a> | 
<a href="https://github.com/vsrsousa/spresso/issues">Report Issue</a>
""")
        about_label.setTextFormat(Qt.RichText)
        about_label.setOpenExternalLinks(True)
        about_label.setWordWrap(True)
        layout.addWidget(about_label)
        
        return sidebar
    
    def _create_pages(self):
        """Create all page widgets."""
        # Create pages in order matching navigation list
        self.pages = [
            MachineConfigPage(self.session_state),
            CodesConfigPage(self.session_state),
            PseudopotentialsConfigPage(self.session_state),
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
        
        file_menu.addSeparator()
        
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # View menu
        view_menu = menubar.addMenu("&View")
        
        # Navigate menu items
        nav_actions = [
            ("Machine Configuration", "Ctrl+1"),
            ("Codes Configuration", "Ctrl+2"),
            ("Pseudopotentials", "Ctrl+3"),
            ("Structure Viewer", "Ctrl+4"),
            ("Calculation Setup", "Ctrl+5"),
            ("Workflow Builder", "Ctrl+6"),
            ("Job Submission", "Ctrl+7"),
            ("Results", "Ctrl+8"),
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
    
    def _setup_statusbar(self):
        """Setup the status bar."""
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        self.statusbar.showMessage("Ready")
    
    def _on_config_item_clicked(self, item):
        """Handle configuration item click."""
        index = self.config_list.currentRow()
        self.nav_list.setCurrentRow(index)
    
    def _on_nav_changed(self, index):
        """Handle navigation change."""
        self.content_stack.setCurrentIndex(index)
        
        # Update status bar
        page_names = [
            "Machine Configuration",
            "Codes Configuration",
            "Pseudopotentials Configuration",
            "Structure Viewer",
            "Calculation Setup",
            "Workflow Builder",
            "Job Submission",
            "Results & Post-Processing"
        ]
        if 0 <= index < len(page_names):
            self.statusbar.showMessage(f"Current page: {page_names[index]}")
    
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
        """Start a new session (reset state)."""
        reply = QMessageBox.question(
            self,
            "New Session",
            "Are you sure you want to start a new session? All unsaved data will be lost.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            self.session_state._initialize_defaults()
            # Refresh all pages
            for page in self.pages:
                if hasattr(page, 'refresh'):
                    page.refresh()
            self.statusbar.showMessage("New session started")
    
    def _show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About xespresso GUI",
            """<h2>xespresso GUI</h2>
<p><b>Version:</b> 1.0.0</p>
<p>PyQt5 interface for Quantum ESPRESSO calculations.</p>
<p>This application provides a user-friendly interface for:
<ul>
<li>Configuring machines (local/remote execution environments)</li>
<li>Setting up Quantum ESPRESSO codes</li>
<li>Viewing and selecting molecular structures</li>
<li>Configuring calculations and workflows</li>
<li>Submitting computational jobs</li>
</ul></p>
<p><a href="https://github.com/vsrsousa/spresso">GitHub Repository</a></p>
"""
        )


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

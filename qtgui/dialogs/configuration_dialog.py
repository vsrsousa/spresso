"""
Configuration Dialog for xespresso PySide6 GUI.

This module provides a non-blocking dialog window that contains configuration
options for Machines, Codes, and Pseudopotentials in a tabbed interface.
The dialog doesn't block the main application and allows users to configure
while still interacting with the main window.
"""

from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget, QPushButton,
    QLabel, QWidget
)
from PySide6.QtCore import Qt, Signal


class ConfigurationDialog(QDialog):
    """
    Non-blocking configuration dialog with tabs for Machine, Codes, and Pseudopotentials.
    
    This dialog allows configuration without blocking the main application window.
    Changes are applied immediately and reflected in the session state.
    
    Signals:
        configuration_changed: Emitted when any configuration is saved.
    """
    
    configuration_changed = Signal()
    
    def __init__(self, session_state, parent=None):
        """
        Initialize the configuration dialog.
        
        Args:
            session_state: Application session state object
            parent: Parent widget (main window)
        """
        super().__init__(parent)
        self.session_state = session_state
        self.setWindowTitle("⚙️ Configuration")
        self.setMinimumSize(900, 700)
        
        # Make the dialog non-modal (doesn't block main window)
        self.setModal(False)
        
        # Set window flags to keep dialog on top but allow interaction with main window
        self.setWindowFlags(
            Qt.Window | Qt.WindowCloseButtonHint | Qt.WindowMinimizeButtonHint
        )
        
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout(self)
        
        # Header
        header = QLabel("""
<h2>⚙️ Configuration</h2>
<p>Configure your computational environment: machines, QE codes, and pseudopotentials.</p>
<p><b>Note:</b> This window doesn't block the main application. You can switch between windows freely.</p>
""")
        header.setTextFormat(Qt.RichText)
        header.setWordWrap(True)
        layout.addWidget(header)
        
        # Tab widget for different configuration sections
        self.tab_widget = QTabWidget()
        
        # Import page widgets (delayed import to avoid circular imports)
        from qtgui.pages import (
            MachineConfigPage,
            CodesConfigPage,
            PseudopotentialsConfigPage
        )
        
        # Create page instances with session state
        self.machine_page = MachineConfigPage(self.session_state)
        self.codes_page = CodesConfigPage(self.session_state)
        self.pseudopotentials_page = PseudopotentialsConfigPage(self.session_state)
        
        # Add tabs
        self.tab_widget.addTab(self.machine_page, "🖥️ Machines")
        self.tab_widget.addTab(self.codes_page, "⚙️ Codes")
        self.tab_widget.addTab(self.pseudopotentials_page, "🧪 Pseudopotentials")
        
        layout.addWidget(self.tab_widget, 1)  # Stretch factor 1 to take available space

        # Provenance configuration tab
        prov_widget = QWidget()
        prov_layout = QVBoxLayout(prov_widget)
        prov_layout.addWidget(QLabel("Provenance database path:"))
        from PySide6.QtWidgets import QLineEdit, QHBoxLayout, QFileDialog

        self._prov_path_edit = QLineEdit()
        self._prov_browse = QPushButton("Browse...")
        hl = QHBoxLayout()
        hl.addWidget(self._prov_path_edit)
        hl.addWidget(self._prov_browse)
        prov_layout.addLayout(hl)

        # Populate current value from session_state if available
        try:
            current = self.session_state.get('provenance_db_path', '')
            if current:
                self._prov_path_edit.setText(current)
            else:
                # Prefill with the default provenance DB path (respects env var)
                try:
                    from xespresso.provenance import ProvenanceDB

                    default = ProvenanceDB.get_default().path
                    if default:
                        self._prov_path_edit.setText(str(default))
                except Exception:
                    pass
        except Exception:
            pass

        def _browse_prov():
            start = self.session_state.get('provenance_db_path') or self.session_state.get('working_directory') or ''
            path, _ = QFileDialog.getSaveFileName(self, "Provenance DB File", start, "JSON files (*.json);;All Files (*)")
            if path:
                self._prov_path_edit.setText(path)

        self._prov_browse.clicked.connect(_browse_prov)

        # Apply button for provenance path
        apply_btn = QPushButton("Apply Provenance Path")
        prov_layout.addWidget(apply_btn)

        def _apply_prov():
            val = self._prov_path_edit.text().strip()
            if val:
                try:
                    self.session_state['provenance_db_path'] = val
                except Exception:
                    pass
            # also emit configuration_changed so other UI can react
            self.configuration_changed.emit()

        apply_btn.clicked.connect(_apply_prov)

        self.tab_widget.addTab(prov_widget, "📚 Provenance")
        
        # Button row
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        # Refresh button
        refresh_btn = QPushButton("🔄 Refresh All")
        refresh_btn.clicked.connect(self._refresh_all)
        button_layout.addWidget(refresh_btn)
        
        # Close button
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        button_layout.addWidget(close_btn)
        
        layout.addLayout(button_layout)
    
    def _refresh_all(self):
        """Refresh all configuration pages."""
        self.machine_page.refresh()
        self.codes_page.refresh()
        self.pseudopotentials_page.refresh()
        self.configuration_changed.emit()
    
    def show_machine_tab(self):
        """Show the machine configuration tab."""
        self.tab_widget.setCurrentWidget(self.machine_page)
        self.show()
        self.raise_()
        self.activateWindow()
    
    def show_codes_tab(self):
        """Show the codes configuration tab."""
        self.tab_widget.setCurrentWidget(self.codes_page)
        self.show()
        self.raise_()
        self.activateWindow()
    
    def show_pseudopotentials_tab(self):
        """Show the pseudopotentials configuration tab."""
        self.tab_widget.setCurrentWidget(self.pseudopotentials_page)
        self.show()
        self.raise_()
        self.activateWindow()
    
    def refresh(self):
        """Refresh all pages in the dialog."""
        self._refresh_all()

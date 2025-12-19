"""
ASE database utilities for the Qt GUI.

This module provides reusable components for ASE database operations.
"""

import os

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QGroupBox, QFormLayout, QMessageBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QFileDialog
)

try:
    from ase.db import connect as ase_db_connect
    ASE_DB_AVAILABLE = True
except ImportError:
    ASE_DB_AVAILABLE = False

# Default database path
DEFAULT_DB_PATH = os.path.expanduser("~/.xespresso/structures.db")


class DatabaseSaveWidget(QWidget):
    """Widget for saving structures to an ASE database."""
    
    def __init__(self, session_state, parent=None):
        super().__init__(parent)
        self.session_state = session_state
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the user interface."""
        layout = QFormLayout(self)
        
        info = QLabel("Save the current structure to the database for quick access.")
        info.setWordWrap(True)
        info.setStyleSheet("color: gray;")
        layout.addRow(info)
        
        # Database path
        db_path_layout = QHBoxLayout()
        self.db_path_edit = QLineEdit()
        self.db_path_edit.setText(DEFAULT_DB_PATH)
        self.db_path_edit.setToolTip("Path to ASE database file")
        db_path_layout.addWidget(self.db_path_edit)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_db_path)
        db_path_layout.addWidget(browse_btn)
        layout.addRow("Database Path:", db_path_layout)
        
        # Structure name
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("e.g., my_structure, bulk_Fe")
        layout.addRow("Structure Name:", self.name_edit)
        
        # Tags
        self.tags_edit = QLineEdit()
        self.tags_edit.setPlaceholderText("e.g., bulk, metal, test (comma-separated)")
        layout.addRow("Tags:", self.tags_edit)
        
        # Save button
        save_btn = QPushButton("💾 Save Current Structure")
        save_btn.clicked.connect(self._save_to_database)
        layout.addRow(save_btn)
        
        # Status label
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addRow(self.status_label)
    
    def _browse_db_path(self):
        """Browse for database file path."""
        current_path = self.db_path_edit.text() or DEFAULT_DB_PATH
        current_dir = os.path.dirname(current_path)
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Database File",
            current_dir,
            "ASE Database Files (*.db);;All Files (*)"
        )
        
        if file_path:
            self.db_path_edit.setText(file_path)
    
    def _save_to_database(self):
        """Save the current structure to the database."""
        if not ASE_DB_AVAILABLE:
            QMessageBox.warning(self, "Warning", "ASE database module not available")
            return
        
        atoms = self.session_state.get('current_structure')
        if atoms is None:
            QMessageBox.warning(self, "Warning", "No structure loaded. Load or build a structure first.")
            return
        
        db_path = self.db_path_edit.text().strip()
        if not db_path:
            QMessageBox.warning(self, "Warning", "Please enter a database path")
            return
        
        # Normalize the path
        db_path = os.path.abspath(os.path.expanduser(db_path))
        
        try:
            # Create database directory if it doesn't exist
            db_dir = os.path.dirname(db_path)
            if db_dir:
                os.makedirs(db_dir, exist_ok=True)
            
            db = ase_db_connect(db_path)
            
            # Prepare key-value pairs
            key_value_pairs = {}
            
            save_name = self.name_edit.text().strip()
            if save_name:
                key_value_pairs['name'] = save_name
            
            source = self.session_state.get('structure_source', '')
            if source:
                key_value_pairs['source'] = source
            
            tags = self.tags_edit.text().strip()
            if tags:
                for tag in tags.split(','):
                    tag = tag.strip()
                    if tag:
                        key_value_pairs[tag] = True
            
            db.write(atoms, **key_value_pairs)
            
            self.status_label.setText(f"✅ Saved: {atoms.get_chemical_formula()}")
            self.status_label.setStyleSheet("color: green;")
            
            QMessageBox.information(
                self, "Success",
                f"Structure saved to database!\n\n"
                f"Formula: {atoms.get_chemical_formula()}\n"
                f"Database: {db_path}"
            )
            
            self.name_edit.clear()
            self.tags_edit.clear()
            
        except Exception as e:
            self.status_label.setText(f"❌ Error: {e}")
            self.status_label.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error saving to database:\n{e}")


class DatabaseLoadWidget(QWidget):
    """Widget for loading structures from an ASE database."""
    
    # Signal to emit when a structure is loaded
    structure_loaded = None  # Will be set up as a callback
    
    def __init__(self, session_state, on_structure_loaded=None, parent=None):
        super().__init__(parent)
        self.session_state = session_state
        self.on_structure_loaded = on_structure_loaded
        self._setup_ui()
    
    def _setup_ui(self):
        """Setup the user interface."""
        layout = QVBoxLayout(self)
        
        # Database path
        db_path_layout = QHBoxLayout()
        db_path_layout.addWidget(QLabel("Database Path:"))
        self.db_path_edit = QLineEdit()
        self.db_path_edit.setText(DEFAULT_DB_PATH)
        db_path_layout.addWidget(self.db_path_edit, 1)
        
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self._browse_db_path)
        db_path_layout.addWidget(browse_btn)
        layout.addLayout(db_path_layout)
        
        # Refresh button
        refresh_btn = QPushButton("🔄 Refresh Database List")
        refresh_btn.clicked.connect(self.refresh_list)
        layout.addWidget(refresh_btn)
        
        # Status
        self.status_label = QLabel("")
        layout.addWidget(self.status_label)
        
        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["ID", "Formula", "Atoms", "Tags"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setMaximumHeight(250)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        layout.addWidget(self.table)
        
        # Load button
        load_btn = QPushButton("📥 Load Selected Structure")
        load_btn.clicked.connect(self._load_selected)
        layout.addWidget(load_btn)
        
        # Result
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        layout.addWidget(self.result_label)
    
    def _browse_db_path(self):
        """Browse for database file path."""
        current_path = self.db_path_edit.text() or DEFAULT_DB_PATH
        current_dir = os.path.dirname(current_path)
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Database File",
            current_dir,
            "ASE Database Files (*.db);;All Files (*)"
        )
        
        if file_path:
            self.db_path_edit.setText(file_path)
    
    def refresh_list(self):
        """Refresh the database structures list."""
        if not ASE_DB_AVAILABLE:
            return
        
        db_path = self.db_path_edit.text().strip()
        if not db_path:
            self.status_label.setText("⚠️ Please enter a database path")
            self.status_label.setStyleSheet("color: #d97706;")
            return
        
        db_path = os.path.abspath(os.path.expanduser(db_path))
        
        if not os.path.exists(db_path):
            self.status_label.setText("ℹ️ Database does not exist yet.")
            self.status_label.setStyleSheet("color: blue;")
            self.table.setRowCount(0)
            return
        
        try:
            db = ase_db_connect(db_path)
            rows = list(db.select())
            
            if rows:
                self.status_label.setText(f"✅ Found {len(rows)} structure(s)")
                self.status_label.setStyleSheet("color: green;")
                
                self.table.setRowCount(len(rows))
                for i, row in enumerate(rows):
                    self.table.setItem(i, 0, QTableWidgetItem(str(row.id)))
                    self.table.setItem(i, 1, QTableWidgetItem(row.formula))
                    self.table.setItem(i, 2, QTableWidgetItem(str(row.natoms)))
                    tags = ", ".join(row.key_value_pairs.keys()) if row.key_value_pairs else ""
                    self.table.setItem(i, 3, QTableWidgetItem(tags))
            else:
                self.status_label.setText("ℹ️ Database is empty.")
                self.status_label.setStyleSheet("color: blue;")
                self.table.setRowCount(0)
                
        except Exception as e:
            self.status_label.setText(f"❌ Error: {e}")
            self.status_label.setStyleSheet("color: red;")
            self.table.setRowCount(0)
    
    def _load_selected(self):
        """Load the selected structure from the database."""
        if not ASE_DB_AVAILABLE:
            return
        
        if not self.table.selectedItems():
            QMessageBox.warning(self, "Warning", "Please select a structure to load")
            return
        
        selected_row = self.table.currentRow()
        id_item = self.table.item(selected_row, 0)
        if not id_item:
            return
        
        try:
            structure_id = int(id_item.text())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid structure ID")
            return
        
        db_path = os.path.abspath(os.path.expanduser(self.db_path_edit.text().strip()))
        
        try:
            db = ase_db_connect(db_path)
            row = db.get(id=structure_id)
            atoms = row.toatoms()
            
            self.result_label.setText(f"✅ Loaded: {atoms.get_chemical_formula()}")
            self.result_label.setStyleSheet("color: green;")
            
            # Save the database path for session restoration
            self.session_state['structure_db_path'] = db_path
            
            # Call the callback if provided
            if self.on_structure_loaded:
                self.on_structure_loaded(atoms, f"Database: ID {structure_id}")
            
        except Exception as e:
            self.result_label.setText(f"❌ Error: {e}")
            self.result_label.setStyleSheet("color: red;")
            QMessageBox.critical(self, "Error", f"Error loading structure:\n{e}")

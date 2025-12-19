from qtpy.QtWidgets import QWidget, QVBoxLayout, QTabWidget

from .workflow_instance import WorkflowInstancePage


class WorkflowTabsPage(QWidget):
    """Container page that holds independent workflow instance tabs."""

    def __init__(self, session_state, parent=None):
        super().__init__(parent)
        self.session = session_state
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        # allow users to close workflow tabs
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self._on_tab_close_requested)
        # ensure we refresh structure-dependent inputs when switching tabs
        try:
            self.tabs.currentChanged.connect(self._on_tab_changed)
        except Exception:
            pass
        layout.addWidget(self.tabs)

    def add_workflow_tab(self, preset_name: str):
        """Add a workflow tab or focus existing one for the same preset.

        Ensures there is at most one tab per workflow preset. If a tab for
        `preset_name` already exists, focus it and return the page.
        """
        # Check for existing tab with same title
        for idx in range(self.tabs.count()):
            if self.tabs.tabText(idx) == preset_name:
                self.tabs.setCurrentIndex(idx)
                return self.tabs.widget(idx)

        page = WorkflowInstancePage(self.session, preset_name, parent=self)
        title = preset_name
        idx = self.tabs.addTab(page, title)
        self.tabs.setCurrentIndex(idx)
        # ensure the newly created page reflects the current session structure
        try:
            atoms = self.session.get('current_structure')
            try:
                page.config_widget.update_for_structure(atoms)
            except Exception:
                pass
        except Exception:
            pass
        # persist current open tabs into the session state
        try:
            self._save_open_tabs()
        except Exception:
            pass
        return page

    def _on_tab_close_requested(self, index: int):
        widget = self.tabs.widget(index)
        # remove and delete the widget
        self.tabs.removeTab(index)
        try:
            widget.deleteLater()
        except Exception:
            pass
        try:
            self._save_open_tabs()
        except Exception:
            pass

    def _save_open_tabs(self):
        try:
            names = [self.tabs.tabText(i) for i in range(self.tabs.count())]
            try:
                self.session['open_workflow_tabs'] = names
            except Exception:
                try:
                    setattr(self.session, 'open_workflow_tabs', names)
                except Exception:
                    pass
        except Exception:
            pass

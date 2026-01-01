from qtpy.QtWidgets import QDialog, QVBoxLayout, QDialogButtonBox, QTextEdit, QLabel

try:
    from qtgui.utils.pseudopotentials_selector import PseudopotentialsSelectorWidget
    PSEUDO_SELECTOR_AVAILABLE = True
except Exception:
    PSEUDO_SELECTOR_AVAILABLE = False


class PseudopotentialsDialog(QDialog):
    """Dialog to choose pseudopotentials. Uses selector widget when available,
    otherwise provides a simple text area for mapping entries like `Fe=Fe.pseudo`.
    """

    def __init__(self, session=None, parent=None):
        super().__init__(parent)
        self.session = session
        self.setWindowTitle('Select Pseudopotentials')
        self._init_ui()
        # If session has a current structure, pre-populate selector/manual inputs
        try:
            atoms = None
            if getattr(self, 'session', None) is not None:
                try:
                    atoms = self.session.get('current_structure')
                except Exception:
                    atoms = None
            elements = set()
            if atoms is not None:
                try:
                    symbols = []
                    try:
                        symbols = atoms.get_chemical_symbols()
                    except Exception:
                        symbols = getattr(atoms, 'symbols', []) or getattr(atoms, 'get_elements', lambda: [])()
                    elements = set(symbols)
                except Exception:
                    elements = set()
            if elements:
                try:
                    self.set_elements(elements)
                except Exception:
                    pass
        except Exception:
            pass

    def _init_ui(self):
        layout = QVBoxLayout(self)
        if PSEUDO_SELECTOR_AVAILABLE:
            try:
                self.selector = PseudopotentialsSelectorWidget(self.session)
            except Exception:
                self.selector = None
        else:
            self.selector = None
        if self.selector is not None:
            layout.addWidget(self.selector)
        else:
            layout.addWidget(QLabel('Enter pseudopotential mappings (comma-separated, e.g. Fe=Fe.pseudo):'))
            self.text = QTextEdit()
            layout.addWidget(self.text)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def get_config(self):
        if getattr(self, 'selector', None) is not None and hasattr(self.selector, 'get_pseudopotentials'):
            try:
                return {'pseudopotentials': self.selector.get_pseudopotentials() or {}}
            except Exception:
                return {'pseudopotentials': {}}
        else:
            raw = getattr(self, 'text', None)
            if raw is None:
                return {'pseudopotentials': {}}
            txt = raw.toPlainText()
            mapping = {}
            for part in txt.split(','):
                part = part.strip()
                if not part:
                    continue
                if '=' in part:
                    k, v = part.split('=', 1)
                    mapping[k.strip()] = v.strip()
                elif ':' in part:
                    k, v = part.split(':', 1)
                    mapping[k.strip()] = v.strip()
            return {'pseudopotentials': mapping}

    def set_elements(self, elements):
        """Provide element set so selector (or manual input) can populate entries."""
        if getattr(self, 'selector', None) is not None:
            try:
                self.selector.set_elements(elements)
            except Exception:
                pass
        else:
            # populate placeholder text to suggest entries
            try:
                if getattr(self, 'text', None) is not None:
                    hints = ', '.join(f"{el}=" for el in sorted(elements))
                    self.text.setPlainText(hints)
            except Exception:
                pass

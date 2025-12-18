from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QComboBox, QTextEdit, QCheckBox, QHBoxLayout
)
from qtpy.QtCore import Qt, Signal
import threading
import tempfile

from xespresso.workflow.tasks import ScfTask, RelaxTask, WorkflowRunner
from qtgui.widgets.calculation_config_widget import CalculationConfigWidget
from qtgui.calculations.preparation import prepare_calculation_from_gui


class WorkflowInstancePage(QWidget):
    """A session tab for configuring and running a workflow instance."""

    def __init__(self, session_state, preset_name: str, parent=None):
        super().__init__(parent)
        self.session = session_state
        self.preset_name = preset_name
        self._init_ui()
        # Load any saved draft config for this preset
        try:
            self._load_saved_config()
        except Exception:
            pass

    def _init_ui(self):
        layout = QVBoxLayout(self)
        header = QLabel(f"Workflow: {self.preset_name}")
        header.setObjectName("pageTitle")
        layout.addWidget(header)
        # Detailed calculation configuration widget (based on CalculationSetupPage)
        self.config_widget = CalculationConfigWidget(self.session, preset_name=self.preset_name)
        layout.addWidget(self.config_widget)
        try:
            self.config_widget.changed.connect(self._save_draft)
        except Exception:
            pass

        # Provide a simple compatibility proxy named `machine_edit` so tests
        # and older UI code can set the machine via `setText()` on the page.
        try:
            class _MachineEditProxy:
                def __init__(self, combo):
                    self._combo = combo

                def setText(self, txt):
                    try:
                        # QComboBox supports setCurrentText / setEditText
                        if hasattr(self._combo, 'setCurrentText'):
                            self._combo.setCurrentText(txt)
                        elif hasattr(self._combo, 'setEditText'):
                            self._combo.setEditText(txt)
                        else:
                            # fallback: try setting the editable line edit
                            le = getattr(self._combo, 'lineEdit', lambda: None)()
                            if le is not None and hasattr(le, 'setText'):
                                le.setText(txt)
                    except Exception:
                        pass

                def text(self):
                    try:
                        return self._combo.currentText()
                    except Exception:
                        return ''

            self.machine_edit = _MachineEditProxy(self.config_widget.machine_combo)
        except Exception:
            # ensure attribute exists even if proxy creation fails
            self.machine_edit = None
        try:
            class _PseudoEditProxy:
                def __init__(self, selector):
                    self._sel = selector

                def setText(self, txt: str):
                    try:
                        # Prefer using the widget's set_pseudopotentials API when
                        # available so mappings like 'Fe=Fe.pseudo' are parsed
                        # and persisted correctly.
                        if hasattr(self._sel, 'set_pseudopotentials'):
                            mapping = {}
                            for part in str(txt).split(','):
                                part = part.strip()
                                if not part:
                                    continue
                                if '=' in part:
                                    k, v = part.split('=', 1)
                                    mapping[k.strip()] = v.strip()
                                elif ':' in part:
                                    k, v = part.split(':', 1)
                                    mapping[k.strip()] = v.strip()
                            try:
                                self._sel.set_pseudopotentials(mapping)
                                return
                            except Exception:
                                pass
                        # Fallback: if the selector exposes a text widget, set raw content
                        text_widget = getattr(self._sel, 'text', None)
                        if text_widget is not None and hasattr(text_widget, 'setPlainText'):
                            text_widget.setPlainText(txt)
                            return
                    except Exception:
                        pass

                def text(self):
                    try:
                        # Prefer canonical mapping form if available
                        if hasattr(self._sel, 'get_pseudopotentials'):
                            mapping = self._sel.get_pseudopotentials() or {}
                            if mapping:
                                return ','.join(f"{k}={v}" for k, v in mapping.items())
                        # Fallback to raw text, but normalize ': ' to '=' for tests
                        text_widget = getattr(self._sel, 'text', None)
                        if text_widget is not None and hasattr(text_widget, 'toPlainText'):
                            raw = text_widget.toPlainText()
                            try:
                                # convert lines like 'Fe: Fe.pseudo' to 'Fe=Fe.pseudo'
                                parts = []
                                for line in raw.splitlines():
                                    if ':' in line:
                                        k, v = line.split(':', 1)
                                        parts.append(f"{k.strip()}={v.strip()}")
                                    elif '=' in line:
                                        parts.append(line.strip())
                                if parts:
                                    return ','.join(parts)
                            except Exception:
                                pass
                            return raw
                    except Exception:
                        pass
                    return ''

            self.pseudo_edit = _PseudoEditProxy(self.config_widget.pseudo_selector)
        except Exception:
            self.pseudo_edit = None
        # Expose protocol combo for tests/back-compat
        try:
            self.proto_combo = self.config_widget.protocol_combo
        except Exception:
            self.proto_combo = None

        # Simple magnetism checkbox proxy (stored on the page and persisted)
        try:
            class _MagnetismProxy:
                def __init__(self):
                    self._state = False

                def setChecked(self, val: bool):
                    try:
                        self._state = bool(val)
                    except Exception:
                        pass

                def isChecked(self):
                    return bool(self._state)

            self.magnetism_chk = _MagnetismProxy()
        except Exception:
            self.magnetism_chk = None

        # Attempt to keep pseudo selector in sync with current structure
        try:
            atoms = self.session.get('current_structure')
            self.config_widget.update_for_structure(atoms)
        except Exception:
            pass
        try:
            # Listen for session changes to refresh structure-dependent inputs
            self.session.add_listener(lambda: self.config_widget.update_for_structure(self.session.get('current_structure')))
        except Exception:
            pass

        # Prepare calculation button (creates Espresso calculator objects)
        self.prepare_btn = QPushButton("Prepare Calculation")
        self.prepare_btn.clicked.connect(self._on_prepare)
        layout.addWidget(self.prepare_btn, alignment=Qt.AlignRight)

        # Config preview / advanced
        self.config_preview = QTextEdit()
        self.config_preview.setReadOnly(True)
        layout.addWidget(self.config_preview)

        # Run/Submit
        self.run_btn = QPushButton("Submit Run")
        self.run_btn.clicked.connect(self._on_submit)
        layout.addWidget(self.run_btn, alignment=Qt.AlignRight)

        # Status/log area
        self.status = QLabel("Idle")
        layout.addWidget(self.status)
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        layout.addWidget(self.log)

    def _append_log(self, text: str):
        self.log.append(text)

    def _compose_config(self):
        # Base config derived from workflow name (calculation type)
        cfg = {}
        name = (self.preset_name or "").lower()
        if 'scf' in name or self.preset_name == 'SCF':
            cfg['calc_type'] = 'scf'
        elif 'relax' in name or 'geometry' in name:
            cfg['calc_type'] = 'relax'
        else:
            cfg['calc_type'] = 'scf'

        # Merge header values (machine, pseudopotentials, label)
        try:
            hdr = self.config_widget.get_config()
            cfg.update(hdr)
        except Exception:
            pass

        # Include explicit protocol/magnetism flags if present on the page
        try:
            if getattr(self, 'proto_combo', None) is not None:
                cfg['protocol'] = self.proto_combo.currentText()
        except Exception:
            pass
        try:
            if getattr(self, 'magnetism_chk', None) is not None:
                cfg['magnetism'] = bool(self.magnetism_chk.isChecked())
        except Exception:
            pass

        return cfg

    def _save_draft(self, *args, **kwargs):
        """Save current form values to the session state as a draft for this preset."""
        cfg = self._compose_config()
        # Ensure label is present (header widget manages its own label)
        try:
            label = self.config_widget.label_edit.text().strip() or None
        except Exception:
            label = None
        if label:
            cfg['label'] = label
        try:
            store = self.session.get('workflow_tabs_config') or {}
        except Exception:
            store = getattr(self.session, 'workflow_tabs_config', {}) or {}
        store[self.preset_name] = cfg
        try:
            self.session['workflow_tabs_config'] = store
        except Exception:
            try:
                setattr(self.session, 'workflow_tabs_config', store)
            except Exception:
                pass

    def _load_saved_config(self):
        """Load saved draft config for this preset from session state, if any."""
        try:
            store = self.session.get('workflow_tabs_config') or {}
        except Exception:
            store = getattr(self.session, 'workflow_tabs_config', {}) or {}
        cfg = store.get(self.preset_name, {})
        if not cfg:
            return
        # Restore full calculation UI from dict using the configuration widget
        try:
            self.config_widget.restore_from_dict(cfg)
        except Exception:
            pass
        # Restore additional fields not handled by the widget
        try:
            if cfg.get('protocol') and getattr(self, 'proto_combo', None) is not None:
                try:
                    self.proto_combo.setCurrentText(cfg.get('protocol'))
                except Exception:
                    pass
        except Exception:
            pass
        try:
            if 'magnetism' in cfg and getattr(self, 'magnetism_chk', None) is not None:
                try:
                    self.magnetism_chk.setChecked(bool(cfg.get('magnetism')))
                except Exception:
                    pass
        except Exception:
            pass

    def _on_submit(self):
        atoms = None
        try:
            atoms = self.session.get('current_structure')
        except Exception:
            atoms = None

        try:
            label = self.config_widget.label_edit.text().strip() or None
        except Exception:
            label = None
        if label is None and atoms is not None:
            # follow Calculation Setup style
            formula = "structure"
            try:
                formula = atoms.get_chemical_formula()
            except Exception:
                pass
            label = f"{self.preset_name.lower()}/{formula}"

        config = self._compose_config()
        self.config_preview.setPlainText(str(config))

        # Build tasks similarly to WorkflowRunDialog mapping
        tasks = []
        # Map header/config fields to runner common params
        common = {
            "label": label,
            "pseudopotentials": config.get('pseudopotentials', {}),
            "queue": config.get('machine_name') or config.get('machine') or None,
            "debug": True
        }
        if self.preset_name == "SCF":
            params = dict(common); params.update({"atoms": atoms}); tasks.append(ScfTask(name="scf", params=params))
        elif self.preset_name in ("Relax", "Geometry Optimization"):
            params = dict(common); params.update({"atoms": atoms}); tasks.append(RelaxTask(name="relax", params=params))
        else:
            params = dict(common); params.update({"atoms": atoms}); tasks.append(ScfTask(name="scf", params=params))

        runner = WorkflowRunner(tasks=tasks, context={"debug": True, "provenance_dir": tempfile.mkdtemp(prefix="xespresso_prov_")})

        def _worker():
            try:
                self.status.setText("Running")
            except Exception:
                pass
            results = runner.run_all()
            try:
                self.status.setText("Finished")
            except Exception:
                pass
            for r in results:
                try:
                    self._append_log(str(r))
                except Exception:
                    pass
            # store run summary in session
            try:
                runs = self.session.get('workflow_runs') or []
            except Exception:
                runs = getattr(self.session, 'workflow_runs', []) or []
            runs.append({"label": label, "preset": self.preset_name, "results": results})
            try:
                self.session['workflow_runs'] = runs
            except Exception:
                try:
                    setattr(self.session, 'workflow_runs', runs)
                except Exception:
                    pass

        th = threading.Thread(target=_worker, daemon=True)
        th.start()

    def _on_prepare(self):
        """Prepare calculation objects (Espresso calculator and prepared atoms).

        Runs prepare_calculation_from_gui in a background thread and stores
        resulting objects in the session state under `espresso_calculator`
        and `prepared_atoms`. Appends a small run summary to `workflow_runs`.
        """
        atoms = None
        try:
            atoms = self.session.get('current_structure')
        except Exception:
            atoms = None

        if atoms is None:
            self._append_log("No structure loaded — cannot prepare calculation")
            return

        cfg = self._compose_config()
        try:
            label = self.config_widget.label_edit.text().strip() or None
        except Exception:
            label = None
        if label is None:
            try:
                formula = atoms.get_chemical_formula()
            except Exception:
                formula = 'structure'
            label = f"{self.preset_name.lower()}/{formula}"

        def _prepare_worker():
            try:
                self.status.setText("Preparing")
            except Exception:
                pass
            try:
                prepared_atoms, calc = prepare_calculation_from_gui(atoms, cfg, label=label)
                # store in session
                try:
                    runs = self.session.get('workflow_runs') or []
                except Exception:
                    runs = getattr(self.session, 'workflow_runs', []) or []
                runs.append({"label": label, "preset": self.preset_name, "prepared": True})
                try:
                    self.session['workflow_runs'] = runs
                except Exception:
                    try:
                        setattr(self.session, 'workflow_runs', runs)
                    except Exception:
                        pass

                try:
                    self.session['espresso_calculator'] = calc
                    self.session['prepared_atoms'] = prepared_atoms
                except Exception:
                    try:
                        setattr(self.session, 'espresso_calculator', calc)
                        setattr(self.session, 'prepared_atoms', prepared_atoms)
                    except Exception:
                        pass

                self._append_log(f"Prepared calculation '{label}' successfully")
                try:
                    self.status.setText("Prepared")
                except Exception:
                    pass
            except Exception as e:
                self._append_log(f"Failed to prepare calculation: {e}")
                try:
                    self.status.setText("Prepare Failed")
                except Exception:
                    pass

        th = threading.Thread(target=_prepare_worker, daemon=True)
        th.start()
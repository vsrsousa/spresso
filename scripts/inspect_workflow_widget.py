import os, traceback
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

try:
    from qtpy.QtWidgets import QApplication
    app = QApplication([])
except Exception as e:
    print('Failed to create QApplication:', e)
    raise

try:
    from qtgui.session_state import SessionState
    HAVE_ASE = True
    try:
        from ase import Atoms
        atoms = Atoms('Fe2O3', positions=[(0,0,0),(0,0,1),(1,0,0)])
    except Exception:
        HAVE_ASE = False
        class DummyAtoms:
            def get_chemical_symbols(self):
                return ['Fe','O']
            def get_chemical_formula(self):
                return 'FeO'
        atoms = DummyAtoms()

    session = SessionState(isolated=True)
    # set structure in session
    try:
        session['current_structure'] = atoms
    except Exception:
        try:
            setattr(session, 'current_structure', atoms)
        except Exception:
            pass

    # Instantiate WorkflowInstancePage directly to inspect its widget
    import importlib.util
    wf_path = os.path.join(os.path.dirname(__file__), '..', 'qtgui', 'pages', 'workflow_instance.py')
    wf_path = os.path.normpath(wf_path)
    # Provide lightweight dummy xespresso.workflow.tasks to satisfy imports
    import sys, types
    tasks_mod = types.ModuleType('xespresso.workflow.tasks')
    class ScfTask:
        def __init__(self, *a, **k):
            pass
    class RelaxTask:
        def __init__(self, *a, **k):
            pass
    class WorkflowRunner:
        def __init__(self, *a, **k):
            pass
        def run_all(self):
            return []
    tasks_mod.ScfTask = ScfTask
    tasks_mod.RelaxTask = RelaxTask
    tasks_mod.WorkflowRunner = WorkflowRunner
    sys.modules['xespresso.workflow.tasks'] = tasks_mod

    spec = importlib.util.spec_from_file_location('workflow_instance', wf_path)
    wf_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wf_mod)
    WorkflowInstancePage = getattr(wf_mod, 'WorkflowInstancePage')
    page = WorkflowInstancePage(session, preset_name='SCF')
    try:
        page.config_widget.update_for_structure(atoms)
    except Exception:
        pass

    cfg = getattr(page, 'config_widget', None)
    if cfg is None:
        print('No config_widget on WorkflowInstancePage')
    else:
        print('config_widget type:', type(cfg))
        sel = getattr(cfg, 'pseudo_selector', None)
        print('has pseudo_selector:', sel is not None)
        try:
            if sel is not None and hasattr(sel, 'get_pseudopotentials'):
                print('selector mapping:', sel.get_pseudopotentials())
        except Exception as e:
            print('selector get failed:', e)
        edits = getattr(cfg, 'pseudo_edits', None)
        print('has pseudo_edits dict:', isinstance(edits, dict))
        if isinstance(edits, dict):
            print('pseudo_edits keys:', list(edits.keys()))
        try:
            info = cfg.pseudo_info_label.text()
        except Exception:
            info = '<no label>'
        print('pseudo_info_label:', info)

    print('session current_structure present:', session.get('current_structure') is not None)

except Exception:
    traceback.print_exc()
finally:
    try:
        app.quit()
    except Exception:
        pass

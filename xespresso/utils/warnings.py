# xespresso/utils/warnings.py
import warnings

def warn(msg, category=UserWarning):
    warnings.warn(msg, category)

def apply_custom_format():
    warnings.formatwarning = lambda msg, cat, fname, lineno, *_: (
        f"\n⚠️ {cat.__name__} in {fname}:{lineno}\n→ {msg}\n"
    )

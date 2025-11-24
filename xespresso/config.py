# xespresso/config.py

import os

# Global verbosity flag for error handling across the package
VERBOSE_ERRORS = os.getenv("XESPRESSO_VERBOSE_ERRORS", "").lower() in ("1", "true", "yes")

# Optional: add other global flags here later
# FORCE_SCHEDULER = os.getenv("XESPRESSO_FORCE_SCHEDULER") == "1"

"""Statusbar helper."""
from qtpy.QtWidgets import QStatusBar


def setup_statusbar(main_window):
    statusbar = QStatusBar()
    main_window.setStatusBar(statusbar)
    main_window.statusbar = statusbar
    main_window._update_statusbar()

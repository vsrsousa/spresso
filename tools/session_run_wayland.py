#!/usr/bin/env python3
"""Run SessionManagerWindow in the current (Wayland) session and keep it open.

This script will attempt multiple activation methods (`requestActivate`, `raise_`,
`activateWindow`, temporary `WindowStaysOnTopHint`) and log actions to
`/tmp/session_wayland.log`. It will not auto-quit; close the window to exit.
"""

from __future__ import annotations
import sys
import os
import time
from qtpy.QtWidgets import QApplication
from qtpy.QtCore import QTimer
from qtpy.QtGui import QGuiApplication

LOG = '/tmp/session_wayland.log'

def log(msg: str):
    with open(LOG, 'a') as f:
        f.write(msg + '\n')
    print(msg)

def main():
    proj_root = os.getcwd()
    if proj_root not in sys.path:
        sys.path.insert(0, proj_root)

    with open(LOG, 'w') as f:
        f.write('Starting session_run_wayland\n')

    app = QApplication(sys.argv)
    log('Qt platform: ' + QGuiApplication.platformName())

    try:
        from qtgui.main_app import SessionManagerWindow
    except Exception as e:
        log('Import failed: ' + repr(e))
        raise

    w = SessionManagerWindow()
    w.show()
    log('Called show()')

    # Try several activation attempts over the first second
    def try_activate():
        try:
            h = w.windowHandle()
            if h is not None and hasattr(h, 'requestActivate'):
                log('Calling windowHandle().requestActivate()')
                try:
                    h.requestActivate()
                except Exception as e:
                    log('requestActivate failed: ' + repr(e))
            log('Calling raise_() and activateWindow()')
            try:
                w.raise_()
                w.activateWindow()
            except Exception as e:
                log('raise/activate failed: ' + repr(e))
            # temporary stay-on-top
            try:
                w.setWindowFlag(w.windowFlags() | w.windowFlags().WindowStaysOnTopHint)
                w.show()
                log('Set WindowStaysOnTopHint temporarily')
                QTimer.singleShot(300, lambda: (w.setWindowFlag(w.windowFlags() & ~w.windowFlags().WindowStaysOnTopHint), w.show()))
            except Exception as e:
                log('stayonTop attempt failed: ' + repr(e))
        except Exception as e:
            log('Activation routine exception: ' + repr(e))

    # schedule repeated activations
    for i in range(5):
        QTimer.singleShot(100 + i*200, try_activate)

    log('Entering event loop; close window to exit')
    app.exec()
    log('Exited')

if __name__ == '__main__':
    main()

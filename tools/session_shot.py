#!/usr/bin/env python3
"""Launch the real SessionManagerWindow, grab a screenshot, and save it.

Intended to run under Xvfb with XCB platform. Saves to /tmp/session_shot.png
and also to tools/session_shot.png inside the workspace for easy access.
"""

from __future__ import annotations
import sys
import os
from qtpy.QtWidgets import QApplication
from qtpy.QtCore import QTimer

OUT_TMP = '/tmp/session_shot.png'
OUT_WORK = os.path.join(os.getcwd(), 'tools', 'session_shot.png')

def main():
    # ensure project root is on sys.path so `import qtgui` works
    proj_root = os.getcwd()
    if proj_root not in sys.path:
        sys.path.insert(0, proj_root)

    app = QApplication(sys.argv)
    try:
        from qtgui.main_app import SessionManagerWindow
    except Exception as e:
        print('Import failed:', e)
        raise

    w = SessionManagerWindow()
    w.show()

    def grab_and_exit():
        try:
            pix = w.grab()
            pix.save(OUT_TMP, 'PNG')
            # also save a copy into workspace tools/ for easier retrieval
            try:
                pix.save(OUT_WORK, 'PNG')
            except Exception as e2:
                print('Failed save to workspace path:', e2)
            print('Saved', OUT_TMP)
        except Exception as e:
            print('Grab failed:', e)
        QTimer.singleShot(200, app.quit)

    QTimer.singleShot(1200, grab_and_exit)
    app.exec()

if __name__ == '__main__':
    main()

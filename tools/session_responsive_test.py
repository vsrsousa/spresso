#!/usr/bin/env python3
"""Launch SessionManagerWindow and emit a heartbeat periodically to check responsiveness.

Prints timestamp heartbeats to stdout and updates the window title every 300ms.
If heartbeats stop, the Qt main thread is blocked.
"""
from __future__ import annotations
import sys, time, os
from qtpy.QtWidgets import QApplication
from qtpy.QtCore import QTimer

def main():
    proj = os.getcwd()
    if proj not in sys.path:
        sys.path.insert(0, proj)

    app = QApplication(sys.argv)
    try:
        from qtgui.main_app import SessionManagerWindow
    except Exception as e:
        print('Import failed:', e)
        raise

    w = SessionManagerWindow()
    w.show()

    def heartbeat():
        t = time.time()
        print(f'heartbeat {t:.3f}', flush=True)
        try:
            w.setWindowTitle(f"xespresso - heartbeat {int(t*1000)%100000}")
        except Exception:
            pass

    timer = QTimer()
    timer.timeout.connect(heartbeat)
    timer.start(300)

    # also schedule a final timer to quit after 30s so we don't leave it running
    QTimer.singleShot(30000, app.quit)
    app.exec()

if __name__ == '__main__':
    main()

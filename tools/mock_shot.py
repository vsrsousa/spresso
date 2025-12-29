#!/usr/bin/env python3
"""Show a small Qt window, grab it to a PNG, then exit.

Designed to run under Xvfb. Saves to /tmp/mock_shot.png.
"""
from __future__ import annotations
import sys
import time
from qtpy.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel
from qtpy.QtCore import QTimer

OUT = '/tmp/mock_shot.png'

def main():
    app = QApplication(sys.argv)
    w = QWidget()
    w.setWindowTitle('MockShot')
    w.resize(640, 480)
    l = QVBoxLayout(w)
    l.addWidget(QLabel('This is a mock screenshot'))
    w.show()

    def grab_and_exit():
        try:
            pix = w.grab()
            pix.save(OUT, 'PNG')
            print('Saved', OUT)
        except Exception as e:
            print('Grab failed:', e)
        QTimer.singleShot(200, app.quit)

    QTimer.singleShot(800, grab_and_exit)
    app.exec()

if __name__ == '__main__':
    main()

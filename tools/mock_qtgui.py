#!/usr/bin/env python3
"""Minimal mock of the qtgui Session Manager to reproduce startup issues.

Run with:
  python3 tools/mock_qtgui.py [--lazy] [--heavy-init] [--block-eventfilter]

Options:
  --lazy            Create lightweight placeholders and instantiate pages on demand.
  --heavy-init      Simulate heavy imports by sleeping during page initialization.
  --block-eventfilter  Install an eventFilter that does blocking work (to reproduce hangs).

This script logs to /tmp/mock_qtgui.log and prints useful messages.
"""

from __future__ import annotations
import sys
import os
import time
import argparse
from qtpy.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QLabel, QListWidget,
    QStackedWidget, QListWidgetItem, QPushButton, QHBoxLayout
)
from qtpy.QtCore import Qt, QTimer, QObject, QEvent

LOG = '/tmp/mock_qtgui.log'


def log(msg: str):
    with open(LOG, 'a') as f:
        f.write(msg + '\n')
    print(msg)


class BlockingFilter(QObject):
    def eventFilter(self, obj, event):
        # Intentionally blocking/sleeping on window activate to reproduce hang
        if event.type() == QEvent.WindowActivate:
            log('BlockingFilter: WindowActivate received — sleeping 2s')
            time.sleep(2.0)
        return super().eventFilter(obj, event)


class MockPage(QWidget):
    def __init__(self, name: str, heavy: bool = False):
        super().__init__()
        self.name = name
        layout = QVBoxLayout(self)
        label = QLabel(f"Page: {name}")
        layout.addWidget(label)
        if heavy:
            log(f'Initializing heavy page {name} (sleep 1s)')
            time.sleep(1.0)
            log(f'Initialized heavy page {name}')


class MockSessionManager(QMainWindow):
    def __init__(self, lazy: bool = False, heavy_init: bool = False, block_eventfilter: bool = False):
        super().__init__()
        self.setWindowTitle('Mock Session Manager')
        self.resize(800, 600)
        self.lazy = lazy
        self.heavy_init = heavy_init
        self.block_eventfilter = block_eventfilter

        central = QWidget()
        self.setCentralWidget(central)
        layout = QHBoxLayout(central)

        self.sid = QListWidget()
        self.sid.addItem(QListWidgetItem('Home'))
        self.sid.addItem(QListWidgetItem('Workflow'))
        self.sid.addItem(QListWidgetItem('Settings'))
        self.sid.currentRowChanged.connect(self._on_nav)
        layout.addWidget(self.sid)

        self.stack = QStackedWidget()
        layout.addWidget(self.stack, 1)

        self.pages = []
        if self.lazy:
            # placeholders
            for name in ('Home', 'Workflow', 'Settings'):
                ph = QLabel(f'Placeholder for {name}')
                ph._page_name = name
                self.stack.addWidget(ph)
                self.pages.append(None)
        else:
            for name in ('Home', 'Workflow', 'Settings'):
                p = MockPage(name, heavy=self.heavy_init)
                self.stack.addWidget(p)
                self.pages.append(p)

        # select first
        self.sid.setCurrentRow(0)

        if self.block_eventfilter:
            f = BlockingFilter(self)
            QApplication.instance().installEventFilter(f)
            log('Installed BlockingFilter (will sleep on WindowActivate)')

    def _ensure_page(self, index: int):
        if self.pages[index] is None:
            name = self.stack.widget(index)._page_name
            log(f'Lazy instantiating page {name}')
            p = MockPage(name, heavy=self.heavy_init)
            self.stack.insertWidget(index, p)
            old = self.stack.widget(index+1)
            try:
                self.stack.removeWidget(old)
            except Exception:
                pass
            self.pages[index] = p

    def _on_nav(self, index: int):
        if index < 0:
            return
        if self.lazy:
            self._ensure_page(index)
        self.stack.setCurrentIndex(index)


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    p = argparse.ArgumentParser()
    p.add_argument('--lazy', action='store_true')
    p.add_argument('--heavy-init', action='store_true')
    p.add_argument('--block-eventfilter', action='store_true')
    p.add_argument('--stay-open', action='store_true', help='Do not auto-quit; keep app running until user closes')
    args = p.parse_args(argv)

    if os.path.exists(LOG):
        os.remove(LOG)
    log('Starting mock_qtgui with args: ' + str(args))

    app = QApplication(sys.argv)
    w = MockSessionManager(lazy=args.lazy, heavy_init=args.heavy_init, block_eventfilter=args.block_eventfilter)
    w.show()
    log('Shown mock window')

    # allow user to see window for a short time when run in headless tests
    if not args.stay_open:
        QTimer.singleShot(2500, app.quit)
    else:
        log('Running in --stay-open mode; waiting for user to close window')
    app.exec()
    log('Exited')


if __name__ == '__main__':
    main()

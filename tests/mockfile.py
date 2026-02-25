#!/usr/bin/env python3
"""Mock module runnable with: python -m mockfile

Subcommands:
  info           Print environment and repo info
  run            Run a short mock job (optionally wait)
  serve          Start a tiny HTTP server (for quick local testing)

Examples:
  python -m mockfile info
  python -m mockfile run --wait 3
  python -m mockfile serve --port 8001
"""

from __future__ import annotations
import argparse
import sys
import os
import time
import platform


def cmd_info(args: argparse.Namespace) -> int:
    print("mockfile: environment info")
    print("cwd:", os.getcwd())
    print("python:", sys.executable, sys.version.splitlines()[0])
    print("platform:", platform.platform())
    try:
        import git
        print('git available')
    except Exception:
        print('git python library not available')
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    print('mockfile: starting run')
    if args.wait and args.wait > 0:
        print(f'holding for {args.wait} seconds (simulate work)')
        try:
            time.sleep(args.wait)
        except KeyboardInterrupt:
            print('interrupted')
            return 2
    print('mockfile: run finished')
    return 0


def cmd_serve(args: argparse.Namespace) -> int:
    port = args.port or 8000
    bind = args.bind or '127.0.0.1'
    print(f'mockfile: starting HTTP server on {bind}:{port} (ctrl-c to stop)')
    try:
        from http.server import SimpleHTTPRequestHandler
        from socketserver import TCPServer

        class ReuseTCPServer(TCPServer):
            allow_reuse_address = True

        handler = SimpleHTTPRequestHandler
        with ReuseTCPServer((bind, port), handler) as httpd:
            try:
                httpd.serve_forever()
            except KeyboardInterrupt:
                print('\nserver stopped')
                return 0
    except Exception as exc:
        print('failed to start server:', exc)
        return 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog='python -m mockfile')
    sub = p.add_subparsers(dest='cmd')

    s_info = sub.add_parser('info', help='Show environment info')
    s_info.set_defaults(func=cmd_info)

    s_run = sub.add_parser('run', help='Run a short mock job')
    s_run.add_argument('--wait', '-w', type=float, default=0.0, help='Seconds to sleep (simulate work)')
    s_run.set_defaults(func=cmd_run)

    s_serve = sub.add_parser('serve', help='Start a tiny HTTP server')
    s_serve.add_argument('--port', '-p', type=int, default=8000, help='Port to serve on')
    s_serve.add_argument('--bind', '-b', type=str, default='127.0.0.1', help='Bind address')
    s_serve.set_defaults(func=cmd_serve)

    return p


def main(argv: list[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, 'func'):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == '__main__':
    raise SystemExit(main())

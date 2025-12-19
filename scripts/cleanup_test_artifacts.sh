#!/usr/bin/env bash
# Remove known transient test artifact folders created during local test runs
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

rm -rf r1_label t1_label t_save_label Co4Gd2 test

echo "Removed transient test artifact folders (if present)." 

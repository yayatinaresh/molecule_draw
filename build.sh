#!/usr/bin/env bash
# One-folder build of Scheme Molecules (run on the OS you want a build for).
# Needs the python env that has rdkit and pyside6 installed.
set -euo pipefail
cd "$(dirname "$0")"

python -m pip install --quiet pyinstaller
rm -rf build dist SchemeMolecules.spec

python -m PyInstaller --noconfirm --windowed --name SchemeMolecules scheme_gui.py

echo "Smoke test of the built program..."
QT_QPA_PLATFORM=offscreen ./dist/SchemeMolecules/SchemeMolecules --selftest

(cd dist && tar czf SchemeMolecules-linux.tar.gz SchemeMolecules)
echo "Done: dist/SchemeMolecules-linux.tar.gz  (run dist/SchemeMolecules/SchemeMolecules)"

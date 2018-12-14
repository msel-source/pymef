#!/usr/bin/env sh

set -e

for PYVER in "3.6.3"; do
  echo "befor python install"
  /c/ProgramData/chocolatey/bin/choco --version
  /c/ProgramData/chocolatey/bin/choco install python3 --version 3.6.3
  echo "after python install"
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

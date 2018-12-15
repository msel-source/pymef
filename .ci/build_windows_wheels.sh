#!/usr/bin/env bash

set -e

for PYVER in "3.6.3"; do
  choco install python3 --version ${PYVER}
  export PATH="/c/Python36:/c/Python36/Scripts:$PATH"
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

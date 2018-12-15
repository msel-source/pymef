#!/usr/bin/env bash

set -e

for PYVER in "3.6.3"; do
  choco install python3 --version ${PYVER}
  refreshenv
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

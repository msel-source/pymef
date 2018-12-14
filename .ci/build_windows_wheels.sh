#!/usr/bin/env bash

set -e

for PYVER in "3.6.3"; do
  choco install python --version ${PYVER}
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

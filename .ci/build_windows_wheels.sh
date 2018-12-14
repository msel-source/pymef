#!/usr/bin/env sh

set -e

for PYVER in "3.6.3"; do
  echo "befor python install"
  choco install python
  echo "after python install"
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

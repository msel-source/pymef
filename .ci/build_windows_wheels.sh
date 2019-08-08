#!/usr/bin/env bash

set -e

for PYVER in "3.6.3" "3.7.4"; do
  choco install python3 --version ${PYVER}
  if [ $PYVER == '3.6.3']
    export PATH="/c/Python36:/c/Python36/Scripts:$PATH"
  else
    export PATH="/c/Python37:/c/Python37/Scripts:$PATH"
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

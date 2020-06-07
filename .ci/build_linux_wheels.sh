#!/usr/bin/env bash

set -e

for PYVER in "cp36-cp36m" "cp37-cp37m" "cp38-cp38m"; do
  PYBIN="/opt/python/${PYVER}/bin"
  "${PYBIN}/pip" install numpy
  "${PYBIN}/pip" install wheel
  "${PYBIN}/python" setup.py test
  "${PYBIN}/python" setup.py bdist_wheel
done
find dist -name "*.whl" -exec auditwheel repair {} \;

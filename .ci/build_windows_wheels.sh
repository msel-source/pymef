#!/usr/bin/env bash

set -e

declare -a ver_arr=("3.6.3" "3.7.4" "3.8.3")

for PYVER in "${ver_arr[@]}"; do
  choco install python3 --version ${PYVER}
  if [ $PYVER == "3.6.3" ]; then
    export PATH="/c/Python36:/c/Python36/Scripts:$PATH"
  elif [ $PYVER == "3.7.4" ]; then
    export PATH="/c/Python37:/c/Python37/Scripts:$PATH"
  elif [ $PYVER == "3.8.3" ]; then
  	export PATH="/c/Python38:/c/Python38/Scripts:$PATH"
  fi
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

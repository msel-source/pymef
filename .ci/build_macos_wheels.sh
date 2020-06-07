#!/usr/bin/env bash

set -e

brew update
brew list pyenv &>/dev/null || brew install pyenv
brew update && brew upgrade pyenv
export PATH=~/.pyenv/shims:$PATH
for PYVER in "3.6.3" "3.7.4" "3.8.0"; do
  pyenv install ${PYVER}
  pyenv global ${PYVER}
  python -m pip install numpy wheel
  python setup.py test
  python setup.py bdist_wheel
done

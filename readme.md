![Tests](https://github.com/msel-source/pymef/actions/workflows/test_publish.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pymef/badge/?version=latest)](https://pymef.readthedocs.io/en/latest/?badge=latest)

Pymef
====

Pymef is a wrapper library for Multiscale Electrophysiology Format developed by 
[MSEL laboratory](http://msel.mayo.edu/).

Currently available for all major distributions (Linux, Mac OS, Windows). Only python 3 is supported.


Mef v 3.0 basic features
------------------------

-   Support for parallelisation of signal processing
-   Data compression
-   Data encryption
-   Real-time read/write, failure when writing file leaves intact valid files
-   CRC functionality to detect data corruption
-   Support for time discontinuities
-   Support for time series and video channels

Wrapper features
----------------

-   MEF3 files write/read
-   Convenience functions to easily read data and metadata for multiple channels

Installation
------------

To install please use:
```bash
pip install pymef
```

To install from source:
```bash
python setup.py install
```

Usage
------------
```
from pymef.mef_session import MefSession

session_path = '/path/to/session.mefd'
password     = 'mef_password'          // leave blank if no password

# read session metadata
ms = MefSession(session_path, password)

# read data of a single channel from beginning to end
data = ms.read_ts_channels_sample('Ch01', [[None, None]])

# read data of multiple channels from beginning to end
data = ms.read_ts_channels_sample(['Ch01', 'Ch05'], [[None, None]])
```

Documentation
-------------

The MEF3 specification can be found [here](https://osf.io/e3sf9/download).
The PyMef documentation can be found [here](https://pymef.readthedocs.io).

Support
-------

Please report problems to jan.cimbalnik@fnusa.cz.

License
-------

Pymef is licensed under the Apache software license. See [LICENSE.txt](./LICENSE.txt) for details.

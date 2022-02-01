![Tests](https://github.com/msel-source/pymef/actions/workflows/test_publish.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pymef/badge/?version=latest)](https://pymef.readthedocs.io/en/latest/?badge=latest)

Pymef
====

Pymef is a wrapper library for Multiscale Electrophysiology Format developed by 
[MSEL laboratory](http://msel.mayo.edu/).

Currently available for all major distributions (Linux, Mac OS, Windows). Only python 3 is supported.

For smooth usage plese see [documentation](https://pymef.readthedocs.io)

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

Documentation
-------------

Documentation of mef library can be found [here](http://msel.mayo.edu/codes.html).

Support
-------

Please report problems to jan.cimbalnik@fnusa.cz.

License
-------

Pymef is licensed under the Apache software license. See [LICENSE.txt](./LICENSE.txt) for details.

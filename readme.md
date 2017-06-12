Pymef
====

Pymef is a wrapper library for Mayo Electrophysiology Format developed by 
[MSEL laboratory](http://msel.mayo.edu/).

Currently available for linux and Mac OS X.

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
python setup.py install
```
(Note: Will be available through pip in the future.)

Documentation
-------------

Documentation of mef library can be found [here](http://msel.mayo.edu/codes.html).

Support
-------

Please report problems to cimbalnik.jan@mayo.edu.

License
-------

Glue is licensed under the Apache software license. See [LICENSE.txt](./LICENSE.txt) for details.

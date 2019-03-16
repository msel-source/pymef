Welcome to PyMef's documentation!
=================================

Pymef is a wrapper library for Multiscale Electrophysiology Format developed by 
[MSEL laboratory](http://msel.mayo.edu/).

Currently available for all major distributions (Linux, Mac OS, Windows). Only python 3 is supported.

Mef v 3.0 basic features
-------------------------------
-   Support for parallelisation of signal processing
-   Data compression
-   Data encryption
-   Real-time read/write, failure when writing file leaves intact valid files
-   CRC functionality to detect data corruption
-   Support for time discontinuities
-   Support for time series and video channels

Wrapper features
-------------------------
-   MEF3 files write/read
-   Convenience functions to easily read data and metadata for multiple channels
-   Functions to detect and repair corrupt data

License
--------------
Pymef is licensed under the Apache software license.



.. toctree::
   :maxdepth: 3
   
   installation
   getting_started
   api



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

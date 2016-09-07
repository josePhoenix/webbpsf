jwxml
================

*Note: This is probably mostly useful internally at STScI.*

Misc XML support helper code related to some JWST config files, mostly optics and wavefront sensing related. 


So far this provides support for reading in and interacting with: 

  * The so-called SIAF, Science Instrument Aperture Files
  * Wavefront control Segment Update Request files.


More functionality may be added ad hoc as needed; no overall long term development master plan is implied. 




SIAFs
------

Science Instrument Aperture Files contain detailed focal plane and pointing models for the science instruments. 


.. code:: python

    siaf = jwxml.SIAF(filename='/your/path/to/NIRCam_SIAF_2016-01-28.xml')


If you happen to have `webbpsf` also installed, you can just ask for an instrument by name and it
will grab the appropriate file from the `$WEBBPSF_DATA` directory tree. 

.. code:: python

    siaf = jwxml.SIAF('NIRCam')



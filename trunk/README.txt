PyGuide: Python software for telescope guiding by Russell E. Owen
based on algorithms and code by Jim Gunn
and help from Connie Rockosi and Craig Loomis.

To install, use the usual "distutils" procedure:
% python setup.py build
% python setup.py install
or, if that fails due to permission errors:
% sudo python setup.py install

Once the package is installed, you can get help in the usual way:
>>>import PyGuide
>>>help(PyGuide)
Before installation, you can read some information in:
python/PyGuide/__init__.py

To Do as of 2004-08-03:
- Centroid has a mediocre error estimate. It is often too optimistic for large slits, and may sometimes be too pessimistic. I am working on a better estimate.
- Depending on the final centroid error estimate, the various interfaces are likely to change. For instance it may be necessary to perform a star shape fit to get a good estimate of the centroid error. Or it may be necessary to know CCD characteristics such as read noise. Or it may be necessary to provide an estimate of star FWHM.
- The star shape fitter does not fit ellipticity.
- I am not happy having one object for centroid data (CentroidData) and another for star shape data (StarShapeData). At this point I have no idea if this will change.

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

The main routines are:
PyGuide.findStars
PyGuide.centroid
PyGuide.starShape: fits a symmetrical double gaussian
PyGuide.ImUtil: includes skyStats and subFrameCtr utilities
PyGuide.FakeData: construct fake stars with realistic noise

Other files of interest (these are not installed):
- History.txt: version history
- checkPyGuide: runs pychecker on PyGuide.
  requires that pychecker be installed.
- test/...: code to check the various routines and sample output.

Known Limitations:
- The centroid error estimate tends to over-esimate the error
  but not in a linear enough fashion to be trivially fixed.
- The star shape fitter does not fit ellipticity.
- FakeData does not add cosmic rays.
- None of the routines deal with dark current. As long as dark current
  is not a dominant source of noise this should be acceptable.
- The asymmetry reported by centroid is probably not sufficiently well normalized
  to be a useful measure of anything. It tends to get large for bright objects
  that are mostly masked off.

PyGuide: Python software for telescope guiding by Russell E. Owen
based on algorithms and code by Jim Gunn
and help from Connie Rockosi and Craig Loomis.

Required Packages:
- numarray (1.1 or later recommended)
- If you want findStars to display images in ds9 then RO.DS9 is required.
  The RO package is available from http://astro.washington.edu/rowen/

To install, use the usual "distutils" procedure:
% python setup.py install
or, if that fails due to permission errors:
% sudo python setup.py install

Once the package is installed, you can get help in the usual way:
>>> import PyGuide
>>> help(PyGuide)

The main routines are:
PyGuide.findStars
PyGuide.centroid
PyGuide.starShape: fits a symmetrical double gaussian
PyGuide.ImUtil: includes skyStats and subFrameCtr utilities
PyGuide.FakeData: construct fake stars with realistic noise

Other files of interest (these are not installed):
- History.txt: version history
- checkPyGuide: runs pychecker on PyGuide (if pychecker is installed).
- doFindStars.py: scan an image file, display in ds9
  and report star centroids and shape info.
- test/...: code to check the various routines and sample output.

Known Limitations:
- The centroid error estimate tends to over-esimate the error
  but not in a sufficiently predictable fashion to be trivially fixed.
- The star shape fitter does not fit ellipticity.
- FakeData does not add cosmic rays.
- None of the routines deal with dark current. As long as dark current
  is not a dominant source of noise this should be acceptable.
- The asymmetry reported by centroid is probably not sufficiently well normalized
  to be a useful measure of anything. It tends to get large for bright objects
  that are mostly masked off.

Copyright 2004 Russell Owen

This file is part of the PyGuide package.

PyGuide is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

PyGuide is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PyGuide; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

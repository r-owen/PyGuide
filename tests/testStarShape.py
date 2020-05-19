#!/usr/bin/env python
"""Test the star shape code.

Warning; this is a very minimal test.

To do:
- try on generated data
- try on data with pieces missing

History:
2004-05-20 ROwen
2004-07-02 ROwen    Added optional image input and ds9 display of data
2004-08-04 ROwen    Modified to work with 2004-08-03 findStars.
2004-08-06 ROwen    Modified to work with 2004-08-06 findStars.
2005-02-07 ROwen    Modified for PyGuide 1.2.
                    Modified to show shape info for all found stars.
2005-04-19 ROwen    Modified for PyGuide 2.0
2005-10-14 ROwen    Supply null satMask for PyGuide 2.1.
"""
from __future__ import print_function
from __future__ import unicode_literals
import os.path
import sys
import PyGuide
import pyfits

# these values are probably wrong for the given test image
CCDInfo = PyGuide.CCDInfo(
    bias = 200,         # image bias, in ADU
    readNoise = 21.391, # read noise, in e-
    ccdGain = 1.643,    # ccd gain, in e-/pixel
)

UseDS9 = True

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    testDir = os.path.dirname(__file__)
    filename = os.path.join(testDir, "test.fits")

testimg = pyfits.open(filename)
data = testimg[0].data

if UseDS9:
    ds9Win = PyGuide.ImUtil.openDS9Win("testStarShape")
    if ds9Win:
        ds9Win.showArray(data)

mask = None
# use for image /Test\ Images/ecam/2004-04-27/gimg0170.fits
#mask = data < 0
#mask[64:101, 78:92] = 1

print("searching for stars")
ctrDataList, imStats = PyGuide.findStars(
    data = data,
    mask = mask,
    satMask = None,
    ccdInfo = CCDInfo,
    verbosity = 0,
    doDS9 = False,
)[0:2]
for ctrData in ctrDataList:
    xyCtr = ctrData.xyCtr
    rad = ctrData.rad
    print("star xyCtr=%.2f, %.2f, radius=%s" % (xyCtr[0], xyCtr[1], rad))
    
    shapeData = PyGuide.starShape(
        data,
        mask = mask,
        xyCtr = xyCtr,
        rad = rad,
    )
    if not shapeData.isOK:
        print("starShape failed: %s" % (shapeData.msgStr,))
    else:
        print("star ampl=%.1f, fwhm=%.1f, bkgnd=%.1f, chiSq=%.2f" %\
            (shapeData.ampl,shapeData.fwhm, shapeData.bkgnd, shapeData.chiSq))
    print()

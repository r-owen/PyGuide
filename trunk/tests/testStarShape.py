#!/usr/local/bin/python
"""Test the star shape code.

Warning; this is a very minimal test.

To do:
- try on generated data
- try on data with pieces missing

History:
2004-05-20 ROwen
2004-07-02 ROwen	Added optional image input and ds9 display of data
2004-08-04 ROwen	Modified to work with 2004-08-03 findStars.
2004-08-06 ROwen	Modified to work with 2004-08-06 findStars.
"""
import sys
import numarray as num
import PyGuide
import pyfits

# these values are probably wrong for the given test image
Bias = 2176		# image bias, in ADU
ReadNoise = 19	# read noise, in e-
CCDGain = 2.1	# inverse ccd gain, in e-/ADU

# to enable debugging output:
PyGuide.StarShape._StarShapeDebug = True
#PyGuide.StarShape._FitRadProfDebug = True
#PyGuide.StarShape._FitRadProfIterDebug = True
UseDS9 = True

if len(sys.argv) > 1:
	filename = sys.argv[1]
else:
	filename = "test.fits"

testimg = pyfits.open(filename)
data = testimg[0].data

if UseDS9:
	import RO.DS9
	ds9 = RO.DS9.DS9Win("testStarShape")
	ds9.showArray(data)

mask = None
# use for image /Test\ Images/ecam/2004-04-27/gimg0170.fits
#mask = data < 0
#mask[64:101, 78:92] = 1

print "searching for brightest star"
isSat, ctrDataList = PyGuide.findStars(
	data = data,
	mask = mask,
	bias = Bias,
	readNoise = ReadNoise,
	ccdGain = CCDGain,
	verbosity = 0,
	ds9 = False,
)
ctrData = ctrDataList[0]
ijCtr = ctrData.ctr
rad = ctrData.rad
print "found star ijCtr=%.2f, %.2f, radius=%s" % (ijCtr[0], ijCtr[1], rad)

gsData = PyGuide.starShape(
	data,
	mask = mask,
	ijCtr = ijCtr,
	predFWHM = rad/2.0,
)
print "star ampl=%.1f, fwhm=%.1f, bkgnd=%.1f, chiSq=%.2f" %\
	(gsData.ampl, gsData.fwhm, gsData.bkgnd, gsData.chiSq)

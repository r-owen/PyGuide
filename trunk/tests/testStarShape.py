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
2005-02-07 ROwen	Modified for PyGuide 1.2.
					Modified to show shape info for all found stars.
"""
import sys
import numarray as num
import PyGuide
import pyfits

# these values are probably wrong for the given test image
bias = 200		# image bias, in ADU
readNoise = 21.391
ccdGain = 1.643 # e-/pixel

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

print "searching for stars"
isSat, ctrDataList = PyGuide.findStars(
	data = data,
	mask = mask,
	bias = bias,
	readNoise = readNoise,
	ccdGain = ccdGain,
	verbosity = 0,
	ds9 = False,
)
for ctrData in ctrDataList:
	try:
		xyCtr = ctrData.xyCtr
		rad = ctrData.rad
		print "star xyCtr=%.2f, %.2f, radius=%s" % (xyCtr[0], xyCtr[1], rad)
		
		gsData = PyGuide.starShape(
			data,
			mask = mask,
			xyCtr = xyCtr,
			predFWHM = rad/2.0,
		)
		print "star ampl=%.1f, fwhm=%.1f, bkgnd=%.1f, chiSq=%.2f" %\
			(gsData.ampl, gsData.fwhm, gsData.bkgnd, gsData.chiSq)
	except RuntimeError, e:
		print "Failed:", e
	print
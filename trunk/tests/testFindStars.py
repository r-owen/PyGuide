#!/usr/local/bin/python
"""Test the star finder.
At present this is a VERY MINIMAL test.

History:
2004-04-16 ROwen
2004-08-03 ROwen	Modified for FindStars 2004-08-03.
2004-08-06 ROwen	Modified to display results using ds9.
"""
import numarray
import PyGuide
import pyfits

# these values are probably wrong for the given test image
Bias = 2176		# image bias, in ADU
ReadNoise = 19	# read noise, in e-
CCDGain = 2.1	# inverse ccd gain, in e-/ADU

testimg = pyfits.open("test.fits")
data = testimg[0].data
#	data.transpose()
mask = None

isSat, starData = PyGuide.findStars(
	data = data,
	mask = mask,
	bias = Bias,
	readNoise = ReadNoise,
	ccdGain = CCDGain,
	ds9=True,
)
print "%s stars found; isSaturated = %s:" % (len(starData), isSat)
print "xctr	yctr	xerr	yerr	pix	counts	rad"
for sd in starData:
	print "%6.2f	%6.2f	%6.2f	%6.2f	%4d	%7.0f	%3d" % \
		(sd.ctr[1], sd.ctr[0], sd.err[1], sd.err[0], sd.pix, sd.counts, sd.rad)

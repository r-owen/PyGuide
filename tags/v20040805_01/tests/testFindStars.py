#!/usr/local/bin/python
"""Test the star finder.
At present this is a VERY MINIMAL test.

History:
2004-04-16 ROwen
2004-08-03 ROwen	Modified for FindStars 2004-08-03.
"""
import numarray
import PyGuide
import pyfits

testimg = pyfits.open("test.fits")
data = testimg[0].data
#	data.transpose()
mask = None

isSat, starData = PyGuide.findStars(data, mask)
print "%s stars found; isSaturated = %s:" % (len(starData), isSat)
print "xctr	yctr	xerr	yerr	pix	counts	rad"
for sd in starData:
	print "%6.2f	%6.2f	%6.2f	%6.2f	%4d	%7.0f	%3d" % \
		(sd.ctr[1], sd.ctr[0], sd.err[1], sd.err[0], sd.pix, sd.counts, sd.rad)

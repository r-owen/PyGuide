#!/usr/local/bin/python -i
"""Measures stars in a given image file, displaying the image in ds9
and reporting star positions and shape information on stdout.

History:
2004-04-16 ROwen
2004-04-29 ROwen	Modified to use new ds9 on findStars.
2004-05-18 ROwen	Modified to set up ds9Win and to use fewer globals.
2004-08-25 ROwen	Modified for 2004-08-06 PyGuide.
2004-10-14 ROwen	Modified to measure starShape.
2004-12-01 ROwen	Renamed function from starUtil to doFindStars to match module name.
					Replaced arguments with globals to make it easier to change settings.
					Bug fix: if starShape failed, shapeData was mis-set.
"""
import numarray as num
import PyGuide
import pyfits
import RO.DS9

im = None
d = None
mask = None
verbosity = 1
ds9 = True
dataCut = 3.0
satLevel = 2**16
radMult = 1.0
ds9Win = RO.DS9.DS9Win(PyGuide.FindStars._DS9Title)

# new NA2 guider
bias = 1780
readNoise = 21.391
ccdGain = 1.643 # e-/pixel

def doFindStars(
	filename = None,
):
	global im, d, isSat, sd
	if filename:
		im = pyfits.open(filename)
	d = im[0].data
	
	# find stars and centroid
	isSat, posDataList = PyGuide.findStars(
		data = d,
		mask = mask,
		bias = bias,
		readNoise = readNoise,
		ccdGain = ccdGain,
		dataCut = dataCut,
		satLevel = satLevel,
		radMult = radMult,
		verbosity = verbosity,
		ds9 = ds9,
	)

	print "%s stars found; isSaturated = %s:" % (len(posDataList), isSat)
	print "   xctr	   yctr	   xerr	   yerr		 ampl	  bkgnd	   fwhm	 |  rad	    pix	  chiSq"
	for posData in posDataList:
		# measure star shape
		try:
			shapeData = PyGuide.starShape(
				data = d,
				mask = mask,
				ijCtr = posData.ctr,
				predFWHM = posData.rad,
			)
		except RuntimeError, e:
			print "starShape failed: %s" % (e,)
			shapeData = PyGuide.StarShapeData()
		
		# print results
		print "%7.2f	%7.2f	%7.2f	%7.2f	%13.1f	%7.1f	%7.1f	%7d	%7d	%7.1f" % (
			posData.ctr[1], posData.ctr[0],
			posData.err[1], posData.err[0],
			shapeData.ampl, shapeData.bkgnd, shapeData.fwhm,
			posData.rad, posData.pix, shapeData.chiSq,
		)

print "global variables:"
print "im =", im
print "mask =", mask
print "bias =", bias
print "readNoise =", readNoise
print "ccdGain =", ccdGain
print "dataCut =", dataCut
print "satLevel =", satLevel
print "radMult =", radMult
print "verbosity =", verbosity
print "ds9 =", ds9
print
print """ds9Win.showArray(arry) will display an array

Computed values:
im: a fits image; if filename is specified, im is loaded from that file;
    if filename is omitted, the current im is used
d: the data from the image (a numarray array)
sd: star data returned by PyGuide.findStars

Reported values include::
rad: radius used to compute centroid
pix: the number of unmasked pixels used to compute the centroid
chiSq: chi squared for shape fit

Notes:
- For a slitviewer image, be sure to specify a suitable mask.
- For optimal centroiding and a reasonable centroid error estimate
  you must set bias, readNoise and ccdGain correctly for your image.

Function call:
doFindStars(filename=None); if you omit filename then the current im is used
"""


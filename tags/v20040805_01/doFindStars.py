#!/usr/local/bin/python -i
"""Exercise findStars

History:
2004-04-16 ROwen
2004-04-29 ROwen	Modified to use new ds9 on findStars.
2004-05-18 ROwen	Modified to set up ds9Win and to use fewer globals.
"""
import numarray as num
import PyGuide
import pyfits
import RO.DS9

mask = None
verbosity = 1
ds9 = True
dataCut = 4.5
satLevel = 2**16
radMult = 1.0
ds9Win = RO.DS9.DS9Win(PyGuide.FindStars._DS9Title)

def dofind(
	filename,
	mask=mask,
	dataCut=dataCut,
	satLevel=satLevel,
	radMult=radMult,
	verbosity=verbosity,
	ds9=ds9,
):
	global im, d, isSat, sd
	im = pyfits.open(filename)
	d = im[0].data
	isSat, sd = PyGuide.findStars(
		data = d,
		mask = mask,
		dataCut=dataCut,
		satLevel = satLevel,
		radMult=radMult,
		verbosity=verbosity,
		ds9=ds9,
	)

	print "Found %s stars:" % len(sd)
	if isSat:
		print "WARNING: some pixels are saturated!"
	print "x ctr\ty ctr\t    counts\t  rad\tpoints"
	for counts, ctr, rad, totPts in sd:
		print "%.1f\t%.1f\t%10.0f\t%5.1f\t%6d" % (ctr[1], ctr[0], counts, rad, totPts)

print "Defaults for dofind:"
print "mask =", mask
print "dataCut =", dataCut
print "satLevel =", satLevel
print "radMult =", radMult
print "verbosity =", verbosity
print "ds9 =", ds9
print
print "ds9Win.showArray(arry) will display an array"
print
print "Computed values:"
print "im: the fits image (including header and data)"
print "d: the data from the image (a numarray array)"
print "sd: star data returned by PyGuide.findStars"
print
print "Function call:"
print "dofind(filename, ...)"
print

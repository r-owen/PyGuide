#!/usr/local/bin/python -i
"""Exercise findStars

Warning: fr best results, set bias, readNoise and ccdGain
correctly for your image.

History:
2004-04-16 ROwen
2004-04-29 ROwen	Modified to use new ds9 on findStars.
2004-05-18 ROwen	Modified to set up ds9Win and to use fewer globals.
2004-08-25 ROwen	Modified for 2004-08-06 PyGuide.
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

# new NA2 guider
bias = 1780
readNoise = 21.391
ccdGain = 1.643 # e-/pixel

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
	isSat, starData = PyGuide.findStars(
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

	print "%s stars found; isSaturated = %s:" % (len(starData), isSat)
	print "xctr	yctr	xerr	yerr	pix	counts	rad"
	for sd in starData:
		print "%6.2f	%6.2f	%6.2f	%6.2f	%4d	%7.0f	%3d" % \
			(sd.ctr[1], sd.ctr[0], sd.err[1], sd.err[0], sd.pix, sd.counts, sd.rad)

print "Defaults for dofind:"
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

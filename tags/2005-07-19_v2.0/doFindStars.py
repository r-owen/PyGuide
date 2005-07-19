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
2005-02-07 ROwen	Modified for findStars 1.2.
2005-03-31 ROwen	Allowed one to specify a mask name to doFindStars.
					Modified doFindStars to allow specifying parameters as keyword arguments.
					Added showDef() method to display current defaults.
					Modified global vars so im, not d, is the image array.
2005-04-01 ROwen	Updated for modified FindStars and StarShape.
2005-04-11 ROwen	Modified to use PyGuide.Constants.DS9Title.
2005-04-22 ROwen	Added support for the rad argument.
2005-05-16 ROwen	Modified for overhauled findStars.
2005-05-20 ROwen	Added doCentroid.
2005-06-17 ROwen	Added an invertMask flag.
"""
import numarray as num
import PyGuide
import pyfits
import RO.DS9

im = None
imfits = None
mask = None
maskfits = None

# Default Parameters
# these settings are for the new NA2 guider
bias = 1780
readNoise = 21.391
ccdGain = 1.643 # e-/pixel
# these are general settings
thresh = 3.0
radMult = 1.0
rad = None
satLevel = (2**16)-1
verbosity = 1
doDS9 = True

# set up a ds9 window
ds9Win = RO.DS9.DS9Win(PyGuide.Constants.DS9Title)

CCDInfoNames = ("bias", "readNoise", "ccdGain", "satLevel")
FindParamNames = CCDInfoNames + ("thresh", "radMult", "rad", "verbosity", "doDS9")
CentroidParamNames = CCDInfoNames + ("thresh", "rad", "verbosity", "doDS9")

def doFindStars(
	imName = None,
	maskName = None,
	invertMask = False,
	**kargs
):
	global isSat, sd
	im, mask = loadFiles(imName, maskName, invertMask)
	
	# check keyword arguments
	for paramName in kargs:
		if paramName not in FindParamNames:
			raise RuntimeError("Invalid argument: %s" % (paramName,))
	
	# fill in defaults
	globalDict = globals()
	for paramName in FindParamNames:
		if paramName not in kargs:
			kargs[paramName] = globalDict[paramName]
	
	# split off ccd info
	ccdInfoDict = {}
	for paramName in CCDInfoNames:
		ccdInfoDict[paramName] = kargs.pop(paramName)
	ccdInfo = PyGuide.CCDInfo(**ccdInfoDict)
	
	# find stars and centroid
	ctrDataList, imStats = PyGuide.findStars(
		data = im,
		mask = mask,
		ccdInfo = ccdInfo,
	**kargs)

	print "%s stars found:" % (len(ctrDataList),)
	print "   xctr    yctr    xerr    yerr         ampl   bkgnd    fwhm  |  rad     pix    nSat  chiSq"
	for ctrData in ctrDataList:
		# measure star shape
		try:
			shapeData = PyGuide.starShape(
				data = im,
				mask = mask,
				xyCtr = ctrData.xyCtr,
				rad = ctrData.rad,
			)
		except RuntimeError, e:
			print "starShape failed: %s" % (e,)
			shapeData = PyGuide.StarShapeData()
		
		# print results
		print "%7.2f %7.2f %7.2f %7.2f %13.1f %7.1f %7.1f %7d %7d %7d %7.1f" % (
			ctrData.xyCtr[0], ctrData.xyCtr[1],
			ctrData.xyErr[0], ctrData.xyErr[1],
			shapeData.ampl, shapeData.bkgnd, shapeData.fwhm,
			ctrData.rad, ctrData.pix, ctrData.nSat, shapeData.chiSq,
		)

def doCentroid(
	imName = None,
	maskName = None,
	xyGuess = None,
	invertMask = False,
	**kargs
):
	global im, imfits, mask, maskfits, isSat, sd
	im, mask = loadFiles(imName, maskName, invertMask)
	if xyGuess == None:
		print "xyGuess is required"
		return
	
	# check keyword arguments
	for paramName in kargs:
		if paramName not in CentroidParamNames:
			raise RuntimeError("Invalid argument: %s" % (paramName,))
	
	# fill in defaults
	globalDict = globals()
	for paramName in CentroidParamNames:
		if paramName not in kargs:
			kargs[paramName] = globalDict[paramName]
	
	# split off ccd info
	ccdInfoDict = {}
	for paramName in CCDInfoNames:
		ccdInfoDict[paramName] = kargs.pop(paramName)
	ccdInfo = PyGuide.CCDInfo(**ccdInfoDict)
	
	# centroid
	ctrData = PyGuide.centroid(
		data = im,
		mask = mask,
		xyGuess = xyGuess,
		ccdInfo = ccdInfo,
	**kargs)

	if not ctrData.isOK:
		print "centroid failed:", ctrData.msgStr
		return

	print "   xctr    yctr    xerr    yerr         ampl   bkgnd    fwhm  |  rad     pix    nSat  chiSq"
	shapeData = PyGuide.starShape(
		data = im,
		mask = mask,
		xyCtr = ctrData.xyCtr,
		rad = ctrData.rad,
	)
	# print results
	print "%7.2f %7.2f %7.2f %7.2f %13.1f %7.1f %7.1f %7d %7d %7d %7.1f" % (
		ctrData.xyCtr[0], ctrData.xyCtr[1],
		ctrData.xyErr[0], ctrData.xyErr[1],
		shapeData.ampl, shapeData.bkgnd, shapeData.fwhm,
		ctrData.rad, ctrData.pix, ctrData.nSat, shapeData.chiSq,
	)

	if not shapeData.isOK:
		print "starShape failed:", shapeData.msgStr

def loadFiles(
	imName = None,
	maskName = None,
	invertMask = False,
):
	"""Load a new image and/or mask from a fits file.
	invertMask is ignored unless maskName is specified.
	"""
	global im, imfits, mask, maskfits, isSat, sd
	if imName:
		imfits = pyfits.open(imName)
		im = imfits[0].data
	if maskName:
		maskfits = pyfits.open(maskName)
		if invertMask:
			mask = maskfits[0].data < 0.1
		else:
			mask = maskfits[0].data > 0.1
	return im, mask

def showDef():
	"""Show current value of various global variables
	which are used as defaults for parameters of doFindStars.
	"""
	print "Global variables (defaults for doFindStars):"
	globalDict = globals()
	for paramName in FindParamNames:
		print "%s = %s" % (paramName, globalDict[paramName])
	print

showDef()
print """
Computed values:
im: image data array
mask: mask data array, or None if no mask
sd: star data returned by PyGuide.findStars

Reported values include::
rad: radius used to compute centroid
pix: the number of unmasked pixels used to compute the centroid
chiSq: chi squared for shape fit

Notes:
- For a slitviewer image, be sure to specify a suitable mask.
- For optimal centroiding and a reasonable centroid error estimate
  you must set bias, readNoise and ccdGain correctly for your image.
- invertMask is ignored unless maskName is specified.

Function calls:
doFindStars(imName=None, maskName=None, invertMask=False, [, optional_named_params])
doCentroid(imName=None, maskName=None, invertMask=False, xyGuess=(x,y) [, optional_named_params])
loadFiles(imName, maskname, invertMask) loads a new image and/or mask
ds9Win.showArray(arry) displays an array in ds9
showDef() prints the current defaults
"""


"""Find Stars.

WARNING: Python handles images as data[y,x]. This same index order is used
for positions, sizes and such. In an attempt to reduce confusion, the code
uses i,j instead of y,x.

Pixel convention: The point 0,0 is at the corner of first pixel read out.
Hence the center of that pixel is (0.5, 0.5) and the center of a 1024x1024 CCD
is (512.0, 512.0). This does not match the iraf and ds9 convention,
which has 0.5,0.5 as the corner and 1,1 as the center of the first pixel.

Uses an algorithm developed by Jim Gunn with some changes of my own:
- compute median and quartiles (Q1 & Q3) of the data
- stdDev = 0.741 * (Q3-Q1)
- repeat the computation, ignoring data above
  median + 0.235 * stdDev
  (iteration reduces the effects of bright pixels)
- median-filter a copy of the data (filling in masked values with median)
- candidate stars are connected blobs whose pixels
  have value > med + 5 * stdDev. Ignore blobs
  with width or height of only 1.
- centroid each blob.

This is an adaptation of code by Jim Gunn and Connie Rockosi.
The original code does a clever thing to reduce memory requirements:
instead of sorting the data to compute quartiles,
it creates a histogram of the data. This allows it to
allocate a single 65k array (for 16-bit data) and reuse that
for all images. But locating the quartiles is a bit more work.

To Do:
- Make use of the centroid error and minimum asymmetry
  to reject bad stars or sort starts in order of desirability?
- Detect and reject identical stars (one centroid within 1/2 pixel of another)?
  This isn't essential as duplicates will not hurt the guider.
  On the other hand, it's probably not terribly difficult to do, either.

History:
2004-04-16 ROwen	First release. Still pretty basic.
2004-04-20 ROwen	Added the dataCut argument and some debug output.
2004-04-21 ROwen	Added verbosity arg to findStars and removed _FindStarDebug.
2004-04-22 ROwen	Added test for saturation.
2004-04-30 ROwen	Added ds9 flag for visual debugging.
					Bug fix: was converting data to Int16 instead of UInt16.
2004-05-03 ROwen	Modified to auto-open new ds9 as needed.
2004-05-20 ROwen	Modified to return radius used for finding centroid.
					Modified radius calculation to be based on larger box dimension.
					Modified to display masked data in ds9.
2004-06-03 ROwen	Modified the module doc string.
2004-08-03 ROwen	Modified for the Centroid module's new output format.
2004-08-06 ROwen	Removed documentation for unused "rad" argument (thanks, Craig!).
					Stopped importing unused math and os modules.
					Fixed variable name error (only visible when ds9 True).
2004-08-06-v2 ROwen	Updated for new centroid code. Requires additional inputs.
2004-08-25 ROwen	Added __all__.
2004-10-14 ROwen	Changed default dataCut to 3.0 from 4.5.
					No longer displays a 3rd ds9 frame (was smoothed data).
2004-10-15 ROwen	No longer warns if RO.DS9 absent (unless you use the ds9 flag).
2005-02-07 ROwen	Modified to use x,y position convention defined by
					PyGuide.Constants.PosMinusIndex.
2005-03-29 ROwen	Bug fix: mis-displayed masked data in ds9.
"""
__all__ = ['findStars']

import numarray as num
import numarray.nd_image
import numarray.ma
import Centroid
import ImUtil
try:
	import RO.DS9
except ImportError:
	pass
_DS9Title = "FindStars"

def _fmtList(alist):
	"""Return "alist[0], alist[1], ..."
	"""
	return str(alist)[1:-1]

def _reversed(alist):
	"""Return a reversed copy of alist
	"""
	retval = list(alist[:])
	retval.reverse()
	return retval

def findStars(
	data,
	mask,
	bias,
	readNoise,
	ccdGain,
	dataCut = 3.0,
	satLevel = 2**16,
	radMult = 1.0,
	verbosity = 1,
	ds9 = False,
):
	"""Find and centroid stars.
	
	Inputs:
	- data		the image data [i,j]; this is converted to UInt16 if necessary
	- mask		a mask [i,j] of 0's (valid data) or 1's (invalid); None if no mask.
				If mask is specified, it must have the same shape as data.
	- bias		ccd bias (ADU)
	- readNoise	ccd read noise (e-)
	- ccdGain	ccd inverse gain (e-/ADU)
	- dataCut	determines the point above which pixels are considered data;
				cut level = median + dataCut * standard deviation
	- satLevel	The value at or above which a pixel is considered saturated (ADU)
	- radMult	centroid radius = radMult * max(box x rad, box y rad)
	- verbosity	0: no output, 1: print warnings, 2: print information
	- ds9		if True, shows current image and other info in ds9 in current frame.
				For this to work, you must have the RO package installed.
	
	Returns two items:
	- isSaturated	a flag indicating whether any stars were saturated
	- starData		a list of centroid information for unsaturated stars;
					each element is a PyGuide.CentroidData object whose fields include:
		- ctr		the i,j centroid (pixels)
		- err		the predicted i,j 1-sigma error (pixels)
		- asymm		measure of asymmetry (dimensionless)
		(see PyGuide.CentroidData for more information)
	
	Note: the found "stars" are not required to look star-like and so
	are not fit to a stellar profile. This routine was designed to handle
	spectrograph slit viewers and fiber bundles, where only a fragment
	of a star may be visible. It also handles donut-shaped images.
	"""
	if data.type() != num.UInt16:
		data = data.astype(num.UInt16)
	elif not data.iscontiguous():
		data = data.copy()
	maskedData = num.ma.array(data, mask=mask)
	
	if ds9:
		if "RO" not in globals():
			print 'RO.DS9 not available; ignoring ds9 flag'
		else:
			# if not already available, open new DS9 with specified template
			try:
				ds9Win = RO.DS9.DS9Win(_DS9Title)
				ds9Win.xpaset("tile frames")
				ds9Win.xpaset("frame 1")
				if mask != None:
					ds9Win.showArray(data * mask)
				else:
					ds9Win.showArray(data)
				ds9Win.xpaset("frame 2")
				ds9Win.showArray(data)
				ds9Win.xpaset("frame 1")
			except RuntimeError, e:
				print "Error communicating with ds9: %s" % (e,)
				ds9Win = None
	else:
		ds9Win = None
	
	isSaturated = False
	
	med, stdDev = ImUtil.skyStats(maskedData)

	# get a copy with the median used to fill in masked areas
	# and apply a filter to get rid of speckle
	smoothedData = maskedData.filled(med)
	num.nd_image.median_filter(smoothedData, 3, output=smoothedData)
#	if ds9Win:
#		ds9Win.xpaset("frame 3")
#		ds9Win.showArray(smoothedData)
#		ds9Win.xpaset("frame 1")
	
	# look for points larger than median + dataCut * stdDev
	dataCut = med + (dataCut * stdDev)
	shapeArry = num.ones((3,3))
	labels, numElts = num.nd_image.label(smoothedData>dataCut, shapeArry)

		
	smoothedData = None # release the storage
	if verbosity >= 2:
		print "findStars found %s possible stars above dataCut=%s" % (numElts, dataCut)

	# examine the candidate stars and compute centroids
	countsCentroidList = []
	slices = num.nd_image.find_objects(labels)
	for ind in range(len(slices)):
		ijSlc = slices[ind]
		ijSize = [slc.stop - slc.start for slc in ijSlc]
		ijCtrInd = [(slc.stop + slc.start) / 2.0 for slc in ijSlc]
		xyCtrGuess = ImUtil.xyPosFromIJPos(ijCtrInd)
		
		# reject regions only 1 pixel tall or wide
		if 1 in ijSize:
			# object is too small to be of interest
			if verbosity >= 1:
				print "findStars warning: candidate star at %s is too small; size=%s" % (xyCtrGuess, ijSize)
			continue
		
		# reject saturated regions and set isSaturated flag
		if num.nd_image.maximum(data, labels, ind) >= satLevel:
			isSaturated = True
			if verbosity >= 1:
				print "findStars warning: candidate star at %s is saturated" % (xyCtrGuess,)
			continue
		
		# region appears to be valid; centroid it
		rad = max(ijSize[0], ijSize[1]) * radMult / 2.0
		if ds9Win:
			ds9BoxCtr = ImUtil.ds9PosFromXYPos(xyCtrGuess)
			# display box from find_objects
			args = ds9BoxCtr + _reversed(ijSize) + [0]
			ds9Win.xpaset("regions", "image; box %s # group=findbox" % _fmtList(args))
			# display circle showing the centroider input
			args = ds9BoxCtr + [rad]
			ds9Win.xpaset("regions", "image; circle %s # group=ctrcirc" % _fmtList(args))
		try:
			if verbosity >= 2:
				print "findStars centroid at %s with rad=%s" % (xyCtrGuess, rad)
			ctrData = Centroid.centroid(
				data = data,
				mask = mask,
				xyGuess = xyCtrGuess,
				rad = rad,
				bias = bias,
				readNoise = readNoise,
				ccdGain = ccdGain,
			)
		except RuntimeError, e:
			if verbosity >= 1:
				print "findStars warning: centroid at %s with rad=%s failed: %s" % (xyCtrGuess, rad, e)
			continue
		countsCentroidList.append((ctrData.counts, ctrData))
		
		if ds9Win:
			# display x showing centroid
			args = ImUtil.ds9PosFromXYPos(ctrData.xyCtr)
			ds9Win.xpaset("regions", "image; x point %s # group=centroid" % _fmtList(args))
	
	# sort by decreasing counts
	countsCentroidList.sort()
	countsCentroidList.reverse()
	centroidList = [cc[1] for cc in countsCentroidList]
	if verbosity >= 2:
		print "findStars returning data for %s stars:" % len(centroidList)
		if isSaturated:
			print "WARNING: some pixels are saturated!"
		print "x ctr\ty ctr\tx err\ty err\t    counts\tpixels\tradius"
		for cd in centroidList:
			print "%.1f\t%.1f\t%.1f\t%.1f\t%10.0f\t%6d\t%5.1f" % \
				(cd.ctr[1], cd.ctr[0],
				 cd.err[1], cd.err[2],
				 cd.counts, cd.pix, cd.rad)
	return isSaturated, centroidList

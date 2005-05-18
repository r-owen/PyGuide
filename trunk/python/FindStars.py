"""Find Stars.

WARNING: Python handles images as data[y,x]. This same index order is used
for positions, sizes and such. In an attempt to reduce confusion, the code
uses i,j instead of y,x.

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

To Do:
- Make use of the centroid error and minimum asymmetry
  to reject bad stars or sort starts in order of desirability?
- Detect and reject identical stars (one centroid within 1/2 pixel of another)?
  This isn't essential as duplicates will not hurt the guider.
  On the other hand, it's probably not terribly difficult to do, either.
- Consider using a histogram to compute the quartiles (as the original
  code did). Thus one can allocate a single 65k array (for 16-bit data)
  and reuse that for all images. But locating the quartiles is a bit more work.

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
2005-03-31 ROwen	Modified to show smoothed, masked data in ds9 in frame 3
					if verbosity>=2 (and dS9 true), else no frame 3.
					Bug fix: error in star data output when verbosity >= 2.
2005-04-01 ROwen	Modified to return the median of the unmasked data.
2005-04-11 ROwen	Modified to use Constants.DS9Title.
2005-04-22 ROwen	Added rad argument (overrides radMult).
2005-05-17 ROwen	Overhauled findStars:
					- four args were consolidated into ccdInfo
					- renamed dataCut to thresh (to match the APO 3.5m hub)
					- renamed ds9 to doDS9
					- saturated stars are now centroided, and they are likely
					  to show up first in the list, so beware!
					- no longer returns isSaturated; the returned centroid data
					  now includes a count of the # of saturated pixels, instead
					- returns imStats instead of med
					Modified to use PyGuide.Constants.MinThresh
"""
__all__ = ['findStars']

import warnings
import numarray as num
import numarray.nd_image
import numarray.ma
import Centroid
import Constants
import ImUtil

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
	ccdInfo,
	thresh = 3.0,
	radMult = 1.0,
	rad = None,
	verbosity = 1,
	doDS9 = False,
):
	"""Find and centroid stars.
	
	Inputs:
	- data		the image data [i,j]; this is converted to an UInt16 numarray if necessary
	- mask		a mask [i,j] of 0's (valid data) or 1's (invalid); None if no mask.
				If mask is specified, it must have the same shape as data.
	- ccdInfo	bias, read noise, etc: a PyGuide.CCDInfo object.
	- thresh	determines the point above which pixels are considered data;
				valid data >= thresh * standard deviation + median
				values less than PyGuide.Constants.MinThresh are silently increased
	- radMult	centroid radius = radMult * max(box x rad, box y rad);
				ignored if rad specified
	- rad		centroid radius; if specified, overrides radMult
	- verbosity	0: no output, 1: print warnings, 2: print information and
				(if doDS9 true) show smoothed image in ds9 frame 3.
	- doDS9		if True, shows current image and other info in ds9 in current frame.
				For this to work, you must have the RO package installed.
	
	Returns two items:
	- starData		a list of centroid information for unsaturated stars;
					each element is a PyGuide.CentroidData object whose fields include:
		- ctr		the i,j centroid (pixels)
		- err		the predicted i,j 1-sigma error (pixels)
		- asymm		measure of asymmetry (dimensionless)
		(see PyGuide.CentroidData for more information)
	- imStats		median and stdDev masked data; an ImStats object
	
	Note: the found "stars" are not required to look star-like and so
	are not fit to a stellar profile. However, if the object is not
	circularly symmetric then the centroid it returns may not match
	the centroid you think it should return; this is especially a worry
	if the image has a large portion masked out.

	This routine was designed to handle	spectrograph slit viewers
	and fiber bundles, where only a fragment of a star may be visible.
	It also handles donut-shaped images.
	"""
	data = num.array(data)
	if data.type() != num.UInt16:
		data = data.astype(num.UInt16)
	elif not data.iscontiguous():
		data = data.copy()
	maskedData = num.ma.array(data, mask=mask)
	
	thresh = max(Constants.MinThresh, float(thresh))
	
	if doDS9:
		ds9Win = ImUtil.openDS9Win()
	else:
		ds9Win = None
			
	if ds9Win:
		# show masked data in frame 1 and unmasked data in frame 2
		ds9Win.xpaset("tile frames")
		ds9Win.xpaset("frame 1")
		if mask != None:
			ds9Win.showArray(data * (mask==0))
		else:
			ds9Win.showArray(data)
		ds9Win.xpaset("frame 2")
		ds9Win.showArray(data)
		ds9Win.xpaset("frame 1")
		
	med, stdDev = ImUtil.skyStats(maskedData)
	dataCut = med + (thresh * stdDev)

	imStats = Centroid.ImStats(
		thresh = thresh,
		med = med,
		stdDev = stdDev,
		dataCut = dataCut,
	)

	# get a copy with the median used to fill in masked areas
	# and apply a filter to get rid of speckle
	smoothedData = maskedData.filled(med)
	num.nd_image.median_filter(smoothedData, 3, output=smoothedData)
	if ds9Win and verbosity >= 2:
		ds9Win.xpaset("frame 3")
		ds9Win.showArray(smoothedData)
		ds9Win.xpaset("frame 1")
	
	# look for points larger than median + dataCut * stdDev
	shapeArry = num.ones((3,3))
	labels, numElts = num.nd_image.label(smoothedData>dataCut, shapeArry)
	smoothedData = None # release the storage
	if verbosity >= 2:
		print "findStars found %s possible stars above dataCut=%s" % (numElts, dataCut)

	# examine the candidate stars and compute centroids
	countsCentroidList = []
	slices = num.nd_image.find_objects(labels)
	for ijSlice in slices:
		ijSize = [slc.stop - slc.start for slc in ijSlice]
		ijCtrInd = [(slc.stop + slc.start) / 2.0 for slc in ijSlice]
		xyCtrGuess = ImUtil.xyPosFromIJPos(ijCtrInd)
		
		# reject regions only 1 pixel tall or wide
		if 1 in ijSize:
			# object is too small to be of interest
			if verbosity >= 1:
				print "findStars warning: candidate star at %s is too small; size=%s" % (xyCtrGuess, ijSize)
			continue
		
		# region appears to be valid; centroid it
		if rad == None:
			rad = max(ijSize[0], ijSize[1]) * radMult / 2.0
		if ds9Win:
			ds9BoxCtr = ImUtil.ds9PosFromXYPos(xyCtrGuess)
			# display box from find_objects
			args = ds9BoxCtr + _reversed(ijSize) + [0]
			ds9Win.xpaset("regions", "image; box %s # group=findbox" % _fmtList(args))
			# display circle showing the centroider input
			args = ds9BoxCtr + [rad]
			ds9Win.xpaset("regions", "image; circle %s # group=ctrcirc" % _fmtList(args))

		if verbosity >= 2:
			print "findStars centroid at %s with rad=%s" % (xyCtrGuess, rad)
		ctrData = Centroid.basicCentroid(
			data = data,
			mask = mask,
			xyGuess = xyCtrGuess,
			rad = rad,
			ccdInfo = ccdInfo,
			verbosity = verbosity,
		)
		if not ctrData.isOK:
			if verbosity >= 1:
				print "findStars warning: centroid at %s with rad=%s failed: %s" % (xyCtrGuess, rad, ctrData.msgStr)
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
		print "x ctr\ty ctr\tx err\ty err\t    pixels\tcounts\tradius"
		for cd in centroidList:
			print "%5.1f\t%5.1f\t%5.1f\t%5.1f\t%10.0f\t%6d\t%6d" % \
				(cd.xyCtr[0], cd.xyCtr[1],
				 cd.xyErr[0], cd.xyErr[1],
				 cd.pix, cd.counts, cd.rad)
	return centroidList, imStats
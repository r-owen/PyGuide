"""Measure centroids.

To do:
- Improve the estimate of i,j centroid error.
- Smooth the data before centroiding it to handle cosmic rays.
  Consider either a gaussian smoothing or a median filter.
  In either case, make sure to handle masked pixels correctly.

WARNINGS:
- Python handles images as data[y,x]. This same index order is used
for positions, sizes and such. In an attempt to reduce confusion, the code
uses i,j instead of y,x (i=y, j=x).
- Will be thrown off by hot pixels. This could perhaps
be improved by centroiding median-filtered data. The question
is whether the median filtering would adversely affect
centroids, especially for faint objects. This is especially
a concern because at present I have no code to do a proper
median filter of masked data.

Pixel convention: The point 0,0 is at the corner of first pixel read out.
Hence the center of that pixel is (0.5, 0.5) and the center of a 1024x1024 CCD
is (512.0, 512.0)

The centroid is the point of mimimum radial asymmetry:
  sum over rad of var(rad)^2 / weight(rad)
where weight is the expected sigma of var(rad) due to pixel noise:
  weight(rad) = pixNoise(rad) * sqrt(2(numPix(rad) - 1))/numPix(rad)
  pixNoise(rad) = sqrt((readNoise/ccdGain)^2 + (meanVal(rad)-bias)/ccdGain)
	
Warning: asymm is supposed to be normalized, but it gets large
for bright objects with lots of masked pixels. This may be
simply because the value is only computed at the nearest integer pixel
or because the noise is assumed gaussian, or some error.

The minimum is found in two stages:
1) Find the pixel with the minimum radAsymm.
The direction to walk is determined by measuring radAsymm at 9 points.
Each step is one pixel along x and/or y.

2) Find the true centroid (to better than one pixel) by applying
a quadratic fit to the 3x3 radAsymm matrix centered on the
pixel of minimum radAsymm. Only the points along +-x and +-y
are used for this fit; the diagonals are ignored.

Acknowledgements:
- The centroiding algorithm was invented by Jim Gunn
- The code uses a new asymmetry weighting function
  developed with help from Connie Rockosi
- This code is adapted from the SDSS centroiding code,
  which was written by Jim Gunn and cleaned up by Connie Rockosi.
  
History:
2004-03-22 ROwen	First release.
2004-04-07 ROwen	Packaged as part of PyGuide and moved test code elsewhere.
					Also changed array data types to match changes in radProf.
2004-04-12 ROwen	Modified centroid to return totCounts.
2004-04-16 ROwen	Modified centroid to not return minAsymm.
2004-04-30 ROwen	Modified to truncate initGuess (i.e. anywhere within a pixel
					selects that pixel) and round radius to the nearest integer.
					Bug fix: was converting to Int16 instead of UInt16.
2004-06-03 ROwen	Modified to use the initial guess without modification.
2004-08-03 ROwen	Finally added a measure of centroiding error.
2004-08-06 ROwen	Weight asymmetry calculation by radial noise.
2004-08-25 ROwen	Added _MinRad, to more reliably centroid small stars.
					Added __all__.
2004-10-14 ROwen	Stopped computing several unused variables. Improved import of radProf.
"""
__all__ = ['centroid']

import math
import numarray as num
import numarray.nd_image as nd_im
import radProf

# minimum radius
_MinRad = 3.0

# max # of iterations
_MaxIter = 40

# debugging flags
_CTRDEBUG = False
_CTRITERDEBUG = False

class CentroidData:
	"""Centroid data, including the following fields:
	- ctr		the i,j centroid (pixels)
	- err		the predicted i,j 1-sigma error (pixels)

	note: the following three values are computed for that radial profile
	centered on the pixel nearest the centroid (NOT the true centroid):

	- asymm		measure of asymmetry:
				  sum over rad of var(rad)^2 / weight(rad)
				where weight is the expected sigma of var(rad) due to pixel noise:
				  weight(rad) = pixNoise(rad) * sqrt(2(numPix(rad) - 1))/numPix(rad)
				  pixNoise(rad) = sqrt((readNoise/ccdGain)^2 + (meanVal(rad)-bias)/ccdGain)
	- pix		the total number of unmasked pixels (ADU)
	- counts	the total number of counts (ADU)
	
	Warning: asymm is supposed to be normalized, but it gets large
	for bright objects with lots of masked pixels. This may be
	simply because the value is only computed at the nearest integer pixel
	or because the noise is assumed gaussian, or some error.
	
	other items of possible interest
	- rad		radius used to find centroid (pixels)
	"""
	def __init__(self,
		ctr,
		err,
		asymm,
		pix,
		counts,
		rad,
	):
		self.ctr = ctr
		self.err = err
		
		self.asymm = asymm
		self.pix = pix
		self.counts = counts
		
		self.rad = rad
		

def centroid(
	data,
	mask,
	initGuess,
	rad,
	bias,
	readNoise,
	ccdGain,
):
	"""Compute a centroid.

	Inputs:
	- data		image data [i,j]
	- mask		a mask [i,j] of 0's (valid data) or 1's (invalid); None if no mask.
				If mask is specified, it must have the same shape as data.
	- initGuess	initial i,j guess for centroid
	- rad		radius of search (pixels);
				values less than _MinRad are treated as _MinRad
	- bias		ccd bias (ADU)
	- readNoise	ccd read noise (e-)
	- ccdGain	ccd inverse gain (e-/ADU)
		
	Returns a CentroidData object (which see)
	"""
	if _CTRDEBUG:
		print "centroid(data, %r)" % (initGuess,)
	
	# convert input data to UInt16 and make contiguous, if necessary, to speed radProf call
	if data.type() != num.UInt16:
		if _CTRDEBUG:
			print "centroid: converting data to UInt16"
		data = data.astype(num.UInt16)
	elif not data.iscontiguous():
		if _CTRDEBUG:
			print "centroid: copying data to make contiguous"
		data = data.copy()

	# round the initial guess and radius to the nearest integer
	if len(initGuess) != 2:
		raise ValueError("initial guess=%r must have 2 elements" % (initGuess,))
	initGuess = [int(val + 0.5) for val in initGuess]
	rad = int(max(rad, _MinRad) + 0.5)
	
	# OK, use this as first guess at maximum. Extract radial profiles in
	# a 3x3 gridlet about this, and walk to find minimum fitting error
	maxi, maxj = initGuess
	radSq = rad**2
	asymmArr = num.zeros([3,3], num.Float64)
	totPtsArr = num.zeros([3,3], num.Int32)
	totCountsArr = num.zeros([3,3], num.Float64)
	
	niter = 0
	while True:
		niter += 1
		if niter > _MaxIter:
			raise RuntimeError("could not find a star in %s iterations" % (niter,))
		
		for i in range(3):
			ii = maxi + i - 1
			for j in range(3):
				jj = maxj + j - 1
				if totPtsArr[i, j] != 0:
					continue
				asymmArr[i, j], totCountsArr[i, j], totPtsArr[i, j] = radProf.radAsymmWeighted(
					data, mask, (ii, jj), rad, bias, readNoise, ccdGain)

				if _CTRDEBUG and _CTRITERDEBUG:
					print "centroid: asymm = %10.1f, totPts = %s, totCounts = %s" % \
						(asymmArr[i, j], totPtsArr[i, j], totCountsArr[i, j])

		# have error matrix. Find minimum
		ii, jj = nd_im.minimum_position(asymmArr)
		ii -= 1
		jj -= 1

		if _CTRDEBUG:
			print "centroid: error matrix min ii=%d, jj=%d, errmin=%5.1f" % (ii, jj, asymmArr[ii,jj])
			if _CTRITERDEBUG:
				print "centroid: asymm matrix =\n", asymmArr

		if (ii != 0 or jj != 0):
			# minimum error not in center; walk and try again
			if _CTRDEBUG:
				print "shift bj", -ii, -jj
			
			# shift asymmArr and totPtsArr to minimum is in center again
			asymmArr = nd_im.shift(asymmArr, (-ii, -jj))
			totCountsArr = nd_im.shift(totCountsArr, (-ii, -jj))
			totPtsArr = nd_im.shift(totPtsArr, (-ii, -jj))

			maxi += ii
			maxj += jj
			if ((maxi - initGuess[0])**2 + (maxj - initGuess[1])**2) >= radSq:
				raise RuntimeError("could not find star within %r pixels" % (rad,))
		else:
			# Have minimum. Get out and go home.
			break

	if _CTRDEBUG:
		print "centroid: after %r iterations computing final quantities" % (niter,)
	
	# perform a parabolic fit to find true centroid
	# and compute the error estimate
	# y(x) = ymin + a(x-xmin)^2
	# a = (y0 - 2y1 + y2) / 2
	# xmin = b/2a where b = (y2-y0)/2
	# ymin = y1 - b^2/4a  but this is tricky in 2 dimensions so we punt
	# for a given delta-y, delta-x = sqrt(delta-y / a)
	ai = 0.5 * (asymmArr[2, 1] - 2.0*asymmArr[1, 1] + asymmArr[0, 1])
	bi = 0.5 * (asymmArr[2, 1] - asymmArr[0, 1])
	aj = 0.5 * (asymmArr[1, 2] - 2.0*asymmArr[1, 1] + asymmArr[1, 0])
	bj = 0.5 * (asymmArr[1, 2] - asymmArr[1, 0])

#	print "asymmArr[1-3]=", asymmArr[0,1], asymmArr[1,1], asymmArr[2,1], "ai, aj=", ai, aj
	
	di = -0.5*bi/ai
	dj = -0.5*bj/aj
	iCentroid = maxi + di
	jCentroid = maxj + dj
	
	# crude error estime, based on measured asymmetry
	# note: I also tried using the minimum along i,j but that sometimes is negative
	# and this is already so crude that it's not likely to help
	radAsymmSigma = asymmArr[1,1]
	iErr = math.sqrt(radAsymmSigma / ai)
	jErr = math.sqrt(radAsymmSigma / aj)

	return CentroidData(
		ctr = (iCentroid, jCentroid),
		err = (iErr, jErr),
		counts = totCountsArr[1,1],
		pix = totPtsArr[1,1],
		asymm = asymmArr[1,1],
		rad = rad,
	)

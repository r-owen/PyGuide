#!/usr/local/bin/python
"""StarShape

Fit a a star to a symmetrical double gaussian.

Uses an algorithm developed by Jim Gunn:
- Model the star as a double gaussian: a main gaussian
plus a small contribution from a gaussian of 1/10 the amplitude
and twice the sigma. (When we speak of sigma or fwhm we mean
the sigma or fwhm of the main gaussian.)
- Try different widths by walking along a table that spans
the space fwhm = [1, 1.5, ... 20] pixels. Note that the table
actually contains width parameter values, where wp = 1 / sigma**2.

Note: the gaussian function is:
C * e**-(x-xo)**2/(2 * sigma**2)
where:
C = 1 / (sigma * sqrt(2 pi))

The full width at half maximum is given by:
fwhm = 2 * sqrt(2 * ln(2)) * sigma ~= 2.35 * sigma

This code is based on an algorithm and code by Jim Gunn,
with refinements by Connie Rockosi.

Note: Jim's original code uses a clever trick to avoid recomputing
the model profile. It computes the model profile once at a large
number of points, then picks which of those points to use
based on the trial star width. I omitted that trick because I felt
it was easier (and fast enough) to compute a new model profile
for each trial width.

Refinements include:
- The final amplitude, background and chiSq are computed
based on the final width. The original code computed
ampl, bkgnd and chiSq at various width parameters,
then used a parbolic fit to compute a final width
and a final chiSq (but not a final ampl and bkgnd).
As a result, chiSq could be negative in extreme cases.

To do:
- Normalize the chiSq function and possibly refine it.
I tried various alternative weightings including:
nPts(rad)**2 / var(rad)
nPts(rad) / var(rad)
nPts(rad) > 1
but none did any better than nPts.

History:
2004-05-20 ROwen
2004-06-03 ROwen	Modified the module doc string.
2004-07-02 ROwen	Improved the error estimate.
					Changed GStarFit.err to chiSq
2004-08-03 ROwen	Renamed GStarFit->StarShapeData.
					Modified to use Constants.FWHMPerSigma.
2004-08-04 ROwen	Improved calculation of ij pixel index and position.
					Simplified final computation of minimum width parameter.
					If shape computation fails, converts ArithmeticError into RuntimeError
2004-08-06 ROwen	Fixed invalid variable reference when _FitRadProfIterDebug true. 
2004-12-01 ROwen	Modified StarShapeData to use NaN as the default for each argument.
					Added __all__.
2005-02-07 ROwen	Changed starShape argument ctr (i,) to xyCtr.
2005-04-01 ROwen	Added required argument rad and optional argument bkgnd.
					No longer iterate the fit with an updated predFWHM because
					it doesn't seem to help when the data radius is fixed.
					Added constant _MinRad to constrain the minimum radius.
2005-04-04 ROwen	Bug fix: mis-handled case of bkgnd not specified.
2005-04-15 ROwen	Temporarily hacked the weighting function to see if it makes things better.
					Added pylab (matplotlib) debugging graphs.
2005-04-22 ROwen	Modified to use nPts as the weighting function.
					This seems to work slightly better than nPts > 1
					and just as well as a combination of nPts and a very crude estimate of S/N.
2005-04-25 ROwen	Updated doc string to state that nPts is the weighting function.
					Removed givenBkgnd argument; it only causes trouble.
2005-04-28 ROwen	Modified to fit using Scientific python's implementation of the
					Levenberg-Marquardt method. Still to do: rewrite doc strings,
					redo pylab plot of results, remove superfluous debugging constants,
					consider using scipy instead of Scientific Python
					(especially if it gives error estimates on parameters!),
					and test this on the usual batch of test code
					(which should be modified to stress the problem the old code showed)
"""
__all__ = ["StarShapeData", "starShape"]

import math
import numarray as num
import numarray.ma
import radProf as RP
from Constants import FWHMPerSigma, NaN
import ImUtil
from Scientific.Functions.LeastSquares import leastSquaresFit
import Numeric

# minimum radius
_MinRad = 3.0

# range of FWHM that is explored
_FWHMMin = 1.0
_FWHMMax = 30.0
_FWHMDelta = 0.25

# constants that may want to be ditched
_DMax = 4096

# debugging flags
_StarShapeDebug = False
_FitRadProfDebug = False
_FitRadProfIterDebug = False
_StarShapePyLab = False

class StarShapeData:
	"""Guide star fit data
	
	Attributes:
	- ampl		profile amplitude (ADUs)
	- bkgnd		background level (ADUs)
	- fwhm		FWHM (pixels)
	- chiSq		chi squared of fit
	"""
	def __init__(self,
		ampl = NaN,
		fwhm = NaN,
		bkgnd = NaN,
		chiSq = NaN,
	):
		self.ampl = float(ampl)
		self.bkgnd = float(bkgnd)
		self.fwhm = float(fwhm)
		self.chiSq = float(chiSq)


def starShape(
	data,
	mask,
	xyCtr,
	rad,
	predFWHM = None,
):
	"""Fit a double gaussian profile to a star
	
	Inputs:
	- data		a numarray array of signed integer data
	- mask		a numarray boolean array, or None if no mask (all data valid).
				If supplied, mask must be the same shape as data
				and elements are True for masked (invalid data).
	- med		median (used as the background)
	- xyCtr		x,y center of star; use the convention specified by
				PyGuide.Constants.PosMinusIndex
	- rad		radius of data to fit (pixels);
				values less than _MinRad are treated as _MinRad
	- predFWHM	predicted FWHM; if omitted then rad is used.
				You can usually omit this because the final results are not very sensitive
				to predFWHM. However, if the predicted FWHM is much too small
				then starShape may fail or give bad results.
	"""
	if _StarShapeDebug:
		print "starShape: data[%s,%s]; xyCtr=%.2f, %.2f; rad=%.1f" % \
			(data.shape[0], data.shape[1], xyCtr[0], xyCtr[1], rad)

	# compute index of nearest pixel center (pixel whose center is nearest xyCtr)
	ijCtrInd = ImUtil.ijIndFromXYPos(xyCtr)
	
	# compute offset of position from nearest pixel center
	ijCtrFloat = ImUtil.ijPosFromXYPos(xyCtr)
	ijOff = [abs(round(pos) - pos) for pos in ijCtrFloat]
	offSq = ijOff[0]**2 + ijOff[1]**2

	# adjust radius as required
	rad = int(round(max(rad, _MinRad)))

	# compute radial profile and associated data
	radIndArrLen = rad + 2 # radial index arrays need two extra points
	radProf = num.zeros([radIndArrLen], num.Float32)
	var = num.zeros([radIndArrLen], num.Float32)
	nPts = num.zeros([radIndArrLen], num.Long)
	RP.radProf(data, mask, ijCtrInd, rad, radProf, var, nPts)
	
	if _StarShapePyLab:
		global pylab
		import pylab
		pylab.close()
		pylab.subplot(3,1,1)
		pylab.plot(radProf)
		pylab.subplot(3,1,3)
		pylab.plot(nPts)
		pylab.subplot(3,1,1)
	
	# fit data
	if predFWHM == None:
		predFWHM = float(rad)
		
	# generate data array for leastSquaresFit function
	# each point is a tuple of:
	# - radSq
	# - meas(radSq)
	# - variance of data point
	radSq = RP.radSqByRadInd(radIndArrLen)
	radSq = [float(rs) for rs in radSq]
	modVar = numarray.where(var<0.5, 9.9e99, var)
	dataArr = zip(radSq, radProf, modVar)
	if _StarShapeDebug:
		print "starShape dataArr=", dataArr
	for ii in range(radIndArrLen):
		if nPts[ii] > 0:
			predAmpl = radProf[ii]
			break
	else:
		raise ValueError("no data")
	
	for ii in range(radIndArrLen-1, 0-1, -1):
		if nPts[ii] > 0:
			predBkgnd = radProf[ii]
			break
	else:
		raise ValueError("no data")

	predBkgnd = radProf[-1]
	predWP = _wpFromFWHM(predFWHM)
	
	predAmpl = radProf[0] - predBkgnd
	initGuess = (predBkgnd, predAmpl, predWP)
	if _StarShapeDebug:
		print "starShape initGuess=", initGuess
	(bkgnd, ampl, wp), chiSq = leastSquaresFit(_sciModel, initGuess, dataArr)

	if _StarShapePyLab:
		# plot fit profile on top of radial profile
		pylab.subplot(3,1,1)
		# plot fit profile on top of radial profile
		pylab.subplot(3,1,2)
		# plot error graph here
		pylab.subplot(3,1,1)

	gsData = StarShapeData(
		ampl = ampl,
		bkgnd = bkgnd,
		fwhm = _fwhmFromWP(wp),
		chiSq = chiSq
	)
	return gsData

def _fwhmFromWP(wp):
	"""Converts width parameter to fwhm in pixels.
	wp is the width parameter: 1/sig**2
	"""
	return FWHMPerSigma / math.sqrt(wp)


def _wpFromFWHM(fwhm):
	"""Converts fwhm in pixels to width parameter in 1/pix^2 (???).
	wp is the width parameter 1/sig^2.
	"""
	return (FWHMPerSigma / fwhm)**2


def _sciModel(params, radSq):
	"""Computes the predicted star profile for the set of parameters.
	
	Inputs:
	- params	background, amplitude and width parameter
	- radSq		radius squared
	"""
	bkgnd, ampl, wp = params
	
	return bkgnd + (ampl * (Numeric.exp(-0.5 * wp * radSq) + 0.1 * Numeric.exp(-0.125 * wp * radSq)) / 1.1)

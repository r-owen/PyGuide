#!/usr/local/bin/python
"""Time PyGuide.centroid and related routines.

History:
2004-04-07 ROwen	First release
2004-04-12 ROwen	Modified to mask of all 0s (sense of mask changed).
2004-08-03 ROwen	Modified to use fake data.
2004-08-06 ROwen	Modified for new centroid.
"""
import time
import numarray as num
import PyGuide

# settings for fake data
ImWidth = 1024	# image width
Sky = 1000		# sky level, in ADU
ReadNoise = 19	# read noise, in e-
CCDGain = 2.1	# inverse ccd gain, in e-/ADU
Bias = 2176		# image bias, in ADU

def timeCentroid(data, mask, initGuess, niter, rad=20):
	print "timeCentroid: initGuess=%3.0f, %3.0f; niter=%2d; rad=%3d;" % \
		(initGuess[0], initGuess[1], niter, rad),
	begTime = time.time()
	for ii in range(niter):
		PyGuide.centroid(
			data = data,
			mask = mask,
			initGuess = initGuess,
			rad = rad,
			bias = Bias,
			readNoise = ReadNoise,
			ccdGain = CCDGain,
		)
	dTime = time.time() - begTime
	print "time/iter=%.3f" % (dTime/niter,)


def timeRadAsymmWeighted(data, mask, niter, rad=20):
	shape = data.getshape()
	xc = shape[0]/2
	yc = shape[1]/2
	print "timeRadAsymmWeighted: niter=%2d; rad=%3d;" % (niter, rad),
	
	begTime = time.time()
	for ii in range(niter):
		PyGuide.radProf.radAsymmWeighted(
			data, mask, (xc, yc), rad,
			Bias, ReadNoise, CCDGain,
		)
	dTime = time.time() - begTime
	print "time/iter=%.3f" % (dTime/niter,)


def timeRadProf(data, mask, niter, rad):
	"""Time radProf and radSqProf"""
	shape = data.getshape()
	xc = shape[0]/2
	yc = shape[1]/2
	print "timeRadProf: niter=%2d; rad=%3d;" % (niter, rad,),
	
	radSqLen = rad**2 + 1
	radSqMean = num.zeros([radSqLen], num.Float64)
	radSqVar = num.zeros([radSqLen], num.Float64)
	radSqNPts = num.zeros([radSqLen], num.Int32)

	radLen = rad + 2
	radMean = num.zeros([radLen], num.Float64)
	radVar = num.zeros([radLen], num.Float64)
	radNPts = num.zeros([radLen], num.Int32)
	
	begTime = time.time()
	for ii in range(niter):
		PyGuide.radProf.radSqProf(data, mask, (xc, yc), rad, radSqMean, radSqVar, radSqNPts)
	dTime = time.time() - begTime
	print "radSqProf time/iter=%.4f;" % (dTime/niter,),
	
	begTime = time.time()
	for ii in range(niter):
		PyGuide.radProf.radProf(data, mask, (xc, yc), rad, radMean, radVar, radNPts)
	dTime = time.time() - begTime
	print "radProf time/iter=%.4f" % (dTime/niter,),
	
	print


def runTests():
	# settings for fake data
	imWidth = 1024
	fwhm = 2.5
	ampl = 1000
	
	# generate fake data
	imShape = (ImWidth, ImWidth)
	ctr = num.divide(imShape, 2.0)
	sigma = fwhm / PyGuide.FWHMPerSigma
	cleanData = PyGuide.FakeData.fakeStar(imShape, ctr, sigma, ampl)
	data = PyGuide.FakeData.addNoise(
		data = cleanData,
		sky = Sky,
		bias = Bias,
		readNoise = ReadNoise,
		ccdGain = CCDGain,
	)
	data = data.astype(num.Int16)
	
	# let centroiding walk a bit to find the center
	initGuess = num.add(ctr, (2, -2))
	
	print "Time various parts of centroiding as a function of radius"
	print
	print "Settings:"
	print "CCD Size      =", ImWidth, "x", ImWidth, "pix"
	print "Star center   = %d, %d pix" % (ctr[0], ctr[1])
	print "Initial guess = %d, %d pix" % (initGuess[0], initGuess[1])
	print
	print "Amplitude  =", ampl, "ADU"
	print "FWHM       =", fwhm, "pix"
	print "Sky        =", Sky, "ADU"
	print "Bias       =", Bias, "ADU"
	print "Read Noise =", ReadNoise, "e-"
	print "CCD Gain   =", CCDGain, "e-/ADU"
	
	allZerosMask = data.astype(num.Bool)
	allZerosMask[:] = 0
	
	maskData = (
		(None, "No mask"),
		(allZerosMask, "Mask of all zeros"),
	)
	
	for mask, expl in maskData:
		print
		print expl
		
		radNiterList = (
			( 10, 20),
			( 20, 20),
			( 40, 20),
			( 80, 20),
			(160, 20),
			(320, 20),
			(640, 10),
		)
	
		print
		for rad, niter in radNiterList:
			try:
				timeRadProf(data, mask, niter, rad)
			except Exception, e:
				print "timeRadProf(niter=%s, rad=%s) failed: %s" % (niter, rad, e)
		
		print
		for rad, niter in radNiterList:
			try:
				timeRadAsymmWeighted(data, mask, niter, rad)
			except Exception, e:
				print "timeRadAsymmWeighted(niter=%s, rad=%s) failed: %s" % (niter, rad, e)

		radNiterList = (
			( 10, 10),
			( 20, 10),
			( 40, 10),
			( 80, 10),
			(160, 5),
			(320, 2),
			(640, 1),
		)
		
		print
		for rad, niter in radNiterList:
			try:
				timeCentroid(data, mask, initGuess, niter, rad)
			except Exception, e:
				print "timeCentroid(niter=%s, rad=%s) failed: %s" % (niter, rad, e)


if __name__ == "__main__":
	runTests()

#!/usr/local/bin/python
"""Time PyGuide.centroid and related routines.

History:
2004-04-07 ROwen	First release
2004-04-12 ROwen	Modified to mask of all 0s (sense of mask changed).
2004-08-03 ROwen	Modified to use fake data.
"""
import time
import numarray as num
import PyGuide

def timeCentroid(data, mask, initGuess, niter, rad=20):
	print "timeCentroid: initGuess=%3.0f, %3.0f; niter=%2d; rad=%3d;" % \
		(initGuess[0], initGuess[1], niter, rad),
	begTime = time.time()
	for ii in range(niter):
		PyGuide.centroid(data, mask, initGuess, rad)
	dTime = time.time() - begTime
	print "time/iter=%.3f" % (dTime/niter,)


def timeRadAsymm(data, mask, niter, rad=20):
	shape = data.getshape()
	xc = shape[0]/2
	yc = shape[1]/2
	print "timeRadAsymm: niter=%2d; rad=%3d;" % (niter, rad),
	
	begTime = time.time()
	for ii in range(niter):
		PyGuide.radProf.radAsymm(data, mask, (xc, yc), rad)
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
	sky = 1000		# sky level, in ADU
	readNoise = 19	# read noise, in e-
	ccdGain = 2.1	# inverse ccd gain, in e-/ADU
	bias = 2176		# image bias, in ADU
	
	# generate fake data
	imShape = (imWidth, imWidth)
	ctr = num.divide(imShape, 2.0)
	sigma = fwhm / PyGuide.FWHMPerSigma
	cleanData = PyGuide.FakeData.fakeStar(imShape, ctr, sigma, ampl)
	data = PyGuide.FakeData.addNoise(
		data = cleanData,
		sky = sky,
		readNoise = readNoise,
		ccdGain = ccdGain,
		bias = bias,
	)
	data = data.astype(num.Int16)
	
	# let centroiding walk a bit to find the center
	initGuess = num.add(ctr, (2, -2))
	
	print "Time various parts of centroiding as a function of radius"
	print
	print "Settings:"
	print "CCD Size      =", imWidth, "x", imWidth, "pix"
	print "Star center   = %d, %d pix" % (ctr[0], ctr[1])
	print "Initial guess = %d, %d pix" % (initGuess[0], initGuess[1])
	print
	print "Amplitude  =", ampl, "ADU"
	print "FWHM       =", fwhm, "pix"
	print "Sky        =", sky, "ADU"
	print "Read Noise =", readNoise, "e-"
	print "CCD Gain   =", ccdGain, "e-/ADU"
	print "Bias       =", bias, "ADU"
	
	allZerosMask = data.astype(num.Bool)
	allZerosMask[:] = 0
	
	maskData = (
		(None, "No mask"),
		(allZerosMask, "Mask of all zeros"),
	)
	
	for mask, expl in maskData:
		print
		print expl
	
		print
		timeRadProf(data, mask, 20, rad=10)
		timeRadProf(data, mask, 20, rad=20)
		timeRadProf(data, mask, 20, rad=40)
		timeRadProf(data, mask, 20, rad=80)
		timeRadProf(data, mask, 20, rad=160)
		timeRadProf(data, mask, 20, rad=320)
		timeRadProf(data, mask, 10, rad=640)
		
		print
		timeRadAsymm(data, mask, 20, rad=10)
		timeRadAsymm(data, mask, 20, rad=20)
		timeRadAsymm(data, mask, 20, rad=40)
		timeRadAsymm(data, mask, 20, rad=80)
		timeRadAsymm(data, mask, 20, rad=160)
		timeRadAsymm(data, mask, 20, rad=320)
		timeRadAsymm(data, mask, 10, rad=640)
		
		print
		timeCentroid(data, mask, initGuess, 10, rad=10)
		timeCentroid(data, mask, initGuess, 10, rad=20)
		timeCentroid(data, mask, initGuess, 10, rad=40)
		timeCentroid(data, mask, initGuess, 10, rad=80)
		timeCentroid(data, mask, initGuess, 5, rad=160)
		timeCentroid(data, mask, initGuess, 2, rad=320)


if __name__ == "__main__":
	runTests()

#!/usr/local/bin/python
"""Test PyGuide.centroid with lots of fake data.

Limitations:
- Only one star per image
- No cosmic rays

History:
2004-08-04 ROwen
2004-08-06 ROwen	Modified for new centroid.
2005-02-07 ROwen	Modified for PyGuide 1.2.
"""
import numarray as num
import numarray.random_array as num_random
import PyGuide

# image data info
ImWidth = 64
Sky = 1000		# sky level, in ADU
ReadNoise = 19	# read noise, in e-
CCDGain = 2.1	# inverse ccd gain, in e-/ADU
Bias = 2176		# image bias, in ADU

imShape = (ImWidth, ImWidth)
nomCtr = (ImWidth // 2, ImWidth // 2)
mask = num.zeros(imShape, num.Bool)
NaN = float("NaN")

# settings
AmplValues = (100, 1000, 10000)
FWHMValues = (2.0, 3.0, 4.0)
MaskWidthsPerFWHM = (0.0, 0.5, 1.0, 1.5, 2.0) # fractions of a FWHM
NumTries = 20

print "Compare centroid actual error vs estimated error"
print "over a range of fake data"
print
print "Settings:"
print "Sky         =", Sky, "ADU"
print "Read Noise  =", ReadNoise, "e-"
print "CCD Gain    =", CCDGain, "e-/ADU"
print "Bias        =", Bias, "ADU"
print "Amplitudes  =", AmplValues, "ADU"
print "FWHMs       =", FWHMValues, "pixels"
print "Mask Widths =", MaskWidthsPerFWHM, "fractions of a FWHM"
print "Num Ctrs    =", NumTries, "number of centroids per star shape"
print "CCD Size    =", ImWidth, "x", ImWidth, "pixels"
print
print "the equation being minimized is:"
print "   sum(var(r) * numPts(r))"
print "though the final number reported is a normalized version:"
print "   asymm = sum(var(r) * numPts(r)) / (pixNoise^2 * totPts)"
print "where"
print "   pixNoise is the noise/pixel due to read noise and sky"
print "   totPts = sum(numPts(r))"
print
print "The centroids are randomly distributed over"
print "a range of -FWHM/2, FWHM/2 with respect to the"
print "center of the CCD = the center of the slit."
print
print "The slit is along y so the error along y should be smaller than x."
print
print "fwhm	ampl	maskWid	xErr	yErr	xUncert	yUncert	asymm	totPix	totCts	rad"

nBad = 0
for ampl in AmplValues:
	for fwhm in FWHMValues:
		sigma = fwhm / PyGuide.FWHMPerSigma
		for maskMult in MaskWidthsPerFWHM:
			maskWidth = maskMult * fwhm
			maskRad = int(maskWidth / 2.0)
			mask[:,:] = 0
			if maskRad > 0:
				mask[nomCtr[0] - maskRad: nomCtr[0] + maskRad + 1, :] = 1

			num_random.seed(1, 1000)
			for ii in range(NumTries):
				actCtr = num_random.uniform(-fwhm/2.0, fwhm/2.0, shape=(2,)) + nomCtr

				cleanData = PyGuide.FakeData.fakeStar(imShape, actCtr, sigma, ampl)
				data = PyGuide.FakeData.addNoise(
					data = cleanData,
					sky = Sky,
					bias = Bias,
					readNoise = ReadNoise,
					ccdGain = CCDGain,
				)

				try:
					ctrData = PyGuide.centroid(
						data = data,
						mask = mask,
						xyGuess = nomCtr,
						rad = fwhm * 3.0,
						bias = Bias,
						readNoise = ReadNoise,
						ccdGain = CCDGain,
					)
					xyMeasErr = [ctrData.xyCtr[ii] - actCtr[ii] for ii in (0,1)]
					print "%s	%s	%s	%.3f	%.3f	%.3f	%.3f	%.3f	%s	%s	%s" % (
						fwhm, ampl, maskWidth,
						xyMeasErr[0], xyMeasErr[1],
						ctrData.xyErr[0], ctrData.xyErr[1],
						ctrData.asymm, ctrData.pix, ctrData.counts, ctrData.rad,
					)
					
				except RuntimeError, e:
					nBad += 1

if nBad > 0:
	print "number of failures =", nBad

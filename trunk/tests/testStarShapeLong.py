#!/usr/local/bin/python
"""Test PyGuide.starShape with lots of fake data.

Limitations:
- Only one star per image
- No cosmic rays

History:
2004-08-04 ROwen
2005-02-07 ROwen	Modified for PyGuide 1.2.
2005-04-01 ROwen	Modified for new StarShape; modified to print statistics.
2005-04-04 ROwen	Modified for fixed StarShape.
2005-04-05 ROwen	Modified to show failed cases (with NaN for the shape data).
2005-04-06 ROwen	Removed unnecessary float("NaN").
2005-04-25 ROwen	Modified for new starShape (now always fit background).
"""
import numarray as num
import numarray.random_array as num_random
import PyGuide
from Stats import Stats

# image data info
ImWidth = 512
Sky = 1000		# sky level, in ADU
ReadNoise = 19	# read noise, in e-
CCDGain = 2.1	# inverse ccd gain, in e-/ADU
Bias = 2176		# image bias, in ADU

imShape = (ImWidth, ImWidth)
nomCtr = (ImWidth // 2, ImWidth // 2)
mask = num.zeros(imShape, num.Bool)

# settings
AmplValues = (100, 1000, 10000)
FWHMValues = (2.0, 3.0, 4.0)
MaskWidthsPerFWHM = (0.0, 0.5, 1.0, 1.5, 2.0) # fractions of a FWHM
NumTries = 10

print "Compare star shape fit values to correct values"
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
print "Num Tries   =", NumTries, "number of cases per star shape"
print "CCD Size    =", ImWidth, "x", ImWidth, "pixels"
print
print "Each try is a star whose center is randomly distributed over"
print "a range of -FWHM/2, FWHM/2 with respect to the"
print "center of the CCD = the center of the slit."
print
print "Reported errors and statistics on these errors are in percent:"
print "reported error (%) = 100 * (meas value - act value) / act value"
print
print "The slit is along y."

def pctErr(meas, act):
	return (meas - act) * 100.0 / float(act)

fwhmStats = Stats()
amplStats = Stats()
bkgndStats = Stats()
nBad = 0

print
print "fwhm	ampl	bg	xCtr	yCtr	maskWid	fitFWHM	fitAmpl	fitBg	chiSq	fwhmErr	amplErr	bgErr"
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
				xyCtr = num_random.uniform(-fwhm/2.0, fwhm/2.0, shape=(2,)) + nomCtr

				cleanData = PyGuide.FakeData.fakeStar(imShape, xyCtr, sigma, ampl)
				data = PyGuide.FakeData.addNoise(
					data = cleanData,
					sky = Sky,
					readNoise = ReadNoise,
					ccdGain = CCDGain,
					bias = Bias,
				)
				
				maskedData = num.ma.array(data, mask=mask)

				try:
					shapeData = PyGuide.starShape(
						data = data,
						mask = mask,
						xyCtr = xyCtr,
						rad = fwhm * 2,
					)
					bkgnd = Sky + Bias
					fwhmErr = pctErr(shapeData.fwhm, fwhm)
					amplErr = pctErr(shapeData.ampl, ampl)
					bkgndErr = pctErr(shapeData.bkgnd, bkgnd)
					print "%.1f	%.1f	%.1f	%.2f	%.2f	%.2f	%.1f	%.1f	%.1f	%.2f	%.1f	%.1f	%.1f" % (
						fwhm, ampl, bkgnd,
						xyCtr[0], xyCtr[1], maskWidth,
						shapeData.fwhm, shapeData.ampl, shapeData.bkgnd, shapeData.chiSq,
						fwhmErr, amplErr, bkgndErr,
					)
					fwhmStats.append(fwhmErr)
					amplStats.append(amplErr)
					bkgndStats.append(bkgndErr)
					
				except RuntimeError, e:
					print "%.1f	%.1f	%.1f	%.2f	%.2f	%.2f	NaN	NaN	NaN	NaN	NaN	NaN	NaN" % (
						fwhm, ampl, bkgnd,
						xyCtr[0], xyCtr[1], maskWidth,
					)
					nBad += 1

print
print "Error statistics (for %d points)" % fwhmStats.nPoints()
print "            min      max    mean   stdDev"
print "fwhm  %8.1f %8.1f %8.1f %8.1f" % (fwhmStats.min(), fwhmStats.max(), fwhmStats.mean(), fwhmStats.stdDev())
print "ampl  %8.1f %8.1f %8.1f %8.1f" % (amplStats.min(), amplStats.max(), amplStats.mean(), amplStats.stdDev())
print "bkgnd %8.1f %8.1f %8.1f %8.1f" % (bkgndStats.min(), bkgndStats.max(), bkgndStats.mean(), bkgndStats.stdDev())

if nBad > 0:
	print
	print "number of failures =", nBad

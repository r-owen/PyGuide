#!/usr/bin/env python
"""Test PyGuide.centroid with lots of fake data.

Limitations:
- Only one star per image
- No cosmic rays

History:
2004-08-04 ROwen
2004-08-06 ROwen    Modified for new centroid.
2005-02-07 ROwen    Modified for PyGuide 1.2.
2005-05-19 ROwen    Modified for PyGuide 2.0.
                    Reports centroid error stats.
                    Shows inputs when centroid fails.
                    Shows centroid warnings.
2005-10-14 ROwen    Supply null satMask for PyGuide 2.1.
"""
import numarray as num
import numarray.random_array as num_random
import PyGuide
from Stats import Stats

# settings
ImWidth = 64
Sky = 1000,     # sky level, in ADU
CCDInfo = PyGuide.CCDInfo(
    bias = 2176,    # image bias, in ADU
    readNoise = 19, # read noise, in e-
    ccdGain = 2.1,  # inverse ccd gain, in e-/ADU
)
Thresh = 2.5

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
print "Thresh      =", Thresh
print "Sky         =", Sky, "ADU"
print "Read Noise  =", CCDInfo.readNoise, "e-"
print "CCD Gain    =", CCDInfo.ccdGain, "e-/ADU"
print "Bias        =", CCDInfo.bias, "ADU"
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
print "fwhm ampl    maskWid xErr    yErr    xUncert yUncert asymm   totPix  totCts  rad msgs"

nBad = 0
ctrXStats = Stats()
ctrYStats = Stats()
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
                    ccdInfo = CCDInfo,
                )

                ctrData = PyGuide.centroid(
                    data = data,
                    mask = mask,
                    satMask = None,
                    xyGuess = nomCtr,
                    rad = fwhm * 3.0,
                    ccdInfo = CCDInfo,
                    thresh = Thresh,
                )
                if not ctrData.isOK:
                    print "%s   %s  %s  NaN NaN NaN NaN NaN NaN NaN %s  %r" % (
                        fwhm, ampl, maskWidth, ctrData.rad, ctrData.msgStr,
                    )
                    nBad += 1
                    continue
                
                xyMeasErr = [ctrData.xyCtr[ii] - actCtr[ii] for ii in (0,1)]
                print "%s   %s  %s  %.3f    %.3f    %.3f    %.3f    %.3f    %s  %s  %s  %r" % (
                    fwhm, ampl, maskWidth,
                    xyMeasErr[0], xyMeasErr[1],
                    ctrData.xyErr[0], ctrData.xyErr[1],
                    ctrData.asymm, ctrData.pix, ctrData.counts,
                    ctrData.rad, ctrData.msgStr,
                )
                ctrXStats.append(xyMeasErr[0])
                ctrYStats.append(xyMeasErr[1])

print
print "Error statistics (for %d points)" % ctrXStats.nPoints()
print "            min      max    mean   stdDev"
print "xErr  %8.1f %8.1f %8.1f %8.1f" % (ctrXStats.min(), ctrXStats.max(), ctrXStats.mean(), ctrXStats.stdDev())
print "yErr  %8.1f %8.1f %8.1f %8.1f" % (ctrYStats.min(), ctrYStats.max(), ctrYStats.mean(), ctrYStats.stdDev())

if nBad > 0:
    print
    print "number of failures =", nBad

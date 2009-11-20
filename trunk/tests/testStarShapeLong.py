#!/usr/bin/env python
"""Test PyGuide.starShape with lots of fake data.

Limitations:
- Only one star per image
- No cosmic rays

History:
2004-08-04 ROwen
2005-02-07 ROwen    Modified for PyGuide 1.2.
2005-04-01 ROwen    Modified for new StarShape; modified to print statistics.
2005-04-04 ROwen    Modified for fixed StarShape.
2005-04-05 ROwen    Modified to show failed cases (with NaN for the shape data).
2005-04-06 ROwen    Removed unnecessary float("NaN").
2005-04-25 ROwen    Modified for new starShape (now always fit background).
2005-04-19 ROwen    Modified for PyGuide 2.0.
                    Modified to optionally centroid before shape fitting.
                    This should tell us more about how well the shape fitter
                    does on realistic data.
                    Shows centroid and star shape warning messages.
2005-10-14 ROwen    Supply null satMask for PyGuide 2.1.
"""
import sys
import traceback
import numarray as num
import numpy.random as num_random
import PyGuide
from Stats import Stats

# image data info
ImWidth = 512
Sky = 1000      # sky level, in ADU
CCDInfo = PyGuide.CCDInfo(
    bias = 2176,    # image bias, in ADU
    readNoise = 19, # read noise, in e-
    ccdGain = 2.1,  # inverse ccd gain, in e-/ADU
)

# other settings
DoCentroid = True
Thresh = 2.0 # allow marginal cases through to exercise starShape harder
            # though as of 2005-05-18 a threshold of 2.0 does not give significantly
            # worse star fits than a threshold of 2.5.
AmplValues = (100, 1000, 10000)
FWHMValues = (2.0, 3.0, 4.0)
MaskWidthsPerFWHM = (0.0, 0.5, 1.0, 1.5, 2.0) # fractions of a FWHM
NumTries = 10

imShape = (ImWidth, ImWidth)
nomCtr = (ImWidth // 2, ImWidth // 2)
mask = numpy.zeros(imShape, numpy.bool)

print "Compare star shape fit values to correct values"
print "over a range of fake data"
print
print "Settings:"
print "Sky         =", Sky, "ADU"
print "Read Noise  =", CCDInfo.readNoise, "e-"
print "CCD Gain    =", CCDInfo.ccdGain, "e-/ADU"
print "Bias        =", CCDInfo.bias, "ADU"
print "DoCentroid  =", DoCentroid
if DoCentroid:
    print "Thresh      =", Thresh
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
nBadCtr = 0
nBad = 0

print
print "fwhm ampl    bg  xCtr    yCtr    maskWid xCtMeas yCtMeas fitFWHM fitAmpl fitBg   chiSq   fwhmErr amplErr bgErr   msgs"
bkgnd = Sky + CCDInfo.bias
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
                    ccdInfo = CCDInfo,
                )
                
                maskedData = numpy.ma.masked_array(data, mask=mask)
                
                if DoCentroid:
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
                        print "%.1f %.1f    %.1f    %.2f    %.2f    %.2f    NaN NaN NaN NaN NaN NaN NaN NaN NaN %r" % (
                            fwhm, ampl, bkgnd,
                            xyCtr[0], xyCtr[1], maskWidth,
                            ctrData.msgStr,
                        )
                        nBadCtr += 1
                        continue
                else:
                    ctrData = PyGuide.CentroidData(
                        isOK = True,
                        xyCtr = xyCtr,
                    )

                shapeData = PyGuide.starShape(
                    data = data,
                    mask = mask,
                    xyCtr = ctrData.xyCtr,
                    rad = fwhm * 3,
                )
                if not shapeData.isOK:
                    print "%.1f %.1f    %.1f    %.2f    %.2f    %.2f    %.2f    %.2f    NaN NaN NaN NaN NaN NaN NaN %r" % (
                        fwhm, ampl, bkgnd,
                        xyCtr[0], xyCtr[1], maskWidth, ctrData.xyCtr[0], ctrData.xyCtr[1],
                        shapeData.msgStr,
                    )
                    nBad += 1
                    continue
                    
                fwhmErr = pctErr(shapeData.fwhm, fwhm)
                amplErr = pctErr(shapeData.ampl, ampl)
                bkgndErr = pctErr(shapeData.bkgnd, bkgnd)
                print "%.1f %.1f    %.1f    %.2f    %.2f    %.2f    %.2f    %.2f    %.1f    %.1f    %.1f    %.2f    %.1f    %.1f    %.1f    %r" % (
                    fwhm, ampl, bkgnd,
                    xyCtr[0], xyCtr[1], maskWidth, ctrData.xyCtr[0], ctrData.xyCtr[1],
                    shapeData.fwhm, shapeData.ampl, shapeData.bkgnd, shapeData.chiSq,
                    fwhmErr, amplErr, bkgndErr,
                    shapeData.msgStr,
                )
                fwhmStats.append(fwhmErr)
                amplStats.append(amplErr)
                bkgndStats.append(bkgndErr)

print
print "Error statistics (for %d points)" % fwhmStats.nPoints()
print "            min      max    mean   stdDev"
print "fwhm  %8.1f %8.1f %8.1f %8.1f" % (fwhmStats.min(), fwhmStats.max(), fwhmStats.mean(), fwhmStats.stdDev())
print "ampl  %8.1f %8.1f %8.1f %8.1f" % (amplStats.min(), amplStats.max(), amplStats.mean(), amplStats.stdDev())
print "bkgnd %8.1f %8.1f %8.1f %8.1f" % (bkgndStats.min(), bkgndStats.max(), bkgndStats.mean(), bkgndStats.stdDev())

if (nBad > 0) or (nBadCtr > 0):
    print
    print "number of shape fit failures =", nBad
    if DoCentroid:
        print "number of centroid failures  =", nBadCtr

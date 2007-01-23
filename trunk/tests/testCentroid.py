#!/usr/local/bin/python
"""Test PyGuide.centroid with masked data.

History:
2004-04-12 ROwen    First version with history. Modified for centroid 2004-04-12.
2004-04-16 ROwen    Modified for centroid 2004-04-16.
2004-07-07 ROwen    Added noise to the simulation. Added ds9 display.
2004-08-03 ROwen    Modified for centroid 2004-08-03.
2004-08-06 ROwen    Modified for new centroid.
                    Explicitly specifies values for noise components.
2005-02-07 ROwen    Modified for PyGuide 1.2.
                    Modified to show estimated error.
2005-05-19 ROwen    Updated for PyGuide 2.0.
2005-10-14 ROwen    Supply null satMask for PyGuide 2.1.
"""
import sys
import numarray as num
import numarray.random_array as num_random
import PyGuide
import RO.DS9

Sky = 1000,     # sky level, in ADU
CCDInfo = PyGuide.CCDInfo(
    bias = 2176,    # image bias, in ADU
    readNoise = 19, # read noise, in e-
    ccdGain = 2.1,  # inverse ccd gain, in e-/ADU
)

# test data format:
# arrShape, actual center, sigma, ampl, scanRad factor, maskLim
testData = (
    [(51, 51), (27.0, 28.6), 2, 1000, 3, (23,29)],
    [(51, 51), (27.1, 28.6), 2, 1000, 3, (23,29)],
    [(51, 51), (27.1, 28.6), 2, 1000, 4, (23,29)],
    [(51, 51), (27.2, 28.7), 2, 1000, 3, (23,29)],
    [(51, 51), (27.3, 28.8), 2, 1000, 3, (23,29)],
    [(51, 51), (27.4, 28.9), 2, 1000, 3, (23,29)],
    [(51, 51), (27.5, 28.9), 2, 1000, 3, (23,29)],
    [(51, 51), (27.6, 28.9), 2, 1000, 3, (23,29)],
    [(51, 51), (27.7, 28.9), 2, 1000, 3, (23,29)],
    [(51, 51), (27.8, 28.9), 2, 1000, 3, (23,29)],
    [(51, 51), (27.9, 28.9), 2, 1000, 3, (23,29)],
    [(51, 51), (28.0, 28.9), 2, 1000, 3, (23,29)],
)

#Centroid._CTRDEBUG = True
#Centroid._CTRITERDEBUG = True

ds9Win = RO.DS9.DS9Win("testCentroid")

for arrShape, actCtr, sigma, ampl, scanRadFactor, maskLim in testData:
    xyGuess = (
        actCtr[0] - 1.0,
        actCtr[1] + 1.0,
    )
    scanRad = scanRadFactor * sigma

    cleanData = PyGuide.FakeData.fakeStar(arrShape, actCtr, sigma, ampl)
    num_random.seed(1, 1000)
    data = PyGuide.FakeData.addNoise(
        cleanData,
        sky = Sky,
        ccdInfo = CCDInfo,
    )

    print "\nactual center = %6.2f, %6.2f, sigma = %.2f, scanRad = %d" % (actCtr[0], actCtr[1], sigma, scanRad)
    mask = None
    ctrData = PyGuide.centroid(
        data = data,
        mask = mask,
        satMask = None,
        xyGuess = xyGuess,
        rad = scanRad,
        ccdInfo = CCDInfo,
    )
    
    if not ctrData.isOK:
        print "centroid failed: %s" % (ctrData.msgStr,)
        continue
    
    measCtr = ctrData.xyCtr
    nCounts = ctrData.counts
    nPts = ctrData.pix
    print "meas err   = %6.2f, %6.2f; est err = %.2f, %.2f; nCounts = %.0f; nPts = %d" % \
        (measCtr[0] - actCtr[0], measCtr[1] - actCtr[1], ctrData.xyErr[0], ctrData.xyErr[1], nCounts, nPts)
    mask = num.zeros(arrShape, num.Bool)
    for row in range(maskLim[0], maskLim[1]+1):
        mask[row,:] = 1
    ctrData = PyGuide.centroid(
        data = data,
        mask = mask,
        satMask = None,
        xyGuess = xyGuess,
        rad = scanRad,
        ccdInfo = CCDInfo,
    )
    measCtr = ctrData.xyCtr
    nCounts = ctrData.counts
    nPts = ctrData.pix
    print "masked err = %6.2f, %6.2f; est err = %.2f, %.2f; nCounts = %.0f; nPts = %d" % \
        (measCtr[0] - actCtr[0], measCtr[1] - actCtr[1], ctrData.xyErr[0], ctrData.xyErr[1], nCounts, nPts)


ds9Win.xpaset("tile frames")
ds9Win.xpaset("frame 1")
ds9Win.showArray(data - cleanData)

ds9Win.xpaset("frame 2")
ds9Win.showArray(data)
if mask != None:
    ds9Win.xpaset("frame 3")
    ds9Win.showArray(data * (1-mask))

#!/usr/bin/env python
"""Time PyGuide.centroid and related routines.

History:
2004-04-07 ROwen    First release
2004-04-12 ROwen    Modified to mask of all 0s (sense of mask changed).
2004-08-03 ROwen    Modified to use fake data.
2004-08-06 ROwen    Modified for new centroid.
2005-02-07 ROwen    Modified for PyGuide 1.2.
2005-05-19 ROwen    Modified for PyGuide 2.0.
2005-10-14 ROwen    Supply null satMask for PyGuide 2.1.
"""
import time
import numarray as num
import PyGuide

# settings
ImWidth = 1024  # image width
Sky = 1000      # sky level, in ADU
CCDInfo = PyGuide.CCDInfo(
    bias = 2176,    # image bias, in ADU
    readNoise = 19, # read noise, in e-
    ccdGain = 2.1,  # inverse ccd gain, in e-/ADU
)
Thresh = 3.0
FWHM = 2.5
Ampl = 5000

def timeCentroid(data, mask, xyGuess, niter, rad=20):
    print "timeCentroid: xyGuess=%3.0f, %3.0f; niter=%2d; rad=%3d;" % \
        (xyGuess[0], xyGuess[1], niter, rad),
    begTime = time.time()
    for ii in range(niter):
        ctrData = PyGuide.centroid(
            data = data,
            mask = mask,
            satMask = None,
            xyGuess = xyGuess,
            rad = rad,
            thresh = Thresh,
            ccdInfo = CCDInfo,
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
            CCDInfo.bias, CCDInfo.readNoise, CCDInfo.ccdGain,
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
    radSqMean = numpy.zeros([radSqLen], numpy.Float64)
    radSqVar = numpy.zeros([radSqLen], numpy.Float64)
    radSqNPts = numpy.zeros([radSqLen], numpy.Int32)

    radLen = rad + 2
    radMean = numpy.zeros([radLen], numpy.Float64)
    radVar = numpy.zeros([radLen], numpy.Float64)
    radNPts = numpy.zeros([radLen], numpy.Int32)
    
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
    # generate fake data
    imShape = (ImWidth, ImWidth)
    xyCtr = numpy.divide(imShape, 2.0)
    sigma = FWHM / PyGuide.FWHMPerSigma
    cleanData = PyGuide.FakeData.fakeStar(imShape, xyCtr, sigma, Ampl)
    data = PyGuide.FakeData.addNoise(
        data = cleanData,
        sky = Sky,
        ccdInfo = CCDInfo,
    )
    data = data.astype(numpy.Int16)
    
    # let centroiding walk a bit to find the center
    xyGuess = numpy.add(xyCtr, (2, -2))

    print "Time various parts of centroiding as a function of radius"
    print
    print "Settings:"
    print "CCD Size      =", ImWidth, "x", ImWidth, "pix"
    print "Star center   = %d, %d pix" % (xyCtr[0], xyCtr[1])
    print "Initial guess = %d, %d pix" % (xyGuess[0], xyGuess[1])
    print
    print "Amplitude  =", Ampl, "ADU"
    print "FWHM       =", FWHM, "pix"
    print "Sky        =", Sky, "ADU"
    print "Bias       =", CCDInfo.bias, "ADU"
    print "Read Noise =", CCDInfo.readNoise, "e-"
    print "CCD Gain   =", CCDInfo.ccdGain, "e-/ADU"

    allZerosMask = data.astype(numpy.bool)
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
                timeCentroid(data, mask, xyGuess, niter, rad)
            except Exception, e:
                raise
#               print "timeCentroid(niter=%s, rad=%s) failed: %s" % (niter, rad, e)


if __name__ == "__main__":
    runTests()

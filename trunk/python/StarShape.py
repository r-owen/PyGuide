#!/usr/local/bin/python
"""StarShape

Fit a a star to a symmetrical double gaussian.

Uses an algorithm developed by Jim Gunn:
- Model the star as a double gaussian: a main gaussian
plus a small contribution from a gaussian of 1/10 the amplitude
and twice the sigma. (When we speak of sigma or fwhm we mean
the sigma or fwhm of the main gaussian.)
- Try different fwhms.
  - Fit bkgnd and ampl at each one using linear least squares (I think).
  - Compute a chiSq based on the actual profile - the model profile.
  - Iterate to find the fwhm that minimizes chiSq.

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
2004-06-03 ROwen    Modified the module doc string.
2004-07-02 ROwen    Improved the error estimate.
                    Changed GStarFit.err to chiSq
2004-08-03 ROwen    Renamed GStarFit->StarShapeData.
                    Modified to use Constants.FWHMPerSigma.
2004-08-04 ROwen    Improved calculation of ij pixel index and position.
                    Simplified final computation of minimum width parameter.
                    If shape computation fails, converts ArithmeticError into RuntimeError
2004-08-06 ROwen    Fixed invalid variable reference when _StarShapeIterDebug true. 
2004-12-01 ROwen    Modified StarShapeData to use NaN as the default for each argument.
                    Added __all__.
2005-02-07 ROwen    Changed starShape argument ctr (i,) to xyCtr.
2005-04-01 ROwen    Added required argument rad and optional argument bkgnd.
                    No longer iterate the fit with an updated predFWHM because
                    it doesn't seem to help when the data radius is fixed.
                    Added constant _MinRad to constrain the minimum radius.
2005-04-04 ROwen    Bug fix: mis-handled case of bkgnd not specified.
2005-04-15 ROwen    Temporarily hacked the weighting function to see if it makes things better.
                    Added pylab (matplotlib) debugging graphs.
2005-04-22 ROwen    Modified to use nPts as the weighting function.
                    This seems to work slightly better than nPts > 1
                    and just as well as a combination of nPts and a very crude estimate of S/N.
2005-04-25 ROwen    Updated doc string to state that nPts is the weighting function.
                    Removed givenBkgnd argument; it only causes trouble.
2005-05-03 ROwen    Modified to use Brent's method to minimize chiSq.
                    This requires a somewhat messy first pass to bracket fwhm.
                    Also ditched all use of width parameter.
2005-05-18 ROwen    Replaced debuggin flags with verbosity and doPlot.
                    Modified to return an isOK flag and msgStr in StarShapeData
                    instead of raising an exception when fitting fails.
2006-04-17 ROwen    Bug fix: _fitRadProfile had bogus diagnostic print (thanks to pychecker).
                    Removed unused constants _FWHMMin/Max/Delta (thanks to pychecker).
"""
__all__ = ["StarShapeData", "starShape"]

import math
import sys
import traceback
import warnings
import numarray as num
import numarray.ma
import radProf as RP
from Constants import FWHMPerSigma, NaN
import ImUtil
import FitUtil

# minimum radius
_MinRad = 3.0

class StarShapeData:
    """Guide star fit data
    
    Attributes:
    - isOK      if False the fit failed; see msgStr for more info
    - msgStr    warning or error message, if any
    - ampl      profile amplitude (ADUs)
    - bkgnd     background level (ADUs)
    - fwhm      FWHM (pixels)
    - chiSq     chi squared of fit
    """
    def __init__(self,
        isOK = True,
        msgStr = "",
        ampl = NaN,
        fwhm = NaN,
        bkgnd = NaN,
        chiSq = NaN,
    ):
        self.isOK = bool(isOK)
        self.msgStr = msgStr
        self.ampl = float(ampl)
        self.bkgnd = float(bkgnd)
        self.fwhm = float(fwhm)
        self.chiSq = float(chiSq)
    
    def __repr__(self):
        dataList = []
        for arg in ("isOK", "msgStr", "ampl", "fwhm", "bkgnd", "chiSq"):
            val = getattr(self, arg)
            if val not in (None, ""):
                dataList.append("%s=%s" % (arg, val))
        return "%s(%s)" % ", ".join(self.__class__.__name__, dataList)


def starShape(
    data,
    mask,
    xyCtr,
    rad,
    predFWHM = None,
    verbosity = 0,
    doPlot = False,
):
    """Fit a double gaussian profile to a star
    
    Inputs:
    - data      a numarray array of signed integer data
    - mask      a numarray boolean array, or None if no mask (all data valid).
                If supplied, mask must be the same shape as data
                and elements are True for masked (invalid data).
    - med       median (used as the background)
    - xyCtr     x,y center of star; use the convention specified by
                PyGuide.Constants.PosMinusIndex
    - rad       radius of data to fit (pixels);
                values less than _MinRad are treated as _MinRad
    - predFWHM  predicted FWHM. Ignored!
    - verbosity 0: no output, 1: print warnings, 2: print information, 3: print iteration info.
                Note: there are no warnings at this time
    - doPlot    if True, output diagnostics using matplotlib
    """
    if verbosity >= 2:
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
    
    # fit data
    try:
        gsData = _fitRadProfile(radProf, var, nPts, rad, verbosity=verbosity, doPlot=doPlot)
        if verbosity >= 2:
            print "starShape: predFWHM=%.1f; ampl=%.1f; fwhm=%.1f; bkgnd=%.1f; chiSq=%.2f" % \
                (predFWHM, gsData.ampl, gsData.fwhm, gsData.bkgnd, gsData.chiSq)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception, e:
        return StarShapeData(isOK = False, msgStr = str(e))
    
    """Adjust the width for the fact that the centroid
    is not exactly on the center of a pixel
    
    The equivalent sigma^2 of a profile displaced by d from its center
    is sig^2 + d^2/2, so we need to subtract d^2/2 from the sigma^2
    of an offcenter extracted profile to get the true sigma^2.
    Note that this correction is negligable for anything except
    extremely compact stars.
    """
    rawFWHM = gsData.fwhm
    rawSigSq = (FWHMPerSigma * rawFWHM)**2
    corrSigSq = rawSigSq - (0.5 * offSq)
    gsData.fwhm = math.sqrt(corrSigSq) / FWHMPerSigma
    
    if verbosity >= 2:
        print "starShape: ijOff=%.2f, %.2f; offSq=%.2f; rawFWHM=%.3f; corrFWHM=%.3f" % \
            (ijOff[0], ijOff[1], offSq, rawFWHM, gsData.fwhm)
        
    return gsData


def _fitRadProfile(radProf, var, nPts, rad, verbosity=0, doPlot=False):
    """Fit in profile space to determine the width, amplitude, and background.
    
    Inputs:
    - radProf   radial profile around center pixel by radial index
    - var       variance as a function of radius
    - nPts      number of points contributing to profile by radial index
    - rad       radius of data to fit (pixels)
    - verbosity 0: no output, 1: print warnings, 2: print information, 3: print iteration info.
                Note: there are no warnings at this time because warnings are returned in the msgStr field.
    - doPlot    if True, output diagnostics using matplotlib
    
    Returns a StarShapeData object
    """
    if verbosity >= 2:
        print "_fitRadProfile(radProf[%s]=%s\n, var[%s]=%s\n, nPts=%s, rad=%s)" % \
            (len(radProf), radProf, len(var), var, nPts, rad)

    radSq = RP.radSqByRadInd(len(radProf))
    totPnts = num.sum(nPts)
    totCounts = num.sum(nPts*radProf)
    
    # This radial weight is the one used by Jim Gunn and it seems to do as well
    # as anything else I tried. however, it results in a chiSq that is not normalized.
#   radWeight = nPts

    # try a simple normalization
    meanVar = num.sum(var) / float(num.sum(nPts > 1))
    radWeight = nPts / meanVar

    if doPlot:
        try:
            import pylab
        except ImportError:
            warnings.warn("StarShape: cannot import matplotlib.pylab; ignoring doPlot argument")
            pylab = None
    else:
        pylab = None

    if pylab:
        pylab.close()
        pylab.subplot(4,1,1)
        pylab.plot(radProf)
        pylab.subplot(4,1,2)
        pylab.plot(nPts)
        pylab.subplot(4,1,3)
        pylab.plot(radWeight)
    
    def myfunc(fwhm):
        ampl, bkgnd, chiSq, seeProf = _fitIter(radProf, nPts, radWeight, radSq, totPnts, totCounts, fwhm, verbosity=verbosity)
        return chiSq
        
    # brute-force check a lot of values to find a good starting place
    fwhmArr = []
    fwhm = 1.0
    while fwhm < rad*1.5:
        fwhmArr.append(fwhm)
        fwhm += fwhm * 0.1
    nTrials = len(fwhmArr)
#   amplArr = num.zeros([nTrials], num.Float32)
#   bkgndArr = num.zeros([nTrials], num.Float32)
    chiSqArr = num.zeros([nTrials], num.Float32)
    
    minInd = 0
    BadChiSq = 9.9e99
    minChiSq = BadChiSq
    
    for ind in range(nTrials):
        fwhm = fwhmArr[ind]
        ampl, bkgnd, chiSq, seeProf = _fitIter(radProf, nPts, radWeight, radSq, totPnts, totCounts, fwhm)
#       amplArr[ind] = ampl
#       bkgndArr[ind] = bkgnd
        chiSqArr[ind] = chiSq
        
        if ampl > 0 and chiSq < minChiSq:
            minChiSq = chiSq
            minInd = ind

    if pylab:
        pylab.subplot(4,1,4)
        try:
            pylab.plot(fwhmArr, chiSqArr)
        except Exception:
            pass

    if minInd in (0, nTrials-1):
        raise RuntimeError("Could not find bracketing values for fwhm")
            
    firstInd = max(0, minInd - 2)
    lastInd = min(nTrials-1, minInd + 2)
    
    fwhmFirst = fwhmArr[firstInd]
    fwhmMin = fwhmArr[minInd]
    fwhmLast = fwhmArr[lastInd]
    fwhmBracket = (fwhmFirst, fwhmMin, fwhmLast)    

    fwhmMin = FitUtil.brent(myfunc, brack=fwhmBracket)

    # compute final answers at fwhmMin
    ampl, bkgnd, chiSq, seeProf = _fitIter(radProf, nPts, radWeight, radSq, totPnts, totCounts, fwhmMin)

    if doPlot:
        pylab.subplot(4,1,1)
        fitCurve = (seeProf * ampl) + bkgnd
        pylab.plot(fitCurve)
            
    # return StarShapeData containing fit data
    return StarShapeData(
        isOK = True,
        ampl  = ampl,
        fwhm = fwhmMin,
        bkgnd = bkgnd,
        chiSq = chiSq,
    )

def _fitIter(radProf, nPts, radWeight, radSq, totPnts, totCounts, fwhm, verbosity=0):
    # compute the seeing profile for the specified width parameter
    seeProf = _seeProf(radSq, fwhm)
    
    # compute sums
    nPtsSeeProf = nPts*seeProf # temporary array
    sumSeeProf = num.sum(nPtsSeeProf)
    sumSeeProfSq = num.sum(nPtsSeeProf*seeProf)
    sumSeeProfRadProf = num.sum(nPtsSeeProf*radProf)

    if verbosity >= 3:
        print "_fitIter sumSeeProf=%s, sumSeeProfSq=%s, totCounts=%s, sumSeeProfRadProf=%s, totPnts=%s" % \
            (sumSeeProf, sumSeeProfSq, totCounts, sumSeeProfRadProf, totPnts)

    # compute amplitude and background
    # using standard linear least squares fit equations
    # (predicted value = bkgnd + ampl * seeProf)
    try:
        disc = (totPnts * sumSeeProfSq) - sumSeeProf**2
        ampl  = ((totPnts * sumSeeProfRadProf) - (totCounts * sumSeeProf)) / disc
        bkgnd = ((sumSeeProfSq * totCounts) - (sumSeeProf * sumSeeProfRadProf)) / disc
        # diff is the weighted difference between the data and the model
        diff = radProf - (ampl * seeProf) - bkgnd
        chiSq = num.sum(radWeight * diff**2) / totPnts
    except ArithmeticError, e:
        sys.stderr.write("_fitIter failed on fwhm=%s\n" % fwhm)
        traceback.print_exc(file=sys.stderr)
        raise RuntimeError("Could not compute shape: %s" % e)

    if verbosity >= 3:
        print "_fitIter: ampl=%s; bkgnd=%s; fwhm=%s; chiSq=%.2f" % \
            (ampl, bkgnd, fwhm, chiSq)
    
    return ampl, bkgnd, chiSq, seeProf

def _seeProf(radSq, fwhm):
    """Computes the predicted star profile for the given width parameter.
    
    Inputs:
    - radSq     array of radius squared values
    - fwhm      desired fwhm
    """
    norm = 1.0/1.1
    x = radSq * (-0.5 * (FWHMPerSigma / fwhm)**2)
    return (num.exp(x) + (0.1*num.exp(0.25*x))) * norm

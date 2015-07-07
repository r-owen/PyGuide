"""Image processing utilities.

History:
2004-04-05 ROwen    First release.
2004-04-16 ROwen    Added getQuartile, skyStats.
2004-04-21 ROwen    Added verbosity arg to skyStats and removed _SkyStatsDebug.
2004-08-04 ROwen    Removed documentation for nonexistent verbosity argument for getQuartile.
2004-12-01 ROwen    Added __all__.
2005-02-08 ROwen    Added ijIndFromXYPos, ijPosFromXYPos, xyPosFromIJPos, ds9PosFromXYPos, xyPosFromDS9Pos.
                    Changed subFrameCtr arguments ctr and size (i,j) to xyCtr and xySize.
2005-05-16 ROwen    Rewrote subFrameCtr to return a (new) SubFrame object.
                    Added openDS9Win.
2005-05-18 ROwen    Modified skyStats to  accept any kind of array, to take a new argument "thresh",
                    to return an ImStats object and to compute dataCut.
                    Modified getQuartile to reject too-short data.
                    and to handle short data more gracefully.
                    Bug fix: skyStats was using ma.compressed().raw_data() instead of just ma.compressed()
                    (though I doubt this change has any real effect).
2006-04-17 ROwen    Renamed subIJPosOK->subIJOK to match subXYOK.
                    Renamed subIJPosOK argument subIJPos -> subIJ to match other args.
                    Renamed subXYOK argumeng subXYPos -> subXY to match other args.
                    Bug fix: subIJPosOK had an argument name issue.
                    Bug fix: test code broken.
                    Note: thanks to pychecker for catching most of these problems.
2009-11-20 ROwen    Modified to use numpy.
"""
__all__ = ["ImStats", "getQuartile", "skyStats", "subFrameCtr",
    "ijIndFromXYPos", "ijPosFromXYPos", "xyPosFromIJPos",
    "ds9PosFromXYPos", "xyPosFromDS9Pos",
]

import math
import warnings
import numpy
import Constants

_QuartileResidRatios = (
    (1.0, 0.0),
    (1.0, 1.0/3.0),
    (1.0, 1.0),
    (1.0/3.0, 1.0),
)
def getQuartile(sortedData, qnum):
    """Returns a quartile.
    Inputs:
    - sortedData    data sorted in ascending order
    - qnum  which quartile? 1=Q1, 2=Q2=median, 3=Q3
    
    If the position of the quartile is between two values,
    uses linear interpolation to compute the value.

    If the input data is not sorted, returns a meaningless number.
    If qnum not in (1, 2, 3) or len(sortedData) < 3 raises ValueError.
    """
    if qnum not in (1, 2, 3):
        raise ValueError("qnum=%r must be 1, 2 or 3" % qnum)
    dataLen = len(sortedData)
    if dataLen < 3:
        raise ValueError("sortedData too short; len = %s < 3" % len(sortedData))
    ratios = _QuartileResidRatios[((dataLen-1) * qnum) % 4]
    ind0 = (dataLen-1) * qnum // 4
    return ((sortedData[ind0] * ratios[0]) + (sortedData[ind0+1] * ratios[1])) / (ratios[0] + ratios[1])


class ImStats:
    """Information about an image
    (including the settings use to obtain that info).
    
    Values are None if unknown.
    
    - med       median
    - stdDev    std dev
    - nPts      number of points used to compute med and stdDev
    - thresh    threshold used to detect signal
    - dataCut   data cut level
    
    If outerRad != None then med and stdDev are for pixels
    outside a circle of radius "rad" and inside a square
    of size outerRad*2 on a side.
    Otherwise the region used to determine the stats is unknown.
    """
    def __init__(self,
        med = None,
        stdDev = None,
        nPts = None,
        thresh = None,
        dataCut = None,
    ):
        self.med = med
        self.stdDev = stdDev
        self.nPts = nPts
        self.thresh = thresh
        self.dataCut = dataCut
    
    def __repr__(self):
        dataList = []
        for arg in ("med", "stdDev", "nPts", "thresh", "dataCut"):
            val = getattr(self, arg)
            if val not in (None, ""):
                dataList.append("%s=%s" % (arg, val))
        return "%s(%s)" % (self.__class__.__name__, ", ".join(dataList))


def skyStats(
    dataArr,
    thresh = Constants.DefThresh,
    verbosity = 0,
):
    """Computes sky statistics.
    
    Inputs:
    - dataArr: an n-dimensional numpy array or numpy.ma masked array
    - thresh: a threshold for valid data: dataCut = med + (stdDev * thresh);
        values less than PyGuide.Constants.MinThresh are silently increased
    - verbosity 0: no output, 1: print warnings, 2: print information

    Returns an ImStats object containing median, std dev, etc.
    
    The statistics are computed using quartiles. Pixel values above
    2.35 * stdDev are ignored, and this computation is iterated a few times
    to refine the cutoff point.
    
    Standard deviation is computed as stdDev = 0.741 * (Q3 - Q1)
    """
    # creating sorted data
    if isinstance(dataArr, numpy.ma.masked_array):
        sortedData = dataArr.compressed()
    else:
        sortedData = dataArr.flatten()
    dataLen = len(sortedData)
    if verbosity >= 2:
        print "skyStats sorting %d elements" % (dataLen)
    sortedData = numpy.sort(sortedData)
    
    # find sky stats; the iteration improves the values slightly
    MaxIter = 3
    for ii in range(1, MaxIter+1):
        q1, med, q3 = [getQuartile(sortedData[0:dataLen], qnum) for qnum in (1, 2, 3)]
        stdDev = 0.741 * (q3 - q1)
        cutVal = med + (2.35 * stdDev)
        if verbosity >= 2:
            print "skyStats med=%s, q1=%s, q4=%s, stdDev=%s, cutVal=%s" % (med, q1, q3, stdDev, cutVal)
        if ii == MaxIter:
            break
        cutInd = numpy.searchsorted(sortedData, [cutVal])[0]
        if verbosity >= 2:
            print "skStats cutInd=%d, sortedData[cutInd]=%d" % (cutInd, sortedData[cutInd])
        if cutInd < 3:
            if verbosity >= 1:
                print "skStats aborting iteration at step %s; not enough data to cut further" % (ii,)
            break
        dataLen = cutInd
    
    thresh = max(Constants.MinThresh, float(thresh))
    dataCut = med + (stdDev * thresh)
        
    return ImStats(
        med = med,
        stdDev = stdDev,
        nPts = dataLen,
        thresh = thresh,
        dataCut = dataCut,
    )


class SubFrame:
    """Create a subframe and provide useful utility methods.
    
    Inputs:
    - dataArr       data array
    - desBegInd     desired starting i,j index (inclusive)
    - desEndInd     desired ending i,j index (exclusive)
    
    The input indices are called "desired" because they need not
    actually be valid indices. The actual indices used are
    silently shrunk to fit if required.
    """
    def __init__(
        self,
        dataArr,
        desBegInd,
        desEndInd,
    ):
        self.dataArr = numpy.array(dataArr)
        #print "SubFrame(data%s, desBegInd=%s, desEndInd=%s)" % (self.dataArr.shape, desBegInd, desEndInd)
        
        # round desired i,j index (just in case)
        self.desBegInd = [int(round(val)) for val in desBegInd]
        self.desEndInd = [int(round(val)) for val in desEndInd]
        
        # truncated desired i,j index to get actual i,j index
        self.begInd = [max(self.desBegInd[ii], 0) for ii in (0, 1)]
        self.endInd = [min(self.desEndInd[ii], self.dataArr.shape[ii]) for ii in (0, 1)]
        
        # compute amount truncated at beginning and end
        self.begCut = [self.begInd[ii] - self.desBegInd[ii] for ii in (0,1)]
        self.endCut = [self.desEndInd[ii] - self.endInd[ii] for ii in (0,1)]
        #print "SubFrame: begInd=%s; endInd=%s; begCutIJ=%s, endCutIJ=%s" % (self.begInd, self.endInd, self.begCut, self.endCut)

    def getIJLim(self):
        """Return (min i, min j, max i, max j) of subframe in full frame coords"""
        return tuple(self.begInd) + tuple(self.endInd)
    
    def getSubFrame(self):
        """Return the subframe as a numpy array.
        Warning: this is a pointer to the data in the full frame, not a copy.
        """
        return self.dataArr[self.begInd[0]:self.endInd[0], self.begInd[1]:self.endInd[1]]
    
    def fullIJFromSubIJ(self, subIJ):
        """Convert ij position from sub frame to full frame coords.
        Does not check range.
        """
        return [subIJ[ii] + self.begInd[ii] for ii in (0,1)]
    
    def subIJFromFullIJ(self, fullIJ):
        """Convert ij position from full frame to sub frame coords.
        Does not check range.
        """
        return [fullIJ[ii] - self.begInd[ii] for ii in (0,1)]

    def fullXYFromSubXY(self, subXY):
        """Convert xy position from sub frame to full frame coords.
        Does not check range.
        """
        return [subXY[ii] + self.begInd[1-ii] for ii in (0,1)]
        
    def subXYFromFullXY(self, fullXY):
        """Convert xy position from full frame to sub frame coords
        Does not check range.
        """
        return [fullXY[ii] - self.begInd[1-ii] for ii in (0,1)]
    
    def subIJOK(self, subIJ):
        """Return True if subIJ is in bounds.
        """
        for ii in (0,1):
            subIJInd = int(round(subIJ[ii]))
            if self.begInd[ii] > subIJInd[ii] or self.endInd[ii] < subIJInd[ii]:
                return False
        return True

    def subXYOK(self, subXY):
        """Return True if subXY is in bounds.
        """
        return self.subIJOK(ijPosFromXYPos(subXY))


def subFrameCtr(data, xyCtr, xySize):
    """Extract a subframe from a 2d array given a center and size.
    
    Return a SubFrame object.

    Inputs:
    - data      2-d array of data [i,j]
    - xyCtr     desired x,y center of subframe (may be float)
    - xySize    desired x,y size of subframe (may be float);
                e.g. (5,7) returns a 5x7 subframe
    
    Returns the following:
    - subFrame  a SubFrame object; see getSubFrame to get the data as a numpy array.
    
    With ctr in the middle of a pixel:
    - a size of 0.0 to 1.999... returns 1 pixel
    - a size of 2.0 to 3.999... returns 3 pixels
    
    With ctr at the boundary between two pixels:
    - A size of 0.0 to 0.999... returns 0 pixels
    - A size of 1.0 to 2.999... returns 2 pixels
    
    If the requested subframe extends outside the boundaries of data,
    the subframe is truncated without warning. You can tell if it
    has been truncated by comparing its size to the requested size.
    
    If size is odd and the entire subframe fits in data without truncation,
    then xyCtr is truly centered.
    """
    #print "subFrameCtr(data%s, xyCtr=%s, xySize=%s)" % (data.shape, xyCtr, xySize)
    ijCtr = ijPosFromXYPos(xyCtr)
    ijRad = [xySize[ii] / 2.0 for ii in (1, 0)]
    
    desBegInd = [int(math.ceil(ijCtr[ii] - ijRad[ii])) for ii in (0, 1)]
    desEndInd = [int(math.floor(ijCtr[ii] + ijRad[ii])) + 1 for ii in (0, 1)]
    return SubFrame(data, desBegInd, desEndInd)

def ijIndFromXYPos(xyPos):
    """Return the integer index of the pixel whose center is nearest the specified position.
    In other words, the same as ijPosFromXYPos but rounded to the nearest int.
    
    x,y position convention is defined by PyGuide.Constants.PosMinusIndex (which see).
    """
    # use floor(0.5 + x) instead of round(x) because round(-0.5) is -1, not 0,
    # making (0.0, 0,0) convert to (-1,-1) instead of (0,0)
    return [int(round(xyPos[ii] - Constants.PosMinusIndex)) for ii in (1, 0)]
    
def ijPosFromXYPos(xyPos):
    """Convert from x,y position to i,j position.
    
    x,y position convention is defined by PyGuide.Constants.PosMinusIndex (which see).
    i,j position has the axes swapped and has (0,0) as the center of the (0,0) pixel.
    """
    return [float(xyPos[ii] - Constants.PosMinusIndex) for ii in (1, 0)]
    
def xyPosFromIJPos(ijInd):
    """Return the x,y position corresponding to the specified i,j position.
    
    x,y position is defined by PyGuide.Constants.PosMinusIndex (which see).
    i,j position has the axes swapped and has (0,0) as the center of the (0,0) pixel.
    """
    return [float(ijInd[ii] + Constants.PosMinusIndex) for ii in (1, 0)]
    
def ds9PosFromXYPos(xyPos):
    """Convert from PyGuide's x,y position to ds9 x,y position.
    """
    return [float(pos - Constants.PosMinusIndex + 1.0) for pos in xyPos]

def xyPosFromDS9Pos(ds9Pos):
    """Convert from ds9 x,y position to PyGuide's x,y position.
    """
    return [float(pos + Constants.PosMinusIndex - 1.0) for pos in ds9Pos]

def openDS9Win(title=Constants.DS9Title, doRaise=False):
    """Open a ds9 window with default title Constants.DS9Title.
    Issue a warnings.UserWarning and return None on error.
    """
    try:
        import RO.DS9
        return RO.DS9.DS9Win(title, doRaise=doRaise)
    except Exception as e:
        warnings.warn("Could not open ds9 window: %s" % (e,))
    return None

if __name__ == "__main__":
    a = numpy.arange(121, shape=[11, 11])
    print a
    
    for xyCtr, xySize in [
        ((5.5, 5.5), (5.9999, 6.0)),
        ((5.0, 5.0), (4.9999, 5.0)),
        ((1.5, 1.5), (5.9999, 6.0)),
        ((1.0, 1.0), (4.9999, 5.0)),
        ((9.5, 9.5), (5.9999, 6.0)),
        ((10.0, 10.0), (4.9999, 5.0)),
    ]:
        sub = subFrameCtr(a, xyCtr, xySize)
        print "shape=%s, xyCtr=%s, xySubCtr=%s" % (sub.getSubFrame().shape, xyCtr, sub.subXYFromFullXY(xyCtr))
        print

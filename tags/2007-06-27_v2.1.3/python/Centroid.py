"""Measure centroids.

Note: as with all PyGuide routines, the coordinate system origin
is specified by PosMinusIndex.

Warnings:
- Will be thrown off by hot pixels. This could perhaps
be improved by centroiding median-filtered data. The question
is whether the median filtering would adversely affect
centroids, especially for faint objects. This is especially
a concern because at present I have no code to do a proper
median filter of masked data.
- The measure of asymmetry is supposed to be normalized,
but it gets large for bright objects with lots of masked pixels.
This may be simply because the value is only computed at the nearest
integer pixel or because the noise is assumed gaussian, or some error.

How it works:
1) Verify that there is usable signal (optional):
  - Measure median and standard deviation of background
  - Look for connected regions with value > med + (std dev * thresh);
    if no such regions at least 2x2 in size are found
    within a circle of radius rad centered at xyGuess
    then reject the field with "No stars found".

2) Centroid the object (done in basicCentroid):
  - The centroid is the point of mimimum radial asymmetry:
        sum over rad of var(rad)^2 / weight(rad)
    where weight is the expected sigma of var(rad) due to pixel noise:
        weight(rad) = pixNoise(rad) * sqrt(2(numPix(rad) - 1))/numPix(rad)
        pixNoise(rad) = sqrt((readNoise/ccdGain)^2 + (meanVal(rad)-bias)/ccdGain)

  - The minimum is found in two stages:
    a) Find the pixel with the minimum radAsymm.
    The direction to walk is determined by measuring radAsymm at 9 points.
    Each step is one pixel along x and/or y.

    b) Find the true centroid (to better than one pixel) by applying
    a quadratic fit to the 3x3 radAsymm matrix centered on the
    pixel of minimum radAsymm. Only the points along +/-x and +/-y
    are used for this fit; the diagonals are ignored.

3) Conform that there is usable signal (optional):
  - This is just like step 1 but is performed at the final centroid position.
    This makes sure the centroider has not wandered off into a region of bad signal.
    It also allows us to return background statistics for the region about the star.

To do:
- Improve the estimate of centroid error.
- Cosider smoothing the data before centroiding it to handle cosmic rays;
    perhaps gaussian smoothing or a median filter.
  I have at least two major concerns about this:
    - Masked pixels must be handled correctly. This is tricky to do well,
      especially because the brightest pixels may be against the edge
      of a slit mask.
    - Make sure the smoothing does not degrade the accuracy of the centroid
      by too much.

Acknowledgements:
- The centroiding algorithm was invented by Jim Gunn
- The code uses a new asymmetry weighting function
    developed with help from Connie Rockosi
- This code is adapted from the SDSS centroiding code,
    which was written by Jim Gunn and cleaned up by Connie Rockosi.
    
History:
2004-03-22 ROwen    First release.
2004-04-07 ROwen    Packaged as part of PyGuide and moved test code elsewhere.
                    Also changed array data types to match changes in radProf.
2004-04-12 ROwen    Modified centroid to return totCounts.
2004-04-16 ROwen    Modified centroid to not return minAsymm.
2004-04-30 ROwen    Modified to truncate initGuess (i.e. anywhere within a pixel
                    selects that pixel) and round radius to the nearest integer.
                    Bug fix: was converting to Int16 instead of UInt16.
2004-06-03 ROwen    Modified to use the initial guess without modification.
2004-08-03 ROwen    Finally added a measure of centroiding error.
2004-08-06 ROwen    Weight asymmetry calculation by radial noise.
2004-08-25 ROwen    Added _MinRad, to more reliably centroid small stars.
                    Added __all__.
2004-10-14 ROwen    Stopped computing several unused variables. Improved import of radProf.
2005-02-07 ROwen    Changed centroid initGuess (i,j) argument to xyGuess.
                    Changed returned Centroid data object fields ctr (i,j) to xyCtr, err (i,j) to xyErr.
2005-03-31 ROwen    Improved debug output and the efficiency of the "walked too far" test.
                    Noted that rad in CentroidData is integer.
2005-04-01 ROwen    Modified to use round to round the radius instead of adding 0.5 and truncating.
2005-04-11 CLoomis  Added ds9 flag (as per FindStars).
2005-05-18 ROwen    Major overhaul to make sure centroid has usable signal:
                    - Modified centroid to measure usable signal.
                    - Renamed old centroid to basicCentroid
                      and modified it to count saturated pixels.
                    - Replaced bias, etc. with ccdInfo, a PyGuide.CCDInfo object.
                    - Renamed ds9 argument to doDS9.
                    - Default thresh is now set by Constants.DefThresh.
                    - Added verbosity argument, which replaces _CTRDEBUG and _CTRITERDBUG.
                    - Added isOK,  msgStr and imStats fields to CentroidData.
                    - Both basicCentroid and centroid now always return normally
                      unless some serious internal error occurs;
                      if centroiding fails then centroidData.isOK is False.
2005-05-20 ROwen    Added checkSignal method.
                    Rewrite centroid to check for usable signal first, then centroid,
                    then check signal again. Both signal checks are optional
                    The second check avoids trouble when centroid walks into a region of pure noise.
                    Bug fix: conditionMask was mis-handling mask=None.
                    Modified centroid to use conditionData and conditionMask.
                    Stopped auto-tiling frames for doDS9.
2005-10-14 ROwen    Added satMask argument to centroid and basicCentroid;
                    they now ignore ccdInfo.satLevel.
                    Modified to use Float32 omage data instead of UInt16.
                    
2006-04-06 ROwen    CentroidData: stores null ImUtil.ImStats instance if imStats=None.
                    checkSignal:
                    - bug fix: sometimes returned the wrong thing
                    - improved compatibility when the subregion has no pixels
2006-04-17 ROwen    Ditch unused "import warnings" (thanks to pychecker).
"""
__all__ = ['CentroidData', 'centroid',]

import math
import sys
import traceback
# import warnings
import numarray as num
import numarray.nd_image
import numarray.ma
import radProf
import Constants
import ImUtil

def _fmtList(alist):
    """Return "alist[0], alist[1], ..."
    """
    return str(alist)[1:-1]

# parameters
_MinRad = 3.0       # minimum radius
_OuterRadAdd = 10   # amount to add to rad to get outerRad
_MaxIter = 40       # max # of iterations
_MinPixForStats = 20    # minimum # of pixels needed to measure med and std dev

class CentroidData:
    """Centroid data, including the following fields:
    
    flags; check before paying attention to the remaining data:
    - isOK      if False then centroiding failed; see msgStr for more info
    - msgStr    warning or error message, if any
    - nSat      number of saturated pixels; None if unknown

    basic info:
    - rad       radius for centroid search (pix)
    - imStats   image statistics such as median and std. dev. (if known); an ImUtil.ImStats object.
    
    star data:
    - xyCtr     the x,y centroid (pixels)
    - xyErr     the predicted 1-sigma uncertainty in xyCtr (pixels)

    note: the following three values are computed for that radial profile
    centered on the pixel nearest the centroid (NOT the true centroid):

    - asymm     measure of asymmetry:
                    sum over rad of var(rad)^2 / weight(rad)
                where weight is the expected sigma of var(rad) due to pixel noise:
                    weight(rad) = pixNoise(rad) * sqrt(2(numPix(rad) - 1))/numPix(rad)
                    pixNoise(rad) = sqrt((readNoise/ccdGain)^2 + (meanVal(rad)-bias)/ccdGain)
    - pix       the total number of unmasked pixels (ADU)
    - counts    the total number of counts (ADU)
    
    Warning: asymm is supposed to be normalized, but it gets large
    for bright objects with lots of masked pixels. This may be
    simply because the value is only computed at the nearest integer pixel
    or because the noise is assumed gaussian, or some error.
    
    Suggested use:
    - check isOK; if False do not use the data
    - check nSat(); if not None and more than a few then be cautious in using the data
        (I don't know how sensitive centroid accuracy is to # of saturated pixels)
    """
    def __init__(self,
        isOK = True,
        msgStr = "",
        nSat = None,
        rad = None,
        imStats = None,
        xyCtr = None,
        xyErr = None,
        asymm = None,
        pix = None,
        counts = None,
    ):
        self.isOK = isOK
        self.msgStr = msgStr
        self.nSat = nSat
        
        self.rad = rad
        if imStats == None:
            imStats = ImUtil.ImStats()
        self.imStats = imStats

        self.xyCtr = xyCtr
        self.xyErr = xyErr
        
        self.asymm = asymm
        self.pix = pix
        self.counts = counts
    
    def __repr__(self):
        dataList = []
        for arg in ("isOK", "msgStr", "nSat", "xyCtr", "xyErr", "asymm", "pix", "counts", "rad", "imStats"):
            val = getattr(self, arg)
            if val not in (None, ""):
                dataList.append("%s=%s" % (arg, val))
        return "%s(%s)" % (self.__class__.__name__, ", ".join(dataList))


def basicCentroid(
    data,
    mask,
    satMask,
    xyGuess,
    rad,
    ccdInfo,
    verbosity = 0,
    doDS9 = False,
):
    """Compute a centroid.

    Inputs:
    - data      image data [i,j]
    - mask      a mask of invalid data (1 if invalid, 0 if valid); None if no mask.
    - satMask   a maks of of saturated pixels (1 if saturated, 0 if not); None if no mask.
    - xyGuess   initial x,y guess for centroid
    - rad       radius of search (pixels);
                values less than _MinRad are treated as _MinRad
    - ccdInfo   ccd bias, gain, etc.; a PyGuide.CCDInfo object
    - verbosity 0: no output, 1: print warnings, 2: print information,
                3: print basic iteration info, 4: print detailed iteration info.
                Note: there are no warnings at this time because the relevant info is returned.
    - doDS9     if True, diagnostic images are displayed in ds9
    
    Masks are optional. If specified, they must be the same shape as "data"
    and should be of type Bool. None means no mask (all data is OK).
        
    Returns a CentroidData object (which see for more info),
    but with no imStats info.
    """
    if verbosity > 1:
        print "basicCentroid(xyGuess=%s, rad=%s, ccdInfo=%s)" % (xyGuess, rad, ccdInfo)
    # condition and check inputs
    data = conditionData(data)
    mask = conditionMask(mask)
    satMask = conditionMask(satMask)
    if len(xyGuess) != 2:
        raise ValueError("initial guess=%r must have 2 elements" % (xyGuess,))
    rad = int(round(max(rad, _MinRad)))

    # compute index of pixel closest to initial guess
    ijIndGuess = ImUtil.ijIndFromXYPos(xyGuess)

    if doDS9:
        ds9Win = ImUtil.openDS9Win()
    else:
        ds9Win = None

    if ds9Win:
        # show masked data in frame 1 and unmasked data in frame 2
        ds9Win.xpaset("frame 1")
        if mask != None:
            ds9Win.showArray(data * (mask==0))
        else:
            ds9Win.showArray(data)
        ds9Win.xpaset("frame 2")
        ds9Win.showArray(data)
        ds9Win.xpaset("frame 1")

        # display circle showing the centroider input in frame 1
        args = list(xyGuess) + [rad]
        ds9Win.xpaset("regions", "image; circle %s # group=ctrcirc" % _fmtList(args))
    
    try:
        # OK, use this as first guess at maximum. Extract radial profiles in
        # a 3x3 gridlet about this, and walk to find minimum fitting error
        maxi, maxj = ijIndGuess
        asymmArr = num.zeros([3,3], num.Float64)
        totPtsArr = num.zeros([3,3], num.Int32)
        totCountsArr = num.zeros([3,3], num.Float64)
        
        niter = 0
        while True:
            niter += 1
            if niter > _MaxIter:
                raise RuntimeError("could not find a star in %s iterations" % (niter,))
            
            for i in range(3):
                ii = maxi + i - 1
                for j in range(3):
                    jj = maxj + j - 1
                    if totPtsArr[i, j] != 0:
                        continue
                    asymmArr[i, j], totCountsArr[i, j], totPtsArr[i, j] = radProf.radAsymmWeighted(
                        data, mask, (ii, jj), rad, ccdInfo.bias, ccdInfo.readNoise, ccdInfo.ccdGain)
# this version omits noise-based weighting
# (warning: the error estimate will be invalid and chiSq will not be normalized)
#                   asymmArr[i, j], totCountsArr[i, j], totPtsArr[i, j] = radProf.radAsymm(
#                       data, mask, (ii, jj), rad)
    
                    if verbosity > 3:
                        print "basicCentroid: asymm = %10.1f, totPts = %s, totCounts = %s" % \
                            (asymmArr[i, j], totPtsArr[i, j], totCountsArr[i, j])
    
            # have error matrix. Find minimum
            ii, jj = num.nd_image.minimum_position(asymmArr)
            ii -= 1
            jj -= 1
    
            if verbosity > 2:
                print "basicCentroid: error matrix min ii=%d, jj=%d, errmin=%5.1f" % (ii, jj, asymmArr[ii,jj])
                if verbosity > 3:
                    print "basicCentroid: asymm matrix =\n", asymmArr
    
            if (ii != 0 or jj != 0):
                # minimum error not in center; walk and try again
                maxi += ii
                maxj += jj
                if verbosity > 2:
                    print "shift by", -ii, -jj, "to", maxi, maxj
    
                if ((maxi - ijIndGuess[0])**2 + (maxj - ijIndGuess[1])**2) >= rad**2:
                    raise RuntimeError("could not find star within %r pixels" % (rad,))
                
                # shift asymmArr and totPtsArr to minimum is in center again
                asymmArr = num.nd_image.shift(asymmArr, (-ii, -jj))
                totCountsArr = num.nd_image.shift(totCountsArr, (-ii, -jj))
                totPtsArr = num.nd_image.shift(totPtsArr, (-ii, -jj))
            else:
                # Have minimum. Get out and go home.
                break
    
        if verbosity > 2:
            print "basicCentroid: found ijMax=%s after %r iterations" % ((maxi, maxj), niter,)
        
        # perform a parabolic fit to find true centroid
        # and compute the error estimate
        # y(x) = ymin + a(x-xmin)^2
        # a = (y0 - 2y1 + y2) / 2
        # xmin = b/2a where b = (y2-y0)/2
        # ymin = y1 - b^2/4a  but this is tricky in 2 dimensions so we punt
        # for a given delta-y, delta-x = sqrt(delta-y / a)
        ai = 0.5 * (asymmArr[2, 1] - 2.0*asymmArr[1, 1] + asymmArr[0, 1])
        bi = 0.5 * (asymmArr[2, 1] - asymmArr[0, 1])
        aj = 0.5 * (asymmArr[1, 2] - 2.0*asymmArr[1, 1] + asymmArr[1, 0])
        bj = 0.5 * (asymmArr[1, 2] - asymmArr[1, 0])
        
        di = -0.5*bi/ai
        dj = -0.5*bj/aj
        ijCtr = (
            maxi + di,
            maxj + dj,
        )
        xyCtr = ImUtil.xyPosFromIJPos(ijCtr)
    
        if verbosity > 2:
            print "basicCentroid: asymmArr[0:3, 1]=%s, ai=%s, bi=%s, di=%s, iCtr=%s" % (asymmArr[0:3, 1], ai, bi, di, ijCtr[0])
            print "basicCentroid: asymmArr[1, 0:3]=%s, aj=%s, bj=%s, dj=%s, jCtr=%s" % (asymmArr[1, 0:3], aj, bj, dj, ijCtr[1])
        
        # crude error estimate, based on measured asymmetry
        # note: I also tried using the minimum along i,j but that sometimes is negative
        # and this is already so crude that it's not likely to help
        radAsymmSigma = asymmArr[1,1]
        iErr = math.sqrt(radAsymmSigma / ai)
        jErr = math.sqrt(radAsymmSigma / aj)
        xyErr = (jErr, iErr)
    
        if ds9Win:
            # display x at centroid
            ds9Win.xpaset("regions", "image; x point %s # group=centroid" % \
                        _fmtList(xyCtr))
    

        # count # saturated pixels, if satMask available
        if satMask == None:
            nSat = None
        else:
            ctrPixIJ = (maxi, maxj)
            ctrPixXY = ImUtil.xyPosFromIJPos(ctrPixIJ)
            subSize = (rad*2) + 1
            subSatMaskObj = ImUtil.subFrameCtr(
                satMask,
                xyCtr = ctrPixXY,
                xySize = (subSize, subSize),
            )
            subSatMask = subSatMaskObj.getSubFrame()
            subCtrIJ = subSatMaskObj.subIJFromFullIJ(ctrPixIJ)

            def makeDisk(i, j):
                return ((i-subCtrIJ[0])**2 + (j-subCtrIJ[1])**2) <= rad**2
            maybeSatPixel = num.fromfunction(makeDisk, subSatMask.shape)

            if mask != None:
                subMaskObj = ImUtil.subFrameCtr(
                    mask,
                    xyCtr = ctrPixXY,
                    xySize = (subSize, subSize),
                )
                subMask = subMaskObj.getSubFrame()
                num.logical_and(maybeSatPixel, num.logical_not(subMask), maybeSatPixel)
            
            num.logical_and(subSatMask, maybeSatPixel, maybeSatPixel)
            nSat = num.nd_image.sum(maybeSatPixel)
    
        ctrData = CentroidData(
            isOK = True,
            rad = rad,
            nSat = nSat,
            xyCtr = xyCtr,
            xyErr = xyErr,
            counts = totCountsArr[1,1],
            pix = totPtsArr[1,1],
            asymm = asymmArr[1,1],
        )
        if verbosity > 2:
            print "basicCentroid: %s" % (ctrData,)
        return ctrData

    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception, e:
        if verbosity > 1:
            traceback.print_exc(file=sys.stderr)
        elif verbosity > 0:
            print "basicCentroid failed: %s" % (e,)
        return CentroidData(
            isOK = False,
            msgStr = str(e),
            rad = rad,
        )


def centroid(
    data,
    mask,
    satMask,
    xyGuess,
    rad,
    ccdInfo,
    thresh = Constants.DefThresh,
    verbosity = 0,
    doDS9 = False,
    checkSig = (True, True),
):
    """Centroid and then confirm that there is usable signal at the location.
    
    Inputs:
    - data      image data [i,j]
    - mask      a mask of invalid data (1 if invalid, 0 if valid); None if no mask.
    - satMask   a maks of of saturated pixels (1 if saturated, 0 if not); None if no mask.
    - xyGuess   initial x,y guess for centroid
    - rad       radius of search (pixels);
                values less than _MinRad are treated as _MinRad
    - ccdInfo   ccd bias, gain, etc.; a PyGuide.CCDInfo object
    - thresh    determines the point above which pixels are considered data;
                valid data >= thresh * standard deviation + median
                values less than PyGuide.Constants.MinThresh are silently increased
    - verbosity 0: no output, 1: print warnings, 2: print information, 3: print iteration info.
                Note: there are no warnings at this time
    - doDS9     if True, display diagnostic images in ds9
    - checkSig  Verify usable signal for circle at (xyGuess, xy centroid)?
                If both are false then imStats is not computed.
    
    Returns a CentroidData object (which see for more info).
    """
    if verbosity > 2:
        print "centroid(xyGuess=%s, rad=%s, ccdInfo=%s, thresh=%s)" % (xyGuess, rad, ccdInfo, thresh)
    
    if checkSig[0]:
        signalOK, imStats = checkSignal(
            data = data,
            mask = mask,
            xyCtr = xyGuess,
            rad = rad,
            thresh = thresh,
            verbosity = verbosity,
        )
        if not signalOK:
            return CentroidData(
                isOK = False,
                msgStr = "No star found",
            )

    ctrData = basicCentroid(
        data = data,
        mask = mask,
        satMask = satMask,
        xyGuess = xyGuess,
        rad = rad,
        ccdInfo = ccdInfo,
        verbosity = verbosity,
        doDS9 = doDS9,
    )

    if ctrData.isOK and checkSig[1]:
        signalOK, imStats = checkSignal(
            data = data,
            mask = mask,
            xyCtr = ctrData.xyCtr,
            rad = rad,
            thresh = thresh,
            verbosity = verbosity,
        )
        ctrData.imStats = imStats
        if not signalOK:
            ctrData.isOK = False
            ctrData.msgStr = "No star found"
    return ctrData


def checkSignal(
    data,
    mask,
    xyCtr,
    rad,
    thresh = Constants.DefThresh,
    verbosity = 0,
):
    """Check that there is usable signal in a given circle.
    
    Inputs:
    - data      image data [i,j]
    - mask      a mask [i,j] of 0's (valid data) or 1's (invalid); None if no mask.
                If mask is specified, it must have the same shape as data.
    - xyCtr     center of circule (pixels)
    - rad       radius of circle (pixels);
                values less than _MinRad are treated as _MinRad
    - thresh    determines the point above which pixels are considered data;
                valid data >= thresh * standard deviation + median
                values less than PyGuide.Constants.MinThresh are silently increased
    - verbosity 0: no output, 1: print warnings, 2: print information, 3: print iteration info.
                Note: there are no warnings at this time
    
    Return:
    - signalOK  True if usable signal is present
    - imStats   background statistics: a PyGuide.ImStats object
    
    Details of usable signal:
    - Computes median and stdDev in a region extending from a circle of radius "rad"
      to a box of size (rad+_OuterRadAdd)*2 on a side.
    - median-smooths the data inside a circle of radius "rad" and makes sure
      there is usable signal: max(data) >= thresh*stdDev + median
    """
    if verbosity > 2:
        print "checkSignal(xyCtr=%s, rad=%s, thresh=%s)" % (xyCtr, rad, thresh)

    # check check inputs
    if len(xyCtr) != 2:
        raise ValueError("initial guess=%r must have 2 elements" % (xyCtr,))
    rad = int(round(max(rad, _MinRad)))

    outerRad = rad + _OuterRadAdd
    subDataObj = ImUtil.subFrameCtr(
        data,
        xyCtr = xyCtr,
        xySize = (outerRad, outerRad),
    )
    subData = subDataObj.getSubFrame().astype(num.Float32) # force type and copy
    if subData.size() < _MinPixForStats:
        return False, ImUtil.ImStats(
            nPts = subData.size(),
        )
    subCtrIJ = subDataObj.subIJFromFullIJ(ImUtil.ijPosFromXYPos(xyCtr))
    
    if mask != None:
        subMaskObj = ImUtil.subFrameCtr(
            mask,
            xyCtr = xyCtr,
            xySize = (outerRad, outerRad),
        )
        subMask = subMaskObj.getSubFrame().astype(num.Bool) # force type and copy
    else:
        subMask = num.zeros(subData.shape, type=num.Bool)

    # create circleMask; a centered circle of radius rad
    # with 0s in the middle and 1s outside
    def makeCircle(i, j):
        return ((i-subCtrIJ[0])**2 + (j-subCtrIJ[1])**2) > rad**2
    circleMask = num.fromfunction(makeCircle, subData.shape)
    
    # make a copy of the data outside a circle of radius "rad";
    # use this to compute background stats
    bkgndPixels = num.ma.array(
        subData,
        mask = num.logical_or(subMask, num.logical_not(circleMask)),
    )
    if bkgndPixels.count() < _OuterRadAdd**2:
        # too few unmasked pixels in outer region; try not masking off the star
        bkgndPixels = num.ma.array(
            subData,
            mask = subMask,
        )
        if bkgndPixels.count() < _MinPixForStats:
            return False, ImUtil.ImStats(
                nPts = bkgndPixels.count(),
            )
    
    imStats = ImUtil.skyStats(bkgndPixels, thresh)
    del(bkgndPixels)

    # median filter the inner data and look for signal > dataCut
    dataPixels = num.ma.array(
        subData,
        mask = num.logical_or(subMask, circleMask),
    )
    smoothedData = dataPixels.filled(imStats.med)
    num.nd_image.median_filter(smoothedData, 3, output=smoothedData)
    del(dataPixels)
    
    # look for a blob of at least 2x2 adjacent pixels with smoothed value >= dataCut
    # note: it'd be much simpler but less safe to simply test:
    #    if max(smoothedData) < dataCut: # have signal
    shapeArry = num.ones((3,3))
    labels, numElts = num.nd_image.label(smoothedData>imStats.dataCut, shapeArry)
    del(smoothedData)
    slices = num.nd_image.find_objects(labels)
    for ijSlice in slices:
        minSize = min([slc.stop - slc.start for slc in ijSlice])
        if minSize >= 2:
            return True, imStats
    return False, imStats


def conditionData(data):
    """Convert dataArr to the correct type
    such that basicCentroid can operate most efficiently on it.
    
    Warning: does not copy the data unless necessary.
    """
    return conditionArr(data, desType=num.Float32)

def conditionMask(mask):
    """Convert mask to the correct type
    such that basicCentroid can operate most efficiently on it.
    
    Mask is optional, so a value of None returns None.

    Warning: does not copy the data unless necessary.
    """
    if mask == None:
        return None
    return conditionArr(mask, num.Bool)

def conditionArr(arr, desType):
    """Convert a sequence to a numarray array of the desired type.

    Warning: does not copy the data unless necessary.
    """
    arr = num.asarray(arr, type=desType)
    if not arr.iscontiguous():
        arr = arr.copy()
    return arr

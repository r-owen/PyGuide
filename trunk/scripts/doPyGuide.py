#!/usr/bin/env python -i
"""Measures stars in a given image file, displaying the image in ds9
and reporting star positions and shape information on stdout.

History:
2004-04-16 ROwen
2004-04-29 ROwen    Modified to use new ds9 on findStars.
2004-05-18 ROwen    Modified to set up ds9Win and to use fewer globals.
2004-08-25 ROwen    Modified for 2004-08-06 PyGuide.
2004-10-14 ROwen    Modified to measure starShape.
2004-12-01 ROwen    Renamed function from starUtil to doFindStars to match module name.
                    Replaced arguments with globals to make it easier to change settings.
                    Bug fix: if starShape failed, shapeData was mis-set.
2005-02-07 ROwen    Modified for findStars 1.2.
2005-03-31 ROwen    Allowed one to specify a mask name to doFindStars.
                    Modified doFindStars to allow specifying parameters as keyword arguments.
                    Added showDef() method to display current defaults.
                    Modified global vars so im, not d, is the image array.
2005-04-01 ROwen    Updated for modified FindStars and StarShape.
2005-04-11 ROwen    Modified to use PyGuide.Constants.DS9Title.
2005-04-22 ROwen    Added support for the rad argument.
2005-05-16 ROwen    Modified for overhauled findStars.
2005-05-20 ROwen    Added doCentroid.
2005-06-17 ROwen    Added an invertMask flag.
2005-10-14 ROwen    Added support for satMask.
                    Bug fix: modified to handle nSat=None on output.
2006-04-06 ROwen    Bug fix: doCentroid mis-handled return values from loadFiles.
                    Renamed doFindStars.py -> doPyGuide.py.
2006-06-02 ROwen    Bug fix: doCentroid output would fail without a saturation mask.
                    Star data output is more robust. If formatting fails, the data
                    is printed unformatted (along with the error message).
                    One subroutine is now used to output all star data.
2007-01-23 ROwen    Changed #!/usr/local/bin/python -i to #!/usr/bin/env python -i
2008-10-01 ROwen    Added support for the DATASEC keyword.
                    Removed invertMask support.
                    Added help strings for doFindStars and doCentroid
                    and improved other help strings and various output.
                    Set NUMERIX to make PyFits use numarray
2009-11-20 ROwen    Modified to use numpy.
                    Stop setting NUMERIX.
"""
import os
import sys
import numpy
import pyfits
import PyGuide
import RO.DS9

im = None
imFits = None
mask = None
maskFits = None
satMask = None
satMaskFits = None

# Default Parameters
# these settings are for the new NA2 guider
bias = 1000
readNoise = 21.391
ccdGain = 1.643 # e-/pixel
# these are general settings
thresh = 3.0
radMult = 1.0
rad = None
satLevel = (2**16)-1
verbosity = 0
doDS9 = True

# set up a ds9 window
try:
    ds9Win = RO.DS9.DS9Win(PyGuide.Constants.DS9Title)
except Exception, e:
    print "Cannot use ds9; error = %s" % (e,)
    doDS9 = False

CCDInfoNames = ("bias", "readNoise", "ccdGain", "satLevel")
FindParamNames = CCDInfoNames + ("thresh", "radMult", "rad", "verbosity", "doDS9")
CentroidParamNames = CCDInfoNames + ("thresh", "rad", "verbosity", "doDS9")

def doFindStars(
    imName = None,
    maskName = None,
    satMaskName = None,
    invertMask = False,
    **kargs
):
    """Find stars and centroid and shape-fit them.
    
    Inputs:
    - all the arguments for loadFiles plus values shown by showDef
    """
    global isSat, sd
    im, mask, satMask = loadFiles(imName, maskName, satMaskName, invertMask)
    
    # check keyword arguments
    for paramName in kargs:
        if paramName not in FindParamNames:
            raise RuntimeError("Invalid argument: %s" % (paramName,))
    
    # fill in defaults
    globalDict = globals()
    for paramName in FindParamNames:
        if paramName not in kargs:
            kargs[paramName] = globalDict[paramName]
    
    # split off ccd info
    ccdInfoDict = {}
    for paramName in CCDInfoNames:
        ccdInfoDict[paramName] = kargs.pop(paramName)
    ccdInfo = PyGuide.CCDInfo(**ccdInfoDict)
    
    # find stars and centroid
    print "Calling PyGuide.findStars"
    ctrDataList, imStats = PyGuide.findStars(
        data = im,
        mask = mask,
        satMask = satMask,
        ccdInfo = ccdInfo,
    **kargs)

    print "%s stars found:" % (len(ctrDataList),)
    printStarHeader()
    for ctrData in ctrDataList:
        # measure star shape
        try:
            shapeData = PyGuide.starShape(
                data = im,
                mask = mask,
                xyCtr = ctrData.xyCtr,
                rad = ctrData.rad,
            )
        except RuntimeError, e:
            print "starShape failed: %s" % (e,)
            shapeData = PyGuide.StarShapeData()

        printStarData(ctrData, shapeData)

def doCentroid(
    imName = None,
    maskName = None,
    satMaskName = None,
    xyGuess = None,
    invertMask = False,
    **kargs
):
    """Centroid and shape-fit a star
    
    Inputs:
    - all the arguments for loadFiles plus most values shown by showDef
    """
    global im, imFits, mask, maskFits, satMask, satMaskFits, isSat, sd
    im, mask, satMask = loadFiles(imName, maskName, satMaskName, invertMask)
    if xyGuess == None:
        print "xyGuess is required"
        return
    
    # check keyword arguments
    for paramName in kargs:
        if paramName not in CentroidParamNames:
            raise RuntimeError("Invalid argument: %s" % (paramName,))
    
    # fill in defaults
    globalDict = globals()
    for paramName in CentroidParamNames:
        if paramName not in kargs:
            kargs[paramName] = globalDict[paramName]
    
    # split off ccd info
    ccdInfoDict = {}
    for paramName in CCDInfoNames:
        ccdInfoDict[paramName] = kargs.pop(paramName)
    ccdInfo = PyGuide.CCDInfo(**ccdInfoDict)
    
    # centroid
    ctrData = PyGuide.centroid(
        data = im,
        mask = mask,
        satMask = satMask,
        xyGuess = xyGuess,
        ccdInfo = ccdInfo,
    **kargs)

    if not ctrData.isOK:
        print "centroid failed:", ctrData.msgStr
        return

    shapeData = PyGuide.starShape(
        data = im,
        mask = mask,
        xyCtr = ctrData.xyCtr,
        rad = ctrData.rad,
    )
    # print results
    printStarHeader()
    printStarData(ctrData, shapeData)

    if not shapeData.isOK:
        print "starShape failed:", shapeData.msgStr

def loadFiles(
    imName = None,
    maskName = None,
    satMaskName = None,
    invertMask = False,
):
    """Load a new image and/or mask and/or satMask from a fits file.

    Inputs:
    - imName: path to image FITS file; None to use current image
    - maskName: path to bad pixel mask; 0=good unless invertMask is true;
            None to use current mask, if any
    - satMaskName: path to saturated pixel mask; 0=good regardless of invertMask;
            None to use current mask, if any
    """
    global im, imFits, mask, maskFits, satMask, satMaskFits, isSat, sd
    if imName:
        imFits = pyfits.open(imName)
        print "Loading image %s into imFits and im" % (imName,)
        dataSec = parseDataSec(imFits[0].header.get("DATASEC"))
        dataShape = imFits[0].data.shape
        if dataSec == None:
            dataSec = [0, dataShape[0], 0, dataShape[1]]
        im = imFits[0].data[dataSec[0]:dataSec[1], dataSec[2]:dataSec[3]]
    if maskName:
        print "Loading bad pixel mask %s into maskFits and mask" % (maskName,)
        maskFits = pyfits.open(maskName)
        mask = maskFits[0].data[dataSec[0]:dataSec[1], dataSec[2]:dataSec[3]] > 0.1
    if satMaskName:
        print "Loading saturated pixel mask %s into satMaskFits and satMask" % (satMaskName,)
        satMaskFits = pyfits.open(satMaskName)
        satMask = satMaskFits[0].data[dataSec[0]:dataSec[1], dataSec[2]:dataSec[3]] > 0.1
    return im, mask, satMask

def parseDataSec(dataSecStr):
    """Parse DATASEC and return (beg x, end+1 x, beg y, end+1 y)
    
    DATASEC has an origin of 1 and the end is inclusive (FITS standard) and has x and y swapped
    The return value has an origin of 0 and the end is exclusive (numpy/C++ standard)
    
    Return None if DATASEC is None or cannot be parsed.
    
    Input:
    - dataSecStr: a DATASEC in the form [begx:endx,begy:endy]; if None then None is returned
    
    On error prints a message to stderr and returns None
    """
    if dataSecStr == None:
        return None
    try:
        trimStr = dataSecStr[1:-1]
        xyStrList = trimStr.split(",")
        if len(xyStrList) != 2:
            raise RuntimeError("Could not split %s" % (trimStr,))
        xyStrList.reverse()
        retValList = []
        for strList in xyStrList:
            begEndStrList = strList.split(":")
            if len(begEndStrList) != 2:
                raise RuntimeError("Could not split %s" % (begEndStrList))
            retValList += [int(begEndStrList[0]) - 1, int(begEndStrList[1])]
        return retValList
    except Exception, e:
        sys.stderr.write("Could not parse %r; error=%s" % (dataSecStr, e))
        return None

def printStarHeader():
    """Print star position data header"""
    print "   xctr    yctr    xerr    yerr          ampl     bkgnd    fwhm     rad     pix    nSat   chiSq"

def printStarData(ctrData, shapeData):
    """Print star position data"""
    # print results
    try:
        print "%7.2f %7.2f %7.2f %7.2f %13.1f %9.1f %7.1f %7d %7d %7s %7.1f" % (
            ctrData.xyCtr[0], ctrData.xyCtr[1],
            ctrData.xyErr[0], ctrData.xyErr[1],
            shapeData.ampl, shapeData.bkgnd, shapeData.fwhm,
            ctrData.rad, ctrData.pix, ctrData.nSat, shapeData.chiSq,
        )
    except (ValueError, TypeError), e:
        print "(printing free-form due to format error: %s)" % (e,)
        print ctrData.xyCtr[0], ctrData.xyCtr[1], \
            ctrData.xyErr[0], ctrData.xyErr[1], \
            shapeData.ampl, shapeData.bkgnd, shapeData.fwhm, \
            ctrData.rad, ctrData.pix, ctrData.nSat, shapeData.chiSq

def showDef():
    """Show current value of various global variables
    which are used as defaults for parameters of doFindStars.
    """
    print "Global variables (defaults for doFindStars):"
    globalDict = globals()
    for paramName in FindParamNames:
        print "%s = %s" % (paramName, globalDict[paramName])
    print

showDef()
print """The following variables are available:

Default aguments for doFindStars and doCentroid:
bias        bias remaining in the data, if any (ADU)
readNoise   ccd read noise (e-)
ccdGain     ccd inverse gain (e-/ADU)
satLevel    saturation level (ADU); ignored by PyFits
            but potentially useful for generating saturated pixel masks.
rad         centroid radius (pix); if None then doFindStars computes rad using radMult
radMult     for doFindStars: if rad = None then centroid radius is computed as follows:
            centroid radius = radMult * max(rad * blob size x, rad * blob size y);
thresh      for doFindStars: determines the point above which pixels are considered data;
            valid data >= thresh * standard deviation + median
            values less than PyGuide.Constants.MinThresh are silently increased
verbosity   0: no output, 1: print warnings, 2: print information and
            (if doDS9 true) show smoothed image in ds9 frame 3.
doDS9       if True, shows current image and other info in ds9 in current frame.
            For this to work, you must have the RO package installed.

Computed data:
im          image data array (set by loadFiles)
mask        mask data array, or None if no mask (set by loadFiles)
satMask     saturated mask data array, or None of no saturated mask (set by loadFiles)
sd          star data returned by PyGuide.findStars

Reported values include::
rad         radius used to compute centroid
pix         the number of unmasked pixels used to compute the centroid
chiSq       chi squared for shape fit

Notes:
- For a slitviewer image, be sure to specify a suitable mask.
- For optimal centroiding and a reasonable centroid error estimate
  you must set bias, readNoise and ccdGain correctly for your image.

Function calls:
doFindStars(imName=None, maskName=None, [, optional_named_params])
doCentroid(imName=None, maskName=None, xyGuess=(x,y) [, optional_named_params])
loadFiles(imName, maskname, satMaskName) loads a new image, mask and/or satMask
ds9Win.showArray(arry) displays an array in ds9
showDef() prints the current defaults
"""


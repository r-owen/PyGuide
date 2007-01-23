#!/usr/local/bin/python
"""Generates fake ccd data

History:
2004-07-07 ROwen
2004-07-29 ROwen    Added sky noise to addNoise; thanks to Connie Rockosi
                    for catching this important omission.
2004-08-06 ROwen    Reordered the bias arg in addNoise for consistency.
2004-12-01 ROwen    Added __all__.
2005-01-31 ROwen    Bug fix: used == instead of = for __all__.
2005-02-07 ROwen    Changed fakeStar to accept xyCtr instead of (i,j) ctr.
2005-05-16 ROwen    Modified addNoise to take ccdInfo instead of 3 args.
"""
__all__ = ["fakeStar", "addNoise"]

import numarray as num
import numarray.random_array as rand
import ImUtil

_MaxValUInt16 = 2**16 - 1

def fakeStar(
    arrShape,
    xyCtr,
    sigma,
    ampl,
):
    """Return a 2-d array containing a noise-free double gaussian
    and whose values are truncated at 2**16-1
    
    Inputs:
    - arrShape  desired array shape (2 integers)
    - xyCtr     desired x,y center
    - sigma     desired sigma (float)
    - ampl      desired amplitude (float)
    """
    if len(arrShape) != 2:
        raise ValueError("arrShape=%r must have 2 elements" % (arrShape,))
    if len(xyCtr) != 2:
        raise ValueError("xyCtr=%r must have 2 elements" % (xyCtr,))
    sigma = float(sigma)
    ampl = float(ampl)
    ijCtr = ImUtil.ijPosFromXYPos(xyCtr)

    def peakFunc(i, j):
        radSq = (i - ijCtr[0])**2 + (j - ijCtr[1])**2
        expArg = - radSq / (2.0 * sigma**2)
        gauss = ampl * (num.exp(expArg)  + 0.1*num.exp(0.25*expArg))
        gauss = num.where(gauss <= _MaxValUInt16, gauss, _MaxValUInt16)
        return gauss.astype(num.UInt16)
    return num.fromfunction(peakFunc, arrShape)

def addNoise(
    data,
    sky,
    ccdInfo,
):
    """Adds poisson noise and gaussian read noise to the specified data.
    
    Inputs:
    - data      noiseless image, in ADU
    - sky       sky level, in ADU
    - ccdInfo   a PyGuide.CCDInfo object
    """
    outData = num.add(data, sky).astype(num.Int32)
    outData = rand.poisson(mean = outData * ccdInfo.ccdGain) / ccdInfo.ccdGain
    outData += rand.normal(mean = ccdInfo.bias, std = ccdInfo.readNoise/float(ccdInfo.ccdGain), shape = data.shape)
    # truncate data and return as UInt16
    outData = num.where(outData >= 0, outData, 0)
    outData = num.where(outData <= _MaxValUInt16, outData, _MaxValUInt16)
    return outData.astype(num.UInt16)

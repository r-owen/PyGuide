#!/usr/local/bin/python
"""Generates fake ccd data

History:
2004-07-07 ROwen
2004-07-29 ROwen	Added sky noise to addNoise; thanks to Connie Rockosi
					for catching this important omission.
"""
import numarray as num
import numarray.random_array as rand

_MaxValUInt16 = 2**16 - 1


def fakeStar(
	arrShape,
	ctr,
	sigma,
	ampl,
):
	"""Return a 2-d array containing a noise-free double gaussian
	and whose values are truncated at 2**16-1
	
	Inputs:
	- arrShape	desired array shape (2 integers)
	- ctr		desired center (2 floats); uses the iraf/ds9 convention:
				the center of the 0,0 pixel is (0.0, 0.0)
	- sigma		desired sigma (float)
	- ampl		desired amplitude (float)
	"""
	if len(arrShape) != 2:
		raise ValueError("arrShape=%r must have 2 elements" % (arrShape,))
	if len(ctr) != 2:
		raise ValueError("ctr=%r must have 2 elements" % (ctr,))
	sigma = float(sigma)
	ampl = float(ampl)

	def peakFunc(i, j):
		radSq = (i - ctr[0])**2 + (j - ctr[1])**2
		expArg = - radSq / (2.0 * sigma**2)
		gauss = ampl * (num.exp(expArg)  + 0.1*num.exp(0.25*expArg))
		gauss = num.where(gauss <= _MaxValUInt16, gauss, _MaxValUInt16)
		return gauss.astype(num.UInt16)
	return num.fromfunction(peakFunc, arrShape)

def addNoise(
	data,
	sky = 1000,
	readNoise = 13,
	ccdGain = 5,
	bias = 1000,
):
	"""Adds poisson noise and gaussian read noise to the specified data.
	
	Inputs:
	- data		noiseless image, in ADU
	- sky		sky level, in ADU
	- readNoise	ccd read noise, in e-
	- ccdGain	ccd inverse gain, in e-/ADU
	- bias		image bias, in ADU
	"""
	outData = num.add(data, sky).astype(num.Int32)
	outData = rand.poisson(mean = outData * ccdGain) / ccdGain
	outData += rand.normal(mean = bias, std = readNoise/float(ccdGain), shape = data.shape)
	# truncate data and return as UInt16
	outData = num.where(outData >= 0, outData, 0)
	outData = num.where(outData <= _MaxValUInt16, outData, _MaxValUInt16)
	return outData.astype(num.UInt16)

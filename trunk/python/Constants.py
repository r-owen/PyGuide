#!/usr/local/bin/python
"""Constants and parameters

History:
2004-08-03 ROwen
2004-12-01 ROwen	Added NaN.
2005-02-07 ROwen	Added PosMinusIndex.
2005-05-18 ROwen	Added CCDInfo and DefThresh.
"""
import math

# Parameters

# PosMinusIndex defines the position convention.
# It is the position of the center of the pixel with index 0.
# 0.5 (default) means the center of the 0,0 pixel is (0.5, 0.5).
#     It also means the borders of an NxM image are (0.0, 0.0) and (N, M)
#     and the center is (N/2.0, M/2.0).
# 1.0 is the convention used by IRAF and DS9. It means the center
#     of the 0,0 pixel is (1.0, 1.0) and the borders of an NxM image
#     are (0.5, 0.5) and (N+0.5, M+0.5).
PosMinusIndex = 0.5

# default and minimum allowed threshold, where:
# valid data = stdDev * threshold + median
DefThresh = 3.0
MinThresh = 1.5

# title for ds9 diagnostic windows
DS9Title = "PyGuide"

# Constants (do not touch)
FWHMPerSigma = 2 * math.sqrt(2 * math.log(2))
NaN = float("nan")

# Misc
class CCDInfo:
	"""Info about the CCD
	
	- bias		ccd bias (ADU)
	- readNoise	ccd read noise (e-)
	- ccdGain	ccd inverse gain (e-/ADU)
	- satLevel	saturation level (ADU); data >= satLevel is saturated; None means unknown
	"""
	def __init__(self,
		bias,
		readNoise,
		ccdGain,
		satLevel = (2**16)-1,
	):
		self.bias = bias
		self.readNoise = readNoise
		self.ccdGain = ccdGain
		self.satLevel = satLevel
	
	def __repr__(self):
		dataList = []
		for arg in ("bias", "readNoise", "ccdGain", "satLevel"):
			val = getattr(self, arg)
			if val not in (None, ""):
				dataList.append("%s=%s" % (arg, val))
		return "CCDInfo(%s)" % ", ".join(dataList)

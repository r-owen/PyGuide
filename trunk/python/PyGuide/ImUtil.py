"""Image processing utilities.

History:
2004-04-05 ROwen	First release.
2004-04-16 ROwen	Added getQuartile, skyStats.
2004-04-21 ROwen	Added verbosity arg to skyStats and removed _SkyStatsDebug.
2004-08-04 ROwen	Removed documentation for nonexistent verbosity argument for getQuartile.
"""
import numarray as num

_QuartileResidRatios = (
	(1.0, 0.0),
	(1.0, 1.0/3.0),
	(1.0, 1.0),
	(1.0/3.0, 1.0),
)
def getQuartile(sortedData, qnum):
	"""Returns a quartile.
	Inputs:
	- sortedData	data sorted in ascending order
	- qnum	which quartile? 1=Q1, 2=Q2=median, 3=Q3
	
	If the position of the quartile is between two values,
	uses linear interpolation to compute the value.

	If the input data is not sorted, returns a meaningless number.
	If qnum not in (1, 2, 3) raises ValueError.
	"""
	if qnum not in (1, 2, 3):
		raise ValueError("qnum=%r must be 1, 2 or 3" % qnum)
	dataLen = len(sortedData)
	ratios = _QuartileResidRatios[((dataLen-1) * qnum) % 4]
	ind0 = (dataLen-1) * qnum // 4
	return ((sortedData[ind0] * ratios[0]) + (sortedData[ind0+1] * ratios[1])) / (ratios[0] + ratios[1])

def skyStats(maskedData, verbosity=1):
	"""Computes sky statistics.
	
	Inputs:
	- maskedData: a numarray.ma masked array.
	- verbosity	0: no output, 1: print warnings, 2: print information

	Returns median and standard deviation.
	
	The statistics are computed using quartiles. Pixel values above
	2.35 * stdDev are ignored, and this computation is iterated a few times
	to refine the cutoff point.
	
	Standard deviation is computed as stdDev = 0.741 * (Q3 - Q1)
	"""
	# creating sorted data
	sortedData = maskedData.compressed().raw_data()
	dataLen = len(sortedData)
	if verbosity >= 2:
		print "skyStats sorting %d elements" % (dataLen)
	sortedData = num.sort(sortedData)
	
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
		cutInd = num.searchsorted(sortedData, [cutVal])[0]
		if verbosity >= 2:
			print "skStats cutInd=%d, sortedData[cutInd]=%d" % (cutInd, sortedData[cutInd])
		dataLen = cutInd
	
	return med, stdDev

def subFrameCtr(data, ctr, size):
	"""Extract a subframe from a 2d array given a center and size.
	
	Return a pointer (not a copy).

	Inputs:
	- data	2-d array of data [i,j]
	- ctr	desired i,j center of subframe (may be float)
	- size	desired i,j size of subframe (may be float);
			e.g. (5,7) returns a 5x7 subframe
	
	Returns two numarray arrays:
	- offset	i,j index of start of subframe (int);
				location subloc in the subframe
				matches location subloc + offset in data
	- subframe	the subframe as a pointer into data (NOT a copy)
	
	If the requested subframe extends outside the boundaries of data,
	the subframe is truncated without warning. You can tell if it
	has been truncated by comparing its size to the requested size.
	
	If size is odd and the entire subframe fits in data without truncation,
	then ctr is truly centered.
	"""
	rad = num.array(size) / 2.0
	begInd = num.maximum(ctr - rad, 0.0).astype(num.Int)
	endInd = (ctr + rad).astype(num.Int)
	
	subframe = data[begInd[0]:endInd[0], begInd[1]:endInd[1]]
	return begInd, subframe

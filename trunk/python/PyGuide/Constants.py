#!/usr/local/bin/python
"""Constants and parameters

History:
2004-08-03 ROwen
2004-12-01 ROwen	Added NaN.
2005-02-07 ROwen	Added PosMinusIndex.
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

# Constants (do not touch)
FWHMPerSigma = 2 * math.sqrt(2 * math.log(2))
NaN = float("nan")

DS9Title = "PyGuide"
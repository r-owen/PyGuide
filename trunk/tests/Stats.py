import math

class Stats:
    def __init__(self):
        self.clear()
    
    def append(self, val):
        """Append a value to the data set.
        """
        val = float(val)
        self._n += 1
        self._sum += val
        self._sumSq += val*val
        if self._minVal == None or val < self._minVal:
            self._minVal = val
        if self._maxVal == None or val > self._maxVal:
            self._maxVal = val
    
    def clear(self):
        """Clear the data.
        """
        self._n = 0.0
        self._sum = 0.0
        self._sumSq = 0.0
        self._minVal = None
        self._maxVal = None
    
    def min(self):
        """Return the minimum value.
        """
        self._checkPoints()
        return self._minVal
        
    def max(self):
        """Return the maximum value.
        """
        self._checkPoints()
        return self._maxVal
        
    def mean(self):
        """Return the mean value.
        """
        self._checkPoints()
        return self._sum / float(self._n)
    
    def nPoints(self):
        """Return the number of points.
        """
        return self._n
    
    def var(self):
        """Return the variance:
        (sum(val^2) - mean^2) / N-1
        """
        self._checkPoints(2)
        return (self._sumSq - (self.mean()**2)) / float(self._n - 1)

    def stdDev(self):
        """Return the stanard deviation: sqrt(var).
        """
        return math.sqrt(self.var())
        
    def rms(self):
        """Return the RMS value.
        """
        self._checkPoints()
        return math.sqrt(self._sumSq / float(self._n))
    
    def _checkPoints(self, minReq=1):
        """Check that we have enough data points.
        Raise a RuntimeError if we do not.
        """
        if self._n < minReq:
            if minReq <= 1:
                raise RuntimeError("No data points")
            raise RuntimeError(
                "Need at least %d data points, but only have %d" % (minReq, self._n)
            )

if __name__ == "__main__":
    s = Stats()
    NaN = float("NaN")
    
    def printRes():
        try:
            minVal = s.min()
        except RuntimeError:
            minVal = NaN

        try:
            maxVal = s.max()
        except RuntimeError:
            maxVal = NaN

        try:
            mean = s.mean()
        except RuntimeError:
            mean = NaN

        try:
            stdDev = s.stdDev()
        except RuntimeError:
            stdDev = NaN

        try:
            var = s.var()
        except RuntimeError:
            var = NaN

        try:
            rms = s.rms()
        except RuntimeError:
            rms = NaN
        
        print ("%2d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f" % \
            (s.nPoints(), minVal, maxVal, mean, stdDev, var, rms))
    
    print (" N     min     max    mean  stdDev     var     RMS")
    testVals = (1, 2, 3, 4, 5, 6, 7, 8, 9)
    printRes()
    for val in testVals:
        s.append(val)
        printRes()

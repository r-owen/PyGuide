#!/usr/bin/env python
"""Check PyGuide using pychecker

History:
2004-08-06 ROwen
2007-01-23 ROwen    Changed #!/usr/local/bin/python to #!/usr/bin/env python
"""
# import modules that I don't want checked
import numpy
import scipy.ndimage

# start checking
import pychecker.checker
import PyGuide

print "Done"

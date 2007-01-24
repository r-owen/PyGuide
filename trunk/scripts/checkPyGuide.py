#!/usr/bin/env python
"""Check PyGuide using pychecker

History:
2004-08-06 ROwen
2007-01-23 ROwen    Changed #!/usr/local/bin/python to #!/usr/bin/env python
"""
# import modules that I don't want checked
import numarray
import numarray.nd_image
import numarray.random_array
import numarray.ma

# start checking
import pychecker.checker
import PyGuide

print "Done"

#!/usr/local/bin/python
"""Check PyGuide using pychecker

History:
2004-08-06 ROwen
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

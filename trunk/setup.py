#!/usr/bin/env python
from distutils.core import setup, Extension
from numarray.numarrayext import NumarrayExtension
import sys
import os

if not hasattr(sys, 'version_info') or sys.version_info[0:2] < (2,2):
	raise SystemExit("Python 2.2 or later required to build this module.")
    
PkgName = "PyGuide"
SrcDir = "src"
PyDir = "python"

radProfExt = NumarrayExtension(
	"radProf",
	sources = [os.path.join(SrcDir, "RadProfModule.c")],
)

setup(
	name = PkgName,
	version = "2.0",
	description = "support for telescope guiding",
	author = "Russell Owen",
#	url = "http://astro.washington.edu/owen/",
	ext_package = PkgName,
	ext_modules = [radProfExt],
	package_dir = {'PyGuide': PyDir},
	packages = [PkgName],
)

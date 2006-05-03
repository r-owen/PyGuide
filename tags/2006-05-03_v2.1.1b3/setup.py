#!/usr/local/bin/python
from distutils.core import setup, Extension
import distutils.sysconfig
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
	description = "Find stars for telescope guiding",
	author = "Russell Owen",
	url = "http://www.astro.washington.edu/rowen/",
	ext_package = PkgName,
	ext_modules = [radProfExt],
	package_dir = {'PyGuide': PyDir},
	packages = [PkgName],
	scripts = ["doPyGuide.py"],
)

# create ups file for Fermi/Princeton packaging system
sitePkgDir = distutils.sysconfig.get_python_lib()
if not os.path.exists("ups"):
	os.mkdir("ups")
upsfile = file("ups/PyGuide.table", "w")
try:
	upsfile.write("""File=Table
Product=PyGuide
Group:
Flavor=ANY
Common:
  Action=setup
    proddir()
    setupenv()
    pathAppend(PATH, ${UPS_PROD_DIR}/bin)
    pathAppend(PYTHONPATH, ${UPS_PROD_DIR}%s)
End:
""" % (sitePkgDir,))
finally:
	upsfile.close()
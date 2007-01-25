#!/usr/bin/env python
from distutils.core import setup, Extension
import distutils.sysconfig
from numarray.numarrayext import NumarrayExtension
import sys
import os
sys.path.append("python")
import Version

if not hasattr(sys, 'version_info') or sys.version_info[0:2] < (2,2):
    raise SystemExit("Python 2.2 or later required to build this module.")

PkgName = "PyGuide"
SrcDir = "src"
PyDir = "python"
UPSArg = "--ups"

dataFiles = []
if UPSArg in sys.argv and "install" in sys.argv:
    # create data file for the Fermi/Princeton UPS runtime loader
    upsFileName = "ups/%s.table" % (PkgName,)
    sys.argv.pop(sys.argv.index(UPSArg))
    sitePkgDir = distutils.sysconfig.get_python_lib(prefix="")
    if not os.path.exists("ups"):
        os.mkdir("ups")
    upsFile = file(upsFileName, "w")
    try:
        upsFile.write("""File=Table
    Product=%s
    Group:
    Flavor=ANY
    Common:
    Action=setup
    proddir()
    setupenv()
    pathAppend(PATH, ${UPS_PROD_DIR}/bin)
    pathAppend(PYTHONPATH, ${UPS_PROD_DIR}/%s)
    End:
    """ % (PkgName, sitePkgDir,))
    finally:
        upsFile.close()
    dataFiles.append(["ups", [upsFileName]])

radProfExt = NumarrayExtension(
    "radProf",
    sources = [os.path.join(SrcDir, "RadProfModule.c")],
)

setup(
    name = PkgName,
    version = Version.__version__,
    description = "Find stars for telescope guiding",
    author = "Russell Owen",
    url = "http://www.astro.washington.edu/rowen/",
    ext_package = PkgName,
    ext_modules = [radProfExt],
    package_dir = {'PyGuide': PyDir},
    packages = [PkgName],
    data_files = dataFiles,
    scripts = ["scripts/doPyGuide.py", "scripts/checkPyGuide.py"],
)

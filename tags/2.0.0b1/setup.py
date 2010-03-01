#!/usr/bin/env python
import distutils.sysconfig
import numpy.distutils.core
import numpy.distutils.misc_util
import sys
import os
rootDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(rootDir, "python"))
import Version

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

def configuration(parent_package = '', top_path = None):
    config = numpy.distutils.misc_util.Configuration(
        name = PkgName,
        parent_package = parent_package,
        top_path = top_path,
        version = Version.__version__,
        description = "Find stars for telescope guiding",
        author = "Russell Owen",
        url = "http://www.astro.washington.edu/rowen/",
        package_dir = {'PyGuide': PyDir},
        packages = [PkgName],
        data_files = dataFiles,
        scripts = ["scripts/doPyGuide.py"],
    )
    config.add_extension(
        "radProf",
        sources = [os.path.join(SrcDir, "RadProfModule.c")],
    )
    return config

numpy.distutils.core.setup(configuration=configuration)

#!/usr/bin/env python

from distutils.core import setup, Extension
import sys
import os

import numpy as np

rootDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(rootDir, "python"))
import Version

PkgName = "PyGuide"
SrcDir = "src"
PyDir = "python"

dataFiles = []

radProfExt = Extension(
    "PyGuide.radProf",
    ["src/RadProfModule.c"],
    include_dirs = ["src", np.get_include()],
)

setup(
    name = PkgName,
    version = Version.__version__,
    description = "Find stars for telescope guiding",
    author = "Russell Owen",
    url = "http://www.astro.washington.edu/rowen/",
    package_dir = {'PyGuide': PyDir},
    packages = [PkgName],
    ext_modules = [radProfExt],
    data_files = dataFiles,
    scripts = ["scripts/doPyGuide.py"],
)

#!/usr/bin/env python
from distutils.core import setup, Extension
import os

PkgName = "PyGuide"
SrcDir = "src"
PyDir = "python"

radProfExt = Extension(
	"radProf",
	sources = [os.path.join(SrcDir, "radProfModule.c")],
)

setup(name = PkgName,
	version = "1.0",
	description = "support for telescope guiding",
	author = "Russell Owen",
#	url = "http://astro.washington.edu/owen/",
	ext_package = PkgName,
	ext_modules = [radProfExt],
	package_dir = {'': PyDir},
	packages = [PkgName],
)

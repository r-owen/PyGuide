#!/usr/bin/env python
"""A script to release a new version of PyGuide and upload it to PyPI

To use:
    ./releaseNewVersion.py
"""
from __future__ import print_function
from __future__ import unicode_literals

from __future__ import with_statement
import os
import re
import sys
import subprocess

UploadToPyPI = True    # wait until package uses numpy before uploading to PyPI

PkgName = "PyGuide"
PythonDir = "python"
sys.path.insert(0, PythonDir)
import Version  # noqa  must come after modifying sys.path

queryStr = "Version from %s.Version = %s; is this OK? (y/[n]) " % (PkgName, Version.__version__,)
versOK = input(queryStr)
if not versOK.lower() == "y":
    sys.exit(0)

versRegEx = re.compile(r"<h3>(\d.*?)\s+\d\d\d\d-\d\d-\d\d</h3>")
with file(os.path.join("docs", "VersionHistory.html")) as vhist:
    for line in vhist:
        versMatch = versRegEx.match(line)
        if versMatch:
            histVersStr = versMatch.groups()[0]
            if histVersStr == Version.__version__:
                print("Version in VersionHistory.html matches")
                break
            else:
                print("Error: version in VersionHistory.html = %s != %s" % (histVersStr, Version.__version__))
                sys.exit(0)

print("Status of git repository:")

subprocess.call(["git", "status"])

versOK = input("Is the git repository up to date? (y/[n]) ")
if not versOK.lower() == "y":
    sys.exit(0)

print("Git repository OK")

if UploadToPyPI:
    print("Uploading to PyPI")
    status = subprocess.call(["python", "setup.py", "sdist", "upload"])
    if status != 0:
        print("Build and upload failed!")

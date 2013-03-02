#!/usr/bin/env python
"""A script to release a new version of PyGuide and upload it to PyPI

To use:
    ./releaseNewVersion.py
"""
from __future__ import with_statement
import os
import re
import shutil
import sys
import subprocess

UploadToPyPI = False    # wait until package uses numpy before uploading to PyPI

PkgName = "PyGuide"
PythonDir = "python"
sys.path.insert(0, PythonDir)
import Version
queryStr = "Version from %s.Version = %s; is this OK? (y/[n]) " % (PkgName, Version.__version__,)
versOK = raw_input(queryStr)
if not versOK.lower() == "y":
    sys.exit(0)

versRegEx = re.compile(r"<h3>(\d.*?)\s+\d\d\d\d-\d\d-\d\d</h3>")
with file(os.path.join("docs", "VersionHistory.html")) as vhist:
    for line in vhist:
        versMatch = versRegEx.match(line)
        if versMatch:
            histVersStr = versMatch.groups()[0]
            if histVersStr == Version.__version__:
                print "Version in VersionHistory.html matches"
                break
            else:
                print "Error: version in VersionHistory.html = %s != %s" % (histVersStr, Version.__version__)
                sys.exit(0)

print "Status of subversion repository:"

subprocess.call(["svn", "status"])

versOK = raw_input("Is the subversion repository up to date? (y/[n]) ")
if not versOK.lower() == "y":
    sys.exit(0)

print "Subversion repository OK"

exportRoot = os.environ["HOME"]
exportFileName = "%s_%s" % (PkgName, Version.__version__)
exportPath = os.path.abspath(os.path.join(exportRoot, exportFileName))
zipFileName = "%s.zip" % (exportFileName,)
zipFilePath = os.path.abspath(os.path.join(exportRoot, zipFileName))
if os.path.exists(exportPath):
    print "Export directory %r already exists" % (exportPath)
    versOK = raw_input("Should I delete the old %r? (yes/[n]) " % (exportPath,))
    if not versOK.lower() == "yes":
        sys.exit(0)
    print "Deleting %r" % (exportPath,)
    shutil.rmtree(exportPath)
if os.path.exists(zipFilePath):
    getOK = raw_input("File %r already exists! Should I delete it? (yes/[n]) " % (zipFilePath,))
    if not getOK.lower() == "yes":
        sys.exit(0)
    print "Deleting %r" % (zipFilePath,)
    os.remove(zipFilePath)

print "Exporting subversion repository to %r" % (exportPath,)

status = subprocess.call(["svn", "export", ".", exportPath])
if status != 0:
    print "Svn export failed!"
    sys.exit(0)

print "Zipping %r" % (exportPath,)
status = subprocess.call(["zip", "-r", "-q", zipFileName, exportFileName], cwd=exportRoot)
if status != 0:
    print "Zip failed!"
else:
    print "Unix package zipped"
    status = subprocess.call(["open", exportRoot])
    
if UploadToPyPI:
    print "Uploading to PyPI"
    status = subprocess.call(["python", "setup.py", "sdist", "upload"], cwd=exportPath)
    if status != 0:
        print "Build and upload failed!"

delOK = raw_input("OK to delete %r? (y/[n]) " % (exportPath,))
if not delOK.lower() == "y":
    sys.exit(0)

print "Deleting %r" % (exportPath,)
shutil.rmtree(exportPath)

This is the pure python code for PyGuide.

It is placed down one level so that after you use setup.py to install
the software you can run python from the same directory to test the code.
If the python package(s) were at the main level, then when you imported them
you would be importing the local copy (which does not include the compiled
C/C++ extensions and so would fail) instead of the installed copy.

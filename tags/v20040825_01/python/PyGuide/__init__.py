"""Python code for guiding telescopes.

The main functions are:
- centroid measures the centroid of a star.
- findStars finds stars.

This code is written to handle stellar images with portions missing, such as
one might find in a spectrograph slit viewer or a coherent fiber bundle guide
probe. It uses star finding and centroiding algorithms that do not assume a
full gaussian, yet the centroider does a good job with what information it has,
trying to find the true center of the star, not just the center of mass of what
of it is visible.

WARNING: Python handles images as data[y,x]. This same index order is used
for positions, sizes and such. In an attempt to reduce confusion, the code
uses i,j instead of y,x.

Pixel convention: The point 0,0 is at the corner of first pixel read out.
Hence the center of that pixel is (0.5, 0.5) and the center of a 1024x1024 CCD
is (512.0, 512.0)

Acknowledgements:
- Connie Rockosi supplied the star finding algorithm. She got it from
  Jim Gunn and Robert Lupton.
- Jim Gunn supplied the centroiding algorithm. Connie Rockosi explained
  some of the subtleties.
"""
from Constants import *
from Centroid import *
from FindStars import *
from StarShape import *
import FakeData

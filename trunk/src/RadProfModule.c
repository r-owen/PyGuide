// include RadProfModule.h *first* because it includes Python.h,
// which must be done before importing standard C libraries.
#include "RadProfModule.h"

#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>


/*
Python extension module to extract radial profiles
and measure radial symmetry. Useful for centroiding.

History:
2004-03-22 ROwen    First release.
2004-04-02 ROwen    Fix memory leak: did not DECREF maskArry.
2004-04-07 ROwen    Flipped meaning of mask to match numarray.ma.
                    Modified data array types for efficiency.
2004-04-12 ROwen    Modified routines to return totCounts.
2004-04-16 ROwen    Changed x,y notation to i,j since i=y, j=x in common cases.
2004-04-30 ROwen    Modified to use UInt16, not Int16
2004-05-13 ROwen    Added radIndByRadSq and radSqByRadInd
2004-07-06 ROwen    Minor cleanups in radAsymm.
2004-08-06 ROwen    Added radAsymmWeighted.
                    Added doc string warning: positional args only
2004-10-13 ROwen    Modified to include RadProfModule.c first,
                    and so use its include of Python.h and numarray.h.
2005-03-31 ROwen    Fixed doc string for Py_radAsymmWeighted_doc
                    and tried to clarify the documentation for ctr.
2005-10-14 ROwen    Modified to use Float32 image data, not UInt16.
2006-07-11 ROwen    Use PyMODINIT_FUNC as the type of return for the init function
                    instead of void (apparently recc. for python 2.3 and later).
*/

// global working arrays for radProf
static long *g_radProf_radIndByRadSq;
static long g_radProf_nElt = 0;

// global working arrays for radAsymm
static Float64 *g_radAsymm_mean;
static Float64 *g_radAsymm_var;
static Int32 *g_radAsymm_nPts;
static long g_radAsymm_nElt = 0;

char radProfModule_doc [] =
"Code to obtain radial profiles of 2-d arrays\n"
"\n"
"Warning: these routines only take positional arguments, not named arguments.\n"
;

// note: MAX and MIN are defined in nummacro.h, imported by libnumarray.h

/* Py_radAsymm ============================================================
*/
char Py_radAsymm_doc [] =
"Compute a measure of radial asymmetry:\n"
"sum over rad of var(rad) * nPts(rad).\n"
"\n"
"Input (by position only):\n"
"- data         a 2-d array [i,j] (Float32)\n"
"- mask         mask array [i,j] (Bool); True for values to mask out (ignore).\n"
"               None if no mask array.\n"
"- ijCtr        i,j center of scan (int)\n"
"- rad          radius of scan (int)\n"
"Returns:\n"
"- asymm        radial asymmetry (see above)\n"
"- totCounts    the total # of counts (float to avoid overflow)\n"
"- totPts       the total # of points (int)\n"
"\n"
"Points off the data array are ignored.\n"
"Thus the center need not be on the array.\n"
"\n"
"If mask is not None then it must have the same shape as data,\n"
"else raises ValueError;\n"
"\n"
"The code is more efficient if the arrays have the suggested type\n"
"and are continuous.\n"
;
static PyObject *Py_radAsymm(PyObject *self, PyObject *args) {
    PyObject *dataObj, *maskObj;
    PyArrayObject *dataArry, *maskArry;
    long iCtr, jCtr, rad, totPts;
    Float64 asymm, totCounts;
    char ModName[] = "radAsymm";

    if (!PyArg_ParseTuple(args, "OO(ll)l",
            &dataObj, &maskObj, &iCtr, &jCtr, &rad))
        return NULL;
    
    // Convert arrays to well-behaved arrays of correct type and verify
    // These arrays MUST be decrefed before return.
    dataArry = NA_InputArray(dataObj, tFloat32, NUM_C_ARRAY);
    if (maskObj == Py_None) {
        maskArry = NULL;
    } else {
        maskArry = NA_InputArray(maskObj, tBool, NUM_C_ARRAY);
    }

    // Check the input arrays
    if (!dataArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting data input", ModName);
        goto errorExit;
    }
    if (!maskArry && maskObj != Py_None) {
        PyErr_Format(PyExc_ValueError, "%s: error converting mask input", ModName);
        goto errorExit;
    }
    if (dataArry->nd != 2) {
        PyErr_Format(PyExc_ValueError, "%s: data must be 2-dimensional", ModName);
        goto errorExit;
    }
    if (maskArry && !NA_ShapeEqual(dataArry, maskArry)) {
        PyErr_Format(PyExc_ValueError, "%s: mask must be the same shape as data", ModName);
        goto errorExit;
    }
    
    // Call the C code
    totPts = radAsymm(
        dataArry->dimensions[0], dataArry->dimensions[1],
        NA_OFFSETDATA(dataArry),
        maskArry? NA_OFFSETDATA(maskArry): NULL,
        iCtr, jCtr,
        rad,
        &asymm,
        &totCounts
    );
    if (totPts < 0) {
        PyErr_Format(PyExc_ValueError, "radAsymm failed");
        goto errorExit;
    }

    // Done with all arrays, decref them
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);

    return Py_BuildValue("ddl", asymm, totCounts, totPts);

errorExit:
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);
    return NULL;
}

/* Py_radAsymmWeighted ========================================================
*/
char Py_radAsymmWeighted_doc [] =
"Compute a weighted measure of radial asymmetry:\n"
"  sum over rad of var(rad)^2 / weight(rad)\n"
"where weight is the expected sigma of var(rad) due to pixel noise:\n"
"  weight(rad) = pixNoise(rad) * sqrt(2(numPix(rad) - 1))/numPix(rad)\n"
"  pixNoise(rad) = sqrt((readNoise/ccdGain)^2 + (meanVal(rad)-bias)/ccdGain)\n"
"\n"
"Inputs (by position only):\n"
"- data              data array [i,j]\n"
"- mask              mask array [i,j] (NULL if none);\n"
"                    0 for valid values, 1 for values to ignore\n"
"- ijCtr             i,j center of scan\n"
"- rad               radius of scan\n"
"- readNoise         read noise in e-\n"
"- ccdGain           ccd inverse gain in e-/ADU\n"
"- bias              ccd bias in ADU\n"
"\n"
"Returns:\n"
"- asymm        radial asymmetry (see above)\n"
"- totCounts    the total # of counts (float to avoid overflow)\n"
"- totPts       the total # of points (int)\n"
"\n"
"Points off the data array are ignored.\n"
"Thus the center need not be on the array.\n"
"\n"
"If mask is not None then it must have the same shape as data,\n"
"else raises ValueError;\n"
"\n"
"The code is more efficient if the arrays have the suggested type\n"
"and are continuous.\n"
;
static PyObject *Py_radAsymmWeighted(PyObject *self, PyObject *args) {
    PyObject *dataObj, *maskObj;
    PyArrayObject *dataArry, *maskArry;
    long iCtr, jCtr, rad, bias, totPts;
    Float64 readNoise, ccdGain, asymm, totCounts;
    char ModName[] = "radAsymm";

    if (!PyArg_ParseTuple(args, "OO(ll)lldd",
            &dataObj, &maskObj, &iCtr, &jCtr, &rad, &bias, &readNoise, &ccdGain))
        return NULL;
    
    // Convert arrays to well-behaved arrays of correct type and verify
    // These arrays MUST be decrefed before return.
    dataArry = NA_InputArray(dataObj, tFloat32, NUM_C_ARRAY);
    if (maskObj == Py_None) {
        maskArry = NULL;
    } else {
        maskArry = NA_InputArray(maskObj, tBool, NUM_C_ARRAY);
    }

    // Check the input arrays
    if (!dataArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting data input", ModName);
        goto errorExit;
    }
    if (!maskArry && maskObj != Py_None) {
        PyErr_Format(PyExc_ValueError, "%s: error converting mask input", ModName);
        goto errorExit;
    }
    if (dataArry->nd != 2) {
        PyErr_Format(PyExc_ValueError, "%s: data must be 2-dimensional", ModName);
        goto errorExit;
    }
    if (maskArry && !NA_ShapeEqual(dataArry, maskArry)) {
        PyErr_Format(PyExc_ValueError, "%s: mask must be the same shape as data", ModName);
        goto errorExit;
    }
    
    // Call the C code
    totPts = radAsymmWeighted(
        dataArry->dimensions[0], dataArry->dimensions[1],
        NA_OFFSETDATA(dataArry),
        maskArry? NA_OFFSETDATA(maskArry): NULL,
        iCtr, jCtr,
        rad,
        bias,
        readNoise,
        ccdGain,
        &asymm,
        &totCounts
    );
    if (totPts < 0) {
        PyErr_Format(PyExc_ValueError, "radAsymm failed");
        goto errorExit;
    }

    // Done with all arrays, decref them
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);

    return Py_BuildValue("ddl", asymm, totCounts, totPts);

errorExit:
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);
    return NULL;
}


/* Py_radProf ============================================================
*/
char Py_radProf_doc [] =
"Generate a radial profile as a function of radial index\n"
"(an approximation of radius; see below for details)\n"
"\n"
"Inputs (by position only):\n"
"- data         a 2-d array [i,j] (Float32)\n"
"- mask         mask array [i,j] (Bool); True for values to mask out (ignore).\n"
"               None if no mask array.\n"
"- ijCtr        i,j center of profile (int)\n"
"- rad          desired radius of profile (int)\n"
"\n"
"Outputs (by position only):\n"
"- mean         the mean at each radial index; 0 if npts=0 (Float64)\n"
"- var          the variance (stdDev^2) at each radial index; 0 if npts=0 (Float64)\n"
"- nPts         the # of points at each radial index (Int32)\n"
"\n"
"Returns:\n"
"- totCounts    the total # of counts (sum of mean*nPts); float to avoid overflow\n"
"- totPts       the total # of points (sum of nPts)\n"
"\n"
"Radial Index:\n"
"radProf uses the Mirage convention for radial profiles;\n"
"this means that profiles are a function of radial index,\n"
"an approximation to radius that handles the central pixels better.\n"
"By definition:\n"
"radial index[rad**2] = 0, 1, 2, int(sqrt(rad**2)+1.5) for rad**2>2\n"
"As a result, radial index[rad] is rad+1 for rad>1\n"
"Thus your output arrays must have at least rad+2 points\n"
"\n"
"(See also radIndByRadSq and radSqByRadInd, offering mappings\n"
"between radial index and radius squared.)\n"
"\n"
"Points off the data array are ignored.\n"
"Thus the center need not be on the array.\n"
"\n"
"If mask is not None then it must have the same shape as data,\n"
"else raises ValueError;\n"
"\n"
"Each output array must have the same length,\n"
"and that length must be at least rad + 2,\n"
"else raises ValueError\n"
"\n"
"All arguments are coerced to the correct data type,\n"
"but the code is more efficient if the arrays have the suggested type.\n"
;
static PyObject *Py_radProf(PyObject *self, PyObject *args) {
    PyObject *dataObj, *maskObj, *meanObj, *varObj, *nPtsObj;
    PyArrayObject *dataArry, *maskArry, *meanArry, *varArry, *nPtsArry;
    long iCtr, jCtr, rad, outLen, totPts;
    Float64 totCounts;
    char ModName[] = "radProf";
    
    if (!PyArg_ParseTuple(args, "OO(ll)lOOO",
            &dataObj, &maskObj, &iCtr, &jCtr, &rad,
            &meanObj, &varObj, &nPtsObj))
        return NULL;
    
    // Convert arrays to well-behaved arrays of correct type and verify
    // These arrays MUST be decrefed before return.
    dataArry = NA_InputArray(dataObj, tFloat32, NUM_C_ARRAY);
    if (maskObj == Py_None) {
        maskArry = NULL;
    } else {
        maskArry = NA_InputArray(maskObj, tBool, NUM_C_ARRAY);
    }
    meanArry = NA_OutputArray(meanObj, tFloat64, NUM_C_ARRAY);
    varArry = NA_OutputArray(varObj, tFloat64, NUM_C_ARRAY);
    nPtsArry = NA_OutputArray(nPtsObj, tInt32, NUM_C_ARRAY);

    // Check the input arrays
    if (!dataArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting data input", ModName);
        goto errorExit;
    }
    if (!maskArry && maskObj != Py_None) {
        PyErr_Format(PyExc_ValueError, "%s: error converting mask input", ModName);
        goto errorExit;
    }
    if (dataArry->nd != 2) {
        PyErr_Format(PyExc_ValueError, "%s: data must be 2-dimensional", ModName);
        goto errorExit;
    }
    if (maskArry && !NA_ShapeEqual(dataArry, maskArry)) {
        PyErr_Format(PyExc_ValueError, "%s: mask must be the same shape as data", ModName);
        goto errorExit;
    }
    
    // Check output arrays and compute outLen
    if (!meanArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting mean input", ModName);
        goto errorExit;
    }
    if (!varArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting variance input", ModName);
        goto errorExit;
    }
    if (!nPtsArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting nPts input", ModName);
        goto errorExit;
    }
    if (meanArry->nd != 1) {
        PyErr_Format(PyExc_ValueError, "%s: mean must be 1-dimensional", ModName);
        goto errorExit;
    }
    if (varArry->nd != 1) {
        PyErr_Format(PyExc_ValueError, "%s: variance must be 1-dimensional", ModName);
        goto errorExit;
    }
    if (nPtsArry->nd != 1) {
        PyErr_Format(PyExc_ValueError, "%s: nPts must be 1-dimensional", ModName);
        goto errorExit;
    }
    outLen = meanArry->dimensions[0];
    if (outLen != varArry->dimensions[0]) {
        PyErr_Format(PyExc_ValueError, "%s: var array length != mean array length", ModName);
        goto errorExit;
    }
    if (outLen != nPtsArry->dimensions[0]) {
        PyErr_Format(PyExc_ValueError, "%s: nPts array length != mean array length", ModName);
        goto errorExit;
    }
    if (outLen < rad + 2) {
        PyErr_Format(PyExc_ValueError, "%s: output arrays are too short", ModName);
        goto errorExit;
    }
    
    // Call the C code
    totPts = radProf(
        dataArry->dimensions[0], dataArry->dimensions[1],
        NA_OFFSETDATA(dataArry),
        maskArry? NA_OFFSETDATA(maskArry): NULL,
        iCtr, jCtr,
        rad,
        outLen,
        NA_OFFSETDATA(meanArry),
        NA_OFFSETDATA(varArry),
        NA_OFFSETDATA(nPtsArry),
        &totCounts
    );
    if (totPts < 0) {
        PyErr_Format(PyExc_ValueError, "radProf failed");
        goto errorExit;
    }

    // Done with all arrays, decref them
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);
    Py_XDECREF(meanArry);
    Py_XDECREF(varArry);
    Py_XDECREF(nPtsArry);

    return Py_BuildValue("dl", totCounts, totPts);

errorExit:
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);
    Py_XDECREF(meanArry);
    Py_XDECREF(varArry);
    Py_XDECREF(nPtsArry);
    return NULL;
}


/* Py_radIndByRadSq ============================================================
*/
char Py_radIndByRadSq_doc [] =
"Return radial index, indexed by radius squared.\n"
"Radial index is explained under radProf.\n"
"\n"
"Inputs (by position only):\n"
"- nElt     the desired number of elements in the returned array\n"
"\n"
"Returns one array:\n"
"- radInd[nElt] radial index, indexed by radius squared (int)\n"
"\n"
"Raises ValueError if nElt < 0\n"
;
static PyObject *Py_radIndByRadSq(PyObject *self, PyObject *args) {
    long nElt;
    char ModName[] = "radIndByRadSq";

    if (!PyArg_ParseTuple(args, "l", &nElt))
        return NULL;
    
    if (nElt < 0) {
        PyErr_Format(PyExc_ValueError, "%s: nPts < 0", ModName);
        return NULL;
    }

    /* compute radSqByRadInd; use the archived global data since it is handy
    and will usually already be large enough to suit the user */
    if (g_radProf_setup(nElt) != 1) {
        PyErr_Format(PyExc_MemoryError, "%s: insufficient memory", ModName);
        return NULL;
    }
    
    return (PyObject *) NA_NewArray((void *)g_radProf_radIndByRadSq, tLong, 1, nElt);
}


/* Py_radSqByRadInd ============================================================
*/
char Py_radSqByRadInd_doc [] =
"Return radius squared, indexed by radial index.\n"
"Radial index is explained under radProf.\n"
"\n"
"Input (by position only):\n"
"- nElt         the desired number of elements in the output array\n"
"\n"
"Returns one array:\n"
"- radSq[nElt]  radius squared, indexed by radial index (int)\n"
"\n"
"Raises ValueError if nElt < 0\n"
;
static PyObject *Py_radSqByRadInd(PyObject *self, PyObject *args) {
    long nElt;
    Long radInd;
    char ModName[] = "radSqByRadInd";
    Long *radSqByRadInd;

    if (!PyArg_ParseTuple(args, "l", &nElt))
        return NULL;
    
    if (nElt < 0) {
        PyErr_Format(PyExc_ValueError, "%s: nPts < 0", ModName);
        return NULL;
    }
    
    // allocate max(3, nElt) elements (simplifies the computation
    // and avoids calling calloc on 0 elements)
    radSqByRadInd = calloc(MAX(nElt, 3), sizeof *radSqByRadInd);
    if (radSqByRadInd == NULL) {
        PyErr_Format(PyExc_MemoryError, "%s: insufficient memory", ModName);
        return NULL;
    }
    
    for (radInd=0; radInd<3; radInd++) {
        radSqByRadInd[radInd] = radInd;
    }
    for (radInd=3; radInd<nElt; radInd++) {
        radSqByRadInd[radInd] = (radInd - 1) * (radInd - 1);
    }
    
    return (PyObject *) NA_NewArray((void *)radSqByRadInd, tLong, 1, nElt);
}


/* Py_radSqProf ============================================================
*/
char Py_radSqProf_doc [] =
"Generate a radial profile as a function of radius squared\n"
"\n"
"Input (by position only):\n"
"- data         a 2-d array [i,j] (Float32)\n"
"- mask         mask array [i,j] (Bool); True for values to mask out (ignore).\n"
"               None if no mask array.\n"
"- ijCtr        i,j center of profile (int)\n"
"- rad          radius of profile (int)\n"
"Outputs (by position only):\n"
"- mean         the mean at each radius squared; 0 if nPts=0 (float)\n"
"- var          the variance (stdDev^2) at each radius squared; 0 if npts=0 (float)\n"
"- nPts         the # of points at each radius squared (Int32)\n"
"Returns\n"
"- totCounts    the total # of counts (float to avoid overflow)\n"
"- totPts       the total # of points (int)\n"
"\n"
"Points off the data array are ignored.\n"
"Thus the center need not be on the array.\n"
"\n"
"If mask is not None then it must have the same shape as data,\n"
"else raises ValueError;\n"
"\n"
"Each output array must have the same length,\n"
"and that length must be at least rad**2 + 1,\n"
"else raises ValueError;\n"
"\n"
"The code is more efficient if the arrays have the suggested type\n"
"and are continuous.\n"
;
static PyObject *Py_radSqProf(PyObject *self, PyObject *args) {
    PyObject *dataObj, *maskObj, *meanObj, *varObj, *nPtsObj;
    PyArrayObject *dataArry, *maskArry, *meanArry, *varArry, *nPtsArry;
    long iCtr, jCtr, rad, radSq, outLen, totPts;
    Float64 totCounts;
    char ModName[] = "radSqProf";
    
    if (!PyArg_ParseTuple(args, "OO(ll)lOOO",
            &dataObj, &maskObj, &iCtr, &jCtr, &rad,
            &meanObj, &varObj, &nPtsObj))
        return NULL;
    
    // Convert arrays to well-behaved arrays of correct type and verify
    // These arrays MUST be decrefed before return.
    dataArry = NA_InputArray(dataObj, tFloat32, NUM_C_ARRAY);
    if (maskObj == Py_None) {
        maskArry = NULL;
    } else {
        maskArry = NA_InputArray(maskObj, tBool, NUM_C_ARRAY);
    }
    meanArry = NA_OutputArray(meanObj, tFloat64, NUM_C_ARRAY);
    varArry = NA_OutputArray(varObj, tFloat64, NUM_C_ARRAY);
    nPtsArry = NA_OutputArray(nPtsObj, tInt32, NUM_C_ARRAY);
    
    radSq = rad*rad;

    // Check the input arrays
    if (!dataArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting data input", ModName);
        goto errorExit;
    }
    if (!maskArry && maskObj != Py_None) {
        PyErr_Format(PyExc_ValueError, "%s: error converting mask input", ModName);
        goto errorExit;
    }
    if (dataArry->nd != 2) {
        PyErr_Format(PyExc_ValueError, "%s: data must be 2-dimensional", ModName);
        goto errorExit;
    }
    if (maskArry && !NA_ShapeEqual(dataArry, maskArry)) {
        PyErr_Format(PyExc_ValueError, "%s: mask must be the same shape as data", ModName);
        goto errorExit;
    }

    // Check output arrays and compute outLen
    if (!meanArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting mean input", ModName);
        goto errorExit;
    }
    if (!varArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting variance input", ModName);
        goto errorExit;
    }
    if (!nPtsArry) {
        PyErr_Format(PyExc_ValueError, "%s: error converting nPts input", ModName);
        goto errorExit;
    }
    if (meanArry->nd != 1) {
        PyErr_Format(PyExc_ValueError, "%s: mean must be 1-dimensional", ModName);
        goto errorExit;
    }
    if (varArry->nd != 1) {
        PyErr_Format(PyExc_ValueError, "%s: variance must be 1-dimensional", ModName);
        goto errorExit;
    }
    if (nPtsArry->nd != 1) {
        PyErr_Format(PyExc_ValueError, "%s: nPts must be 1-dimensional", ModName);
        goto errorExit;
    }
    outLen = meanArry->dimensions[0];
    if (outLen != varArry->dimensions[0]) {
        PyErr_Format(PyExc_ValueError, "%s: variance length != mean", ModName);
        goto errorExit;
    }
    if (outLen != nPtsArry->dimensions[0]) {
        PyErr_Format(PyExc_ValueError, "%s: nPts length != mean", ModName);
        goto errorExit;
    }
    if (outLen < radSq + 1) {
        PyErr_Format(PyExc_ValueError, "%s: output arrays are too short", ModName);
        goto errorExit;
    }
    
    // Call the C code
    totPts = radSqProf(
        dataArry->dimensions[0], dataArry->dimensions[1],
        NA_OFFSETDATA(dataArry),
        maskArry? NA_OFFSETDATA(maskArry): NULL,
        iCtr, jCtr,
        rad,
        outLen,
        NA_OFFSETDATA(meanArry),
        NA_OFFSETDATA(varArry),
        NA_OFFSETDATA(nPtsArry),
        &totCounts
    );
    if (totPts < 0) {
        PyErr_Format(PyExc_ValueError, "radSqProf failed");
        goto errorExit;
    }

    // Done with all arrays, decref them
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);
    Py_XDECREF(meanArry);
    Py_XDECREF(varArry);
    Py_XDECREF(nPtsArry);

    return Py_BuildValue("dl", totCounts, totPts);

errorExit:
    Py_XDECREF(dataArry);
    Py_XDECREF(maskArry);
    Py_XDECREF(meanArry);
    Py_XDECREF(varArry);
    Py_XDECREF(nPtsArry);
    return NULL;
}


/* g_radProf_setup ============================================================

Set up the global arrays used by radProf.

Inputs:
- nElt  the numer of points for the arrays = maxRad^2 + 1

If the arrays already have enough elements, leaves them alone.
Else deallocates them, allocates anew and fills with new values.
If there is not sufficient memory, deallocates all of them.

Returns 1 on success, 0 on failure (insufficient memory).
*/
long g_radProf_setup(
    long rad
) {
    long nElt, radSq;
    
    // compute nElt and make sure it is long enough for the initialization code
    nElt = MAX(rad*rad + 1, 3);

    if (g_radProf_nElt >= nElt) {
        // array is already long enough; bail out.
        return 1;
    }
    
    g_radProf_free();
    
    g_radProf_radIndByRadSq = calloc(nElt, sizeof *g_radProf_radIndByRadSq);
    if (g_radProf_radIndByRadSq == NULL) {
        g_radProf_free();
        return 0;
    }
    g_radProf_nElt = nElt;
    
    for (radSq=0; radSq<3; radSq++) {
        g_radProf_radIndByRadSq[radSq] = radSq;
    }
    for (radSq=3; radSq<nElt; radSq++) {
        g_radProf_radIndByRadSq[radSq] = (int)(sqrt((Float64)(radSq)) + 1.5);
    }
    
    return 1;
}

/* g_radProf_free ============================================================

Free the global arrays for radProf
*/
void g_radProf_free() {
    free(g_radProf_radIndByRadSq);
    g_radProf_nElt = 0;
}

/* g_radAsymm_alloc ============================================================

Allocate or reallocate the global arrays for radAsymm

Inputs:
- rad   the desired radius (each array must have rad+2 elements)

If the arrays already have enough elements, leaves them alone.
Else deallocates them and allocates anew.
If there is not sufficient memory, deallocates all of them.

Returns 1 on success, 0 on failure.
*/
long g_radAsymm_alloc(
    long rad
) {
    long nElt = rad + 2;

    if (g_radAsymm_nElt >= nElt) {
        return 1;
    }
    
    g_radAsymm_free();
    
    g_radAsymm_mean = calloc(nElt, sizeof *g_radAsymm_mean);
    g_radAsymm_var = calloc(nElt, sizeof *g_radAsymm_var);
    g_radAsymm_nPts = calloc(nElt, sizeof *g_radAsymm_nPts);
    if (g_radAsymm_mean == NULL || g_radAsymm_var == NULL || g_radAsymm_nPts == NULL) {
        g_radAsymm_free();
        return 0;
    }
    g_radAsymm_nElt = nElt;
    return 1;
}

/* g_radAsymm_free ============================================================

Free the global arrays for radSqProf.
*/
void g_radAsymm_free() {
    free(g_radAsymm_mean);
    free(g_radAsymm_var);
    free(g_radAsymm_nPts);
    g_radAsymm_nElt = 0;
}

/* radAsymm ============================================================

Compute a measure of radial asymmetry:
sum over rad of var(rad)^2 * nPts(rad).

Inputs:
- inLenI, inLenJ    dimensions of data and mask
- data              data array [i,j]
- mask              mask array [i,j] (NULL if none);
                    0 for valid values, 1 for values to ignore
- iCtr, jCtr        i,j center of profile
- rad               radius of profile

Outputs:
- asymm             radial asymmetry (see above)
- totCounts         the total # of counts (floating point to avoid overflow)

Returns:
- totPts            the total # of points (sum of nPts); <0 on error

Error Conditions:
- If insufficient memory to generate a working array, returns -2.
- Any other negative return value indicates a bug.

Points off the data array are ignored.
Thus the center need not be on the array.
*/
long radAsymm(
    long inLenI, long inLenJ,
    Float32 data[inLenI][inLenJ],
    Bool mask[inLenI][inLenJ],
    long iCtr, long jCtr,
    long rad,
    Float64 *asymmPtr,
    Float64 *totCountsPtr
) {
    long nElt = rad + 2;
    long ind;
    long totPts;

    // initialize outputs
    *asymmPtr = 0.0;
    *totCountsPtr = 0.0;
    
    // reallocate working arrays if necessary
    if (!g_radAsymm_alloc(nElt)) {
        return -2;
    }
    
    // compute radial profile stats
    totPts = radProf (
        inLenI, inLenJ,
        data,
        mask,
        iCtr, jCtr,
        rad,
        nElt,
        g_radAsymm_mean,
        g_radAsymm_var,
        g_radAsymm_nPts,
        totCountsPtr
    );
    if (totPts <= 0) {
        // <0 indicates a problem, 0 indicates no valid points
        return totPts;
    }
    
    // asymm = sum(std dev^2)
    for(ind=0; ind<nElt; ind++){
        *asymmPtr += g_radAsymm_var[ind] * (Float64) g_radAsymm_nPts[ind];
    }
        
    return totPts;
}


/* radAsymmWeighted =====================================================

Compute a weighted measure of radial asymmetry:
  sum over rad of var(rad)^2 / weight(rad)
where weight is the expected sigma of var(rad) due to pixel noise:
  weight(rad) = pixNoise(rad) * sqrt(2(numPix(rad) - 1))/numPix(rad)
  pixNoise(rad) = sqrt((readNoise/ccdGain)^2 + (meanVal(rad)-bias)/ccdGain)

Inputs:
- inLenI, inLenJ    dimensions of data and mask
- data              data array [i,j]
- mask              mask array [i,j] (NULL if none);
                    0 for valid values, 1 for values to ignore
- iCtr, jCtr        i,j center of profile
- rad               radius of profile
- readNoise         read noise in e-
- ccdGain           ccd inverse gain in e-/ADU
- bias              ccd bias in ADU

Outputs:
- asymm             radial asymmetry (see above)
- totCounts         the total # of counts (floating point to avoid overflow)

Returns:
- totPts            the total # of points (sum of nPts); <0 on error

Note:
- asymm does not include contributions where nPts(rad) < 1
  but totPts and totCounts *do* include such data.

Error Conditions:
- If insufficient memory to generate a working array, returns -2.
- Any other negative return value indicates a bug.

Points off the data array are ignored.
Thus the center need not be on the array.
*/
long radAsymmWeighted(
    long inLenI, long inLenJ,
    Float32 data[inLenI][inLenJ],
    Bool mask[inLenI][inLenJ],
    long iCtr, long jCtr,
    long rad,
    long bias,
    Float64 readNoise,
    Float64 ccdGain,
    Float64 *asymmPtr,
    Float64 *totCountsPtr
) {
    long nElt = rad + 2;
    long ind;
    long totPts;
    long nPts;
    Float64 readNoiseSqADU = (readNoise * readNoise) / (ccdGain * ccdGain);
    Float64 pixNoiseSq;
    Float64 weight;

    // initialize outputs
    *asymmPtr = 0.0;
    *totCountsPtr = 0.0;
    
    // reallocate working arrays if necessary
    if (!g_radAsymm_alloc(nElt)) {
        return -2;
    }
    
    // compute radial profile stats
    totPts = radProf (
        inLenI, inLenJ,
        data,
        mask,
        iCtr, jCtr,
        rad,
        nElt,
        g_radAsymm_mean,
        g_radAsymm_var,
        g_radAsymm_nPts,
        totCountsPtr
    );
    if (totPts <= 0) {
        // <0 indicates a problem, 0 indicates no valid points
        return totPts;
    }
    
    // asymm = sum(std dev^2)
    for(ind=0; ind<nElt; ind++){
        nPts = g_radAsymm_nPts[ind];
        if (nPts > 1) {
            pixNoiseSq = readNoiseSqADU + ((g_radAsymm_mean[ind] - (Float64) bias) / ccdGain);
            weight = sqrt(2.0 * (Float64) (nPts - 1)) * pixNoiseSq / (Float64) nPts;
            *asymmPtr += g_radAsymm_var[ind] / weight;
        }
    }
        
    return totPts;
}


/* radProf ============================================================

Generate a radial profile as a function of radial index
(an approximation of radius; see Py_radProf_doc for details).

Inputs:
- inLenI, inLenJ    dimensions of data and mask
- data              data array [i,j]
- mask              mask array [i,j] (NULL if none);
                    0 for valid values, 1 for values to ignore
- iCtr, jCtr        i,j center of profile
- rad               radius of profile
- outLen            length of output arrays

Outputs:
- mean              the mean at each radius squared; 0 if npts=0
- var               the variance (stdDev^2) at each radius squared; 0 if npts=0
- nPts              the # of points at each radius squared
- totCounts         the total # of counts (floating point to avoid overflow)

Returns:
- totPts            the total # of points (sum of nPts); <0 on error

Error Conditions:
- If any output arrays have fewer than outLen points,
  writes off the end of an array.

- If outLen < rad + 2, returns -1.
- If insufficient memory to generate a working array, returns -2.
- If any value in g_radProf_radIndByRadSq[0:rad^2] > rad, returns -3. 

Points off the data array are ignored.
Thus the center need not be on the array.
*/
long radProf(
    long inLenI, long inLenJ,
    Float32 data[inLenI][inLenJ],
    Bool mask[inLenI][inLenJ],
    long iCtr, long jCtr,
    long rad,
    long outLen,
    Float64 *mean,
    Float64 *var,
    Int32 *nPts,
    Float64 *totCountsPtr
) {
    long desOutLen = rad + 2;
    long maxRadSq = rad*rad;
    long jj, ii, currRadSq, outInd;
    long minJJ, maxJJ, minII, maxII;
    long totPts;
    Float64 d;
    char ModName[]="radProf";
    
    // test inputs
    if (outLen < desOutLen) {
        printf("%s: outLen too small\n", ModName);
        return -1;
    }
    
    // set up index array
    if (!g_radProf_setup(rad)) {
        printf("%s: insufficient memory\n", ModName);
        return -2;
    }

    // initialize outputs to 0
    totPts = 0;
    for(outInd=0; outInd<outLen; outInd++){
        nPts[outInd] = 0;
        mean[outInd] = 0.0;
        var[outInd] = 0.0;
    }
    *totCountsPtr = 0;

    // compute sums
    minII = MAX(iCtr - rad, 0);
    minJJ = MAX(jCtr - rad, 0);
    maxII = MIN(iCtr + rad, inLenI - 1);
    maxJJ = MIN(jCtr + rad, inLenJ - 1);
    for (ii=minII; ii<=maxII; ii++) {
        for (jj=minJJ; jj<=maxJJ; jj++) {
            if (mask==NULL || !mask[ii][jj]) {
                currRadSq = (ii - iCtr)*(ii - iCtr) + (jj - jCtr)*(jj - jCtr);
                if (currRadSq > maxRadSq)
                    continue;
                outInd = g_radProf_radIndByRadSq[currRadSq];
                if (outInd >= desOutLen) {
                    printf("radProf failed: outInd=%ld, rad=%ld\n", outInd, rad);
                    return -3;
                }
    
                d = (Float64) data[ii][jj];
                mean[outInd] += d;
                var[outInd] += d*d;
                nPts[outInd]++;
                *totCountsPtr += d;
                totPts++;
            }
        }
    }

    /* normalize outputs */
    for(outInd=0; outInd<desOutLen; outInd++) {
        if (nPts[outInd] != 0) {
            mean[outInd] /= nPts[outInd];
            var[outInd] = (var[outInd]/(Float64)nPts[outInd]) - (mean[outInd]*mean[outInd]);
        }
    }
    
    return totPts;
}


/* radSqProf ============================================================

Generate a radial profile as a function of radius squared.
    
Inputs:
- inLenI, inLenJ    dimensions of data and mask
- data              data array [i,j]
- mask              mask array [i,j] (NULL if none);
                    0 for valid values, 1 for values to ignore
- iCtr, jCtr        i,j center of profile
- rad               radius of profile
- outLen            length of output arrays

Outputs:
- mean              the mean at each radius squared; 0 if npts=0
- var               the variance (stdDev^2) at each radius squared; 0 if npts=0
- nPts              the # of points at each radius squared
- totCounts         the total # of counts (floating point to avoid overflow)

Returns:
- totPts            the total # of points (sum of nPts); <0 on error

If any output arrays have fewer than outLen points,
writes off the end of an array.

If outLen < radSq**2 + 1, returns -1.

Points off the data array are ignored.
Thus the center need not be on the array.
*/
long radSqProf(
    long inLenI, long inLenJ,
    Float32 data[inLenI][inLenJ],
    Bool mask[inLenI][inLenJ],
    long iCtr, long jCtr,
    long rad,
    long outLen,
    Float64 *mean,
    Float64 *var,
    Int32 *nPts,
    Float64 *totCountsPtr
) {
    long desOutLen = rad*rad + 1;
    long jj, ii, outInd;
    long minJJ, maxJJ, minII, maxII;
    Float64 d;
    long totPts;
    
    // test inputs
    if (outLen < desOutLen) {
        return -1;
    }
    
    // initialize outputs to 0
    for(outInd=0; outInd<outLen; outInd++){
        nPts[outInd] = 0;
        mean[outInd] = 0.0;
        var[outInd] = 0.0;
    }
    *totCountsPtr = 0.0;
    totPts = 0;

    // compute sums
    minII = MAX(iCtr - rad, 0);
    minJJ = MAX(jCtr - rad, 0);
    maxII = MIN(iCtr + rad, inLenI - 1);
    maxJJ = MIN(jCtr + rad, inLenJ - 1);
    for (ii=minII; ii<=maxII; ii++) {
        for (jj=minJJ; jj<=maxJJ; jj++) {
            if (mask==NULL || !mask[ii][jj]) {
                outInd = (ii - iCtr)*(ii - iCtr) + (jj - jCtr)*(jj - jCtr);
                if (outInd >= desOutLen)
                    continue;
    
                d = (Float64) data[ii][jj];
                mean[outInd] += d;
                var[outInd] += d*d;
                nPts[outInd]++;
                totPts++;
                *totCountsPtr += d;
            }
        }
    }

    /* normalize outputs */
    for(outInd=0; outInd<desOutLen; outInd++) {
        if (nPts[outInd] != 0) {
            mean[outInd] /= nPts[outInd];
            var[outInd] = (var[outInd]/(Float64)nPts[outInd]) - (mean[outInd]*mean[outInd]);
        }
    }
    return totPts;
}


static PyMethodDef radProfMethods[] = {
    {"radAsymm", Py_radAsymm, METH_VARARGS, Py_radAsymm_doc},
    {"radAsymmWeighted", Py_radAsymmWeighted, METH_VARARGS, Py_radAsymmWeighted_doc},
    {"radProf", Py_radProf, METH_VARARGS, Py_radProf_doc},
    {"radIndByRadSq", Py_radIndByRadSq, METH_VARARGS, Py_radIndByRadSq_doc},
    {"radSqByRadInd", Py_radSqByRadInd, METH_VARARGS, Py_radSqByRadInd_doc},
    {"radSqProf", Py_radSqProf, METH_VARARGS, Py_radSqProf_doc},
    {NULL, NULL} /* Sentinel */
};

// Module initialization function
PyMODINIT_FUNC initradProf(void) {
    PyObject *m;
    m = Py_InitModule3("radProf", radProfMethods, radProfModule_doc);
    
    if (m == NULL)
        return;

    import_libnumarray();
}

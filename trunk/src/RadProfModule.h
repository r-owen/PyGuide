#ifndef __radprofmodule__
#define __radprofmodule__

/*
Header file for RadProfModule.c, which see for information.

History:
2004-10-13 ROwen    modified libnumarray.h include to match numarray 1.1 docs.
2008-10-01 ROwen    radAsymmWeighted: changed bias from int to double.
2009-11-19 ROwen    Modified to use numpy instead of numarray.
*/

#include "Python.h"
#include "numpy/arrayobject.h"

#ifdef __cplusplus
extern "C" {
#endif

// routines visible to Python
static PyObject *Py_radAsymm(PyObject *dumObj, PyObject *args);
static PyObject *Py_radProf(PyObject *dumObj, PyObject *args);
static PyObject *Py_radIndByRadSq(PyObject *dumObj, PyObject *args);
static PyObject *Py_radSqByRadInd(PyObject *dumObj, PyObject *args);
static PyObject *Py_radSqProf(PyObject *dumObj, PyObject *args);

// internal routines
int g_radProf_setup(
    int rad
);
void g_radProf_free(
    void
);
int g_radAsymm_alloc(
    int nElt
);
void g_radAsymm_free(
    void
);
int radAsymm(
    int inLenI, int inLenJ,
    npy_float32 data[inLenI][inLenJ],
    npy_bool mask[inLenI][inLenJ],
    int iCtr, int jCtr,
    int rad,
    double *asymmPtr,
    double *totCountsPtr
);
int radAsymmWeighted(
    int inLenI, int inLenJ,
    npy_float32 data[inLenI][inLenJ],
    npy_bool mask[inLenI][inLenJ],
    int iCtr, int jCtr,
    int rad,
    double bias,
    double readNoise,
    double ccdGain,
    double *asymmPtr,
    double *totCountsPtr
);
int radProf(
    int inLenI, int inLenJ,
    npy_float32 data[inLenI][inLenJ],
    npy_bool mask[inLenI][inLenJ],
    int iCtr, int jCtr,
    int rad,
    int outLen,
    double *mean,
    double *var,
    long *nPts,
    double *totCountsPtr
);
int radSqProf(
    int inLenI, int inLenJ,
    npy_float32 data[inLenI][inLenJ],
    npy_bool mask[inLenI][inLenJ],
    int iCtr, int jCtr,
    int rad,
    int outLen,
    double *mean,
    double *var,
    long *nPts,
    double *totCountsPtr
);

#ifdef __cplusplus
}
#endif

#endif

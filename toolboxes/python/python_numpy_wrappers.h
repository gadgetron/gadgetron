#ifndef GADGETRON_PYTHON_NUMPY_WRAPPERS_H
#define GADGETRON_PYTHON_NUMPY_WRAPPERS_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/ndarraytypes.h"

namespace Gadgetron {

/// Wrappers for NumPy C-API functions. These functions must be used
/// in the same C++ source file as the call to `import_array()`. In this
/// case, that is Python.cpp. The simplest solution is to lightly wrap the API.
EXPORTPYTHON int NumPyArray_NDIM(PyObject* obj);
EXPORTPYTHON npy_intp NumPyArray_DIM(PyObject* obj, int i);
EXPORTPYTHON void *NumPyArray_DATA(PyObject* obj);
EXPORTPYTHON int NumPyArray_ITEMSIZE(PyObject* obj);
EXPORTPYTHON PyObject *NumPyArray_SimpleNew(int nd, npy_intp* dims, int typenum);

}

#endif // GADGETRON_PYTHON_NUMPY_WRAPPERS_H
